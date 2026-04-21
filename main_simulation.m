clear; clc; close all;

busdata0   = readmatrix('busdata.xlsx');
linedata0  = readmatrix('linedata.xlsx');
nBuses     = size(busdata0, 1);
basemva    = 100;
linedata0(:,7) = (1:size(linedata0,1))';  
Vmax       = 1.05;   
Vmin       = 0.95;   
maxCascade = 10;     
slackBusID = 1;      
dt         = 0.5;    
tau       = 1.1; 

json_raw = load_substation_json();  
data     = jsondecode(strjoin(json_raw, ''));
n        = length(data);
lats     = zeros(n,1);
lons     = zeros(n,1);
for i = 1:n
    item = data{i};
    if isfield(item, 'lat')
        lats(i) = item.lat;
        lons(i) = item.lon;
    else
        lats(i) = item.center.lat;
        lons(i) = item.center.lon;
    end
end

json_to_busID = [14,23,7,11,24,9,20,29,1,13,26,27,30,6,17,22,15,5,10,16,18,28,2,4,19,8,12,3,25,21];



                 
target_date    = datetime(2015,12,26);
time_profile   = load_demand_profile('demand2015.xlsx', target_date);
apx_price      = load_apx_prices('marketindexprice2015.xlsx', target_date);
[national_CCGT, national_COAL, national_NUCLEAR] = ...
    load_generation_mix('generationmix2015.xlsx', target_date);
gen_bus_indices = [1, 2, 13, 22, 23, 27];
gen_fuel_types  = {'CCGT','CCGT','COAL','COAL','NUCLEAR','NUCLEAR'};
heat_rate       = struct('CCGT',0.9, 'COAL',0.6, 'NUCLEAR',0.22);
ngen_buses      = length(gen_bus_indices);
nT              = 48;
national_total    = national_CCGT + national_COAL + national_NUCLEAR;
national_total(national_total == 0) = 1; 

gen_marginal_costs_ts = zeros(nT, ngen_buses);
gen_Pmax_ts           = zeros(nT, ngen_buses);
for t = 1:nT
    spot = max(apx_price(t), 20);   % floor price at £20/MWh
    frac = struct('CCGT',  national_CCGT(t)   / national_total(t), ...
                  'COAL',  national_COAL(t)   / national_total(t), ...
                  'NUCLEAR', national_NUCLEAR(t) / national_total(t));
    for g = 1:ngen_buses
        fuel = gen_fuel_types{g};
        gen_marginal_costs_ts(t,g) = spot * heat_rate.(fuel);
        gen_Pmax_ts(t,g)           = (frac.(fuel) * local_system_peak) / 2;
    end
end




sector_emergency   = [1, 12, 22];
sector_commercial  = [2, 3, 6, 10, 15, 18, 19, 24, 25];
sector_residential = [4, 5, 7, 11, 17, 20, 23, 27, 28, 29, 30];
sector_public      = [8, 9, 13, 14, 16, 21, 26];

VOLL_EMERGENCY   = 50000;  
VOLL_COMMERCIAL  = 15000;
VOLL_RESIDENTIAL = 8000;
VOLL_PUBLIC      = 2000;

VOLL_static = zeros(nBuses, 1);
for i = 1:nBuses
    bid = busdata0(i,1);
    if     ismember(bid, sector_emergency),   VOLL_static(i) = VOLL_EMERGENCY;
    elseif ismember(bid, sector_commercial),  VOLL_static(i) = VOLL_COMMERCIAL;
    elseif ismember(bid, sector_residential), VOLL_static(i) = VOLL_RESIDENTIAL;
    elseif ismember(bid, sector_public),      VOLL_static(i) = VOLL_PUBLIC;
    else
        error('Bus ID %d not assigned to a VoLL sector.', bid);
    end
end




depth_matrix = readmatrix('York_Flood_Depth.xlsx', 'Range', 'B36');
d0     = 0.10678;   
d50    = 1.1441;    
beta_f = 0.4381;    

Pf        = zeros(size(depth_matrix));
mask      = depth_matrix > d0;
Pf(mask)  = normcdf((log(depth_matrix(mask)) - log(d50)) / beta_f);

Nmc = 1000;
rng(42);  
failure_tensor = rand(size(depth_matrix,1), nT, Nmc) < Pf;



[critical_proximity, density_norm] = compute_criticality_spatial(lats, lons, n);
alpha1 = 0.4;   
alpha2 = 0.3;   
alpha3 = 0.3;   

C_matrix   = zeros(n, nT);
d_max_all  = max(depth_matrix(:));
for tt = 1:nT
    C_matrix(:,tt) = alpha1 * (depth_matrix(:,tt) / d_max_all) +  alpha2 * critical_proximity +  alpha3 * density_norm;
end

C_byBusID = zeros(nBuses, nT);
for i = 1:n
    row = busdata0(:,1) == json_to_busID(i);
    if any(row)
        C_byBusID(row,:) = C_matrix(i,:);
    end
end

nT_sim = length(time_profile);
fields = {'TotalLoad','ENS','CascadeSteps','LinesTripped','Islands', 'TotalCost','GenCost','Rcrit'};


opf = struct();
heu = struct();
for f = 1:length(fields)
    opf.(fields{f}) = zeros(nT_sim, Nmc);
    heu.(fields{f}) = zeros(nT_sim, Nmc);
end



fprintf('Starting simulation: %d MC realisations x %d timesteps\n', Nmc, nT_sim);
tic;

for mc = 1:Nmc
    if mod(mc,100) == 0
        fprintf('  MC iteration %d / %d  (%.1f min elapsed)\n', mc, Nmc, toc/60);
    end

    for t_idx = 1:nT_sim
        scale_factor      = time_profile(t_idx);
        current_gen_costs = gen_marginal_costs_ts(t_idx,:);

        busdata_base      = busdata0;
        busdata_base(:,5) = busdata_base(:,5) * scale_factor;
        busdata_base(:,6) = busdata_base(:,6) * scale_factor;
        for g = 1:ngen_buses
            row_g = busdata_base(:,1) == gen_bus_indices(g);
            if any(row_g)
                busdata_base(row_g,9) = gen_Pmax_ts(t_idx,g);
            end
        end

        linedata_base = linedata0;
        outages = find(failure_tensor(:,t_idx,mc))';
        if ~isempty(outages)
            busdata_base(ismember(busdata_base(:,1),outages),:) = [];
            linedata_base(~ismember(linedata_base(:,1),busdata_base(:,1)) | ...
                          ~ismember(linedata_base(:,2),busdata_base(:,1)),:) = [];
        end

        for scenario = 1:2
            busdata  = busdata_base;
            linedata = linedata_base;

            Pg_final        = zeros(ngen_buses,1);
            stable          = false;
            iter            = 0;
            linesTrippedThis = 0;
            islandsThis     = 0;

            while ~stable && iter < maxCascade
                iter = iter + 1;

                busIDs   = busdata(:,1);
                nbus     = length(busIDs);
                mapGraph = containers.Map(busIDs, 1:nbus);
                Adj      = zeros(nbus);
                for k = 1:size(linedata,1)
                    if isKey(mapGraph,linedata(k,1)) && isKey(mapGraph,linedata(k,2))
                        ii = mapGraph(linedata(k,1));
                        jj = mapGraph(linedata(k,2));
                        Adj(ii,jj) = 1;
                        Adj(jj,ii) = 1;
                    end
                end
                comps = conncomp(graph(Adj));
                numC  = max(comps);

                if isKey(mapGraph, slackBusID)
                    islandsThis = islandsThis + (numC - 1);
                else
                    islandsThis = islandsThis + numC;
                end

                for c = 1:numC
                    members = busIDs(comps == c);
                    if ~ismember(slackBusID, members)
                        busdata(ismember(busdata(:,1),members),:) = [];
                        linedata(~ismember(linedata(:,1),busdata(:,1)) | ...
                                 ~ismember(linedata(:,2),busdata(:,1)),:) = [];
                    end
                end

                originalBusIDs = busdata(:,1);
                nbus_current   = length(originalBusIDs);
                if nbus_current == 0, break; end
                busMap = containers.Map(originalBusIDs, 1:nbus_current);

                C_current    = zeros(nbus_current,1);
                VOLL_current = zeros(nbus_current,1);
                for i = 1:nbus_current
                    bid  = originalBusIDs(i);
                    brow = busdata0(:,1) == bid;
                    C_current(i)    = C_byBusID(brow, t_idx);
                    VOLL_current(i) = VOLL_static(brow);
                end

                busdata_s  = busdata;  busdata_s(:,1)  = (1:nbus_current)';
                linedata_s = linedata;
                for k = 1:size(linedata,1)
                    linedata_s(k,1) = busMap(linedata(k,1));
                    linedata_s(k,2) = busMap(linedata(k,2));
                end

                [busdata_s, V] = solver_NR(busdata_s, linedata_s);
                V = V(:);
                if any(isnan(V)), break; end

                nbranch = size(linedata_s,1);
                Sline   = zeros(nbranch,1);
                for k = 1:nbranch
                    fb  = linedata_s(k,1);
                    tb  = linedata_s(k,2);
                    Z   = linedata_s(k,3) + 1j*linedata_s(k,4);
                    Bsh = 1j*linedata_s(k,5)/2;
                    Iij = (V(fb)-V(tb))/Z + V(fb)*Bsh;
                    Sline(k) = V(fb)*conj(Iij);
                end
                MVA_line = abs(Sline) * basemva;
                RateA    = linedata(:,6) * tau;   
                Vmag     = abs(V);

                overloaded = find(MVA_line > RateA);
                overV      = find(Vmag > Vmax);
                underV     = find(Vmag < Vmin);

                if isempty(overloaded) && isempty(overV) && isempty(underV)
                    stable = true;
                    for g = 1:ngen_buses
                        if isKey(busMap, gen_bus_indices(g))
                            Pg_final(g) = busdata_s(busMap(gen_bus_indices(g)),7);
                        end
                    end

                elseif scenario == 1
                    [Pg_opt, alpha_opt, ~, ~, exitflag] = ...
                        solve_DC_OPF(busdata, linedata, busMap, gen_bus_indices, current_gen_costs, VOLL_current, C_current, basemva, 10);

                    if exitflag == 1
                        pg_ptr = 1;
                        for g = 1:ngen_buses
                            tid   = gen_bus_indices(g);
                            if isKey(busMap, tid)
                                row_g = find(busdata(:,1) == tid);
                                if ~isempty(row_g) && pg_ptr <= length(Pg_opt)
                                    busdata(row_g,7) = Pg_opt(pg_ptr);
                                    Pg_final(g)      = Pg_opt(pg_ptr);
                                    pg_ptr = pg_ptr + 1;
                                end
                            end
                        end
                        for i = 1:nbus_current
                            orig_id  = originalBusIDs(i);
                            row_cur  = busdata(:,1)  == orig_id;
                            row_orig = busdata0(:,1) == orig_id;
                            busdata(row_cur,5) = busdata0(row_orig,5) * scale_factor * alpha_opt(i);
                            busdata(row_cur,6) = busdata0(row_orig,6) * scale_factor * alpha_opt(i);
                        end
                    else
                        if ~isempty(overloaded)
                            [~,worst] = max(MVA_line ./ RateA);
                            linedata(worst,:) = [];
                            linesTrippedThis  = linesTrippedThis + 1;
                        end
                    end

                else
                    % Shed 10% uniformly, re-run NR, trip worst line if needed
                    for i = 1:nbus_current
                        busdata(i,5) = busdata(i,5) * 0.9;
                        busdata(i,6) = busdata(i,6) * 0.9;
                    end

                    busdata_s2  = busdata;  busdata_s2(:,1)  = (1:nbus_current)';
                    linedata_s2 = linedata;
                    for k = 1:size(linedata,1)
                        linedata_s2(k,1) = busMap(linedata(k,1));
                        linedata_s2(k,2) = busMap(linedata(k,2));
                    end
                    [~, V2] = solver_NR(busdata_s2, linedata_s2);
                    V2 = V2(:);

                    if ~any(isnan(V2))
                        Sline2 = zeros(nbranch,1);
                        for k = 1:nbranch
                            fb   = linedata_s2(k,1);
                            tb   = linedata_s2(k,2);
                            Z    = linedata_s2(k,3) + 1j*linedata_s2(k,4);
                            Bsh  = 1j*linedata_s2(k,5)/2;
                            Iij2 = (V2(fb)-V2(tb))/Z + V2(fb)*Bsh;
                            Sline2(k) = V2(fb)*conj(Iij2);
                        end
                        MVA2 = abs(Sline2) * basemva;
                        if ~isempty(find(MVA2 > RateA, 1))
                            [~,worst] = max(MVA2 ./ RateA);
                            linedata(worst,:)  = [];
                            linesTrippedThis   = linesTrippedThis + 1;
                        end
                    else
                        if ~isempty(overloaded)
                            [~,worst] = max(MVA_line ./ RateA);
                            linedata(worst,:)  = [];
                            linesTrippedThis   = linesTrippedThis + 1;
                        end
                    end

                    for g = 1:ngen_buses
                        if isKey(busMap, gen_bus_indices(g))
                            Pg_final(g) = busdata_s(busMap(gen_bus_indices(g)),7);
                        end
                    end
                end % scenario branch
            end % cascade loop

            fullDemand_ts = sum(busdata0(:,5)) * scale_factor;

            gen_cost_ts = 0;
            for g = 1:ngen_buses
                gen_cost_ts = gen_cost_ts + Pg_final(g) * current_gen_costs(g);
            end

            totalServed_ts    = 0;
            timestep_ENS_Cost = 0;
            w_ens_num = 0;
            w_ens_den = 0;

            for i = 1:nBuses
                orig_id     = busdata0(i,1);
                demand_orig = busdata0(busdata0(:,1)==orig_id, 5) * scale_factor;

                cur_row = busdata(:,1) == orig_id;
                if any(cur_row)
                    demand_served = busdata(cur_row,5);
                else
                    demand_served = 0;
                end

                totalServed_ts    = totalServed_ts + demand_served;
                bus_ens           = max(0, demand_orig - demand_served);
                timestep_ENS_Cost = timestep_ENS_Cost + bus_ens * dt * VOLL_static(i);

                C_i = C_byBusID(i, t_idx);
                w_ens_num = w_ens_num + C_i * bus_ens       * dt;
                w_ens_den = w_ens_den + C_i * demand_orig   * dt;


            end

            ens_ts       = fullDemand_ts - totalServed_ts;
            gen_cost_half = gen_cost_ts * dt;
            total_cost    = timestep_ENS_Cost + gen_cost_half;

            if w_ens_den > 1e-9
                rcrit = 1 - w_ens_num / w_ens_den;
            else
                rcrit = 1;
            end

            if scenario == 1, S = opf; else, S = heu; end
            S.TotalLoad(t_idx,mc)    = totalServed_ts;
            S.ENS(t_idx,mc)          = ens_ts;
            S.CascadeSteps(t_idx,mc) = iter;
            S.LinesTripped(t_idx,mc) = linesTrippedThis;
            S.Islands(t_idx,mc)      = islandsThis;
            S.GenCost(t_idx,mc)      = gen_cost_half;
            S.TotalCost(t_idx,mc)    = total_cost;
            S.Rcrit(t_idx,mc)        = rcrit;
            if scenario == 1, opf = S; else, heu = S; end

        end % scenario loop
    end % timestep loop
end % MC loop

fprintf('Simulation complete in %.1f minutes.\n', toc/60);




peakDemand    = sum(busdata0(:,5));
demandProfile = time_profile * peakDemand;
totalDemand   = sum(demandProfile) * dt;

opf_mean = struct();
heu_mean = struct();
for f = 1:length(fields)
    opf_mean.(fields{f}) = mean(opf.(fields{f}), 2);
    heu_mean.(fields{f}) = mean(heu.(fields{f}), 2);
end

opf_res = 1 - opf_mean.ENS ./ demandProfile;
heu_res = 1 - heu_mean.ENS ./ demandProfile;
opf_res(demandProfile == 0) = 1;
heu_res(demandProfile == 0) = 1;

opf_R = 1 - sum(opf_mean.ENS) * dt / totalDemand;
heu_R = 1 - sum(heu_mean.ENS) * dt / totalDemand;

opf_ENS_cost = sum(mean(opf.TotalCost - opf.GenCost, 2));
heu_ENS_cost = sum(mean(heu.TotalCost - heu.GenCost, 2));
opf_gen_cost = sum(opf_mean.GenCost);
heu_gen_cost = sum(heu_mean.GenCost);




fprintf('\n%s\n', repmat('=',1,60));
fprintf('%-22s  %10s  %10s  %10s\n','Metric','OPF','Heuristic','Improvement');
fprintf('%s\n', repmat('-',1,60));
fprintf('%-22s  %10.4f  %10.4f  %+10.4f\n','R (standard)',  opf_R, heu_R, opf_R-heu_R);
fprintf('%-22s  %10.4f  %10.4f  %+10.4f\n','R_crit (mean)',mean(opf_mean.Rcrit), mean(heu_mean.Rcrit), ...
    mean(opf_mean.Rcrit)-mean(heu_mean.Rcrit));
fprintf('%-22s  %10.0f  %10.0f  %+10.0f\n','ENS cost (£)',   opf_ENS_cost, heu_ENS_cost, heu_ENS_cost-opf_ENS_cost);
fprintf('%-22s  %10.0f  %10.0f  %+10.0f\n','Gen cost (£)',  opf_gen_cost, heu_gen_cost, heu_gen_cost-opf_gen_cost);
fprintf('%-22s  %10.0f  %10.0f  %+10.0f\n','Total cost (£)',opf_ENS_cost+opf_gen_cost, heu_ENS_cost+heu_gen_cost, (heu_ENS_cost+heu_gen_cost)-(opf_ENS_cost+opf_gen_cost));
fprintf('%s\n', repmat('=',1,60));




plot_results(opf_mean, heu_mean, opf_res, heu_res, demandProfile, time_profile, opf, heu);
plot_pareto(busdata0, linedata0, C_byBusID, gen_bus_indices, gen_marginal_costs_ts, VOLL_static, basemva, time_profile);




function json_raw = load_substation_json()
json_raw = [ "["
  "{""type"":""node"",""lat"":53.946011,""lon"":-1.054898,""tags"":{""name"":""Ridgeway Yor""}},"
  "{""type"":""node"",""lat"":53.94815,""lon"":-1.139708,""tags"":{""name"":""583738755""}},"
  "{""type"":""node"",""lat"":53.944681,""lon"":-1.032498,""tags"":{""name"":""228170123""}},"
  "{""type"":""way"",""center"":{""lat"":53.973215,""lon"":-1.070525},""tags"":{""name"":""464845082""}},"
  "{""type"":""node"",""lat"":53.9752084,""lon"":-1.1316600,""tags"":{""name"":""11413046791""}},"
  "{""type"":""node"",""lat"":53.9738289,""lon"":-1.0486766,""tags"":{""name"":""6833442632""}},"
  "{""type"":""node"",""lat"":53.943189,""lon"":-1.088544,""tags"":{""name"":""533841593""}},"
  "{""type"":""way"",""center"":{""lat"":53.985308,""lon"":-1.113994},""tags"":{""name"":""828602195""}},"
  "{""type"":""node"",""lat"":53.9581094,""lon"":-1.0198141,""tags"":{""name"":""Osbaldwick Substation""}},"
  "{""type"":""node"",""lat"":53.9532039,""lon"":-1.0747363,""tags"":{""name"":""Barbican Centre""}},"
  "{""type"":""way"",""center"":{""lat"":53.959319,""lon"":-1.0827282},""tags"":{""name"":""9656108099""}},"
  "{""type"":""way"",""center"":{""lat"":53.964886,""lon"":-1.1041897},""tags"":{""name"":""10026009297""}},"
  "{""type"":""way"",""center"":{""lat"":53.99021,""lon"":-1.065846},""tags"":{""name"":""Huntington school""}},"
  "{""type"":""way"",""center"":{""lat"":53.947454,""lon"":-1.032762},""tags"":{""name"":""TFTV substation""}},"
  "{""type"":""way"",""center"":{""lat"":53.9246974,""lon"":-1.0946668},""tags"":{""name"":""6845278898""}},"
  "{""type"":""node"",""lat"":53.9536498,""lon"":-1.0818034,""tags"":{""name"":""Bishops Wharf""}},"
  "{""type"":""way"",""center"":{""lat"":53.9392038,""lon"":-1.0615282},""tags"":{""name"":""6845249758""}},"
  "{""type"":""way"",""center"":{""lat"":53.9651497,""lon"":-1.041832},""tags"":{""name"":""6845230898""}},"
  "{""type"":""way"",""center"":{""lat"":53.977455,""lon"":-1.069089},""tags"":{""name"":""10007239969""}},"
  "{""type"":""way"",""center"":{""lat"":53.9368354,""lon"":-1.0673886},""tags"":{""name"":""6845249759""}},"
  "{""type"":""way"",""center"":{""lat"":53.936055,""lon"":-1.120184},""tags"":{""name"":""573761380""}},"
  "{""type"":""way"",""center"":{""lat"":53.99047,""lon"":-1.096466},""tags"":{""name"":""762926050""}},"
  "{""type"":""way"",""center"":{""lat"":53.9565544,""lon"":-1.0209307},""tags"":{""name"":""9051665121""}},"
  "{""type"":""way"",""center"":{""lat"":53.9675889,""lon"":-1.0426743},""tags"":{""name"":""6845230899""}},"
  "{""type"":""way"",""center"":{""lat"":53.9408454,""lon"":-1.091714},""tags"":{""name"":""6266581976""}},"
  "{""type"":""way"",""center"":{""lat"":53.9746775,""lon"":-1.0425614},""tags"":{""name"":""3593578896""}},"
  "{""type"":""way"",""center"":{""lat"":53.9573467,""lon"":-1.0741447},""tags"":{""name"":""Student Castle""}},"
  "{""type"":""way"",""center"":{""lat"":53.951138,""lon"":-1.018246},""tags"":{""name"":""228170156""}},"
  "{""type"":""way"",""center"":{""lat"":53.9673365,""lon"":-1.128289},""tags"":{""name"":""6833356670""}},"
  "{""type"":""way"",""center"":{""lat"":53.949748,""lon"":-1.108713},""tags"":{""name"":""636376649""}}"
  "]" ];
end
