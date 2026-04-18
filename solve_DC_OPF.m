function [Pg_opt, alpha_opt, f_econ_val, f_crit_val, exitflag, ...
          pareto_econ, pareto_crit] = ...
    solve_DC_OPF(busdata, linedata, busMap, gen_bus_indices, ...
                 gen_costs, VOLL_vec, C_vec, basemva, n_pareto)
%SOLVE_DC_OPF  Bi-objective DC Optimal Power Flow via epsilon-constraint.
%
%   Minimises two competing objectives over generation dispatch (Pg),
%   load retention fraction (alpha), and voltage angles (theta):
%
%   f_econ = sum_g(C_g * Pg) + sum_i(VoLL_i * (1-alpha_i) * Pd_i)   [eq.12]
%   f_crit = sum_i(VoLL_i * C_i * (1-alpha_i)^2 * Pd_i)              [eq.17]
%
%   Both objectives are in GBP (£). The epsilon-constraint method sweeps
%   f_crit from its minimum to its value at the economic optimum, solving
%   a sequence of single-objective subproblems to trace the Pareto front.
%   The knee-point (minimum normalised distance to ideal) is returned.
%
%   INPUTS:
%     busdata          Current surviving bus data matrix
%     linedata         Current surviving line data matrix
%     busMap           containers.Map: original bus ID -> local index
%     gen_bus_indices  Generator bus IDs [1 x ngen]
%     gen_costs        Marginal generation costs [1 x ngen], £/MWh
%     VOLL_vec         Value of Lost Load per bus [nbus x 1], £/MWh
%     C_vec            Criticality index per bus [nbus x 1], in [0,1]
%     basemva          System base MVA (default 100)
%     n_pareto         Number of Pareto points (default 20)
%
%   OUTPUTS:
%     Pg_opt           Optimal generator dispatch [ngen x 1], MW
%     alpha_opt        Optimal load retention fractions [nbus x 1]
%     f_econ_val       Economic objective at knee-point, £
%     f_crit_val       Criticality objective at knee-point, £
%     exitflag         1 = success, -1 = infeasible
%     pareto_econ      Economic objective values on Pareto front
%     pareto_crit      Criticality objective values on Pareto front

if nargin < 9, n_pareto = 20; end


originalBusIDs   = busdata(:,1);
nbus             = length(originalBusIDs);
nlines           = size(linedata,1);
active_gen_mask  = ismember(gen_bus_indices, originalBusIDs);
active_gen_ids   = gen_bus_indices(active_gen_mask);
active_gen_costs = gen_costs(active_gen_mask);
ngen             = length(active_gen_ids);

% Decision variable layout: [Pg(ngen) | alpha(nbus) | theta(nbus)]
theta_offset = ngen + nbus;
nvars        = ngen + nbus + nbus;
Bbus = zeros(nbus);
for k = 1:nlines
    ii = busMap(linedata(k,1));
    jj = busMap(linedata(k,2));
    b  = 1 / linedata(k,4);   % susceptance = 1/reactance
    Bbus(ii,ii) = Bbus(ii,ii) + b;
    Bbus(jj,jj) = Bbus(jj,jj) + b;
    Bbus(ii,jj) = Bbus(ii,jj) - b;
    Bbus(jj,ii) = Bbus(jj,ii) - b;
end

%% --- Equality constraints: power balance + slack reference ---
Aeq  = zeros(nbus+1, nvars);
beq  = zeros(nbus+1, 1);
Pd_all = busdata(:,5);

for i = 1:nbus
    bus_id = originalBusIDs(i);
    g_idx  = find(active_gen_ids == bus_id);
    if ~isempty(g_idx)
        Aeq(i, g_idx) = 1 / basemva;
    end
    Aeq(i, ngen+i)                     = -Pd_all(i) / basemva;
    Aeq(i, theta_offset+(1:nbus))      = -Bbus(i,:);
end
Aeq(nbus+1, theta_offset+1) = 1;   % slack bus angle = 0

%% --- Inequality constraints: thermal line limits ---
% -Fmax <= b_ij*(theta_i - theta_j) <= Fmax
Aineq = zeros(2*nlines, nvars);
bineq = zeros(2*nlines, 1);
for k = 1:nlines
    ii   = busMap(linedata(k,1));
    jj   = busMap(linedata(k,2));
    b    = 1 / linedata(k,4);
    Fmax = linedata(k,6) / basemva;
    Aineq(2*k-1, theta_offset+ii) =  b;
    Aineq(2*k-1, theta_offset+jj) = -b;
    Aineq(2*k,   theta_offset+ii) = -b;
    Aineq(2*k,   theta_offset+jj) =  b;
    bineq(2*k-1:2*k) = Fmax;
end




lb = -inf(nvars,1);
ub =  inf(nvars,1);
for g = 1:ngen
    row_idx = find(busdata(:,1) == active_gen_ids(g));
    lb(g)   = 0;
    ub(g)   = busdata(row_idx, 9); 
end
lb(ngen+1:ngen+nbus)    = 0;    
ub(ngen+1:ngen+nbus)    = 1;
lb(theta_offset+1:end)  = -pi;
ub(theta_offset+1:end)  =  pi;



% Economic objective (eq. 12): linear in x, all buses included
f_econ_vec = zeros(nvars,1);
f_econ_vec(1:ngen) = active_gen_costs;
for i = 1:nbus
    f_econ_vec(ngen+i) = -VOLL_vec(i) * Pd_all(i);
end
f_econ_const = sum(VOLL_vec .* Pd_all);   % shifts ENS term to be positive
f_econ_fh    = @(x) f_econ_vec'*x + f_econ_const;

% Criticality objective (eq. 17): nonlinear, VOLL-scaled
% Combined weight: high-VoLL + high-criticality buses penalised most
VOLL_C    = VOLL_vec .* C_vec;
f_crit_fh = @(x) sum(VOLL_C .* (1 - x(ngen+1:ngen+nbus)).^2 .* Pd_all);

lp_opts  = optimoptions('linprog',  'Display','none');
fmc_opts = optimoptions('fmincon',  'Display','none', ...
    'Algorithm','interior-point', ...
    'MaxIterations',500, ...
    'OptimalityTol',1e-6, ...
    'ConstraintTol',1e-6, ...
    'SpecifyObjectiveGradient',  true, ...
    'SpecifyConstraintGradient', true);

% Initial point: generators at half capacity, full load retention
x0 = zeros(nvars,1);
for g = 1:ngen, x0(g) = ub(g) * 0.5; end
x0(ngen+1:ngen+nbus) = 1;
x0 = max(lb, min(ub, x0));

%% --- Step 1: Economic anchor (minimise f_econ alone) ---
[x_econ, ~, flag_econ] = linprog(f_econ_vec, Aineq, bineq, Aeq, beq, lb, ub, lp_opts);
if flag_econ ~= 1
    % Infeasible problem — return empty
    Pg_opt=[];  alpha_opt=[];  f_econ_val=NaN;  f_crit_val=NaN;
    exitflag=-1;  pareto_econ=[];  pareto_crit=[];
    return;
end
f_econ_at_minecon = f_econ_fh(x_econ);
f_crit_at_minecon = f_crit_fh(x_econ);
eps_hi            = f_crit_at_minecon;

%% --- Step 2: Criticality anchor (minimise f_crit alone) ---
obj_crit = @(x) fcrit_with_grad(x, ngen, nbus, VOLL_C, Pd_all);
[x_crit, ~, flag_crit] = fmincon(obj_crit, x0, Aineq, bineq, Aeq, beq, lb, ub, [], fmc_opts);
if flag_crit > 0
    f_crit_at_mincrit = f_crit_fh(x_crit);
else
    f_crit_at_mincrit = f_crit_at_minecon;
    x_crit = x_econ;
end
eps_lo = f_crit_at_mincrit;

if abs(eps_hi - eps_lo) < 1e-6
    Pg_opt    = x_econ(1:ngen);
    alpha_opt = x_econ(ngen+1:ngen+nbus);
    f_econ_val = f_econ_at_minecon;
    f_crit_val = f_crit_at_minecon;
    exitflag   = 1;
    pareto_econ = f_econ_val;
    pareto_crit = f_crit_val;
    return;
end

%% --- Step 3: Epsilon-constraint Pareto sweep (eq. 23) ---
eps_values  = linspace(eps_lo, eps_hi, n_pareto);
pareto_econ = nan(n_pareto,1);
pareto_crit = nan(n_pareto,1);
pareto_x    = cell(n_pareto,1);
obj_econ    = @(x) fecon_with_grad(x, f_econ_vec, f_econ_const);

for ep = 1:n_pareto
    epsilon = eps_values(ep);
    nl_con  = @(x) fcrit_constraint(x, ngen, nbus, VOLL_C, Pd_all, epsilon);

    frac   = (ep-1) / max(n_pareto-1, 1);
    x_init = max(lb, min(ub, (1-frac)*x_econ + frac*x_crit));

    [x_ep, ~, flag_ep] = fmincon(obj_econ, x_init, Aineq, bineq, Aeq, beq, lb, ub, nl_con, fmc_opts);
    if flag_ep > 0
        pareto_econ(ep) = f_econ_fh(x_ep);
        pareto_crit(ep) = f_crit_fh(x_ep);
        pareto_x{ep}    = x_ep;
    end
end

%% --- Step 4: Select knee-point (eq. 24) ---
valid = ~isnan(pareto_econ);
if ~any(valid)
    Pg_opt    = x_econ(1:ngen);
    alpha_opt = x_econ(ngen+1:ngen+nbus);
    f_econ_val = f_econ_at_minecon;
    f_crit_val = f_crit_at_minecon;
    exitflag   = 1;
    return;
end

pe   = pareto_econ(valid);
pc   = pareto_crit(valid);
pe_n = (pe - min(pe)) / max(max(pe)-min(pe), 1e-9);  
pc_n = (pc - min(pc)) / max(max(pc)-min(pc), 1e-9);
[~,ki]  = min(pe_n.^2 + pc_n.^2);
vi      = find(valid);
x_knee  = pareto_x{vi(ki)};

Pg_opt     = x_knee(1:ngen);
alpha_opt  = x_knee(ngen+1:ngen+nbus);
f_econ_val = pareto_econ(vi(ki));
f_crit_val = pareto_crit(vi(ki));
exitflag   = 1;
end 



function [val, grad] = fcrit_with_grad(x, ngen, nbus, VOLL_C, Pd_all)
%   f_crit = sum_i(VOLL_C_i * (1-alpha_i)^2 * Pd_i)
%   df/d(alpha_i) = -2 * VOLL_C_i * (1-alpha_i) * Pd_i
    alpha = x(ngen+1:ngen+nbus);
    shed  = 1 - alpha;
    val   = sum(VOLL_C .* shed.^2 .* Pd_all);
    grad  = zeros(size(x));
    grad(ngen+1:ngen+nbus) = -2 * VOLL_C .* shed .* Pd_all;
end



function [val, grad] = fecon_with_grad(x, f_econ_vec, f_econ_const)
%   f_econ = f_econ_vec' * x + f_econ_const  (linear)
    val  = f_econ_vec'*x + f_econ_const;
    grad = f_econ_vec;
end



function [c, ceq, dc, dceq] = fcrit_constraint(x, ngen, nbus, VOLL_C, Pd_all, epsilon)
%   c(x) = f_crit(x) - epsilon <= 0
    alpha = x(ngen+1:ngen+nbus);
    shed  = 1 - alpha;
    c     = sum(VOLL_C .* shed.^2 .* Pd_all) - epsilon;
    ceq   = [];
    dc    = zeros(size(x));
    dc(ngen+1:ngen+nbus) = -2 * VOLL_C .* shed .* Pd_all;
    dceq  = [];
end
