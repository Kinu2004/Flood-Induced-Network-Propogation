function time_profile = load_demand_profile(filename, target_date)

    T  = readtable(filename);
    T.SETTLEMENT_DATE = datetime(T.SETTLEMENT_DATE);
    day_data = T(T.SETTLEMENT_DATE == target_date, :);

    if height(day_data) ~= 48
        error('Demand data: expected 48 half-hour periods for %s, got %d.', datestr(target_date), height(day_data));
    end

    annual_peak  = max(T.ND);
    time_profile = day_data.ND / annual_peak;
end


function apx_price = load_apx_prices(filename, target_date)

    T = readtable(filename);
    T.SettlementDate = datetime(T.SettlementDate, 'InputFormat','dd-MMM-yy');

    apx     = T(strcmp(strtrim(T.MarketIndexDataProviderId),'APXMIDP'), :);
    apx     = sortrows(apx, {'SettlementDate','SettlementPeriod'});
    apx_day = sortrows(apx(apx.SettlementDate == target_date, :), 'SettlementPeriod');

    if height(apx_day) ~= 48
        error('APX price data: expected 48 periods for %s, got %d.', datestr(target_date), height(apx_day));
    end

    apx_price = apx_day.MarketIndexPrice___MWh_;
end


function [CCGT, COAL, NUCLEAR] = load_generation_mix(filename, target_date)
    T = readtable(filename);
    T.SettlementDate = datetime(T.x_SettlementDate, 'InputFormat','dd/MM/yyyy');
    T       = sortrows(T, {'SettlementDate','SettlementPeriod'});
    gen_day = sortrows(T(T.SettlementDate == target_date, :), 'SettlementPeriod');

    if height(gen_day) ~= 48
        error('Generation mix: expected 48 periods for %s, got %d.', datestr(target_date), height(gen_day));
    end

    CCGT    = gen_day.CCGT;
    COAL    = gen_day.COAL;
    NUCLEAR = gen_day.NUCLEAR;
end
