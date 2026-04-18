function [critical_proximity, density_norm] =  compute_criticality_spatial(lats, lons, n)
rd_deg = 0.005;   
rc_deg = 0.002; 
[b_lats, b_lons] = load_geojson_points('export.geojson');
N_buildings      = length(b_lats);
building_density = zeros(n,1);

for i = 1:n
    dist_deg = sqrt((b_lats - lats(i)).^2 + (b_lons - lons(i)).^2);
    building_density(i) = sum(dist_deg <= rd_deg);
end

rho_max      = max(building_density);
density_norm = building_density / rho_max;

[c_lats, c_lons] = load_geojson_points('critical.geojson');
critical_proximity = zeros(n,1);

for i = 1:n
    dist_deg = sqrt((c_lats - lats(i)).^2 + (c_lons - lons(i)).^2);
    critical_proximity(i) = double(min(dist_deg) <= rc_deg);
end

slack_json_idx = find(ismember(1:n, 9));   
if ~isempty(slack_json_idx)
    critical_proximity(slack_json_idx) = 1;
end

end 

function [point_lats, point_lons] = load_geojson_points(filename)

    fid  = fopen(filename, 'r');
    raw  = fread(fid, inf, 'uint8=>char')';
    fclose(fid);
    geojson   = jsondecode(raw);
    features  = geojson.features;
    N         = length(features);

    point_lons = zeros(N,1);
    point_lats = zeros(N,1);
    valid      = true(N,1);

    for k = 1:N
        try
            coords        = features(k).geometry.coordinates;
            point_lons(k) = coords(1);   % GeoJSON: [lon, lat]
            point_lats(k) = coords(2);
        catch
            valid(k) = false;
        end
    end

    point_lats = point_lats(valid);
    point_lons = point_lons(valid);
end
