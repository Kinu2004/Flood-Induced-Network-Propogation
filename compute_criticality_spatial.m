function [critical_proximity, density_norm] =  compute_criticality_spatial(lats, lons, n)
lat_to_m = 111320;
lon_to_m = 111320 * cos(deg2rad(53.96));

rd_metres = 700; 
[b_lats, b_lons] = load_geojson_points('export.geojson');
N_buildings      = length(b_lats);
building_density = zeros(n,1);

for i = 1:n
    dlat_m = (b_lats - lats(i)) * lat_to_m;
    dlon_m = (b_lons - lons(i)) * lon_to_m;
    dist_m = sqrt(dlat_m.^2 + dlon_m.^2);
    building_density(i) = sum(dist_m <= rd_metres);
end

rho_max      = max(building_density);
density_norm = building_density / rho_max;

rc_metres = 500;   
[c_lats, c_lons] = load_geojson_points('critical.geojson');
critical_proximity = zeros(n,1);

for i = 1:n
    dlat_c = (c_lats - lats(i)) * lat_to_m;
    dlon_c = (c_lons - lons(i)) * lon_to_m;
    dist_c = sqrt(dlat_c.^2 + dlon_c.^2);
    critical_proximity(i) = double(min(dist_c) <= rc_metres);
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
