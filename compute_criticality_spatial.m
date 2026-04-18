function [critical_proximity, density_norm] = ...
    compute_criticality_spatial(lats, lons, n)
%COMPUTE_CRITICALITY_SPATIAL  Compute spatial components of criticality index.
%
%   Computes two of the three components of the composite criticality
%   index defined in eq. (16):
%
%     Y_i  = binary proximity indicator (eq. 14)
%     rho_i/rho_max = normalised building density (eq. 15)
%
%   Building and critical infrastructure coordinates are loaded from
%   GeoJSON files exported from OpenStreetMap for the York study area.
%
%   INPUTS:
%     lats   Substation latitudes  [n x 1]
%     lons   Substation longitudes [n x 1]
%     n      Number of substations
%
%   OUTPUTS:
%     critical_proximity   Binary proximity indicator Y_i [n x 1]
%     density_norm         Normalised building density rho_i/rho_max [n x 1]
%
%   NOTE: rc = 500m (~0.002 deg) and rd = 700m (~0.005 deg) at York latitude.
%   These correspond to the values reported in Table II of the paper.

% Local projection factors for York (~53.96 N)
lat_to_m = 111320;
lon_to_m = 111320 * cos(deg2rad(53.96));

%% --- Building density (eq. 15) ---
rd_metres = 700;   % building density radius

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

%% --- Critical infrastructure proximity (eq. 14) ---
rc_metres = 500;   % proximity radius (hospitals, police, fire)

[c_lats, c_lons] = load_geojson_points('critical.geojson');
critical_proximity = zeros(n,1);

for i = 1:n
    dlat_c = (c_lats - lats(i)) * lat_to_m;
    dlon_c = (c_lons - lons(i)) * lon_to_m;
    dist_c = sqrt(dlat_c.^2 + dlon_c.^2);
    critical_proximity(i) = double(min(dist_c) <= rc_metres);
end

% Manual override: Bus 1 (slack/reference) always marked critical
% This ensures the reference bus is protected under all conditions
slack_json_idx = find(ismember(1:n, 9));   % Bus 1 is OSM index 9
if ~isempty(slack_json_idx)
    critical_proximity(slack_json_idx) = 1;
end

end % compute_criticality_spatial

%% =========================================================
%  LOCAL HELPER
% =========================================================

function [point_lats, point_lons] = load_geojson_points(filename)
%LOAD_GEOJSON_POINTS  Extract point coordinates from a GeoJSON file.
%   Expects Point geometry features with [longitude, latitude] coordinate order.

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
