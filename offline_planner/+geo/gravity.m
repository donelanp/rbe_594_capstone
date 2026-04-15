function [g_anom] = gravity(x, y, origin, region)
% Compute gravity anomaly from IGPP Earth free-air anomaly data with
% automatic download and caching.
%
% Inputs:
%   x      - x positions [m]
%   y      - y positions [m]
%   origin - geographic origin [lat0, lon0] [deg]
%   region - region bounds [east_min, east_max, north_min, north_max] [m]
%
% Outputs:
%   g_anom - gravity anomaly (1 x n) [mGal]

persistent Fg cached_file

lat0       = origin(1);
lon0       = origin(2);
cache_name = sprintf('gravity_lat%0.2f_lon%0.2f.mat', lat0, lon0);

% load or download data on first call or when origin changes
if isempty(cached_file) || ~strcmp(cached_file, cache_name)
    if isfile(cache_name)
        fprintf('Loading cached gravity anomaly data from %s ...\n', cache_name);
        data = load(cache_name);
    else
        fprintf('Downloading gravity anomaly data for origin [%0.2f, %0.2f] ...\n', lat0, lon0);
        data = download_gravity_data(origin, region);
        save(cache_name, '-struct', 'data', '-v7.3');
        fprintf('Saved %s\n', cache_name);
    end

    Fg          = griddedInterpolant({data.lon, data.lat}, data.z, 'linear', 'nearest');
    cached_file = cache_name;
end

% convert from local coordinates
lla    = enu2lla([x(:) y(:) zeros(numel(x), 1)], [lat0 lon0 0], 'flat');
lat    = lla(:,1)';
lon    = wrapTo360(lla(:,2))';
g_anom = Fg(lon, lat);
end

function [data] = download_gravity_data(origin, region)
lat0 = origin(1);
lon0 = origin(2);

% convert from local coordinates
lla0    = [lat0, lon0, 0];
sw      = enu2lla([region(1), region(3), 0], lla0, 'flat');
ne      = enu2lla([region(2), region(4), 0], lla0, 'flat');
lat_min = sw(1);
lat_max = ne(1);
lon_min = wrapTo360(sw(2));
lon_max = wrapTo360(ne(2));

% extract regional grid via GMT
nc_file = [tempname() '.nc'];
gmt_cmd = sprintf('gmt grdcut @earth_faa_01m -R%.6f/%.6f/%.6f/%.6f -G%s', ...
    lon_min, lon_max, lat_min, lat_max, nc_file);

[status, output] = system(gmt_cmd);
assert(status == 0, 'GMT grdcut failed:\n%s', output);

% read grid
try
    lon = double(ncread(nc_file, 'lon'));
    lat = double(ncread(nc_file, 'lat'));
catch
    lon = double(ncread(nc_file, 'x'));
    lat = double(ncread(nc_file, 'y'));
end
z = double(ncread(nc_file, 'z'));

delete(nc_file);

data.lon = lon;
data.lat = lat;
data.z   = z;
end
