function [v_c] = hycom(t, x, y, origin, region)
% Compute ocean current from HYCOM data with automatic download and caching.
% Downloads the latest 7-day cycle from HYCOM.
%
% Inputs:
%   t      - time [s]
%   x      - x positions [m]
%   y      - y positions [m]
%   origin - geographic origin [lat0, lon0] [deg]
%   region - region of interest [east_min, east_max, north_min, north_max] [m]
%
% Outputs:
%   v_c - current velocity (2 x n) [m/s]
%
% Reference:
%   http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z/2024.html 

persistent Fu Fv t_max cached_file

lat0       = origin(1);
lon0       = origin(2);
cache_name = sprintf('hycom_lat%0.2f_lon%0.2f.mat', lat0, lon0);

% load or download data on first call or when origin changes
if isempty(cached_file) || ~strcmp(cached_file, cache_name)
    if isfile(cache_name)
        fprintf('Loading cached HYCOM data from %s ...\n', cache_name);
        data = load(cache_name);
    else
        fprintf('Downloading HYCOM data for origin [%0.2f, %0.2f] ...\n', lat0, lon0);
        data = download_hycom_data(origin, region);
        save(cache_name, '-struct', 'data', '-v7.3');
        fprintf('Saved %s\n', cache_name);
    end

    Fu          = griddedInterpolant({data.lon, data.lat, data.time}, data.u_east, 'linear', 'nearest');
    Fv          = griddedInterpolant({data.lon, data.lat, data.time}, data.v_north, 'linear', 'nearest');
    t_max       = data.time(end);
    cached_file = cache_name;
end

% convert from local coordinates
lla = enu2lla([x(:) y(:) zeros(numel(x), 1)], [lat0 lon0 0], 'flat');
lat = lla(:,1)';
lon = wrapTo360(lla(:,2))';

% wrap queried time in case it extends beyond maximum time
t_wrapped = mod(t, t_max);
v_c       = [Fu(lon, lat, t_wrapped); Fv(lon, lat, t_wrapped)];
end

function data = download_hycom_data(origin, region)
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

% url for dataset
url = 'http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0/uv3z/2024';

% read coordinate vectors
lat_full  = ncread(url, 'lat');
lon_full  = ncread(url, 'lon');
time_full = ncread(url, 'time');

% find index ranges for spatial subset
lat_idx = find(lat_full >= lat_min & lat_full <= lat_max);
lon_idx = find(lon_full >= lon_min & lon_full <= lon_max);
lat_sub = lat_full(lat_idx);
lon_sub = lon_full(lon_idx);
n_lat   = numel(lat_idx);
n_lon   = numel(lon_idx);

% find index ranges for temporal subset (latest 7-day cycle)
n_cycle  = 7 * 24 / 3;
n_avail  = numel(time_full);
i_start  = max(1, n_avail - n_cycle + 1);
time_idx = i_start:n_avail;
time_sub = time_full(time_idx);
n_time   = numel(time_idx);

% read surface currents (ncread appears to reverse dimension order)
start   = [lon_idx(1), lat_idx(1), 1, time_idx(1)];
count   = [n_lon, n_lat, 1, n_time];
stride  = [1, 1, 1, 1];
u_east  = squeeze(ncread(url, 'water_u', start, count, stride));
v_north = squeeze(ncread(url, 'water_v', start, count, stride));

% replace NaN (land) with zero
u_east(isnan(u_east))   = 0;
v_north(isnan(v_north)) = 0;

% pack output
data.lat     = lat_sub;
data.lon     = lon_sub;
data.time    = (time_sub - time_sub(1)) * 3600;
data.u_east  = u_east;
data.v_north = v_north;
end
