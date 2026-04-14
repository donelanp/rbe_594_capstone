enable_current = true;

maxVelocity = 2; % m/s
%% populate reference positions over time
if ~exist("x_i",'var')

    % Define sample time and total steps (simple case
    
    total_time = 30; % seconds
    

    positions_x = [0 1 4];
    positions_y = [0 2 5];
    positions_z = [0 0 0];
    ref_times = [0 10 total_time];
    water_current_params = [];

% use ref states from offline_planner
else
    % x_i is a 3xN list of positions
    % Ensure x_i is defined and has the correct dimensions
    if size(x_i, 1) ~= 3
        error('x_i must be a 3xN matrix of positions.');
    end
    % normalize
    % max_value = (max(abs(x_i),[],"all"))/9;
    m_to_km = 1000;
    positions_x = x_i(1,:)' ./ m_to_km;
    positions_y = x_i(2,:)' ./ m_to_km;
    positions_z = x_i(3,:)' ./ m_to_km;
    % set starting point to 0,0,0
    % start_value = x_i(:,1) ./ m_to_km;
    % positions_x = positions_x - start_value(1);
    % positions_y = positions_y - start_value(2);
    % positions_z = positions_z - start_value(3);
    start_value = [0 0 0];
    positions_x = positions_x - start_value(1);
    positions_y = positions_y - start_value(2);
    positions_z = positions_z - start_value(3);


    % assign data to x_i
    start_pos = [positions_x(1,1);positions_y(1,1);positions_z(1,1)];
    end_pos = [positions_x(end,1);positions_y(end,1);positions_z(end,1)];
    total_distance = norm(start_pos-end_pos);
    total_time = (total_distance / (0.001*maxVelocity*0.5))/3600;  % in hours

    % create ref_times based on position data
    ref_times = linspace(0, total_time-5, size(x_i, 2));

    % populate simulation map
    water_current_params.field = current_field;
    water_current_params.x1_min = (x1_min + x1_min*0.2);
    water_current_params.x1_max = (x1_max + x1_max*0.2);
    water_current_params.x2_min = (x2_min + x2_min*0.2);
    water_current_params.x2_max = (x2_max + x2_max*0.2);

end





%% Create a time/reference vector for the entire length of the simulation
sampleTime = round(total_time * 0.001);
sampleTime = 0.1;
numSteps = round(total_time / sampleTime);
time = sampleTime * (0:numSteps-1)';
set_param('ON', 'StopTime', num2str(total_time));

% create a ref position for each time step
ref_positions_x = zeros(numSteps,1);
ref_positions_y = zeros(numSteps,1);
ref_positions_z = zeros(numSteps,1);
for i=2:length(ref_times)
    start_ind = round(ref_times(i-1) / sampleTime)+1;
    end_ind = round(ref_times(i) / sampleTime);
    ref_positions_x(start_ind:end_ind) = positions_x(i);
    ref_positions_y(start_ind:end_ind) = positions_y(i);
    ref_positions_z(start_ind:end_ind) = positions_z(i);
    if i == length(ref_times)
        ref_positions_x(end_ind+1:end) = positions_x(i);
        ref_positions_y(end_ind+1:end) = positions_y(i);
        ref_positions_z(end_ind+1:end) = positions_z(i);
    end
end

% create data, which is the 3D reference positions of the uuv
data = [ ref_positions_x, ...
         ref_positions_y, ...
         ref_positions_z ];
% plot3(data(:,1),data(:,2),data(:,3));
total_distance = norm(data(1,:)-data(end,:));

%% Create a timeseries object and assign it to a variable (e.g., 'simin')
% create timeseries object
simin = timeseries(data, time); 
% Initialize the timeseries object with appropriate properties
simin.Name = 'Reference Positions';
simin.TimeInfo.Units = 'hours';
simin.DataInfo.Units = 'kilometers';

%% create obstacle parameters
obstacle_params.num_obstacles = 1;
obstacle_params.obstacle_sizes = [];
obstacle_positions = [];
for i = 1:obstacle_params.num_obstacles
    % define obstacle size
    obstacle_size = [3; 3; 3]; % km
    obstacle_params.obstacle_sizes = [obstacle_params.obstacle_sizes obstacle_size];
    % define obstacle center position over the span of the simulation
    obstacle_i_positions = flip(data, 1);
    % add obstacle positions to list
    obstacle_positions = [obstacle_positions obstacle_i_positions];
    
end
if obstacle_params.num_obstacles > 0
    obstacle_params.obstacle_positions = timeseries(obstacle_positions, time);
end