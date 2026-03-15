% Define sample time and total steps
sampleTime = 0.01;
totalTime = 30; % seconds
numSteps = totalTime / sampleTime;
% Create a time vector
time = sampleTime * (0:numSteps-1)';

% populate reference positions over time
positions_x = [0 1 4];
positions_y = [0 2 5];
positions_z = [0 0 0];
ref_times = [0 10 20];

ref_positions_x = zeros(size(time));
ref_positions_y = zeros(size(time));
ref_positions_z = zeros(size(time));
for i=2:length(ref_times)
    start_ind = (ref_times(i-1) / sampleTime)+1;
    end_ind = (ref_times(i) / sampleTime);
    ref_positions_x(start_ind:end_ind) = positions_x(i);
    ref_positions_y(start_ind:end_ind) = positions_y(i);
    ref_positions_z(start_ind:end_ind) = positions_z(i);
end

data = [ ref_positions_x, ...
         ref_positions_y, ...
         ref_positions_z ];


% Create a timeseries object and assign it to a variable (e.g., 'simin')
simin = timeseries(data, time); 