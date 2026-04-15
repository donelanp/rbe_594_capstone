function waterCurrent(block)
%ASBAUVAX Plots the live position of the AUV.

%   Copyright 2023-2024 The MathWorks, Inc.

    setup(block);
end

function setup(block)
    mdlInitializeSizes(block);
    block.RegBlockMethod('Outputs', @mdlGetTimeOfNextVarHit);
    block.RegBlockMethod('Update', @mdlUpdate);
end

%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
function mdlInitializeSizes(block)
% Register parameters
    block.NumDialogPrms = 1; % Config
    Config = block.DialogPrm(1).Data;

    % Register number of ports
    block.NumInputPorts = 2;
    block.NumOutputPorts = 1;

    % Override input port properties
    block.InputPort(2).DatatypeID = 0;
    block.InputPort(2).Dimensions = 3;
    block.InputPort(2).DirectFeedthrough = false;

    block.InputPort(1).DatatypeID = 0;
    block.InputPort(1).Dimensions = 1;
    block.InputPort(1).DirectFeedthrough = false;

    % Override output port properties
    block.OutputPort(1).DatatypeID = 0;
    block.OutputPort(1).Dimensions = 3;
    block.OutputPort(1).Complexity = 'Real';
    
    %
    % initialize the array of sample times
    %
    block.SampleTimes = [1 0]; % fixed sample time

    %
    % specify that the simState for this s-function is same as the default
    %
    block.SimStateCompliance = 'DefaultSimState';




end

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function mdlUpdate(block)
    time = block.InputPort(1).Data;
    estimated_position = block.InputPort(2).Data;
    Config = block.DialogPrm(1).Data;
    
    x = estimated_position(1)*1000;
    y = estimated_position(2)*1000;
    % disp('x:');
    % disp(x);
    % disp('y:');
    % disp(y);
    vc = [Config.current_field(time, x, y); 0];
    % disp('vc:');
    % disp(vc);
    block.OutputPort(1).Data = vc;

end

function mdlGetTimeOfNextVarHit(block)
% mdlGetTimeOfNextVarHit  Return the time of the next hit for this block.
    t = block.CurrentTime;
    Config = block.DialogPrm(1).Data;
    block.NextTimeHit = t;
end

