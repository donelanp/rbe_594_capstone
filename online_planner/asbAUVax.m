function asbAUVax(block)
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
    block.NumOutputPorts = 0;

    % Override input port properties
    block.InputPort(1).DatatypeID = 0;
    block.InputPort(1).Dimensions = 6;
    block.InputPort(1).DirectFeedthrough = false;

    block.InputPort(2).DatatypeID = 0;
    block.InputPort(2).Dimensions = 3;
    block.InputPort(2).DirectFeedthrough = false;
    
    %
    % initialize the array of sample times
    %
    block.SampleTimes = [Config.update 0]; % fixed sample time

    %
    % specify that the simState for this s-function is same as the default
    %
    block.SimStateCompliance = 'DefaultSimState';

    if ~Config.Animenable
        return
    end

    %
    % Initialize Figure Window
    %

    h_f=findobj('type','figure','Tag','asbAUVFigure');

    if isempty(h_f)
        h_anim=figure;
        delete(uigettool(h_anim,'Exploration.Brushing'));
        delete(uigettool(h_anim,'DataManager.Linking'));
        delete(findall(h_anim,'Tag','figDataManagerBrushTools'));
        delete(findall(h_anim,'Tag','figBrush'));
        delete(findall(h_anim,'Tag','figLinked'));
    else
        h_anim=h_f;
    end

    if ~strcmpi(get(h_anim,'windowStyle'),'docked')
        set(h_anim,'resize','off','position',[253 252 546 449]);
    end

    set(h_anim,'name','Animation Figure', ...
               'clipping','off', ...
               'Tag','asbAUVFigure');

    if ~isempty(h_anim)
        h_del = findobj(h_anim,'type','axes');
        delete(h_del);
        figure(h_anim);
    end

    %
    % Initialize Axes
    %
    handle.axes(1)=axes;

    % default simulation bounds
    x_lims = [-2 12];
    y_lims = [-2 12];
    z_lims = [-2 12];

    if ~isempty(Config.water_current_params)
        x_lims_water = [ Config.water_current_params.x1_min ...
                   Config.water_current_params.x1_max ];
        y_lims_water = [ Config.water_current_params.x2_min ...
                   Config.water_current_params.x2_max ];
        % plot current speed
        current_field = Config.water_current_params.field;
        x1 = linspace(x_lims_water(1), x_lims_water(2), 41);
        x2 = linspace(y_lims_water(1), y_lims_water(2), 41);
        [X1, X2]         = meshgrid(x1, x2);
        C     = current_field(zeros(1, numel(X1)), X1(:)', X2(:)');
        speed = reshape(sqrt(C(1,:).^2 + C(2,:).^2), [41,41]);
        m_to_km = 1000;
      
        % x_lims = [0 (abs(X1(1,1)) + Config.water_current_params.x1_max)/m_to_km];
        % y_lims = [0 (abs(X2(1,1)) + Config.water_current_params.x2_max)/m_to_km];
        x_lims = x_lims_water / m_to_km;
        y_lims = y_lims_water / m_to_km;

        % X1 = X1 - X1(1,1);
        % X2 = X2 - X2(1,1);
        hold on;
        pcolor(handle.axes(1),X1/m_to_km, X2/m_to_km, speed, 'EdgeColor',  'none', 'FaceAlpha', 0.7, 'FaceColor', 'interp', 'HandleVisibility', 'off');
        colorbar;
        % plot current direction
        quiver(handle.axes(1),X1(:)'/m_to_km, X2(:)'/m_to_km, C(1,:), C(2,:), 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end

    
    set(handle.axes(1),'visible','on', ...
       'xLim',x_lims,'yLim',y_lims,'zLim',z_lims,...
       'xGrid','on','yGrid','on','zGrid','on',...
       'units','normal', ...
       'clipping','off');

    if Config.dimension==1
        % 3D
        view(handle.axes(1),[-38 35]);
    else
        % 2D
        view(handle.axes(1),2);
    end
    xlabel(handle.axes(1),'X (km)');
    ylabel(handle.axes(1),'Y (km)');
    zlabel(handle.axes(1),'Z (km)');


    %
    % Initialize Time & Utilities for back-stepping support
    %
    handle.Id = 'asbAUV';
    handle.Model = Config.Model;
    handle.Index = 0;
    handle.snapshot = [];

    %
    % Draw in Body Shape
    %
    [xm,ym,zm] = bodyShape;
    handle.craft = patch(xm,ym,zm,[1 1 1],"facealpha",0);
    vert = get(handle.craft,'vertices');
    set(handle.axes(1),'userData',vert)
    set(handle.craft,'edgeColor',[1 0 0],'clipping','off');
    %
    % Set Handles of graphics in Figure UserData
    %
    set(h_anim,'userData',handle);

    %
    % MCOS support for back-stepping
    %
    data.handle = h_anim;
    data.sixDofData = 0;
    b=snapshotanim(data);

    try
        Simulink.SimulationStepper(bdroot).addSnapshotInterface(b);
    catch ME
        if ~bdIsLibrary(bdroot)
            rethrow(ME);
        end
    end

end

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function mdlUpdate(block)
    u = block.InputPort(1).Data;
    reference_position = block.InputPort(2).Data;
    Config = block.DialogPrm(1).Data;
    
    if ~Config.Animenable
        return;
    end

    %
    % Obtain Handles of Figure Objects
    %
    handle = get(findobj('type','figure','Tag','asbAUVFigure'),'userData');

    if isempty(findobj('type','figure','Tag','asbAUVFigure'))
        %figure has been manually closed
        return
    end

    hold on;
    if Config.dimension==1
        % 3D
        plot3(handle.axes(1),reference_position(1),reference_position(2),reference_position(3),'y*');
    else
        % 2D
        plot(handle.axes(1),reference_position(1),reference_position(2),'y*'); 
    end
    legend('UUV','Local Goal Points','Location','northwest');

    %
    % Form Transformation Matrix
    %
    cph = cos(u(4));            % Roll
    sph = sin(u(4));
    cth = cos(u(5));            % Pitch
    sth = sin(u(5));
    cps = cos(u(6));            % Yaw
    sps = sin(u(6));
    attitude = [cth*cps sph*sth*cps-cph*sps cph*sth*cps+sph*sps
                cth*sps sph*sth*sps+cph*cps cph*sth*sps-sph*cps
                -sth         sph*cth         cph*cth];

    %
    % Update Craft Object
    %
    vert = get(handle.axes(1),'userData');
    a=size(vert,1);
    dum =attitude*vert'+repmat(u(1:3,1),1,a);
    set(handle.craft,'vertices',dum');

    %
    % Force MATLAB to Update Drawing
    %
    drawnow

end

function mdlGetTimeOfNextVarHit(block)
% mdlGetTimeOfNextVarHit  Return the time of the next hit for this block.
    t = block.CurrentTime;
    Config = block.DialogPrm(1).Data;
    block.NextTimeHit = ceil(t/Config.update)*Config.update+Config.update;
end

function [x,y,z]=bodyShape    
% Function to draw shape of body
    % xyz = 2*[0 2 2   0   0 0   2   2   0   0   0   0   2   2   2   2
    %          0 0 0.4 0.4 0 0   0   0.4 0.4 0   0.4 0.4 0.4 0.4 0   0
    %          0 0 0   0   0 0.4 0.4 0.4 0.4 0.4 0.4 0   0   0.4 0.4 0];
    xyz = 20*[0 2 2   0   0 0   2   2   0   0   0   0   2   2   2   2
                0 0 0.4 0.4 0 0   0   0.4 0.4 0   0.4 0.4 0.4 0.4 0   0
                0 0 0   0   0 0.4 0.4 0.4 0.4 0.4 0.4 0   0   0.4 0.4 0];
    x = xyz(1,:);
    y = xyz(2,:);
    z = xyz(3,:);
end
