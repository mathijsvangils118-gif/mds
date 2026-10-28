%% Q3.2
% z = [x vx ix y vy iy Tn iN uC]'
%% Parameters
p.rm = 0.01; % [m]
p.Lm = 0.006; % [H]
p.Rm = 6; % [Omega]
p.Jm = 1.25e-6; % [kg m^2]
p.Km = 0.020;  % [Nm / A]
p.m = 0.15;  % [kg]
p.mg = 0.4;  % [kg]
p.mur = 0.05; % [Ns/m]
p.mulx = 0.1; % [Ns/m]
p.muly = 0.1; % [Ns/m]
p.dx = 0.25; % [m]
p.dy = 0.25; % [m]
p.Mx = p.m + 2*p.Jm/p.rm^2; % [kg]
p.My = p.m + p.mg + 2*p.Jm/p.rm^2; % [kg]
p.bx = 2*p.mulx + 8*p.mur; % [Ns / m]
p.by = 2*p.muly + 4*p.mur; % [Ns / m]

p.Rn = 4; % [Omega]
p.Ln = 2e-3; % [H]
p.Cn = 25e-6; % [F]
p.Ctn = 4; % [J/K]
p.Rtn = 10; % [K/W]
p.Te = 25 + 273.15; % [K]
p.Thigh = 250 + 273.15;  % [K]
p.Uhold = sqrt(p.Rn*(p.Thigh-p.Te)/p.Rtn);  % [V]

p.voltage = 0.8;  % [V]
p.sampleTime = 0.02;  % [s]
p.omegaTol = 1e-3;  % [rad/s]
p.sideAxis = [1 2 1 2];
p.sideSign = [1 1 -1 -1];
p.sideStart = [3 2 3 2];

% Rows: q0, qH, qV, qP
% Columns: H+, H0, V+, V0.
parkTable = [1 2 1 3; 1 2 2 4; 3 4 1 3; 3 4 2 4];

%% Initial state and modes
z = zeros(9, 1);
z([1 4]) = [p.dx/2; p.dy/2];
z(7:9) = [p.Thigh; p.Uhold/p.Rn; p.Uhold];
d.q = 1;
d.park = 1;
d.axis = 2;             
d.wall = [0 0];        

t = 0;
sample = 1;
time = 0;
state = z.';
mode = d.q;
switchTime = 0;
stopTime = [];
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',p.sampleTime);

%% flow, physical events and sampled controller
while d.q ~= 11
    tNext = sample*p.sampleTime;
    if tNext > 60
        error('The square was not completed within 60 s.');
    end

    while t < tNext-1e-12
        u = inputs(d,p);
        for k = 1:2
            if d.wall(k)*u.motor(k) < 0 && d.wall(k)*z(3*k) <= 0
                d.wall(k) = 0;
            end
        end
        u = inputs(d,p);
        options = odeset(options,'Events',@(t,z) guards(z,u,p));
        [tt,zz,te,ze,ie] = ode15s(@(t, z) sys(z, u, p),[t tNext],z,options);

        if isempty(te)
            t = tt(end);
            z = zz(end,:).';
        else
            % first event
            t = te(1);
            z = ze(1,:).';
            keep = tt < t;
            tt = [tt(keep); t];
            zz = [zz(keep,:); z.'];
        end
        time = [time; tt];
        state = [state; zz];
        mode = [mode; d.q*ones(size(tt))];

        if ~isempty(te)
            events = ie(abs(te - t) < 1e-10);
            for event = events(:).'
                k = ceil(event/ 3);
                pos = 3*k - 2;
                edge = mod(event - 1, 3) +1;
                if edge == 3
                    d.wall(k) = 0;
                else
                    limits = [p.dx p.dy];
                    z(pos) = (edge-1)*limits(k);
                    z(pos+1) = 0;
                    d.wall(k) = 2*edge-3;
                end
            end
            time(end+1,1) = t;
            state(end+1,:) = z.';
            mode(end+1,1) = d.q;
        end
    end

    t = tNext;
    omega = [(z(2) + z(5))/p.rm; (z(2) - z(5))/p.rm];
    moving = max(abs(omega)) > p.omegaTol;
    if d.q == 1
        ax = d.axis;
    else
        ax = p.sideAxis(d.q-6);
    end
    symbol = 2*ax-double(moving);
    d.park = parkTable(d.park,symbol);
    if ~moving
        stopTime(end + 1, 1) = t;
    end

    oldMode = d.q;
    if d.q == 1
        if d.park == 4
            d.q = 7;
        elseif d.park == 3
            d.axis = 1;
        end
    elseif d.park == 4
        d.q = d.q+1;
    end

    if d.q ~= oldMode
        if d.q <= 10
            d.park = p.sideStart(d.q-6);
        end
        switchTime(end+1, 1) = t;
        time(end+1, 1) = t;
        state(end+1, :) = z.';
        mode(end+1, 1) = d.q;
    end
    sample = sample + 1;
end

%% Results
modeNames = ["Homing"; "Bottom side"; "Right side"; ...
             "Top side"; "Left side"; "Finished"];
disp(table(switchTime,modeNames,'VariableNames',{'Time_s','Mode'}));

%% path, motor speeds and mode
figure('Color','w','Position',[100 100 1100 560]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile([2 1]); hold on;
grid on;
homing = mode == 1;
printing = mode >= 7 & mode <= 10;
hSquare = plot(100*state(printing,1),100*state(printing,4),'b','LineWidth',1.8);
hHome = plot(100*state(homing,1),100*state(homing,4),'--',...
    'Color',[0.5 0.5 0.5],'LineWidth',1.2);
corners = [0 0; p.dx 0; p.dx p.dy; 0 p.dy];
plot(100*corners(:, 1),100*corners(:, 2),'ko','MarkerFaceColor','k',...
    'HandleVisibility','off');
axis equal;
xlim([-1 100*p.dx+1]);
ylim([-1 100*p.dy+1]);
xlabel('x (cm)');
ylabel('y (cm)');
title('Path of the head');
legend([hHome hSquare],{'Homing','Printed square'},...
    'Location','southoutside','Orientation','horizontal');

nexttile; hold on;
grid on;
plot(time,(state(:,2)+state(:,5))/p.rm,'LineWidth',1.1);
plot(time,(state(:,2)-state(:,5))/p.rm,'LineWidth',1.1);
plot(stopTime,zeros(size(stopTime)),'kd','MarkerFaceColor','k');
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');
title('Motor speeds and stopped samples');
legend('\omega_l','\omega_r','Stopped sample',...
    'Location','south','Orientation','horizontal');
xlim([0 ceil(t)]);

nexttile; hold on;
grid on;
stairs(switchTime,0:5,'k','LineWidth',1.2);
plot(t,5,'ko','MarkerFaceColor','k');
yticks(0:5);
yticklabels({'q_1','q_7','q_8','q_9','q_{10}','q_{11}'});
xlim([0 ceil(t)]);
ylim([-0.3 5.3]);
xlabel('Time (s)');
ylabel('Mode');
title('Square-drawing mode');

%% functions used
function dz = sys(z,u,p)
    dz = zeros(9,1);
    dz(1) = z(2);
    dz(2) = (p.Km/p.rm*z(3) - p.bx*z(2))/p.Mx;
    dz(3) = (u.motor(1) - p.Rm*z(3) - 2*p.Km/p.rm*z(2))/p.Lm;
    dz(4) = z(5);
    dz(5) = (p.Km/p.rm*z(6) - p.by*z(5))/p.My;
    dz(6) = (u.motor(2) - p.Rm*z(6) - 2*p.Km/p.rm*z(5))/p.Lm;
    if u.wall(1) ~= 0, dz(1:2) = 0; end
    if u.wall(2) ~= 0, dz(4:5) = 0; end
    dz(7) = (z(9)^2/p.Rn - (z(7) - p.Te)/p.Rtn)/p.Ctn;
    dz(8) = (u.nozzle - z(9))/p.Ln;
    dz(9) = (z(8) - z(9)/p.Rn)/p.Cn;
end

function u = inputs(d,p)
    u.motor = [0;0];
    u.nozzle = p.Uhold;
    u.wall = d.wall;
    if d.q == 1
        u.motor(d.axis) = -p.voltage;
    elseif d.q >= 7 && d.q <= 10
        k = d.q-6;
        u.motor(p.sideAxis(k)) = p.sideSign(k)*p.voltage;
    end
end

function [value,isterminal,direction] = guards(z,u,p)
    value = ones(6,1);
    limits = [p.dx p.dy];
    for k = 1:2
        pos = 3*k-2;
        if u.wall(k) == 0
            if u.motor(k) < 0
                value(pos) = z(pos);
            elseif u.motor(k) > 0
                value(pos+1) = limits(k) - z(pos);
            end
        elseif u.wall(k)*u.motor(k) < 0
            value(pos + 2) = u.wall(k)*z(pos + 2);
        end
    end
    isterminal = ones(6, 1);
    direction = -ones(6, 1);
end
