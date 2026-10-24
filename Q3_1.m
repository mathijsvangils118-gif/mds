%% Q3.1

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
p.Thigh = 250 + 273.15; % [K]
p.Tlow = 60 + 273.15; % [K]
p.Uheat = 12; % [V]
p.Uhold = sqrt(p.Rn*(p.Thigh-p.Te)/p.Rtn); % [V]

p.xc = p.dx/2; % [m]
p.yc = p.dy/2; % [m]
p.a = 0.10; % [m]
p.b = 0.06; % [m]
p.halfTime = 3; % [s]
p.returnTime = 2; % [s]
p.wipeTime = 1; % [s]
p.loadTime = 3; % [s]
p.sampleTime = 0.02; % [s]
p.homeVoltage = 0.8; % [V]
p.omegaTol = 1e-3;  % [rad/s]
p.positionTol = 1e-5; % [m]
p.speedTol = 1e-4; % [m/s]

% The filament keeps cooling during the replacement
p.N = 800; % []
p.rho = 1040;  % [kg/m^3]
p.C = 1847; % [J/ (kg K)]
p.kappa = 0.18;  % [W /(m K)]
p.af = 0.4e-3;  % [m]
p.kappa_air = 15;  % [W /(m^2 K)]
p.kappa_bed = 100;  % [W /(m^2 K)]
p.Tb = 100 + 273.15;  % [K]
p.alpha = p.kappa/(p.rho*p.C);
p.beta = (3*p.kappa_air+p.kappa_bed)/(p.rho*p.C*p.af);
p.Tinf = (3*p.kappa_air*p.Te+p.kappa_bed*p.Tb)/(3*p.kappa_air+p.kappa_bed);
r = linspace(0,1,4001);
s = 10*r.^3 - 15*r.^4 + 6*r.^5;
phi = [pi*s, pi+pi*s(2:end)];
tp = [p.halfTime*r, p.halfTime*(1+r(2:end))];
arc = cumtrapz(phi,hypot(p.a*sin(phi),p.b*cos(phi)));
p.ds = arc(end)/p.N;
p.birthTime = interp1(arc,tp,(1:p.N)*p.ds);
p.birthTime(p.N/2) = p.halfTime;
p.birthTime(end) = 2*p.halfTime;

%% Initial conditions
z = zeros(11+p.N,1);
z(1) = p.xc+p.a;
z(4) = p.yc;
z(7:9) = [p.Thigh; p.Uhold/p.Rn; p.Uhold];
z(12:end) = p.Tinf;

d.q = 0;
d.changed = false;
d.colour = 0;
d.park = 1;
d.axis = 2;
d.contact = [false false];
d.n = 0;
d.saved = [z(1);z(4)];

% Rows: q0, qH, qV, qP
% Columns: H+, H0, V+, V0.
parkTable = [1 2 1 3; 1 2 2 4; 3 4 1 3; 3 4 2 4];

t = 0;
finished = false;
time = [];
state = [];
mode = [];
colour = [];
heater = [];
switchTime = 0;
switchMode = 0;

nStates = numel(z);
pattern = spdiags(ones(nStates,3),-1:1,nStates,nStates);
pattern(1:11,1:11) = 1;
pattern(12,end) = 1;
pattern(end,12) = 1;
options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1,...
                 'JPattern',pattern);

%% Integrate until an event, update the discrete state and restart
while ~finished
    options = odeset(options,'Events',@(t,z) guards(t,z,d,p));
    [tt,zz,te,ze,ie] = ode15s(@(t,z) sys(z,inputs(z,d,p),p),...
                           [t 300],z,options);
    if isempty(ie)
        error('A required event was not reached before 300 s.');
    end

    % process the first event
    tEvent = te(1);
    zEvent = ze(1,:);
    keep = tt < tEvent;
    tt = [tt(keep); tEvent];
    zz = [zz(keep,:); zEvent];

    time = [time;tt];
    state = [state;zz];
    mode = [mode;d.q*ones(size(tt))];
    colour = [colour;d.colour*ones(size(tt))];
    u = inputs(z,d,p);
    heater = [heater;u.nozzle*ones(size(tt))];
    t = tt(end);
    z = zz(end,:).';

    % New deposited cell enters at the current nozzle temperature
    while d.q == 0 && d.n < p.N && ...
            z(10) >= p.birthTime(d.n+1)-1e-9
        d.n = d.n+1;
        z(11+d.n) = z(7);
    end

    % inelastic contact
    if d.q == 1 && d.axis == 1 && ~d.contact(1) && z(1) <= 1e-9
        z(1:2) = 0;
        d.contact(1) = true;
    elseif d.q == 1 && d.axis == 2 && ~d.contact(2) && z(4) <= 1e-9
        z(4:5) = 0;
        d.contact(2) = true;
    end

    g = guards(t, z, d, p);
    oldMode = d.q;
    if g(1) <= 1e-8
        switch d.q
            case 0
                if ~d.changed
                    d.saved = [z(1);z(4)];
                    d.changed = true;
                    d.park = 1;
                    d.axis = 2;
                    d.q = 1;
                else
                    finished = true;
                end
            case 1
                omega = [(z(2) + z(5))/p.rm; (z(2) - z(5))/p.rm];
                moving = max(abs(omega)) > p.omegaTol;
                symbol = 2*d.axis-double(moving);
                d.park = parkTable(d.park,symbol);
                if d.park == 4
                    d.q = 2;
                elseif d.park == 2
                    d.axis = 2;
                elseif d.park == 3
                    d.axis = 1;
                end
            case 2
                d.q = 3;
            case 3
                d.q = 4;
            case 4
                d.colour = 1;
                d.q = 5;
            case 5
                d.contact = [false false];
                d.q = 6;
            case 6
                d.q = 0;
        end
        z(11) = 0;
    end

    if d.q ~= oldMode
        switchTime(end+1,1) = t;
        switchMode(end+1,1) = d.q;
    end
end

%% Results
names = {'Printing';'Moving to wipe';'Wiping';'Cooling';...
         'Loading';'Heating';'Returning'};
transitionTable = table(switchTime,names(switchMode+1),...
    'VariableNames',{'Time_s','Mode'});
disp(transitionTable);
fprintf('Printing completed at %.3f s.\n',t);

%% Printed path and deposited filament
figure('Color','w','Name','Filament replacement','Position',[100 100 1100 450]);
subplot(1,2,1); hold on;
k = mode == 0 & colour == 0;
plot(100*state(k,1),100*state(k,4),'b','LineWidth',2);
k = mode == 0 & colour == 1;
plot(100*state(k,1),100*state(k,4),'Color',[1 0.5 0],'LineWidth',2);
k = mode == 1;
plot(100*state(k,1),100*state(k,4),'k--','LineWidth',1);
k = mode == 6;
plot(100*state(k,1),100*state(k,4),'Color',[0.5 0.5 0.5],...
    'LineStyle',':','LineWidth',1.5);
plot(100*d.saved(1),100*d.saved(2),'ko','MarkerFaceColor','k',...
    'HandleVisibility','off');
plot(0,0,'ks','MarkerFaceColor','k','HandleVisibility','off');
axis equal; xlim([-1 100*p.dx+1]); ylim([-1 100*p.dy+1]); grid on;
xlabel('x (cm)'); ylabel('y (cm)'); title('Blue to orange');
legend('Blue filament','Orange filament','Moving to wipe','Returning',...
       'Location','northeast');

subplot(1,2,2); hold on;
centres = ((1:p.N)-0.5)*p.ds;
plot(100*centres(1:p.N/2),z(12:11+p.N/2)-273.15,'b','LineWidth',1.5);
plot(100*centres(p.N/2+1:end),z(12+p.N/2:end)-273.15,...
    'Color',[1 0.5 0],'LineWidth',1.5);
yline(p.Tinf-273.15,'k--','Equilibrium');
grid on;
xlabel('Distance along the ellipse (cm)');
ylabel('Filament temperature (deg C)');
title('At printing completion');
legend('Blue filament','Orange filament','Location','northwest');

%% Mode, nozzle temperature and heater input
figure('Color','w','Name','Hybrid simulation','Position',[100 100 1050 680]);
subplot(3,1,1);
stairs([switchTime;t],[switchMode;switchMode(end)],'k','LineWidth',1.3);
yticks(0:6);
yticklabels({'q_0','q_1','q_2','q_3','q_4','q_5','q_6'});
ylim([-0.4 6.4]);
xlim([0 t]);
grid on;
ylabel('Mode');

subplot(3,1,2);
plot(time,state(:,7)-273.15,'b','LineWidth',1.4); hold on;
yline(p.Tlow-273.15,'k--','Low');
yline(p.Thigh-273.15,'k:','Printing');
xlim([0 t]); grid on; ylabel('Nozzle temperature (deg C)');

subplot(3,1,3);
plot(time,heater,'Color',[0.8 0.3 0],'LineWidth',1.4);
xlim([0 t]);
ylim([-1 p.Uheat+1]);
grid on;
xlabel('Time (s)');
ylabel('Nozzle voltage (V)');

%% Functions used
function dz = sys(z,u,p)
    dz = zeros(size(z));
    dz(1) = z(2);
    dz(2) = (p.Km/p.rm*z(3)- p.bx*z(2))/p.Mx;
    dz(3) = (u.motor(1) - p.Rm*z(3) - 2*p.Km/p.rm*z(2))/p.Lm;
    dz(4) = z(5);
    dz(5) = (p.Km/p.rm*z(6) - p.by*z(5))/p.My;
    dz(6) = (u.motor(2) - p.Rm*z(6)- 2*p.Km/p.rm*z(5))/p.Lm;
    if u.contact(1), dz(1:2) = 0; end
    if u.contact(2), dz(4:5) = 0; end
    dz(7) = (z(9)^2/p.Rn - (z(7) - p.Te)/p.Rtn)/p.Ctn;
    dz(8) = (u.nozzle - z(9))/p.Ln;
    dz(9) = (z(8) - z(9)/p.Rn)/p.Cn;
    dz(10) = u.printing;
    dz(11) = 1;

    n = u.n;
    if n > 0
        T = z(12:11+n);
        dT = -p.beta*(T-p.Tinf);
        flux = p.alpha/p.ds^2*(T(2:end) - T(1:end-1));
        dT(1:end-1) = dT(1:end-1)+flux;
        dT(2:end) = dT(2:end)-flux;
        if n == p.N
            flux = p.alpha/p.ds^2*(T(1) - T(end));
            dT(end) = dT(end)+flux;
            dT(1) = dT(1) - flux;
        end
        dz(12:11+n) = dT;
    end
end

function u = inputs(z,d,p)
    u.motor = [0;0];
    u.nozzle = p.Uhold;
    u.printing = double(d.q == 0);
    u.contact = d.contact;
    u.n = d.n;
    if d.q == 0
        [~,v,a,j] = ellipse_reference(z(10),p);
        u.motor = motor_voltage(v,a,j,p);
    elseif d.q == 1
        u.motor(d.axis) = -p.homeVoltage;
    elseif d.q == 6
        [~,v,a,j] = smooth_step(z(11),p.returnTime);
        u.motor = motor_voltage(d.saved*v,d.saved*a,d.saved*j,p);
    end
    if d.q == 3 || d.q == 4
        u.nozzle = 0;
    elseif d.q == 5
        u.nozzle = p.Uheat;
    end
end

function [value,isterminal,direction] = guards(~,z,d,p)
    value = ones(4,1);
    switch d.q
        case 0
            value(1) = p.halfTime*(1 + double(d.changed))- z(10);
        case 1
            value(1) = p.sampleTime-z(11);
            if d.axis == 1 && ~d.contact(1), value(2) = z(1); end
            if d.axis == 2 && ~d.contact(2), value(3) = z(4); end
        case 2
            value(1) = p.wipeTime - z(11);
        case 3
            value(1) = z(7) - p.Tlow;
        case 4
            value(1) = p.loadTime - z(11);
        case 5
            value(1) = p.Thigh - z(7);
        case 6
            error = norm([z(1);z(4)]- d.saved)/p.positionTol-1;
            speed = norm([z(2);z(5)])/p.speedTol-1;
            value(1) = max([1-z(11)/p.returnTime,error,speed]);
    end
    if d.q == 0 && d.n < p.N
        value(4) = p.birthTime(d.n+1)-z(10);
    end
    isterminal = ones(4,1);
    direction = -ones(4,1);
end

function u = motor_voltage(v,a,j,p)
    M = [p.Mx;p.My];
    b = [p.bx;p.by];
    i = p.rm/p.Km*(M.*a+b.*v);
    di = p.rm/p.Km*(M.*j+b.*a);
    u = p.Lm*di+p.Rm*i+2*p.Km/p.rm*v;
end

function [q,v,a,j] = ellipse_reference(t,p)
    if t < p.halfTime
        [f,f1,f2,f3] = smooth_step(t,p.halfTime);
        phi = pi*f;
    else
        [f,f1,f2,f3] = smooth_step(t-p.halfTime,p.halfTime);
        phi = pi+pi*f;
    end
    w = pi*f1;
    dw = pi*f2;
    ddw = pi*f3;
    q = [p.xc + p.a*cos(phi);p.yc + p.b*sin(phi)];
    q1 = [-p.a*sin(phi);p.b*cos(phi)];
    q2 = [-p.a*cos(phi);-p.b*sin(phi)];
    q3 = [p.a*sin(phi);-p.b*cos(phi)];
    v = q1*w;
    a = q2*w^2 + q1*dw;
    j = q3*w^3 + 3*q2*w*dw + q1*ddw;
end
 
function [s,v,a,j] = smooth_step(t,T)
    r = max(0,min(1,t/T));
    s = 10*r^3 - 15*r^4 + 6*r^5;
    v = (30*r^2 - 60*r^3 + 30*r^4)/T;
    a = (60*r - 180*r^2 + 120*r^3)/T^2;
    j = (60 - 360*r + 360*r^2)/T^3;
    if t < 0 || t > T
        v = 0; a = 0; j = 0;
    end
end
