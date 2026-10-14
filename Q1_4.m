%% Q1.4
% State = [theta_l; theta_r; omega_l; omega_r; i_m_l; i_m_r]

%% Parameters
p.rm = 0.01; % [m]
p.Lm = 0.006; % [H]
p.Rm = 6; % [omega]
p.Jm = 1.25e-6; % [kg m^2]
p.Km = 0.020;  % [Nm / A]
p.m = 0.15;  % [kg]
p.mg = 0.4;  % [kg]
p.mur = 0.05; % [Ns/m]
p.mulx = 0.1; % [Ns/m]
p.muly = 0.1; % [Ns/m]
p.dx = 0.25; % [m]
p.dy = 0.25; % [m]
% Q1.3 parameters
p.Mx = p.m + 2*p.Jm/p.rm^2; % [kg]
p.My = p.m + p.mg + 2*p.Jm/p.rm^2; % [kg]
p.bx = 2*p.mulx + 8*p.mur; % [Ns / m]
p.by = 2*p.muly + 4*p.mur; % [Ns / m]

%% Shapes
p.xc = p.dx/2;
p.yc = p.dy/2;
p.a = 0.10;
p.b = 0.06;
p.heartY = 0.14;
p.heartScale = 0.006;
shapeNames = {'Ellipse','Heart'};
drawingTime = [5 6];
rampTime = 1;
restTime = 1;

%% State-space model from Q1.2
a = p.Jm + p.rm^2/4*(2*p.m + p.mg);
c = -p.rm^2*p.mg/4;
Mtheta = [a c; c a];
b0 = p.rm^2/2*(p.mulx + p.muly + 6*p.mur);
b1 = p.rm^2/2*(p.mulx - p.muly + 2*p.mur);
Ntheta = [b0 b1; b1 b0];

A = [zeros(2),  eye(2), zeros(2);
     zeros(2), -(Mtheta\Ntheta), p.Km*(Mtheta\eye(2));
     zeros(2), -p.Km/p.Lm*eye(2),  -p.Rm/p.Lm*eye(2)];
B = [zeros(2); zeros(2); eye(2)/p.Lm];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.01,'Jacobian',A);

%% Simulation
for shape = 1:2
    Tf = drawingTime(shape);
    tspan = linspace(0, Tf+restTime, 3001);
    q0 = reference_path(0, p, shape, Tf, rampTime);
    theta0 = [q0(1)+q0(2); q0(1)-q0(2)]/p.rm;
    z0 = [theta0; zeros(4,1)];

    rhs = @(t,z) A*z + B*motor_voltage(t ,p ,shape , Tf ,rampTime);
    [t,z] = ode15s(rhs, tspan, z0, options);

    x = p.rm/2*(z(:,1)+z(:,2));
    y = p.rm/2*(z(:,1)-z(:,2));
    vx = p.rm/2*(z(:,3)+z(:,4));
    vy = p.rm/2*(z(:,3)-z(:,4));
    speed = hypot(vx,vy);
    qd = zeros(numel(t),2);
    u = zeros(numel(t),2);
    for k = 1:numel(t)
        qd(k,:) = reference_path(t(k),p,shape,Tf,rampTime).';
        u(k,:) = motor_voltage(t(k),p,shape,Tf,rampTime).';
    end
    positionError = hypot(x-qd(:,1),y-qd(:,2));

    fprintf('%s: max error = %.3e m, max |u_l|= %.3f V, max |u_r|= %.3f V\n',...
        shapeNames{shape},max(positionError),...
        max(abs(u(:,1))),max(abs(u(:,2))));

    results(shape).name = shapeNames{shape};
    results(shape).time = t;
    results(shape).state = z;
    results(shape).position = [x y];
    results(shape).reference = qd;
    results(shape).motorVoltage = u;
    results(shape).speed = speed;
    results(shape).positionError = positionError;

    figure('Name',shapeNames{shape},'NumberTitle','off','Color','w',...
                 'WindowStyle','normal','Units','normalized',...
                 'OuterPosition',[0.08 0.10 0.84 0.78]);

    subplot(2,2,[1 3]);
    plot(100*x,100*y,'b','LineWidth',2); hold on;
    plot(100*qd(:,1),100*qd(:,2),'r--','LineWidth',1);
    plot(100*x(1),100*y(1),'ko','MarkerFaceColor','k',...
         'MarkerSize',5,'HandleVisibility','off');
    axis equal; xlim([0 100*p.dx]); ylim([0 100*p.dy]); grid on;
    xlabel('x (cm)'); ylabel('y (cm)');
    title([shapeNames{shape},': head path']);
    legend('Simulated','Reference','Location','northwest');

    subplot(2,2,2);
    plot(t,u(:,1),'LineWidth',1.3); hold on;
    plot(t,u(:,2),'LineWidth',1.3); grid on;
    xlim([0 Tf+restTime]); xlabel('Time (s)'); ylabel('Voltage (V)');
    title('Motor voltages');
    legend('u_l','u_r','Location','northeast');

    subplot(2,2,4);
    plot(t,1000*speed,'b','LineWidth',1.3); grid on;
    xlim([0 Tf+restTime]); xlabel('Time (s)'); ylabel('Speed (mm/s)');
    title('Head speed');

end

%% Functions used
function u = motor_voltage(t,p,shape,Tf,tr)
    [~,vel,acc,jerk] = reference_path(t,p,shape,Tf,tr);

    ix = p.rm/p.Km*(p.Mx*acc(1) + p.bx*vel(1));
    iy = p.rm/p.Km*(p.My*acc(2) + p.by*vel(2));
    dix = p.rm/p.Km*(p.Mx*jerk(1) + p.bx*acc(1));
    diy = p.rm/p.Km*(p.My*jerk(2) + p.by*acc(2));
    ux = p.Lm*dix + p.Rm*ix + 2*p.Km/p.rm*vel(1);
    uy = p.Lm*diy + p.Rm*iy + 2*p.Km/p.rm*vel(2);

    u = [(ux+uy)/2; (ux-uy)/2];
end

function [q,vel,acc,jerk] = reference_path(t,p,shape,Tf,tr)
    [phi,dphi,ddphi,dddphi] = phase_law(t,Tf,tr);

    % Derivatives with respect to phi
    if shape == 1
        q = [p.xc+p.a*cos(phi); p.yc+p.b*sin(phi)];
        q1 = [-p.a*sin(phi); p.b*cos(phi)];
        q2 = [-p.a*cos(phi); -p.b*sin(phi)];
        q3 = [p.a*sin(phi); -p.b*cos(phi)];
    else
        k = p.heartScale;
        q = [p.xc+4*k*(3*sin(phi)-sin(3*phi));
             p.heartY+k*(13*cos(phi)-5*cos(2*phi)-2*cos(3*phi)-cos(4*phi))];
        q1 = [12*k*(cos(phi)-cos(3*phi));
              k*(-13*sin(phi)+10*sin(2*phi)+6*sin(3*phi)+4*sin(4*phi))];
        q2 = [12*k*(-sin(phi)+3*sin(3*phi));
              k*(-13*cos(phi)+20*cos(2*phi)+18*cos(3*phi)+16*cos(4*phi))];
        q3 = [12*k*(-cos(phi)+9*cos(3*phi));
              k*(13*sin(phi)-40*sin(2*phi)-54*sin(3*phi)-64*sin(4*phi))];
    end

    % Derivatives with respect to time
    vel = q1*dphi;
    acc = q2*dphi^2 + q1*ddphi;
    jerk = q3*dphi^3 + 3*q2*dphi*ddphi + q1*dddphi;
end

function [phi,dphi,ddphi,dddphi] = phase_law(t,Tf,tr)
    W = 2*pi/(Tf-tr);

    S = @(s) 10*s^3 - 15*s^4 + 6*s^5;
    S1 = @(s) 30*s^2 - 60*s^3 + 30*s^4;
    S2 = @(s) 60*s - 180*s^2 + 120*s^3;
    Si = @(s) 2.5*s^4 - 3*s^5 + s^6;

    if t <= 0
        phi = 0; dphi = 0; ddphi = 0; dddphi = 0;
    elseif t < tr
        s = t/tr;
        phi = W*tr*Si(s);
        dphi = W*S(s);
        ddphi = W/tr*S1(s);
        dddphi = W/tr^2*S2(s);
    elseif t <= Tf-tr
        phi = W*tr/2 + W*(t-tr);
        dphi = W;
        ddphi = 0;
        dddphi = 0;
    elseif t < Tf
        s = (Tf-t)/tr;
        phi = 2*pi - W*tr*Si(s);
        dphi = W*S(s);
        ddphi = -W/tr*S1(s);
        dddphi = W/tr^2*S2(s);
    else
        phi = 2*pi;
        dphi = 0;
        ddphi = 0;
        dddphi = 0;
    end
end
