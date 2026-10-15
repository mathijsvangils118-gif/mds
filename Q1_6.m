%% Q1.6 

%% Parameters
% Gotten from Hofstaetter(2015)
p.rho = 1040; % [kg/m^3]
p.C = 1847; % [J/(kg K)]
p.kappa = 0.18; % [W/(m K)]

%filament size and heat-transfer coefficients
p.af = 0.4e-3; % [m]         
p.kappa_air = 15; % [W/(m^2 K)]          
p.kappa_bed = 100; % [W/(m^2 K)]         
p.tolerance = 1; % []       
p.Tdep = 250 + 273.15; % [K]
p.Tb = 100 + 273.15; % [K]
p.Te = 25 + 273.15; % [K]   
p.thresholdTemp = 100 + 273.15; % [K]

%% Ellipse from Q1.4
p.a = 0.10; p.b = 0.06;
p.dx = 0.25; p.dy = 0.25;
p.xc = p.dx/2; p.yc = p.dy/2;
p.Tf = 5;              
p.tr = 1;              
p.area = p.af^2;
p.alpha = p.kappa/(p.rho*p.C);
p.beta = (3*p.kappa_air+p.kappa_bed)/(p.rho*p.C*p.af);
p.Tinf = (3*p.kappa_air*p.Te+p.kappa_bed*p.Tb)/(3*p.kappa_air+p.kappa_bed);

% distance along the ellipse
p.phiTable = linspace(0,2*pi,40001).';
p.sTable = cumtrapz(p.phiTable,...
    hypot(p.a*sin(p.phiTable),p.b*cos(p.phiTable)));
p.length = p.sTable(end);

%% Simulation
settings = [400 0.005; 800 0.005; 800 0.0025];
checkResults = zeros(size(settings,1),4);

for run = 1:size(settings, 1)
    N = settings(run, 1);
    dt = settings(run, 2);
    keepHistory = run == size(settings, 1);
    [metric,data] = simulate_cooling(p, N, dt, keepHistory);
    checkResults(run,:) = [N, dt, metric.uniformTime, metric.waitTime];
end

vmax = max(p.a,p.b)*2*pi/(p.Tf-p.tr);
depositionTimeScale = (p.length/N)/vmax;

fprintf('Ellipse length: %.6f m\n',p.length);
fprintf('Equilibrium temperature: %.3f degrees Celsius\n',p.Tinf-273.15);
fprintf('ds = %.4f mm, dt = %.4f ms\n',1000*p.length/N,1000*dt);
fprintf('ds/vmax: %.4f ms\n',1000*depositionTimeScale);
fprintf('\n     N       dt (s)    Uniform time (s)\n');
fprintf('%6d     %8.4f       %10.4f\n',checkResults(:,1:3).');
fprintf('\nNearly uniform at %.4f s; waiting time: %.4f s\n',...
    metric.uniformTime,metric.waitTime);
fprintf('Temperature range: %.3f to %.3f degrees Celcius\n',...
    metric.minFinal-273.15,metric.maxFinal-273.15);

if isnan(metric.thresholdTime)
    fprintf('The %.1f degrees Celcius threshold was not reached.\n',p.thresholdTemp-273.15);
else
    fprintf('\nMaximum temperature reaches %.1f degrees Celcius at %.4f s\n',...
        p.thresholdTemp-273.15,metric.thresholdTime);
    fprintf('Waiting time: %.4f s; temperature spread: %.3f K\n',...
        metric.thresholdWaitTime,metric.spreadAtThreshold);
end

%% Cooling
new_figure('Q1.6 - Cooling');
t = data.stats(:,1);
subplot(1,2,1);
plot(t,data.stats(:,4)-273.15,'r','LineWidth',1.4); hold on;
plot(t,data.stats(:,3)-273.15,'k','LineWidth',1.4);
plot(t,data.stats(:,2)-273.15,'b','LineWidth',1.4);
yline(p.Tinf-273.15,'k:','HandleVisibility','off');
xline(p.Tf,'k:','HandleVisibility','off');
xline(metric.uniformTime,'r--','HandleVisibility','off');
grid on;
xlabel('Time (s)');
ylabel('Temperature (^oC)');
title('Temperature of the deposited filament');
legend('Maximum','Length-weighted mean','Minimum','Location','northeast');
xlim([0,ceil(metric.uniformTime)]);

subplot(1,2,2);
spread = data.stats(:,5);
spread(spread <= 0) = NaN;
semilogy(t,spread,'b','LineWidth',1.5); hold on;
yline(p.tolerance,'k--','1 K');
xline(p.Tf,'k:','Printing complete');
xline(metric.uniformTime,'r--','Nearly uniform');
grid on; xlabel('Time (s)'); ylabel('Maximum - minimum temperature (K)');
title('Uniformity criterion');
xlim([0,ceil(metric.uniformTime)]); ylim([0.5 200]);

%% Temperature map
new_figure('Q1.6 - Temperature map');
surf(data.time,100*data.s,(data.temperature-273.15).',...
    'EdgeColor','none'); view(2); hold on;
plot3(data.time,100*data.front,p.Tdep*ones(size(data.time)),...
    'k--','LineWidth',1.2);
axis tight; set(gca,'YDir','normal');
colormap(turbo);
caxis([p.Tinf-273.15,p.Tdep-273.15]);
cb = colorbar;
cb.Label.String = 'Temperature (^oC)';
xlabel('Time (s)');
ylabel('Distance along the ellipse (cm)');
title('Filament temperature along the ellipse');
xline(p.Tf,'k:','Printing complete');


%% Functions used
function [metric,data] = simulate_cooling(p,N,dt,keepHistory)
    ds = p.length/N;
    left = (0:N-1).'*ds;
    s = left+ds/2;
    capPerLength = p.rho*p.C*p.area;
    lossPerLength = p.af*(3*p.kappa_air+p.kappa_bed);
    tEnd = max(60,p.Tf+12/p.beta);
    tGrid = (0:dt:tEnd).';
    phi = ellipse_phase(tGrid,p.Tf,p.tr);
    front = interp1(p.phiTable,p.sTable,phi,'linear');
    front(tGrid >= p.Tf) = p.length;

    Told = zeros(N,1);
    capOld = zeros(N,1);
    stats = zeros(numel(tGrid)-1,5);
    recordEvery = max(1,round(0.02/dt));
    if keepHistory
        maxRecords = ceil(numel(tGrid)/recordEvery)+5;
        tempHistory = NaN(maxRecords,N);
        timeHistory = zeros(maxRecords,1);
        frontHistory = zeros(maxRecords,1);
        record = 1; 
    end
    energyResidual = 0;
    uniformTime = NaN;

    for step = 1:numel(tGrid)-1
        t = tGrid(step+1);
        cellLength = min(max(front(step+1)-left,0),ds);
        if t >= p.Tf, cellLength(:) = ds; end
        n = find(cellLength > 0,1,'last');
        li = cellLength(1:n);
        cap = capPerLength*li;
        H = lossPerLength*li;

        % Conduction
        i = (1:n-1).'; j = i+1;
        G = p.kappa*p.area./(0.5*(li(1:end-1)+li(2:end)));
        if tGrid(step) >= p.Tf-1e-12
            % closing the loop after printing
            i = [i;n]; j = [j;1];
            G = [G;p.kappa*p.area/ds];
        end
        K = sparse([i;j;i;j],[i;j;j;i],[G;G;-G;-G],n,n);

        % Backward euler
        addedCap = cap-capOld(1:n);
        lhs = spdiags(cap+dt*H,0,n,n)+dt*K;
        rhs = capOld(1:n).*Told(1:n)+addedCap*p.Tdep+dt*H*p.Tinf;
        T = lhs\rhs;

        energyChange = sum(cap.*T)- sum(capOld(1:n).*Told(1:n));
        energyInput = sum(addedCap)*p.Tdep-dt*sum(H.*(T - p.Tinf));
        energyResidual = max(energyResidual,abs(energyChange-energyInput));
        Tmin = min(T);
        Tmax = max(T);
        Tmean = sum(li.*T)/sum(li);
        spread = Tmax-Tmin;
        stats(step,:) = [t,Tmin,Tmean,Tmax,spread];
        if abs(t-p.Tf) < 0.1*dt, Tclosure = T; end

        done = t >= p.Tf && spread <= p.tolerance;
        if done
            uniformTime = t;
            if step > 1 && stats(step - 1, 1) >= p.Tf && ...
                    stats(step - 1,5) > p.tolerance
                prev = stats(step - 1, :);
                uniformTime = prev(1)+(p.tolerance - prev(5))*...
                    (t - prev(1))/(spread - prev(5));
            end
        end
        if keepHistory && (step == 1 || mod(step,recordEvery) == 0 || ...
                abs(t-p.Tf) < 0.1*dt || done)
            record = record + 1;
            tempHistory(record, 1:n) = T.';
            timeHistory(record) = t;
            frontHistory(record) = front(step+1);
        end
        Told(:) = 0;
        Told(1:n) = T;
        capOld(:) = 0;
        capOld(1:n) = cap;
        if done, break; end
    end
    if isnan(uniformTime)
        error('Uniformity tolerance was not reached.');
    end
    stats = stats(1:step, :);

    %temperature threshold after printing
    thresholdTime = NaN;
    spreadAtThreshold = NaN;
    k = find(stats(:,1) >= p.Tf & stats(:,4) <= p.thresholdTemp,1,'first');
    if ~isempty(k)
        thresholdTime = stats(k, 1);
        spreadAtThreshold = stats(k, 5);
        if k > 1 && stats(k - 1, 1) >= p.Tf && ...
                stats(k - 1, 4) > p.thresholdTemp
            fraction = (p.thresholdTemp-stats(k - 1, 4))/...
                (stats(k, 4)-stats(k - 1, 4));
            thresholdTime = stats(k - 1, 1)+fraction*...
                (stats(k, 1)-stats(k - 1, 1));
            spreadAtThreshold = stats(k - 1, 5)+fraction*...
                (stats(k, 5)-stats(k - 1, 5));
        end
    end
    after = stats(:,1) >= p.Tf;
    meanExact = p.Tinf+(mean(Tclosure)-p.Tinf)*...
        exp(-p.beta*(stats(after,1)-p.Tf));

    metric.uniformTime = uniformTime;
    metric.waitTime = uniformTime-p.Tf;
    metric.thresholdTime = thresholdTime;
    metric.thresholdWaitTime = thresholdTime-p.Tf;
    metric.spreadAtThreshold = spreadAtThreshold;
    metric.minAtClosure = min(Tclosure);
    metric.maxAtClosure = max(Tclosure);
    metric.minFinal = min(T);
    metric.maxFinal = max(T);
    metric.meanFinal = mean(T);
    metric.energyResidual = energyResidual;
    metric.meanDecayError = max(abs(stats(after,3)-meanExact));
    data.stats = stats;
    if keepHistory
        data.s = s;
        data.time = timeHistory(1:record);
        data.front = frontHistory(1:record);
        data.temperature = tempHistory(1:record,:);
    end
end

function phi = ellipse_phase(t,Tf,tr)
    W = 2*pi/(Tf-tr);
    SIntegral = @(q) 2.5*q.^4-3*q.^5+q.^6;
    phi = zeros(size(t));
    first = t < tr;
    middle = t >= tr & t <= Tf-tr;
    last = t > Tf-tr & t < Tf;
    phi(first) = W*tr*SIntegral(t(first)/tr);
    phi(middle) = W*(t(middle)-tr/2);
    phi(last) = 2*pi-W*tr*SIntegral((Tf-t(last))/tr);
    phi(t >= Tf) = 2*pi;
end

function new_figure(name)
    figure('Name',name,'NumberTitle','off','Color','w',...
        'WindowStyle','normal','Units','normalized',...
        'OuterPosition',[0.08 0.10 0.84 0.78]);
end
