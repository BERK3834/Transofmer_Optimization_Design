%% ================================================================
% EE3006 Task 4 Optimization
% Maxwell-consistent safe optimization
% Baseline: N1=230, N2=115, Wc=45 mm
%% ================================================================

clc; clear; close all;

%% Given data
S  = 1000;     % VA
V1 = 220;      % V rms
V2 = 110;      % V rms
f  = 50;       % Hz
a  = V1/V2;

I1 = S/V1;
I2 = S/V2;

D    = 80e-3;
Wwin = 40e-3;
Hwin = 70e-3;

rho_cu  = 0.0175;   % ohm mm^2/m
dens_cu = 8960;     % kg/m^3
dens_fe = 7650;     % kg/m^3
kf      = 0.95;

price_cu = 9.0;
price_fe = 2.0;

B_ref   = 1.5;
Pfe_ref = 1.5;

% Maxwell calibration:
% analytical baseline B = 1.260 T, Maxwell baseline B ≈ 1.43 T
Bcorr = 1.43 / 1.260;

B_limit = 1.50;

%% Fixed conductor areas
Aw1 = 2.08;   % AWG14 primary [mm^2]
Aw2 = 5.26;   % AWG10 secondary [mm^2]

%% Baseline design
base.N1  = 230;
base.N2  = 115;
base.Wc  = 45e-3;
base.Aw1 = Aw1;
base.Aw2 = Aw2;

base = evaluateDesign(base,S,V1,V2,f,a,I1,D,Wwin,Hwin, ...
    rho_cu,dens_cu,dens_fe,kf,price_cu,price_fe,B_ref,Pfe_ref,Bcorr);

%% Optimization range
N1_values = 220:1:230;
Wc_values = (46:0.25:47)*1e-3;

best.F = inf;
results = [];

wCost = 0.45;
wLoss = 0.55;

%% Optimization loop
for N1 = N1_values
    for Wc = Wc_values

        d.N1  = N1;
        d.N2  = round(N1/2);
        d.Wc  = Wc;
        d.Aw1 = Aw1;
        d.Aw2 = Aw2;

        d = evaluateDesign(d,S,V1,V2,f,a,I1,D,Wwin,Hwin, ...
            rho_cu,dens_cu,dens_fe,kf,price_cu,price_fe,B_ref,Pfe_ref,Bcorr);

        if d.Bmax > B_limit
            continue;
        end

        if d.Cost > base.Cost
            continue;
        end

        if abs(d.N1/d.N2 - 2) > 0.02
            continue;
        end

        F = wCost*(d.Cost/base.Cost) + wLoss*(d.Ploss/base.Ploss);
        d.F = F;

        results = [results; ...
            d.N1 d.N2 d.Wc*1e3 d.Bmax d.Pcu d.Pfe d.Ploss ...
            d.eff*100 d.m_cu d.m_fe d.Cost d.F];

        if F < best.F
            best = d;
            best.F = F;
        end
    end
end

if isempty(results)
    error('No feasible design found. Relax constraints or check input values.');
end

%% Print summary
fprintf('\n========================================================\n');
fprintf(' TASK 4 OPTIMIZATION SUMMARY\n');
fprintf('========================================================\n');
fprintf('%-14s %12s %12s\n','Quantity','Baseline','Optimized');
fprintf('%-14s %12.3f %12.3f\n','N1',base.N1,best.N1);
fprintf('%-14s %12.3f %12.3f\n','N2',base.N2,best.N2);
fprintf('%-14s %12.3f %12.3f\n','Wc [mm]',base.Wc*1e3,best.Wc*1e3);
fprintf('%-14s %12.3f %12.3f\n','Aw1 [mm2]',base.Aw1,best.Aw1);
fprintf('%-14s %12.3f %12.3f\n','Aw2 [mm2]',base.Aw2,best.Aw2);
fprintf('%-14s %12.3f %12.3f\n','Bmax [T]',base.Bmax,best.Bmax);
fprintf('%-14s %12.3f %12.3f\n','Pcu [W]',base.Pcu,best.Pcu);
fprintf('%-14s %12.3f %12.3f\n','Pfe [W]',base.Pfe,best.Pfe);
fprintf('%-14s %12.3f %12.3f\n','Ploss [W]',base.Ploss,best.Ploss);
fprintf('%-14s %12.3f %12.3f\n','Eff [%]',base.eff*100,best.eff*100);
fprintf('%-14s %12.3f %12.3f\n','m_cu [kg]',base.m_cu,best.m_cu);
fprintf('%-14s %12.3f %12.3f\n','m_fe [kg]',base.m_fe,best.m_fe);
fprintf('%-14s %12.3f %12.3f\n','Cost [USD]',base.Cost,best.Cost);
fprintf('========================================================\n');
fprintf('Cost reduction vs baseline : %6.2f %%\n',100*(1-best.Cost/base.Cost));
fprintf('Loss reduction             : %6.2f %%\n',100*(1-best.Ploss/base.Ploss));
fprintf('Efficiency gain            : %+6.3f points\n',100*(best.eff-base.eff));
fprintf('========================================================\n');

%% Output folder
desktopPath = fullfile(getenv('USERPROFILE'), 'Desktop');
outFolder = fullfile(desktopPath, 'Transformer_Report_Figures');

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%% Plot 1
figure('Color','w');
scatter(results(:,11),results(:,7),45,results(:,4),'filled'); hold on;
plot(base.Cost,base.Ploss,'ks','MarkerSize',10,'LineWidth',2);
plot(best.Cost,best.Ploss,'rp','MarkerSize',14,'LineWidth',2);
grid on;
xlabel('Material cost [USD]');
ylabel('Total loss [W]');
title('Optimization: Cost vs Total Loss');
legend('Feasible designs','Baseline','Selected optimized design','Location','best');
cb = colorbar;
ylabel(cb,'Bmax [T]');
exportgraphics(gcf, fullfile(outFolder, 'opt_cost_loss.png'), 'Resolution', 300);

%% Plot 2
figure('Color','w');
scatter(results(:,11),results(:,8),45,results(:,4),'filled'); hold on;
plot(base.Cost,base.eff*100,'ks','MarkerSize',10,'LineWidth',2);
plot(best.Cost,best.eff*100,'rp','MarkerSize',14,'LineWidth',2);
grid on;
xlabel('Material cost [USD]');
ylabel('Efficiency [%]');
title('Optimization: Efficiency vs Cost');
legend('Feasible designs','Baseline','Selected optimized design','Location','best');
cb = colorbar;
ylabel(cb,'Bmax [T]');
exportgraphics(gcf, fullfile(outFolder, 'opt_efficiency_cost.png'), 'Resolution', 300);

%% Plot 3
figure('Color','w');
scatter(results(:,1),results(:,3),45,results(:,4),'filled'); hold on;
plot(best.N1,best.Wc*1e3,'rp','MarkerSize',14,'LineWidth',2);
grid on;
xlabel('Primary turns N1');
ylabel('Center limb width Wc [mm]');
title('Bmax Map of Feasible Designs');
cb = colorbar;
ylabel(cb,'Bmax [T]');
exportgraphics(gcf, fullfile(outFolder, 'opt_bmax_map.png'), 'Resolution', 300);

fprintf('\nFigures saved to:\n%s\n', outFolder);

%% ================= LOCAL FUNCTION =================
function d = evaluateDesign(d,S,V1,V2,f,a,I1,D,Wwin,Hwin, ...
    rho_cu,dens_cu,dens_fe,kf,price_cu,price_fe,B_ref,Pfe_ref,Bcorr)

    Ac = d.Wc*D*kf;
    d.Bmax = Bcorr * V1/(4.44*f*d.N1*Ac);

    MLTbase = 2*(d.Wc + D);
    MLT1 = MLTbase + 0.144;
    MLT2 = MLTbase + 0.056;

    l1 = d.N1*MLT1;
    l2 = d.N2*MLT2;

    R1 = rho_cu*l1/d.Aw1;
    R2 = rho_cu*l2/d.Aw2;

    d.Rcu_eq = R1 + a^2*R2;
    d.Pcu = I1^2*d.Rcu_eq;

    Vcu = (l1*d.Aw1 + l2*d.Aw2)*1e-6;
    d.m_cu = Vcu*dens_cu;

    Wtot = d.Wc + d.Wc + 2*Wwin;
    Htot = d.Wc + Hwin;

    Airon = Wtot*Htot - 2*(Wwin*Hwin);
    Vfe = Airon*D*kf;

    d.m_fe = Vfe*dens_fe;

    Pspec = Pfe_ref*(d.Bmax/B_ref)^2;
    d.Pfe = d.m_fe*Pspec;

    d.Ploss = d.Pcu + d.Pfe;
    d.eff = S/(S+d.Ploss);

    d.Cost = price_cu*d.m_cu + price_fe*d.m_fe;
end
