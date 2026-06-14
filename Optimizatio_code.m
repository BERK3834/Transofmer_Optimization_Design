%% ================================================================
% EE3006 Task 4 BONUS Optimization
% Differential Evolution + Weighted Multi-Objective Function
% Variables:
%   x(1) = N1 primary turns
%   x(2) = Wc center limb width [mm]
%   x(3) = primary AWG index
%   x(4) = secondary AWG index
%
% Objective:
%   minimize material cost and total loss
%   F = wCost*(Cost/Cost_base) + wLoss*(Ploss/Ploss_base)
%
% Baseline:
%   N1 = 230, N2 = 115, Wc = 45 mm, AWG14 / AWG10
%% ================================================================

clc; clear; close all;
rng(7);   % fixed seed for repeatable result

%% ---------------- GIVEN DATA ----------------
p.S  = 1000;      % VA
p.V1 = 220;       % V rms
p.V2 = 110;       % V rms
p.f  = 50;        % Hz
p.a  = p.V1/p.V2;

p.I1 = p.S/p.V1;
p.I2 = p.S/p.V2;

p.D    = 80e-3;   % stack depth [m]
p.Wwin = 40e-3;   % window width [m]
p.Hwin = 70e-3;   % window height [m]

p.rho_cu  = 0.0175;   % ohm mm^2/m
p.dens_cu = 8960;     % kg/m^3
p.dens_fe = 7650;     % kg/m^3
p.kf      = 0.95;

p.price_cu = 9.0;     % USD/kg
p.price_fe = 2.0;     % USD/kg

p.B_ref   = 1.5;      % T
p.Pfe_ref = 1.5;      % W/kg at B_ref

% Maxwell calibration:
% Analytical baseline B = 1.260 T, Maxwell baseline B ≈ 1.43 T
p.Bcorr = 1.43 / 1.260;

%% ---------------- CONSTRAINTS ----------------
p.B_limit  = 1.50;    % T
p.Jmax     = 3.0;     % A/mm^2 practical limit
p.fillMax  = 0.45;    % maximum total window fill factor
p.ratioTol = 0.02;    % allowed turns-ratio error
p.costHardLimit = true;

%% ---------------- STANDARD AWG TABLE ----------------
% Area values in mm^2
p.awg.name = ["AWG16","AWG15","AWG14","AWG13","AWG12","AWG11","AWG10","AWG9","AWG8"];
p.awg.area = [1.31,   1.65,   2.08,   2.62,   3.31,   4.17,   5.26,   6.63,  8.37];

idx_AWG14 = find(p.awg.name == "AWG14");
idx_AWG10 = find(p.awg.name == "AWG10");

%% ---------------- BASELINE DESIGN ----------------
% x = [N1, Wc_mm, primary_AWG_index, secondary_AWG_index]
x_base = [230, 45, idx_AWG14, idx_AWG10];
base = evaluateX(x_base, p);

%% ---------------- OPTIMIZATION SETTINGS ----------------
% Same Maxwell-safe search range as the report
lb = [220, 46.00, 1, 1];
ub = [230, 47.00, numel(p.awg.area), numel(p.awg.area)];

wCost = 0.45;
wLoss = 0.55;

nVar   = 4;
nPop   = 70;
maxGen = 140;
Fmut   = 0.75;
CR     = 0.85;

%% ---------------- INITIAL POPULATION ----------------
pop = zeros(nPop, nVar);
fit = zeros(nPop, 1);

for i = 1:nPop
    pop(i,:) = lb + rand(1,nVar).*(ub-lb);
    pop(i,:) = repairX(pop(i,:), lb, ub);
    fit(i) = objectiveWeighted(pop(i,:), p, base, wCost, wLoss);
end

bestHistory = zeros(maxGen,1);

%% ---------------- DIFFERENTIAL EVOLUTION LOOP ----------------
for gen = 1:maxGen
    for i = 1:nPop

        candidates = setdiff(1:nPop, i);
        r = candidates(randperm(numel(candidates), 3));

        x1 = pop(r(1),:);
        x2 = pop(r(2),:);
        x3 = pop(r(3),:);

        mutant = x1 + Fmut*(x2 - x3);
        mutant = repairX(mutant, lb, ub);

        trial = pop(i,:);
        jrand = randi(nVar);

        for j = 1:nVar
            if rand <= CR || j == jrand
                trial(j) = mutant(j);
            end
        end

        trial = repairX(trial, lb, ub);
        trialFit = objectiveWeighted(trial, p, base, wCost, wLoss);

        if trialFit <= fit(i)
            pop(i,:) = trial;
            fit(i) = trialFit;
        end
    end

    [bestHistory(gen), bestIndex] = min(fit);
end

[bestF, bestIndex] = min(fit);
bestX = pop(bestIndex,:);
best = evaluateX(bestX, p);
best.F = bestF;

%% ---------------- FEASIBLE DATABASE FOR REPORT PLOTS ----------------
results = [];

for N1 = 220:1:230
    for Wc_mm = 46:0.25:47
        for idx1 = 1:numel(p.awg.area)
            for idx2 = 1:numel(p.awg.area)

                x = [N1, Wc_mm, idx1, idx2];
                d = evaluateX(x, p);

                if isFeasible(d, p, base)
                    Fval = wCost*(d.Cost/base.Cost) + wLoss*(d.Ploss/base.Ploss);

                    results = [results; ...
                        d.N1, d.N2, d.Wc*1e3, d.idx1, d.idx2, ...
                        d.Bmax, d.Pcu, d.Pfe, d.Ploss, d.eff*100, ...
                        d.m_cu, d.m_fe, d.Cost, d.J1, d.J2, d.fillFactor, Fval];
                end
            end
        end
    end
end

if isempty(results)
    error('No feasible design found. Relax constraints or check input values.');
end

[~, checkIdx] = min(results(:,17));
checkBest = results(checkIdx,:);

%% ---------------- PRINT SUMMARY ----------------
fprintf('\n========================================================\n');
fprintf(' TASK 4 BONUS OPTIMIZATION SUMMARY\n');
fprintf(' Algorithm: Differential Evolution + Weighted Objective\n');
fprintf('========================================================\n');
fprintf('%-20s %14s %14s\n','Quantity','Baseline','Optimized');
fprintf('%-20s %14.3f %14.3f\n','N1',base.N1,best.N1);
fprintf('%-20s %14.3f %14.3f\n','N2',base.N2,best.N2);
fprintf('%-20s %14.3f %14.3f\n','Wc [mm]',base.Wc*1e3,best.Wc*1e3);
fprintf('%-20s %14s %14s\n','Primary AWG',char(base.awg1),char(best.awg1));
fprintf('%-20s %14s %14s\n','Secondary AWG',char(base.awg2),char(best.awg2));
fprintf('%-20s %14.3f %14.3f\n','Aw1 [mm2]',base.Aw1,best.Aw1);
fprintf('%-20s %14.3f %14.3f\n','Aw2 [mm2]',base.Aw2,best.Aw2);
fprintf('%-20s %14.3f %14.3f\n','J1 [A/mm2]',base.J1,best.J1);
fprintf('%-20s %14.3f %14.3f\n','J2 [A/mm2]',base.J2,best.J2);
fprintf('%-20s %14.3f %14.3f\n','Fill factor',base.fillFactor,best.fillFactor);
fprintf('%-20s %14.3f %14.3f\n','Bmax [T]',base.Bmax,best.Bmax);
fprintf('%-20s %14.3f %14.3f\n','Pcu [W]',base.Pcu,best.Pcu);
fprintf('%-20s %14.3f %14.3f\n','Pfe [W]',base.Pfe,best.Pfe);
fprintf('%-20s %14.3f %14.3f\n','Ploss [W]',base.Ploss,best.Ploss);
fprintf('%-20s %14.3f %14.3f\n','Eff [%]',base.eff*100,best.eff*100);
fprintf('%-20s %14.3f %14.3f\n','m_cu [kg]',base.m_cu,best.m_cu);
fprintf('%-20s %14.3f %14.3f\n','m_fe [kg]',base.m_fe,best.m_fe);
fprintf('%-20s %14.3f %14.3f\n','Cost [USD]',base.Cost,best.Cost);
fprintf('%-20s %14.3f %14.3f\n','Objective F',1.000,best.F);
fprintf('========================================================\n');
fprintf('Cost reduction vs baseline : %8.3f %%\n',100*(1-best.Cost/base.Cost));
fprintf('Loss reduction             : %8.3f %%\n',100*(1-best.Ploss/base.Ploss));
fprintf('Efficiency gain            : %+8.4f points\n',100*(best.eff-base.eff));
fprintf('========================================================\n');

fprintf('\nBest feasible design from verification database:\n');
fprintf('N1 = %.0f, N2 = %.0f, Wc = %.2f mm, Primary = %s, Secondary = %s, F = %.5f\n', ...
    checkBest(1), checkBest(2), checkBest(3), ...
    char(p.awg.name(checkBest(4))), char(p.awg.name(checkBest(5))), checkBest(17));

%% ---------------- OUTPUT FOLDER ----------------
desktopPath = fullfile(getenv('USERPROFILE'), 'Desktop');
if isempty(desktopPath) || ~exist(desktopPath,'dir')
    desktopPath = pwd;
end

outFolder = fullfile(desktopPath, 'Transformer_Report_Figures_Bonus');

if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%% ---------------- PLOT 1: CONVERGENCE ----------------
figure('Color','w');
plot(bestHistory,'LineWidth',2);
grid on;
xlabel('Generation');
ylabel('Best objective value, F');
title('Differential Evolution Convergence');
exportgraphics(gcf, fullfile(outFolder, 'bonus_DE_convergence.png'), 'Resolution', 300);

%% ---------------- PLOT 2: COST VS LOSS ----------------
figure('Color','w');
scatter(results(:,13), results(:,9), 45, results(:,6), 'filled'); hold on;
plot(base.Cost, base.Ploss, 'ks', 'MarkerSize', 10, 'LineWidth', 2);
plot(best.Cost, best.Ploss, 'rp', 'MarkerSize', 14, 'LineWidth', 2);
grid on;
xlabel('Material cost [USD]');
ylabel('Total analytical loss [W]');
title('Bonus Optimization: Cost vs Total Loss');
legend('Feasible standard-gauge designs','Baseline','DE selected design','Location','best');
cb = colorbar;
ylabel(cb,'Bmax [T]');
exportgraphics(gcf, fullfile(outFolder, 'bonus_cost_loss.png'), 'Resolution', 300);

%% ---------------- PLOT 3: EFFICIENCY VS COST ----------------
figure('Color','w');
scatter(results(:,13), results(:,10), 45, results(:,6), 'filled'); hold on;
plot(base.Cost, base.eff*100, 'ks', 'MarkerSize', 10, 'LineWidth', 2);
plot(best.Cost, best.eff*100, 'rp', 'MarkerSize', 14, 'LineWidth', 2);
grid on;
xlabel('Material cost [USD]');
ylabel('Efficiency [%]');
title('Bonus Optimization: Efficiency vs Cost');
legend('Feasible standard-gauge designs','Baseline','DE selected design','Location','best');
cb = colorbar;
ylabel(cb,'Bmax [T]');
exportgraphics(gcf, fullfile(outFolder, 'bonus_efficiency_cost.png'), 'Resolution', 300);

%% ---------------- PLOT 4: DESIGN SPACE ----------------
figure('Color','w');
scatter(results(:,1), results(:,3), 45, results(:,6), 'filled'); hold on;
plot(best.N1, best.Wc*1e3, 'rp', 'MarkerSize', 14, 'LineWidth', 2);
grid on;
xlabel('Primary turns N1');
ylabel('Center limb width Wc [mm]');
title('Bonus Optimization: Feasible Bmax Map');
cb = colorbar;
ylabel(cb,'Bmax [T]');
exportgraphics(gcf, fullfile(outFolder, 'bonus_bmax_map.png'), 'Resolution', 300);

fprintf('\nFigures saved to:\n%s\n', outFolder);

%% ================================================================
% LOCAL FUNCTIONS
%% ================================================================

function x = repairX(x, lb, ub)
    x = max(x, lb);
    x = min(x, ub);

    % Discrete/integer variables
    x(1) = round(x(1));             % N1
    x(2) = round(x(2)/0.25)*0.25;   % Wc in 0.25 mm steps
    x(3) = round(x(3));             % AWG index primary
    x(4) = round(x(4));             % AWG index secondary

    x = max(x, lb);
    x = min(x, ub);
end

function F = objectiveWeighted(x, p, base, wCost, wLoss)
    d = evaluateX(x, p);

    F0 = wCost*(d.Cost/base.Cost) + wLoss*(d.Ploss/base.Ploss);

    penalty = 0;

    % Magnetic saturation constraint
    penalty = penalty + 1e4*max(0, d.Bmax - p.B_limit)^2;

    % Current density constraints
    penalty = penalty + 1e3*max(0, d.J1 - p.Jmax)^2;
    penalty = penalty + 1e3*max(0, d.J2 - p.Jmax)^2;

    % Window fill constraint
    penalty = penalty + 1e3*max(0, d.fillFactor - p.fillMax)^2;

    % Voltage ratio constraint
    ratioError = abs(d.N1/d.N2 - p.a);
    penalty = penalty + 1e3*max(0, ratioError - p.ratioTol)^2;

    % Cost not higher than baseline
    if p.costHardLimit
        penalty = penalty + 1e3*max(0, d.Cost/base.Cost - 1)^2;
    end

    F = F0 + penalty;
end

function ok = isFeasible(d, p, base)
    ok = true;

    if d.Bmax > p.B_limit
        ok = false;
    end

    if d.J1 > p.Jmax || d.J2 > p.Jmax
        ok = false;
    end

    if d.fillFactor > p.fillMax
        ok = false;
    end

    if abs(d.N1/d.N2 - p.a) > p.ratioTol
        ok = false;
    end

    if p.costHardLimit && d.Cost > base.Cost
        ok = false;
    end
end

function d = evaluateX(x, p)

    d.N1 = round(x(1));
    d.N2 = round(d.N1 / p.a);

    d.Wc = x(2)*1e-3;       % mm to m

    d.idx1 = round(x(3));
    d.idx2 = round(x(4));

    d.Aw1 = p.awg.area(d.idx1);
    d.Aw2 = p.awg.area(d.idx2);

    d.awg1 = p.awg.name(d.idx1);
    d.awg2 = p.awg.name(d.idx2);

    %% Magnetic flux density
    Ac = d.Wc * p.D * p.kf;
    d.Bmax = p.Bcorr * p.V1/(4.44*p.f*d.N1*Ac);

    %% Mean length per turn approximation
    MLTbase = 2*(d.Wc + p.D);
    MLT1 = MLTbase + 0.144;
    MLT2 = MLTbase + 0.056;

    l1 = d.N1*MLT1;
    l2 = d.N2*MLT2;

    %% Copper resistance and copper loss
    R1 = p.rho_cu*l1/d.Aw1;
    R2 = p.rho_cu*l2/d.Aw2;

    d.R1 = R1;
    d.R2 = R2;
    d.Rcu_eq = R1 + p.a^2*R2;
    d.Pcu = p.I1^2*d.Rcu_eq;

    %% Copper mass
    Vcu = (l1*d.Aw1 + l2*d.Aw2)*1e-6;
    d.m_cu = Vcu*p.dens_cu;

    %% Core mass
    Wtot = 2*d.Wc + 2*p.Wwin;
    Htot = d.Wc + p.Hwin;

    Airon = Wtot*Htot - 2*(p.Wwin*p.Hwin);
    Vfe = Airon*p.D*p.kf;

    d.m_fe = Vfe*p.dens_fe;

    %% Core loss approximation
    Pspec = p.Pfe_ref*(d.Bmax/p.B_ref)^2;
    d.Pfe = d.m_fe*Pspec;

    %% Total loss and efficiency
    d.Ploss = d.Pcu + d.Pfe;
    d.eff = p.S/(p.S + d.Ploss);

    %% Cost
    d.Cost = p.price_cu*d.m_cu + p.price_fe*d.m_fe;

    %% Practical checks
    d.J1 = p.I1/d.Aw1;
    d.J2 = p.I2/d.Aw2;

    % Total available winding area is approximated as two windows.
    windowArea_mm2 = 2*p.Wwin*p.Hwin*1e6;
    copperArea_mm2 = d.N1*d.Aw1 + d.N2*d.Aw2;
    d.fillFactor = copperArea_mm2/windowArea_mm2;
end
