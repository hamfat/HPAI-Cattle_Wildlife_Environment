clc; clear; close all;

%% =====================================================================
%                     BASE PARAMETERS
% =====================================================================
params = struct();
params.beta_c = 0.0003034;
params.beta_b = 3.943e-7;
params.beta_bc = 1.367e-5;     
params.beta_env = 1.587e-5;
params.sigma_c = 0.5075;
params.gamma_c = 0.9782;
params.gamma_b = 1.655e-03;
params.mu_c = 0.02;
params.mu_b = 0.02;
params.theta_c = 0.62853;
params.theta_b = 0.273;
params.mu_env = 0.07746;
params.LambdaC = 30;
params.LambdaB = 10;

params.muSc = 0.02;
params.muEc = 0.02;
params.muIc = 0.02;
params.muSb = 0.02;
params.muIb = 0.02;

params.mudc = 0.02;
params.mudb = 0.02;

%% =====================================================================
%                     INITIAL CONDITIONS & TIME
% =====================================================================
y0 = [1500, 0, 0, 499, 1, 0];
tspan = [0 365];
T = 1000;
t_fixed = linspace(tspan(1), tspan(2), T);

%% =====================================================================
%                     MONTE CARLO SETTINGS
% =====================================================================
nSim = 200;
nVars = 6;
param_var_frac = 0.005;

Y_all = zeros(nSim, T, nVars);
R0_values = zeros(nSim, 1);
rng(1);

param_names = fieldnames(params);
sampled_params = zeros(nSim, numel(param_names));

%% =====================================================================
%                     MONTE CARLO LOOP
% =====================================================================
for k = 1:nSim
    p_sample = params;

    for f = 1:numel(param_names)
        val = params.(param_names{f});
        lower = (1 - param_var_frac) * val;
        upper = (1 + param_var_frac) * val;
        sampled_val = lower + (upper - lower) * rand;
        p_sample.(param_names{f}) = sampled_val;
        sampled_params(k, f) = sampled_val;
    end

    % Compute R0
    Rc = (p_sample.beta_c * p_sample.sigma_c * p_sample.LambdaC) / ...
         (p_sample.mudc * (p_sample.sigma_c + p_sample.mu_c + p_sample.mudc) * ...
          (p_sample.gamma_c + p_sample.mu_c + p_sample.mudc));

    Rb = (p_sample.beta_b * p_sample.LambdaB) / ...
         (p_sample.mudb * (p_sample.gamma_b + p_sample.mu_b + p_sample.mudb));

    Renv = (p_sample.LambdaC * p_sample.beta_env * p_sample.theta_c) / ...
           (p_sample.mudc * (p_sample.sigma_c + p_sample.mu_c + p_sample.mudc) * ...
            (p_sample.gamma_c + p_sample.mu_c + p_sample.mudc) * p_sample.mu_env) + ...
           (p_sample.LambdaB * p_sample.beta_env * p_sample.theta_b) / ...
           (p_sample.mudb * (p_sample.gamma_b + p_sample.mu_b + p_sample.mudb) * p_sample.mu_env);

    R0_values(k) = Rc + Rb + Renv;

    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    for v = 1:nVars
        Y_all(k,:,v) = interp1(t, y(:,v), t_fixed, 'linear', 'extrap');
    end
end

%% =====================================================================
%                     FINAL STATES
% =====================================================================
Ic_final = squeeze(Y_all(:,end,3));
Ib_final = squeeze(Y_all(:,end,5));
B_final  = squeeze(Y_all(:,end,6));

%% =====================================================================
%                     SENSITIVITY TARGET (DEFAULT = R0)
% =====================================================================
Y = R0_values;

%% =====================================================================
%                     PRCC + p-values
% =====================================================================

exclude_R0 = {'beta_bc','muSc','muEc','muIc','muSb','muIb'};

if isequal(Y, R0_values)
    idx_rm = ismember(param_names, exclude_R0);
    X = sampled_params(:, ~idx_rm);
    names_used = param_names(~idx_rm);
else
    X = sampled_params;
    names_used = param_names;
end

Xr = tiedrank(X);
Yr = tiedrank(Y);

nParams = size(Xr,2);
PRCC = zeros(nParams,1);
pvals = zeros(nParams,1);

for j = 1:nParams
    others = setdiff(1:nParams, j);
    bj = Xr(:,others) \ Xr(:,j);
    rXj = Xr(:,j) - Xr(:,others)*bj;

    bY = Xr(:,others) \ Yr;
    rY = Yr - Xr(:,others)*bY;

    PRCC(j) = corr(rXj, rY);
    [~, pvals(j)] = corr(rXj, rY);
end



%% =====================================================================
%                     LATEX PARAMETER LABELS
% =====================================================================
latex_names = strrep(names_used,'_','\_');
latex_names = strcat('$', latex_names, '$');

%% =====================================================================
%                     SINGLE-PANEL PRCC PLOT
% =====================================================================
figure;
barh(PRCC); hold on;
set(gca,'YTick',1:nParams,'YTickLabel',latex_names,'YDir','reverse','TickLabelInterpreter','latex');
xlabel('PRCC','Interpreter','latex');
title('PRCC Sensitivity (Filled Red = p < 0.05)','Interpreter','latex');
grid on;

for j = 1:nParams
    if pvals(j) < 0.05
        plot(PRCC(j), j, 'ro', 'MarkerSize', 6, ...
            'MarkerFaceColor','r', 'LineWidth', 1.5);
    end
end


%% =====================================================================
%     MULTI-PANEL PRCC COMPARISON (NO R0)
% =====================================================================

Targets = struct();

Targets(1).name = 'Final $I_c$';
Targets(1).Y = Ic_final;

Targets(2).name = 'Peak $I_c$';
Targets(2).Y = max(squeeze(Y_all(:,:,3)), [], 2);

[~, idx] = max(squeeze(Y_all(:,:,3)), [], 2);
Targets(3).name = 'Time-to-peak $I_c$';
Targets(3).Y = t_fixed(idx)';

figure;
for m = 1:3
    Y = Targets(m).Y;

    X = sampled_params;
    names_m = param_names;

    Xr = tiedrank(X);
    Yr = tiedrank(Y);

    nParams_m = size(Xr,2);
    PRCC_m = zeros(nParams_m,1);
    pvals_m = zeros(nParams_m,1);

    for j = 1:nParams_m
        others = setdiff(1:nParams_m, j);
        bj = Xr(:,others) \ Xr(:,j);
        rXj = Xr(:,j) - Xr(:,others)*bj;

        bY = Xr(:,others) \ Yr;
        rY = Yr - Xr(:,others)*bY;

        PRCC_m(j) = corr(rXj, rY);
        [~, pvals_m(j)] = corr(rXj, rY);
    end

    latex_m = strcat('$', strrep(names_m,'_','\_'), '$');

    subplot(2,2,m);
    barh(PRCC_m); hold on;
    set(gca,'YTick',1:nParams_m,'YTickLabel',latex_m,'YDir','reverse','TickLabelInterpreter','latex');
    xlabel('PRCC','Interpreter','latex');
    title(Targets(m).name,'Interpreter','latex');
    grid on;

    for j = 1:nParams_m
        if pvals_m(j) < 0.05
            plot(PRCC_m(j), j, 'ro', 'MarkerSize', 5, ...
                'MarkerFaceColor','r', 'LineWidth', 1.5);
        end
    end
end

sgtitle('PRCC Comparison (Final, Peak, Time-to-Peak)','Interpreter','latex');

%% =====================================================================
%                     TORNADO-STYLE SORTED PRCC PLOT
% =====================================================================

% Sort by absolute PRCC (descending)
[~, sort_idx] = sort(abs(PRCC), 'descend');

PRCC_sorted = PRCC(sort_idx);
pvals_sorted = pvals(sort_idx);
names_sorted = latex_names(sort_idx);

figure;
barh(PRCC_sorted); hold on;

set(gca,'YTick',1:nParams,'YTickLabel',names_sorted, ...
    'YDir','reverse','TickLabelInterpreter','latex');

xlabel('PRCC','Interpreter','latex');
title('Tornado-Style Sorted PRCC Plot','Interpreter','latex');
grid on;

% Filled red dots for significant parameters
for j = 1:nParams
    if pvals_sorted(j) < 0.05
        plot(PRCC_sorted(j), j, 'ro', 'MarkerSize', 6, ...
            'MarkerFaceColor','r', 'LineWidth', 1.5);
    end
end

legend('PRCC','p < 0.05','Location','best');

%% =====================================================================
%                     ODE SYSTEM
% =====================================================================
function dydt = ode_system(~, y, params)
    S_c = y(1); E_c = y(2); I_c = y(3);
    S_b = y(4); I_b = y(5); B = y(6);

    dS_c = params.LambdaC - params.beta_c*S_c*I_c - params.beta_bc*S_c*I_b ...
           - params.beta_env*S_c*B - params.muSc*S_c;

    dE_c = params.beta_c*S_c*I_c + params.beta_bc*S_c*I_b + params.beta_env*S_c*B ...
           - params.sigma_c*E_c - params.mu_c*E_c - params.muEc*E_c;

    dI_c = params.sigma_c*E_c - params.gamma_c*I_c - params.mu_c*I_c - params.muIc*I_c;

    dS_b = params.LambdaB - params.beta_b*S_b*I_b - params.beta_env*S_b*B - params.muSb*S_b;

    dI_b = params.beta_b*S_b*I_b + params.beta_env*S_b*B ...
           - params.gamma_b*I_b - params.mu_b*I_b - params.muIb*I_b;

    dB = params.theta_c*I_c + params.theta_b*I_b - params.mu_env*B;

    dydt = [dS_c; dE_c; dI_c; dS_b; dI_b; dB];
end