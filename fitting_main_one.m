%% ============================================================
%   CATTLE–WILDBIRD–ENVIRONMENT MODEL: FITTING + MCMC
%   Full integrated version with posterior predictive plots
% =============================================================
rng(1)
clear; clc;

%% ------------------------------------------------------------
%  Load data
% -------------------------------------------------------------
data = readtable('cattle_wildbirds_count1.csv');

t_data = data.time;
Ic_obs = data.cattle_cases;
Ib_obs = data.wildbird_cases;

%% ------------------------------------------------------------
%  Initial conditions
% -------------------------------------------------------------
y0 = [
    1000;          % Sc
    10;            % Ec
    Ic_obs(1);     % Ic
    500;           % Sb
    Ib_obs(1);     % Ib
    10             % B
];

%% ------------------------------------------------------------
%  Parameter vector (10 parameters)
% -------------------------------------------------------------
p0 = [
    8e-5;    % beta_c
    5e-5;    % beta_b
    8e-5;    % beta_bc
    5e-5;    % beta_env
    0.25;    % sigma_c
    0.2;     % gamma_c
    0.2;     % gamma_b
    0.1;     % theta_c
    0.1;     % theta_b
    0.5      % mu_env
];

param_names = { ...
    'beta_c','beta_b','beta_bc','beta_env', ...
    'sigma_c','gamma_c','gamma_b', ...
    'theta_c','theta_b','mu_env'};

% Escape underscores for plotting
param_names_tex = strrep(param_names, '_', '\_');

%% ------------------------------------------------------------
%  Fixed parameters
% -------------------------------------------------------------
par.LambdaC = 30;
par.LambdaB = 10;

par.muSc = 0.02;
par.muEc = 0.02;
par.muIc = 0.02;
par.muSb = 0.02;
par.muIb = 0.02;

%% ------------------------------------------------------------
%  Fit parameters using LSQNONLIN
% -------------------------------------------------------------
opts = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',5000);

lb = zeros(size(p0));
ub = [];

[p_hat,~,residual] = lsqnonlin( ...
    @(p) residual_fun(p,t_data,y0,Ic_obs,Ib_obs,par), ...
    p0, lb, ub, opts);

chi2_min = sum(residual.^2);

%% ------------------------------------------------------------
%  Print MLE estimates
% -------------------------------------------------------------
disp('=== MLE Parameter Estimates ===')
for i = 1:length(p_hat)
    fprintf('%12s : %.6g\n', param_names{i}, p_hat(i));
end

%% ------------------------------------------------------------
%  Plot fitted model
% -------------------------------------------------------------
[t,y_hat] = ode45(@(t,y) cattle_bird_env_ode(t,y,p_hat,par), t_data, y0);

figure;
plot(t_data, Ic_obs, 'o', t_data, y_hat(:,3), '-')
xlabel('Time','Interpreter','none')
ylabel('Infectious cattle','Interpreter','none')
legend({'Observed','Fitted'},'Interpreter','none')

figure;
plot(t_data, Ib_obs, 'o', t_data, y_hat(:,5), '-')
xlabel('Time','Interpreter','none')
xlim([0 21])
ylabel('Infectious wild birds','Interpreter','none')
legend({'Observed','Fitted'},'Interpreter','none')

%% ------------------------------------------------------------
%  Profile likelihood for beta_env
% -------------------------------------------------------------
beta_env_vals = linspace(0.3*p_hat(4), 3*p_hat(4), 30);
chi2_prof = zeros(size(beta_env_vals));

for i = 1:length(beta_env_vals)

    beta_env_fixed = beta_env_vals(i);

    p_free0 = p_hat([1:3,5:end]);
    lb_free = lb([1:3,5:end]);

    p_fit_free = lsqnonlin( ...
        @(p_free) residual_fun(insert_beta_env(p_free,beta_env_fixed), ...
        t_data,y0,Ic_obs,Ib_obs,par), ...
        p_free0, lb_free, [], opts);

    p_full = insert_beta_env(p_fit_free, beta_env_fixed);

    r = residual_fun(p_full,t_data,y0,Ic_obs,Ib_obs,par);
    chi2_prof(i) = sum(r.^2);
end

figure;
plot(beta_env_vals, chi2_prof,'LineWidth',2)
hold on
yline(chi2_min + 1,'--')
yline(chi2_min + 3.84,':')
xlabel('\beta_{env}','Interpreter','tex')
ylabel('\chi^2','Interpreter','none')
title('Profile Likelihood for \beta_{env}','Interpreter','tex')

%% ------------------------------------------------------------
%  MCMC
% -------------------------------------------------------------
n_iter = 20000;
chain = zeros(n_iter,length(p_hat));

chain(1,:) = p_hat';
chi2_curr = chi2_min;

proposal_sd = 0.05*p_hat;

for i = 2:n_iter

    p_prop = chain(i-1,:) + proposal_sd'.*randn(1,length(p_hat));

    if any(p_prop <= 0)
        chain(i,:) = chain(i-1,:);
        continue
    end

    r_prop = residual_fun(p_prop',t_data,y0,Ic_obs,Ib_obs,par);
    chi2_prop = sum(r_prop.^2);

    alpha = exp(-0.5*(chi2_prop - chi2_curr));

    if rand < alpha
        chain(i,:) = p_prop;
        chi2_curr = chi2_prop;
    else
        chain(i,:) = chain(i-1,:);
    end
end

%% ------------------------------------------------------------
%  Posterior summaries
% -------------------------------------------------------------
posterior_mean  = mean(chain);
posterior_median = median(chain);
cred_int = prctile(chain,[2.5 97.5]);

disp('=== Posterior Summary ===')
for i = 1:length(p_hat)
    fprintf('%12s : mean=%.4g, median=%.4g, 95%% CI=[%.4g, %.4g]\n', ...
        param_names{i}, ...
        posterior_mean(i), posterior_median(i), ...
        cred_int(1,i), cred_int(2,i));
end

%% ------------------------------------------------------------
%  Posterior predictive: Infectious Cattle (Ic)
% -------------------------------------------------------------
idx = randsample(size(chain,1), 500);
Ymat_Ic = zeros(length(t_data), length(idx));

for k = 1:length(idx)
    p_k = chain(idx(k),:)';
    [~, yk] = ode45(@(t,y) cattle_bird_env_ode(t,y,p_k,par), t_data, y0);
    Ymat_Ic(:,k) = yk(:,3);
end

Ic_mean = mean(Ymat_Ic, 2);
Ic_CI = prctile(Ymat_Ic, [2.5 97.5], 2);

figure; hold on;
fill([t_data; flipud(t_data)], ...
     [Ic_CI(:,1); flipud(Ic_CI(:,2))], ...
     [0.8 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.4);
plot(t_data, Ic_mean, 'b-', 'LineWidth', 2)
plot(t_data, Ic_obs, 'ko', 'MarkerFaceColor','k')
xlim([0 21])
xlabel('Time','Interpreter','none')
ylabel('Infectious cattle','Interpreter','none')
title('Posterior Predictive: Infectious Cattle (Ic)','Interpreter','none')
legend({'95% CI','Posterior mean','Observed data'},'Interpreter','none')

%% ------------------------------------------------------------
%  Posterior predictive: Infectious Birds (Ib)
% -------------------------------------------------------------
Ymat_Ib = zeros(length(t_data), length(idx));

for k = 1:length(idx)
    p_k = chain(idx(k),:)';
    [~, yk] = ode45(@(t,y) cattle_bird_env_ode(t,y,p_k,par), t_data, y0);
    Ymat_Ib(:,k) = yk(:,5);
end

Ib_mean = mean(Ymat_Ib, 2);
Ib_CI = prctile(Ymat_Ib, [2.5 97.5], 2);

figure; hold on;
fill([t_data; flipud(t_data)], ...
     [Ib_CI(:,1); flipud(Ib_CI(:,2))], ...
     [0.8 1 0.8], 'EdgeColor','none', 'FaceAlpha',0.4);
plot(t_data, Ib_mean, 'g-', 'LineWidth', 2)
plot(t_data, Ib_obs, 'ko', 'MarkerFaceColor','k')
xlabel('Time','Interpreter','none')
xlim([0 21])
ylabel('Infectious wild birds','Interpreter','none')
title('Posterior Predictive: Infectious Wild Birds (Ib)','Interpreter','none')
legend({'95% CI','Posterior mean','Observed data'},'Interpreter','none')

%% ------------------------------------------------------------
%  MCMC Diagnostics: Trace Plots
% ------------------------------------------------------------
figure('Name','MCMC Trace Plots');
for i = 1:length(p_hat)
    subplot(5,2,i);
    plot(chain(:,i));
    ylabel(param_names_tex{i});
    if i >= 9, xlabel('Iteration'); end
    title(['Trace for ', param_names_tex{i}]);
end

%% ------------------------------------------------------------
%  MCMC Diagnostics: Posterior Distributions
% ------------------------------------------------------------
figure('Name','Posterior Histograms');
for i = 1:length(p_hat)
    subplot(5,2,i);
    histogram(chain(:,i), 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8]);
    hold on;
    % Draw a line at the MLE estimate (p_hat)
    xline(p_hat(i), 'r--', 'LineWidth', 1.5);
    ylabel('Density');
    title(param_names_tex{i});
end
legend('Posterior','MLE Estimate');





%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function dydt = cattle_bird_env_ode(~, y, p, par)
    Sc = y(1); Ec = y(2); Ic = y(3);
    Sb = y(4); Ib = y(5); B  = y(6);

    beta_c   = p(1);
    beta_b   = p(2);
    beta_bc  = p(3);
    beta_env = p(4);
    sigma_c  = p(5);
    gamma_c  = p(6);
    gamma_b  = p(7);
    theta_c  = p(8);
    theta_b  = p(9);
    mu_env   = p(10);

    dSc = par.LambdaC - beta_c*Sc*Ic - beta_bc*Sc*Ib - beta_env*Sc*B - par.muSc*Sc;
    dEc = beta_c*Sc*Ic + beta_bc*Sc*Ib + beta_env*Sc*B - sigma_c*Ec - par.muEc*Ec;
    dIc = sigma_c*Ec - gamma_c*Ic - par.muIc*Ic;

    dSb = par.LambdaB - beta_b*Sb*Ib - beta_env*Sb*B - par.muSb*Sb;
    dIb = beta_b*Sb*Ib + beta_env*Sb*B - gamma_b*Ib - par.muIb*Ib;

    dB  = theta_c*Ic + theta_b*Ib - mu_env*B;

    dydt = [dSc; dEc; dIc; dSb; dIb; dB];
end

function res = residual_fun(p, t_data, y0, Ic_obs, Ib_obs, par)
    [~, y] = ode45(@(t,y) cattle_bird_env_ode(t,y,p,par), t_data, y0);
    res = [y(:,3) - Ic_obs; y(:,5) - Ib_obs];
end

function p_full = insert_beta_env(p_free, beta_env)
    p_full = zeros(length(p_free)+1,1);
    p_full(1:3) = p_free(1:3);
    p_full(4) = beta_env;
    p_full(5:end) = p_free(4:end);
end