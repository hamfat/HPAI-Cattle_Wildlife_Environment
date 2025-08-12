clear; clc;

%% === Baseline Parameters ===
params.beta_c   = 0.00008;
params.beta_b   = 0.00005;
%params.beta_bc  = 0.00008;
params.beta_env = 0.00005;

params.sigma_c = 0.25;
params.gamma_c = 0.2;
params.gamma_b = 0.2;
params.mu_c    = 0.02;
params.mu_b    = 0.02;

params.theta_c = 0.1;
params.theta_b = 0.1;
params.mu_env  = 0.5;

params.LambdaC = 30;
params.LambdaB = 10;

params.mudc = 0.02;
params.mudb = 0.02;

%% === Parameters to Vary ===
param_names = {'beta_c', 'beta_b', 'beta_env', 'sigma_c', ...
               'gamma_c', 'gamma_b', 'mu_c', 'mu_b','mudb', 'mudc', ...
               'mu_env', 'theta_c', 'theta_b', ...
               'LambdaC', 'LambdaB'};

base_vals = cellfun(@(f) params.(f), param_names);

%% === Sampling ===
N = 500; % Number of samples
rng(1);  % Reproducibility
sampled_params = randn(N, length(base_vals)) .* (0.1 * base_vals) + base_vals;

%% === Compute R0 for Each Sample ===
R0_values = zeros(N, 1);

for i = 1:N
    for j = 1:length(param_names)
        params.(param_names{j}) = sampled_params(i, j);
    end

    Rc = (params.beta_c * params.sigma_c * params.LambdaC) / ...
         (params.mudc * (params.sigma_c + params.mu_c + params.mudc) * ...
          (params.gamma_c + params.mu_c + params.mudc));

    Rb = (params.beta_b * params.LambdaB) / ...
         (params.mudb * (params.gamma_b + params.mu_b + params.mudb));

    Renv = (params.LambdaC * params.beta_env * params.theta_c) / ...
           (params.mudc * (params.sigma_c + params.mu_c + params.mudc) * ...
            (params.gamma_c + params.mu_c + params.mudc) * params.mu_env) + ...
           (params.LambdaB * params.beta_env * params.theta_b) / ...
           (params.mudb * (params.gamma_b + params.mu_b + params.mudb) * params.mu_env);

    R0_values(i) = Rc + Rb + Renv;
end

%% === Compute PRCC ===
ranked_params = tiedrank(sampled_params);
ranked_R0     = tiedrank(R0_values);

PRCC = zeros(1, size(ranked_params, 2));
for j = 1:size(ranked_params, 2)
    others = setdiff(1:size(ranked_params, 2), j);
    PRCC(j) = partialcorr(ranked_params(:, j), ranked_R0, ranked_params(:, others));
end

%% === LaTeX-formatted Parameter Names ===
latex_param_names = {'$\beta_c$', '$\beta_b$', '$\beta_{env}$', '$\sigma_c$', ...
                     '$\gamma_c$', '$\gamma_b$', '$\mu_c$', '$\mu_b$', '$\mu_{db}$', '$\mu_{dc}$', ...
                     '$\mu_{env}$', '$\theta_c$', '$\theta_b$', ...
                     '$\Lambda_C$', '$\Lambda_B$'};

%% === Plot PRCC as Rectangular Horizontal Boxes ===
figure;
barh(PRCC, 'FaceColor', [0.2 0.6 0.5], 'BarWidth', 0.6);
yticks(1:length(latex_param_names));
yticklabels(latex_param_names);
set(gca, 'TickLabelInterpreter', 'latex');
xlabel('Partial Rank Correlation Coefficient (PRCC)', 'Interpreter', 'latex');
title('Global Sensitivity Analysis using PRCC', 'Interpreter', 'latex');
grid on;
box on;
