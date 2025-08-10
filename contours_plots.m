clc; clear; close all;

% === Define the base parameters ===
params = struct();
params.beta_c = 0.00008;
params.beta_b = 0.00005;
params.beta_bc = 0.00008;
params.beta_env = 0.00005;
params.sigma_c = 0.25;
params.gamma_c = 0.2;
params.gamma_b = 0.2;
params.mu_c = 0.02;
params.mu_b = 0.02;
params.theta_c = 0.1;
params.theta_b = 0.1;
params.mu_env = 0.5;
params.LambdaC = 30;
params.LambdaB = 10;
params.muSc = 0.02;
params.muEc = 0.02;
params.muIc = 0.02;
params.muSb = 0.02;
params.muIb = 0.02;

%% === Parameter ranges for contour plots ===
beta_c_vals = linspace(0.00001, 0.002, 200);
gamma_c_vals = linspace(0.005, 0.5, 200);
beta_b_vals = linspace(0.0001, 0.002, 200);
gamma_b_vals = linspace(0.05, 0.5, 200);
beta_env_vals = linspace(0.00001, 0.001, 200);
theta_c_vals = linspace(0.01, 0.2, 200);

R0_grid = zeros(length(beta_c_vals), length(gamma_c_vals));
R0_bb_grid = zeros(length(beta_b_vals), length(gamma_b_vals));
R0_env_grid = zeros(length(beta_env_vals), length(theta_c_vals));

%% === Compute R0 surfaces ===
for i = 1:length(beta_c_vals)
    for j = 1:length(gamma_c_vals)
        beta_c = beta_c_vals(i); gamma_c = gamma_c_vals(j);
        Rc = (beta_c * params.sigma_c * params.LambdaC) / ...
             (params.muSc * (params.sigma_c + params.mu_c + params.muEc) * ...
              (gamma_c + params.mu_c + params.muIc));
        Rb = (params.beta_b * params.LambdaB) / ...
             (params.muSb * (params.gamma_b + params.mu_b + params.muIb));
        Renv = (params.LambdaC * params.beta_env * params.theta_c) / ...
               (params.muSc * (params.sigma_c + params.mu_c + params.muEc) * ...
                (gamma_c + params.mu_c + params.muIc) * params.mu_env) + ...
               (params.LambdaB * params.beta_env * params.theta_b) / ...
               (params.muSb * (params.gamma_b + params.mu_b + params.muIb) * params.mu_env);
        R0_grid(i,j) = Rc + Rb + Renv;
    end
end

for i = 1:length(beta_b_vals)
    for j = 1:length(gamma_b_vals)
        beta_b = beta_b_vals(i); gamma_b = gamma_b_vals(j);
        Rc = (params.beta_c * params.sigma_c * params.LambdaC) / ...
             (params.muSc * (params.sigma_c + params.mu_c + params.muEc) * ...
              (params.gamma_c + params.mu_c + params.muIc));
        Rb = (beta_b * params.LambdaB) / ...
             (params.muSb * (gamma_b + params.mu_b + params.muIb));
        Renv = (params.LambdaC * params.beta_env * params.theta_c) / ...
               (params.muSc * (params.sigma_c + params.mu_c + params.muEc) * ...
                (params.gamma_c + params.mu_c + params.muIc) * params.mu_env) + ...
               (params.LambdaB * params.beta_env * params.theta_b) / ...
               (params.muSb * (gamma_b + params.mu_b + params.muIb) * params.mu_env);
        R0_bb_grid(i,j) = Rc + Rb + Renv;
    end
end

for i = 1:length(beta_env_vals)
    for j = 1:length(theta_c_vals)
        beta_env = beta_env_vals(i); theta_c = theta_c_vals(j);
        Rc = (params.beta_c * params.sigma_c * params.LambdaC) / ...
             (params.muSc * (params.sigma_c + params.mu_c + params.muEc) * ...
              (params.gamma_c + params.mu_c + params.muIc));
        Rb = (params.beta_b * params.LambdaB) / ...
             (params.muSb * (params.gamma_b + params.mu_b + params.muIb));
        Renv = (params.LambdaC * beta_env * theta_c) / ...
               (params.muSc * (params.sigma_c + params.mu_c + params.muEc) * ...
                (params.gamma_c + params.mu_c + params.muIc) * params.mu_env) + ...
               (params.LambdaB * beta_env * params.theta_b) / ...
               (params.muSb * (params.gamma_b + params.mu_b + params.muIb) * params.mu_env);
        R0_env_grid(i,j) = Rc + Rb + Renv;
    end
end


figure;

subplot(1,3,1);
contourf(gamma_c_vals, beta_c_vals, R0_grid, 40, 'LineColor', 'none');
colormap(turbo);
caxis([0 50]);
colorbar;
xlabel('\gamma_c');
ylabel('\beta_c');
title('R₀: \beta_c vs \gamma_c');
hold on;
contour(gamma_c_vals, beta_c_vals, R0_grid, [1 1], 'k', 'LineWidth', 2);

subplot(1,3,2);
contourf(gamma_b_vals, beta_b_vals, R0_bb_grid, 40, 'LineColor', 'none');
colormap(turbo);
caxis([0 10]);
colorbar;
xlabel('\gamma_b');
ylabel('\beta_b');
title('R₀: \beta_b vs \gamma_b');
hold on;
contour(gamma_b_vals, beta_b_vals, R0_bb_grid, [1 1], 'k', 'LineWidth', 2);

subplot(1,3,3);
contourf(theta_c_vals, beta_env_vals, R0_env_grid, 40, 'LineColor', 'none');
colormap(turbo);
caxis([0 10]);
colorbar;
xlabel('\theta_c');
ylabel('\beta_{env}');
title('R₀: \beta_{env} vs \theta_c');
hold on;
contour(theta_c_vals, beta_env_vals, R0_env_grid, [1 1], 'k', 'LineWidth', 2);
