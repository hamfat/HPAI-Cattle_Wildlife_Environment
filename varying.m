clc; clear; close all;

% === Base Parameters ===
params.beta_c = 0.0008;
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

% === Time Span ===
tspan = [0 365];

% === Single Initial Condition ===
y0 = [1500, 0, 0, 499, 10, 0]; % [S_c, E_c, I_c, S_b, I_b, B]

% === Vary beta_c for I_c ===
beta_c_values = [0.0004, 0.0008, 0.0015];
colors_c = lines(length(beta_c_values));

figure('Name', 'I_c vs beta_c');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$I_c$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\beta_c$ on Infected Cattle ($I_c$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(beta_c_values)
    p_sample = params;
    p_sample.beta_c = beta_c_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,3), 'LineWidth', 2, 'Color', colors_c(b,:));
end
legend(arrayfun(@(x) sprintf('$\\beta_c = %.5f$', x), beta_c_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

% === Vary beta_b for I_b ===
beta_b_values = [0.00002, 0.00005, 0.0001];
colors_b = lines(length(beta_b_values));

figure('Name', 'I_b vs beta_b');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$I_b$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\beta_b$ on Infected Birds ($I_b$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(beta_b_values)
    p_sample = params;
    p_sample.beta_b = beta_b_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,5), 'LineWidth', 2, 'Color', colors_b(b,:));
end
legend(arrayfun(@(x) sprintf('$\\beta_b = %.5f$', x), beta_b_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

% === Vary beta_env for I_c ===
beta_env_values = [0.00002, 0.00005, 0.0001];
colors_env = lines(length(beta_env_values));

figure('Name', 'I_c vs beta_env');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$I_c$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\beta_{env}$ on Infected Cattle ($I_c$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(beta_env_values)
    p_sample = params;
    p_sample.beta_env = beta_env_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,3), 'LineWidth', 2, 'Color', colors_env(b,:));
end
legend(arrayfun(@(x) sprintf('$\\beta_{env} = %.5f$', x), beta_env_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

% === Vary beta_env for I_b ===
figure('Name', 'I_b vs beta_env');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$I_b$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\beta_{env}$ on Infected Birds ($I_b$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(beta_env_values)
    p_sample = params;
    p_sample.beta_env = beta_env_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,5), 'LineWidth', 2, 'Color', colors_env(b,:));
end
legend(arrayfun(@(x) sprintf('$\\beta_{env} = %.5f$', x), beta_env_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');


% === Vary theta_c for B ===
theta_c_values = [0.05, 0.15, 0.4];
colors_b = lines(length(theta_c_values));

figure('Name', 'B vs theta_c');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$B$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\theta_c$ on Environment ($B$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(theta_c_values)
    p_sample = params;
    p_sample.theta_c = theta_c_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,6), 'LineWidth', 2, 'Color', colors_b(b,:));
end
legend(arrayfun(@(x) sprintf('$\\theta_c = %.5f$', x), theta_c_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

% === Vary gamma_c for I_c ===
gamma_c_values = [0.05, 0.15, 0.4];
colors_b = lines(length(gamma_c_values));

figure('Name', 'I_c vs gamma_c');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$I_c$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\gamma_c$ on Infected Cattle ($I_c$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(gamma_c_values)
    p_sample = params;
    p_sample.gamma_c = gamma_c_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,3), 'LineWidth', 2, 'Color', colors_b(b,:));
end
legend(arrayfun(@(x) sprintf('$\\gamma_c = %.5f$', x), theta_c_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');

% === Vary mu_env for B ===
mu_env_values = [0.4, 0.6, 0.8];
colors_b = lines(length(mu_env_values));

figure('Name', 'B vs mu_env');
hold on; grid on;
xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$B$', 'Interpreter', 'latex', 'FontSize', 16);
title('Effect of $\mu_{env}$ on Environment ($B$)', ...
    'Interpreter', 'latex', 'FontSize', 18);

for b = 1:length(mu_env_values)
    p_sample = params;
    p_sample.theta_c = mu_env_values(b);
    
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    plot(t, y(:,6), 'LineWidth', 2, 'Color', colors_b(b,:));
end
legend(arrayfun(@(x) sprintf('$\\mu_{env} = %.5f$', x), mu_env_values, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');


% === ODE System ===
function dydt = ode_system(~, y, p)
    S_c = y(1); E_c = y(2); I_c = y(3);
    S_b = y(4); I_b = y(5); B = y(6);

    dS_c = p.LambdaC - p.beta_c*S_c*I_c - p.beta_bc*S_c*I_b - p.beta_env*S_c*B - p.muSc*S_c;
    dE_c = p.beta_c*S_c*I_c + p.beta_bc*S_c*I_b + p.beta_env*S_c*B - p.sigma_c*E_c - p.mu_c*E_c - p.muEc*E_c;
    dI_c = p.sigma_c*E_c - p.gamma_c*I_c - p.mu_c*I_c - p.muIc*I_c;

    dS_b = p.LambdaB - p.beta_b*S_b*I_b - p.beta_env*S_b*B - p.muSb*S_b;
    dI_b = p.beta_b*S_b*I_b + p.beta_env*S_b*B - p.gamma_b*I_b - p.mu_b*I_b - p.muIb*I_b;

    dB = p.theta_c*I_c + p.theta_b*I_b - p.mu_env*B;

    dydt = [dS_c; dE_c; dI_c; dS_b; dI_b; dB];
end
