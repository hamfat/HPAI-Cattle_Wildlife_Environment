clc; clear; close all;

% === Define the base parameters ===
params = struct();
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

% === Initial conditions ===
y0 = [1500, 0, 0, 499, 1, 0];  % [S_c, E_c, I_c, S_b, I_b, B]

% === Time span ===
tspan = [0 365];
T = 1000;
t_fixed = linspace(tspan(1), tspan(2), T);

% === Monte Carlo Simulation Settings ===
nSim = 200;
nVars = 6;
param_var_frac = 0.005; % 10% bounds for uniform sampling

Y_all = zeros(nSim, T, nVars);
R0_values = zeros(nSim, 1);
rng(1);

param_names = fieldnames(params);
sampled_params = zeros(nSim, numel(param_names));

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
         (p_sample.muSc * (p_sample.sigma_c + p_sample.mu_c + p_sample.muEc) * ...
          (p_sample.gamma_c + p_sample.mu_c + p_sample.muIc));
    
    Rb = (p_sample.beta_b * p_sample.LambdaB) / ...
         (p_sample.muSb * (p_sample.gamma_b + p_sample.mu_b + p_sample.muIb));
    
    Renv = (p_sample.LambdaC * p_sample.beta_env * p_sample.theta_c) / ...
           (p_sample.muSc * (p_sample.sigma_c + p_sample.mu_c + p_sample.muEc) * ...
            (p_sample.gamma_c + p_sample.mu_c + p_sample.muIc) * p_sample.mu_env) + ...
           (p_sample.LambdaB * p_sample.beta_env * p_sample.theta_b) / ...
           (p_sample.muSb * (p_sample.gamma_b + p_sample.mu_b + p_sample.muIb) * p_sample.mu_env);
    
    R0_values(k) = Rc + Rb + Renv;
    mean(R0_values(k))

    % Simulate ODE
    [t, y] = ode45(@(t, y) ode_system(t, y, p_sample), tspan, y0);
    for v = 1:nVars
        Y_all(k,:,v) = interp1(t, y(:,v), t_fixed, 'linear', 'extrap');
    end
end

% === Compute Mean and 95% CI ===
Y_mean = squeeze(mean(Y_all, 1));
Y_std  = squeeze(std(Y_all, 0, 1));
Y_CI   = 1.96 * Y_std / sqrt(nSim);

% === Final Infection Levels for Global Stability ===
Ic_final = squeeze(Y_all(:,end,3));  % I_c
Ib_final = squeeze(Y_all(:,end,5));  % I_b
B_final  = squeeze(Y_all(:,end,6));  % B

tol = 1e-3;
Ic_extinct = sum(Ic_final < tol);
Ib_extinct = sum(Ib_final < tol);
B_extinct = sum(B_final < tol);

fprintf('--- GLOBAL STABILITY NUMERICAL CHECK ---\n');
fprintf('R₀ < 1 in %d out of %d simulations\n', sum(R0_values < 1), nSim);
fprintf('I_c extinct in %d out of %d simulations\n', Ic_extinct, nSim);
fprintf('I_b extinct in %d out of %d simulations\n', Ib_extinct, nSim);
fprintf('B extinct in %d out of %d simulations\n', B_extinct, nSim);

% === Plot mean ± 95% CI ===
% var_names = {'S_c', 'E_c', 'I_c', 'S_b', 'I_b', 'B'};
% colors = lines(nVars);
% 
% for v = 1:nVars
%     figure;
%     fill([t_fixed, fliplr(t_fixed)], ...
%          [Y_mean(:,v)' - Y_CI(:,v)', fliplr(Y_mean(:,v)' + Y_CI(:,v)')], ...
%          colors(v,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none'); hold on;
% 
%     plot(t_fixed, Y_mean(:,v), 'Color', colors(v,:), 'LineWidth', 2);
%     xlabel('Time (days)');
%     ylabel(var_names{v});
%     title(['Mean ± 95% CI for ', var_names{v}]);
%     xlim([0 max(tspan)])
%     legend('95% CI', 'Mean');
%     grid on;
% end

% === Simulations for Different Initial Conditions (Separate Figures) ===
init_conditions = [
    1492, 6, 2,   490, 10,   5;
    1470, 20, 10,  499, 0,   50;
    1500, 0, 0,   499, 10,  20;
    1350, 50, 100,   499, 0,   10;
    1500, 0, 5,   499, 5,   5;
];
nInitCases = size(init_conditions, 1);
colors = lines(nInitCases);

% Infected Cattle (I_c)
figure;
for k = 1:nInitCases
    y0_case = init_conditions(k,:);
    [t, y] = ode45(@(t, y) ode_system(t, y, params), tspan, y0_case);
    plot(t, y(:,3), 'LineWidth', 2, 'Color', colors(k,:)); hold on;
end
%title('Infected Cattle (I_c)');
xlabel('Time (days)');
ylabel('$I_c$');
%legend({'Ic=2,Ib=10', 'Ic=10', 'Ib=10', 'B=10', 'Ic=5,Ib=5,B=5'});
xlim([0 365])
grid on;


% Exposed Cattle (E_c)
figure;
for k = 1:nInitCases
    y0_case = init_conditions(k,:);
    [t, y] = ode45(@(t, y) ode_system(t, y, params), tspan, y0_case);
    plot(t, y(:,2), 'LineWidth', 2, 'Color', colors(k,:)); hold on;
end
%title('Infected Cattle (E_c)');
xlabel('Time (days)');
ylabel('$E_c$');
%legend({'Ic=2,Ib=10', 'Ic=10', 'Ib=10', 'B=10', 'Ic=5,Ib=5,B=5'});
xlim([0 365])
grid on;


% Infected Birds (I_b)
figure;
for k = 1:nInitCases
    y0_case = init_conditions(k,:);
    [t, y] = ode45(@(t, y) ode_system(t, y, params), tspan, y0_case);
    plot(t, y(:,5), 'LineWidth', 2, 'Color', colors(k,:)); hold on;
end
%title('Infected Birds (I_b)');
xlabel('Time (days)');
ylabel('$I_b$');
%legend({'Ic=2,Ib=10', 'Ic=10', 'Ib=10', 'B=10', 'Ic=5,Ib=5,B=5'});
xlim([0 365])
grid on;

% Environmental Contamination (B)
figure;
for k = 1:nInitCases
    y0_case = init_conditions(k,:);
    [t, y] = ode45(@(t, y) ode_system(t, y, params), tspan, y0_case);
    plot(t, y(:,6), 'LineWidth', 2, 'Color', colors(k,:)); hold on;
end
%title('Environmental Contamination (B)');
xlabel('Time (days)');
ylabel('$B$');
%legend({'Ic=2,Ib=10', 'Ic=10', 'Ib=10', 'B=10', 'Ic=5,Ib=5,B=5'});
xlim([0 365])
grid on;

% === ODE System ===
function dydt = ode_system(~, y, params)
    S_c = y(1); E_c = y(2); I_c = y(3);
    S_b = y(4); I_b = y(5); B = y(6);

    beta_c = params.beta_c;     beta_bc = params.beta_bc; beta_env = params.beta_env;
    beta_b = params.beta_b;     sigma_c = params.sigma_c; gamma_c  = params.gamma_c;
    gamma_b = params.gamma_b;   mu_c = params.mu_c;       mu_b = params.mu_b;
    theta_c = params.theta_c;   theta_b = params.theta_b; mu_env = params.mu_env;

    LambdaC = params.LambdaC;   LambdaB = params.LambdaB;
    muSc = params.muSc; muEc = params.muEc; muIc = params.muIc;
    muSb = params.muSb; muIb = params.muIb;

    dS_c = LambdaC - beta_c*S_c*I_c - beta_bc*S_c*I_b - beta_env*S_c*B - muSc*S_c;
    dE_c = beta_c*S_c*I_c + beta_bc*S_c*I_b + beta_env*S_c*B - sigma_c*E_c - mu_c*E_c - muEc*E_c;
    dI_c = sigma_c*E_c - gamma_c*I_c - mu_c*I_c - muIc*I_c;

    dS_b = LambdaB - beta_b*S_b*I_b - beta_env*S_b*B - muSb*S_b;
    dI_b = beta_b*S_b*I_b + beta_env*S_b*B - gamma_b*I_b - mu_b*I_b - muIb*I_b;

    dB = theta_c*I_c + theta_b*I_b - mu_env*B;

    dydt = [dS_c; dE_c; dI_c; dS_b; dI_b; dB];
end
