clc; clear; close all;

% === Define the base parameters ===
params = struct();
% === Parameters ===
params.beta_c = 0.003034;
params.beta_b = 3.943e-7;
params.beta_bc = 1.367e-5;
params.beta_env = 1.587e-5;
params.sigma_c =0.5075;
params.gamma_c = 0.9782;
params.gamma_b = 1.655e-03;
params.mu_c = 0.02;
params.mu_b = 0.02;
params.theta_c =0.2853;
params.theta_b = .0173;
params.mu_env = 0.07746;
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
    plot(t, y(:,3), 'LineWidth', 3, 'Color', colors(k,:)); hold on;
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
    plot(t, y(:,2), 'LineWidth', 3, 'Color', colors(k,:)); hold on;
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
    plot(t, y(:,5), 'LineWidth', 3, 'Color', colors(k,:)); hold on;
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
    plot(t, y(:,6), 'LineWidth', 3, 'Color', colors(k,:)); hold on;
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


% clear; close all; clc;
% 
% %% ============================================================
% %  PARAMETERS
% %  ============================================================
% 
% params = struct();
% params.beta_c   = 0.003034;     % will be overwritten
% params.beta_b   = 3.943e-7;
% params.beta_bc  = 1.367e-5;
% params.beta_env = 1.587e-5;
% 
% params.sigma_c = 0.5075;
% params.gamma_c = 0.9782;
% params.gamma_b = 1.655e-03;
% 
% params.mu_c = 0.02;
% params.mu_b = 0.02;
% 
% params.theta_c = 0.2853;
% params.theta_b = 0.0173;
% 
% params.mu_env = 0.07746;
% 
% params.LambdaC = 30;
% params.LambdaB = 10;
% 
% params.muSc = 0.02;
% params.muEc = 0.02;
% params.muIc = 0.02;
% params.muSb = 0.02;
% params.muIb = 0.02;
% 
% param_names = fieldnames(params);
% 
% %% ============================================================
% %  INITIAL CONDITIONS AND TIME GRID
% %  ============================================================
% 
% y0 = [1500, 0, 0, 499, 1, 0];   % [Sc Ec Ic Sb Ib B]
% tspan = [0 365];
% 
% T = 1000;
% t_fixed = linspace(tspan(1), tspan(2), T);
% 
% %% ============================================================
% %  MONTE CARLO SETTINGS
% %  ============================================================
% 
% nSim = 200;
% nVars = 6;
% param_var_frac = 0.005;   % 0.5% variation
% rng(1);
% 
% beta_c_values = [0.00012, 0.0032];
% nBeta = numel(beta_c_values);
% 
% % Storage
% mean_tpeak = nan(nBeta,3);
% CI_low     = nan(nBeta,3);
% CI_high    = nan(nBeta,3);
% 
% %% ============================================================
% %  MAIN LOOP OVER beta_c VALUES
% %  ============================================================
% 
% for b = 1:nBeta
% 
%     params_local = params;
%     params_local.beta_c = beta_c_values(b);
% 
%     Y_all  = zeros(nSim, T, nVars);
%     R0_all = zeros(nSim,1);
% 
%     for k = 1:nSim
% 
%         % --- Sample parameters ---
%         p = params_local;
%         for f = 1:numel(param_names)
%             val = params_local.(param_names{f});
%             p.(param_names{f}) = val * (1 + param_var_frac*(2*rand-1));
%         end
% 
%         % --- Compute R0 ---
%         Rc = (p.beta_c * p.sigma_c * p.LambdaC) / ...
%              (p.muSc * (p.sigma_c + p.mu_c + p.muEc) * ...
%               (p.gamma_c + p.mu_c + p.muIc));
% 
%         Rb = (p.beta_b * p.LambdaB) / ...
%              (p.muSb * (p.gamma_b + p.mu_b + p.muIb));
% 
%         Renv = (p.LambdaC * p.beta_env * p.theta_c) / ...
%                (p.muSc * (p.sigma_c + p.mu_c + p.muEc) * ...
%                 (p.gamma_c + p.mu_c + p.muIc) * p.mu_env) + ...
%                (p.LambdaB * p.beta_env * p.theta_b) / ...
%                (p.muSb * (p.gamma_b + p.mu_b + p.muIb) * p.mu_env);
% 
%         R0_all(k) = Rc + Rb + Renv;
% 
%         % --- Solve ODE ---
%         [t,y] = ode45(@(t,y) ode_system(t,y,p), tspan, y0);
% 
%         for v = 1:nVars
%             Y_all(k,:,v) = interp1(t, y(:,v), t_fixed, 'linear', 'extrap');
%         end
%     end
% 
%     % --- Time to peak CLC ---
%     t_peak = zeros(nSim,1);
%     for k = 1:nSim
%         [~,idx] = max(squeeze(Y_all(k,:,3)));
%         t_peak(k) = t_fixed(idx);
%     end
% 
%     % --- R0 regimes ---
%     regimes = {
%         R0_all < 1
%         R0_all >= 1 & R0_all < 1.5
%         R0_all >= 1.5
%     };
% 
%     for r = 1:3
%         tp = t_peak(regimes{r});
%         if isempty(tp); continue; end
%         mean_tpeak(b,r) = mean(tp);
%         CI = prctile(tp,[2.5 97.5]);
%         CI_low(b,r)  = CI(1);
%         CI_high(b,r) = CI(2);
%     end
% end
% 
% %% ============================================================
% %  PLOTTING (CLEAR + ROBUST)
% %  ============================================================
% 
% %% ============================================================
% %  PLOTTING (LaTeX AXES, CLEAR + ROBUST)
% %  ============================================================
% 
% figure;
% 
% for b = 1:nBeta
%     subplot(1,nBeta,b); hold on;
% 
%     means = mean_tpeak(b,:);
%     low   = CI_low(b,:);
%     high  = CI_high(b,:);
% 
%     valid = ~isnan(means);
%     x = find(valid);
% 
%     errorbar(x, means(valid), ...
%              means(valid)-low(valid), ...
%              high(valid)-means(valid), ...
%              'o', 'LineWidth', 2, 'MarkerSize', 8);
% 
%     % --- X-axis (LaTeX) ---
%     set(gca,'XTick',1:3,...
%         'XTickLabel',{'$R_0<1$','$1\le R_0<1.5$','$R_0\ge1.5$'},...
%         'TickLabelInterpreter','latex');
% 
%     xlim([0.5 3.5]);
% 
%     % --- Y-axis (LaTeX) ---
%     ylabel('$\text{Time to peak CLC (days)}$', ...
%            'Interpreter','latex');
% 
%     grid on;
% 
%     % --- Adaptive y-axis with minimum span ---
%     ymin = min(low(valid));
%     ymax = max(high(valid));
%     span = ymax - ymin;
% 
%     if span < 1
%         pad = 0.5;
%     else
%         pad = 0.1 * span;
%     end
% 
%     ylim([ymin - pad, ymax + pad]);
% 
%     % --- Panel title (LaTeX) ---
%     title(sprintf('$\\beta_c = %.5f$', beta_c_values(b)), ...
%           'Interpreter','latex');
% end
% 
% % --- Super title (LaTeX) ---
% sgtitle('Time to Peak CLC Aligned with $R_0$ Regimes', ...
%         'Interpreter','latex');
% 
% 
% %% ============================================================
% %  ODE SYSTEM
% %  ============================================================
% 
% function dydt = ode_system(~, y, p)
% 
% S_c = y(1); E_c = y(2); I_c = y(3);
% S_b = y(4); I_b = y(5); B = y(6);
% 
% dS_c = p.LambdaC - p.beta_c*S_c*I_c - p.beta_bc*S_c*I_b ...
%        - p.beta_env*S_c*B - p.muSc*S_c;
% 
% dE_c = p.beta_c*S_c*I_c + p.beta_bc*S_c*I_b ...
%        + p.beta_env*S_c*B ...
%        - (p.sigma_c + p.mu_c + p.muEc)*E_c;
% 
% dI_c = p.sigma_c*E_c - (p.gamma_c + p.mu_c + p.muIc)*I_c;
% 
% dS_b = p.LambdaB - p.beta_b*S_b*I_b ...
%        - p.beta_env*S_b*B - p.muSb*S_b;
% 
% dI_b = p.beta_b*S_b*I_b + p.beta_env*S_b*B ...
%        - (p.gamma_b + p.mu_b + p.muIb)*I_b;
% 
% dB = p.theta_c*I_c + p.theta_b*I_b - p.mu_env*B;
% 
% dydt = [dS_c; dE_c; dI_c; dS_b; dI_b; dB];
% end
