clc; clear; close all;

% === Parameters ===
params.beta_c = 0.0003034;
params.beta_b = 3.943e-7;
params.beta_bc = 1.367e-5;
params.beta_env = 1.587e-5;
params.sigma_c =0.5075;
params.gamma_c = 0.9782;
params.gamma_b = 1.655e-03;
params.mu_c = 0.02;
params.mu_b = 0.02;
params.theta_c =0.2853;
params.theta_b = .0170;
params.mu_env = 0.07746;
params.LambdaC = 30;
params.LambdaB = 10;
params.muSc = 0.02;
params.muEc = 0.02;
params.muIc = 0.02;
params.muSb = 0.02;
params.muIb = 0.02;



% === Time span ===
tspan = [0 365];

% === Initial Conditions (10 different cases) ===
init_conditions = [
    1500, 0, 0, 499, 0, 20;
    1498, 0, 2, 490, 10, 0;
    1500, 0, 10, 499, 0, 0;
    1500, 0, 0, 499, 10, 40;
    1500, 0, 0, 499, 0, 20;
    1500, 0, 5, 499, 5, 5;
    1480, 5, 10, 480, 10, 10;
    1500, 2, 1, 495, 2, 1;
    1505, 0, 0, 499, 1, 0;
    1490, 3, 5, 485, 5, 3;
];

nInit = size(init_conditions, 1);
tol = 1e-3;

% Containers for extinction results
Ic_final = zeros(nInit,1);
Ib_final = zeros(nInit,1);
B_final = zeros(nInit,1);

% Colors for plots
colors = lines(nInit);

% Variables names for plotting
varNames = {'S_c', 'E_c', 'I_c', ...
    'S_b', 'I_b', 'B'};

% Prepare figure handles for all 6 variables
figs = gobjects(6,1);
for i = 1:6
    figs(i) = figure('Name', varNames{i});
    hold on; grid on;
    xlabel('Time (days)', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
    ylabel(['' varNames{i} ''], 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times New Roman');
   % title(['\textbf{Time evolution of ' varNames{i} '}'], 'Interpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman'); % Axes font
end

% Plot all variables in separate figures with line width 2
for k = 1:nInit
    y0 = init_conditions(k,:);
    [t, y] = ode45(@(t,y) ode_system(t,y,params), tspan, y0);

    for varIdx = 1:6
        figure(figs(varIdx));
        plot(t, y(:,varIdx), 'LineWidth', 2, 'Color', colors(k,:));
        xlim([0 max(tspan)])
    end
end

% % legendStrings = arrayfun(@(x) sprintf('IC %d', x), 1:nInit, 'UniformOutput', false);
% % 
% % % Add legends with LaTeX and font formatting
% % for i = 1:6
% %     figure(figs(i));
% %     legend(legendStrings, 'Interpreter', 'latex', 'FontSize', 15, 'FontName', 'Times New Roman', 'Location', 'best');
% % end

% Run simulations and plot each variable in its own figure
for k = 1:nInit
    y0 = init_conditions(k,:);
    [t, y] = ode45(@(t,y) ode_system(t,y,params), tspan, y0);

    % Record final values for extinction check
    Ic_final(k) = y(end,3);
    Ib_final(k) = y(end,5);
    B_final(k) = y(end,6);

    % Plot all variables in separate figures
    for varIdx = 1:6
        figure(figs(varIdx));
        plot(t, y(:,varIdx), 'Color', colors(k,:), 'LineWidth', 1.5);
    end
end

% % Add legends to each figure
% for i = 1:6
%     figure(figs(i));
%     legend(legendStrings, 'Location', 'best');
% end

% Check extinction (below tolerance)
Ic_extinct = sum(Ic_final < tol);
Ib_extinct = sum(Ib_final < tol);
B_extinct = sum(B_final < tol);

fprintf('--- GLOBAL STABILITY OF DFE CHECK ---\n');
fprintf('Number of initial conditions: %d\n', nInit);
fprintf('I_c extinct in %d out of %d initial conditions\n', Ic_extinct, nInit);
fprintf('I_b extinct in %d out of %d initial conditions\n', Ib_extinct, nInit);
fprintf('B extinct in %d out of %d initial conditions\n', B_extinct, nInit);

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
