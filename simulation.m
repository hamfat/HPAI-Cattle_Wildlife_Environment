clc; clear; close all;

% === Define the base parameters ===
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
y0 = [1500, 0, 0, 499, 10, 0];  % [S_c, E_c, I_c, S_b, I_b, B]

% === Time span ===
tspan = [0 365];

% === Solve the system of ODEs ===
[t, y] = ode45(@(t, y) odefun(t, y, params), tspan, y0);

% === Extract variables ===
Sc = y(:,1);
Ec = y(:,2);
Ic = y(:,3);
Sb = y(:,4);
Ib = y(:,5);
B  = y(:,6);

% === Figure 1: Cattle Dynamics in (3,1,3) ===
figure;

subplot(3,1,1);
plot(t, Sc, 'b', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('S_c');
title('Cattle Susceptible (S_c)');
xlim([0, 365])
grid on;

subplot(3,1,2);
plot(t, Ec, 'm', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('E_c');
xlim([0, 365])
title('Cattle Exposed (E_c)');
grid on;

subplot(3,1,3);
plot(t, Ic, 'r', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('I_c');
xlim([0, 365])
title('Cattle Infectious (I_c)');
grid on;

% === Figure 2: Bird & Environment Dynamics in (3,2,3) ===
figure;

subplot(3,1,1);
plot(t, Sb, 'g', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('S_b');
xlim([0, 365])
title('Bird Susceptible (S_b)');
grid on;

subplot(3,1,2);
plot(t, Ib, 'r', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('I_b');
xlim([0, 365])
title('Bird Infectious (I_b)');
grid on;

subplot(3,1,3);
plot(t, B, 'k', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('B');
xlim([0, 365])
title('Environmental Pathogen (B)');

grid on;



% === ODE function ===
function dydt = odefun(~, y, p)
    Sc = y(1);
    Ec = y(2);
    Ic = y(3);
    Sb = y(4);
    Ib = y(5);
    B  = y(6);
    
    dSc = p.LambdaC - p.beta_c*Sc*Ic - p.beta_bc*Sc*Ib - p.beta_env*Sc*B - p.muSc*Sc;
    dEc = p.beta_c*Sc*Ic + p.beta_bc*Sc*Ib + p.beta_env*Sc*B - p.sigma_c*Ec - p.mu_c*Ec - p.muEc*Ec;
    dIc = p.sigma_c*Ec - p.gamma_c*Ic - p.mu_c*Ic - p.muIc*Ic;
    dSb = p.LambdaB - p.beta_b*Sb*Ib - p.beta_env*Sb*B - p.muSb*Sb;
    dIb = p.beta_b*Sb*Ib + p.beta_env*Sb*B - p.gamma_b*Ib - p.mu_b*Ib - p.muIb*Ib;
    dB  = p.theta_c*Ic + p.theta_b*Ib - p.mu_env*B;
    
    dydt = [dSc; dEc; dIc; dSb; dIb; dB];
end
