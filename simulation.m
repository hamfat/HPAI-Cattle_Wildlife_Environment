clc; clear; close all;

% === Define the base parameters ===
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
params.theta_b = .170;
params.mu_env = 0.07746;
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


%% === Sweep beta_c and compute time-to-peak ===

% === Sweep beta_c and compute time-to-peak ===
beta_vals    = linspace(0.0007, 0.005, 50);   % adjust range if needed
time_to_peak = zeros(size(beta_vals));
peak_Ic      = zeros(size(beta_vals));

for k = 1:length(beta_vals)
    p = params;                 % fresh copy each time
    p.beta_c = beta_vals(k);    % update beta_c
    
    [t, y] = ode45(@(t, y) odefun(t, y, p), [0 500], y0);  % extend time
    Ic = y(:,3);
    
    % Ignore trivial “peak” at t=0 by forcing at least some growth
    if max(Ic) < 1e-3
        time_to_peak(k) = NaN;  % essentially no outbreak
        peak_Ic(k)      = max(Ic);
        continue
    end
    
    [peak_Ic(k), idx] = max(Ic);
    time_to_peak(k)   = t(idx);
end

% === Plot time-to-peak vs beta_c ===
figure;
plot(beta_vals, time_to_peak, '-', 'LineWidth', 3)
xlabel('\beta_c')
ylabel('Time to peak I_c (days)')
title('Time to peak of I_c vs \beta_c')
grid on

% % (Optional) peak height vs beta_c
% figure;
% plot(beta_vals, peak_Ic, 'o-', 'LineWidth', 2)
% xlabel('\beta_c')
% ylabel('Peak I_c')
% title('Peak I_c vs \beta_c')
% grid on

%% === Sweep beta_c for several gamma_c values ===
beta_vals    = linspace(0.0007, 0.005, 50); 
gamma_vals = [0.6, 0.8, 0.9, 1.2];   % choose any set you want

time_to_peak = zeros(length(gamma_vals), length(beta_vals));

for g = 1:length(gamma_vals)
    for k = 1:length(beta_vals)

        p = params;   % fresh copy
        p.beta_c  = beta_vals(k);
        p.gamma_c = gamma_vals(g);

        [t, y] = ode45(@(t,y) odefun(t,y,p), [0 500], y0);
        Ic = y(:,3);

        if max(Ic) < 1e-3
            time_to_peak(g,k) = NaN;
            continue
        end

        [~, idx] = max(Ic);
        time_to_peak(g,k) = t(idx);
    end
end

%% === Plot ===
figure; hold on
colors = lines(length(gamma_vals));

for g = 1:length(gamma_vals)
    plot(beta_vals, time_to_peak(g,:), 'LineWidth', 2, 'Color', colors(g,:))
end

xlabel('\beta_c')
ylabel('Time to peak I_c (days)')
title('Time-to-peak vs \beta_c for different \gamma_c')
legend(arrayfun(@(x) sprintf('\\gamma_c = %.2f', x), gamma_vals, 'UniformOutput', false))
grid on

beta_vals    = linspace(0.001, 0.005, 50); 
gamma_vals = linspace(0.9, 1.4, 50);

T = zeros(length(gamma_vals), length(beta_vals));

for i = 1:length(gamma_vals)
    for j = 1:length(beta_vals)

        p = params;
        p.beta_c  = beta_vals(j);
        p.gamma_c = gamma_vals(i);

        [t, y] = ode45(@(t,y) odefun(t,y,p), [0 500], y0);
        Ic = y(:,3);

        if max(Ic) < 1e-3
            T(i,j) = NaN;
            continue
        end

        [~, idx] = max(Ic);
        T(i,j) = t(idx);
    end
end

figure;
imagesc(beta_vals, gamma_vals, T)
set(gca, 'YDir', 'normal')
colorbar
xlabel('\beta_c')
ylabel('\gamma_c')
title('Heatmap of Time-to-Peak of I_c')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

function TTP = compute_TTP(beta_c, gamma_c, params, y0)
    p = params;
    p.beta_c  = beta_c;
    p.gamma_c = gamma_c;

    [t, y] = ode45(@(t,y) odefun(t,y,p), [0 500], y0);
    Ic = y(:,3);

    if max(Ic) < 1e-3
        TTP = NaN;   % no outbreak
        return
    end

    [~, idx] = max(Ic);
    TTP = t(idx);
end

%% === 1D Plot: Gamma_c vs TTP ===
gamma_vals = linspace(0.9, 1.4, 100);
TTP_gamma  = zeros(size(gamma_vals));

fixed_beta = params.beta_c;   % baseline beta_c

for i = 1:length(gamma_vals)
    TTP_gamma(i) = compute_TTP(fixed_beta, gamma_vals(i), params, y0);
end

figure;
plot(gamma_vals, TTP_gamma, 'LineWidth', 2)
xlabel('\gamma_c')
ylabel('Time to Peak (days)')
title('Time-to-Peak vs \gamma_c (for fixed \beta_c)')
grid on

%% === 2D Surface: Beta_c vs Gamma_c vs TTP ===

beta_vals    = linspace(0.001, 0.005, 50); 
gamma_vals = linspace(0.9, 1.4, 50);

TTP_surf = zeros(length(gamma_vals), length(beta_vals));

for i = 1:length(gamma_vals)
    for j = 1:length(beta_vals)
        TTP_surf(i,j) = compute_TTP(beta_vals(j), gamma_vals(i), params, y0);
    end
end

figure;
surf(beta_vals, gamma_vals, TTP_surf, 'EdgeColor', 'none')
xlabel('\beta_c')
ylabel('\gamma_c')
zlabel('Time to Peak (days)')
title('Surface of Time-to-Peak over (\beta_c, \gamma_c)')
colorbar
view(45, 30)