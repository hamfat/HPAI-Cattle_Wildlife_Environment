clc; clear; close all;

%% ================================================================
%   BASE PARAMETERS
% ================================================================
params.beta_c   = 0.0003034;
params.beta_b   = 3.943e-7;
params.beta_bc  = 1.367e-5;
params.beta_env = 1.587e-5;

params.sigma_c  = 0.5075;
params.gamma_c  = 0.9782;
params.gamma_b  = 1.655e-03;

params.mu_c     = 0.02;
params.mu_b     = 0.02;

params.theta_c  = 0.3853;
params.theta_b  = 0.173;
params.mu_env   = 0.07746;

params.LambdaC  = 30;
params.LambdaB  = 10;

params.mudc     = 0.02;
params.mudb     = 0.02;

%% ================================================================
%   R0 FUNCTION
% ================================================================
compute_R0 = @(p) ...
    (p.beta_c * p.sigma_c * p.LambdaC) ./ ...
    (p.mudc .* (p.sigma_c + p.mu_c + p.mudc) .* (p.gamma_c + p.mu_c + p.mudc)) + ...
    (p.beta_b * p.LambdaB) ./ ...
    (p.mudb .* (p.gamma_b + p.mu_b + p.mudb)) + ...
    (p.LambdaC * p.beta_env * p.theta_c) ./ ...
    (p.mudc .* (p.sigma_c + p.mu_c + p.mudc) .* (p.gamma_c + p.mu_c + p.mudc) .* p.mu_env) + ...
    (p.LambdaB * p.beta_env * p.theta_b) ./ ...
    (p.mudb .* (p.gamma_b + p.mu_b + p.mudb) .* p.mu_env);

%% ================================================================
%   PARAMETER PAIRS FOR R0 SURFACES
% ================================================================
pairs = {
    'beta_c',  'theta_c';
    'beta_env','mu_env';
    'theta_b', 'mudb';
};

rangeFrac = 0.25;

%% ================================================================
%   R0 SURF PLOTS
% ================================================================
for k = 1:3
    param1 = pairs{k,1};
    param2 = pairs{k,2};

    p1_base = params.(param1);
    p2_base = params.(param2);

    p1_vals = linspace((1-rangeFrac)*p1_base, (1+rangeFrac)*p1_base, 60);
    p2_vals = linspace((1-rangeFrac)*p2_base, (1+rangeFrac)*p2_base, 60);

    [P1, P2] = meshgrid(p1_vals, p2_vals);

    R0_grid = zeros(size(P1));

    for i = 1:numel(P1)
        p = params;
        p.(param1) = P1(i);
        p.(param2) = P2(i);
        R0_grid(i) = compute_R0(p);
    end

    figure('Position',[100 100 1200 400]);
    surf(P1, P2, R0_grid, 'EdgeColor','none');
    xlabel(['$', strrep(param1,'_','\_'), '$'], 'Interpreter','latex');
    ylabel(['$', strrep(param2,'_','\_'), '$'], 'Interpreter','latex');
    zlabel('$R_0$', 'Interpreter','latex');
    title(['$R_0$ Surface: $', strrep(param1,'_','\_'), '$ vs $', ...
           strrep(param2,'_','\_'), '$'], 'Interpreter','latex');
    colormap turbo; colorbar; shading interp;
end


%% =====================================================================
%   ODE SYSTEM
% =====================================================================
function dydt = ode_system(~, y, params)

    S_c = y(1); 
    E_c = y(2); 
    I_c = y(3);
    S_b = y(4); 
    I_b = y(5); 
    B   = y(6);

    dS_c = params.LambdaC ...
           - params.beta_c*S_c*I_c ...
           - params.beta_bc*S_c*I_b ...
           - params.beta_env*S_c*B ...
           - params.mu_c*S_c;

    dE_c = params.beta_c*S_c*I_c ...
           + params.beta_bc*S_c*I_b ...
           + params.beta_env*S_c*B ...
           - params.sigma_c*E_c ...
           - params.mu_c*E_c;

    dI_c = params.sigma_c*E_c ...
           - params.gamma_c*I_c ...
           - params.mu_c*I_c;

    dS_b = params.LambdaB ...
           - params.beta_b*S_b*I_b ...
           - params.beta_env*S_b*B ...
           - params.mu_b*S_b;

    dI_b = params.beta_b*S_b*I_b ...
           + params.beta_env*S_b*B ...
           - params.gamma_b*I_b ...
           - params.mu_b*I_b;

    dB   = params.theta_c*I_c ...
           + params.theta_b*I_b ...
           - params.mu_env*B;

    dydt = [dS_c; dE_c; dI_c; dS_b; dI_b; dB];
end
