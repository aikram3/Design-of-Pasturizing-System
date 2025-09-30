%% Shell and Tube Optimization
% 12/5/24
% Alejandro Plata, Ryan Smith, and Ben Hong
% Optimize a shelll and tube heat exchanger for cost as a function of
% pressure drop and material cost

%% Clear

clc; clearvars; close all

%% Properties

% [Section 1, Section 2]

cp_cold = [3931	3931];          % J/kg/K
cp_hot = [3931	4187];

rho_cold = [1030 1030];         % kg/m^3
rho_hot = [1030	980.45];    

mu_cold = [0.00213	0.00213];   % Pa*s
mu_hot = [0.00213	0.000434];

Pr_cold = [14.97858676	14.97858676];
Pr_hot = [14.97858676	2.75];

k_cold = [0.559	0.559];         % W/m/K
k_hot = [0.559	0.63];

Tin_cold = [4.44 35.6];         % dC
Tout_cold = [Tin_cold(2) 72.778];

Tin_hot = [72.778 93.3];
% Tout_hot = [2x1];

% Other Properties and Constants

tube_passes = 2;            
thickness = .001;           % m
g = 9.81;                   % m/s^2           
mdot_cold = 4.490744792;    % kg/s

% h_o calcs 

C_coeffs = polyfit([1.25 1.5 2 3],[.386 .278 .254 .22],2);
n_coeffs = polyfit([1.25 1.5 2 3],[0.592 0.62 0.632 0.608],2);
Fo = 0.6;

% Fouling

Rf_o = .00035;  % K/W 
Rf_i = .00035;

% Thermal Resistance (Stainless Steel)

k_ss = 15;  % W/m/K

% Price and Pumping (Stainless Steel)

rho_ss = 8000;      % kg/m^3
price_ss = 1.3;     % $/kg
e_ss = .015*10^-3;  % m

%% Variables
iterations = 1000000;

% Mass flow rate of water
mdot_water_min = 4.5;   % Minimum necessary to satisfy Q requirements
mdot_water_max = 20;    % Arbitrarily set but the optimization usually settles around 7
mdot_hot = (mdot_water_max - mdot_water_min)*rand(iterations,1) + mdot_water_min;
mdot_hot_out = [repmat(mdot_cold,iterations,1) mdot_hot];

% Tube diameter
tube_d_min = .01;   % Min and max are based on heuristics for shell and tube design
tube_d_max = .025;
D_i_tube_out = (tube_d_max - tube_d_min)*rand(iterations,1) + tube_d_min;

% Number of hexagons (0 is a concentric annulus, 1 is a tube surrounded by
% 6 tubes in hexagon pattern, etc.)
num_hex_min = 1;
num_hex_max = 10;
num_hex_out = randi([num_hex_min num_hex_max],iterations,1);

%% Initial Calculations

% Initalizations

L_final = zeros(iterations,2);
cost = zeros(iterations,1);
head = zeros(iterations,1);
head_cost = zeros(iterations,1);
material_cost = zeros(iterations,1);
D_i_shell_out = zeros(iterations,1);

for i = 1:iterations

    % Assigning values for variables

    mdot_hot = mdot_hot_out(i,:);
    D_i_tube = D_i_tube_out(i);
    num_hex = num_hex_out(i);
    num_tubes = 1;
    for j = 1:num_hex
        num_tubes = num_tubes + j* 6;
    end

    % Log Mean Temperature Difference

    Qdot = cp_cold .* mdot_cold .* (Tout_cold - Tin_cold);

    Tout_hot = Tin_hot - Qdot ./ mdot_hot ./ cp_hot;

    dT1 = Tin_hot - Tout_cold;
    dT2 = Tout_hot - Tin_cold;

    dTlm = (dT1 - dT2)./log(dT1./dT2);

    % Correction factor calculations

    R = (Tin_cold - Tout_cold)./(Tin_cold - Tin_hot);
    P = (Tout_hot - Tin_hot)./(Tin_cold - Tin_hot);
    A = 2./P - 1 - R;
    F = sqrt(R.^2+1) ./ (R-1) .* log((1-P)./(1-P.*R)) ./ ...
        log((A+sqrt(R.^2+1))./(A-sqrt(R.^2+1)));

    % h_i calculation
    
    V_i = mdot_hot ./ rho_hot ./ (num_tubes./tube_passes .* D_i_tube .^2 / 4 .* pi);
    Re_i = rho_hot .* V_i .* D_i_tube ./ mu_hot;
    Nu_i = 0.023 .* Re_i .^ 0.8 .* Pr_cold .^ 0.4;

    h_i = Nu_i .* k_hot ./ D_i_tube;

    % h_o calculation
    
        % Resources suggest that the pitch and baflle_spacing should be
        % functions of the tube outer diameter and the shell inner diameter
        % respectively, rather than variables.
    D_o_tube = D_i_tube + thickness*2;
    pitch = 1.25 * D_o_tube;                        
    D_i_shell = (num_hex.*pitch + D_o_tube/2)*1.1*2;
    baffle_spacing = .5 * D_i_shell;
    D_i_shell_out(i) = D_i_shell;

    cross_flow_area = D_i_shell .* baffle_spacing .* (pitch-D_o_tube) ./pitch;
    Vmax = mdot_cold ./ rho_cold ./ cross_flow_area;
    Re_o = D_o_tube .* Vmax .* rho_cold ./ mu_cold;
    
    x = pitch ./ D_o_tube;
    n = n_coeffs(1) .* x.^2 + n_coeffs(2) .* x + n_coeffs(3);
    C = C_coeffs(1) .* x.^2 + C_coeffs(2) .* x + C_coeffs(3);

    h_o = C .* Re_o .^ n .* Pr_hot .^ 1/3 * Fo.* k_cold ./ D_o_tube;


    L_check = zeros(1,2);
    

    for j = 1:2

        L_guess = .01;
        L_check(j) = 10;

        while abs(L_guess-L_check(j)) > 1

            L_guess = L_guess + .9*(L_check(j)-L_guess);
            A_o = D_o_tube .* pi .* L_guess .* 2;
            A_i = D_i_tube .* pi .* L_guess .* 2;
            A_lm = (A_o - A_i) ./ (log(A_o) - log(A_i));

            U = (1./h_o(j) + thickness ./ k_ss .* A_o ./ A_lm + 1./h_i(j) ...
                .* A_o ./ A_i ...
                + Rf_i + Rf_o).^-1;
            A_ht = Qdot(j)./U./dTlm(j)./F(j);
            L_check(j) = A_ht ./ ((D_o_tube .* pi -  D_i_tube .* pi) ./ (log(D_o_tube .* pi) -  log(D_i_tube .* pi)));
        end
    end

    L_tube = L_check ./ (num_tubes);
    L_final(i,:) = L_tube;
        
        % Total volume is the sume of the shell volume and tube volume
    volume = ((D_i_shell + thickness).^2 - D_i_shell.^2) ./ 4 .* pi .* sum(L_final(i,:))...
        + (D_o_tube.^2-D_i_tube.^2) ./ 4 .* pi .* sum(L_final(i,:)) .* num_tubes;

    % Final material cost

    material_cost(i) = volume .* rho_ss .* price_ss;

    % Pumping

    % Tube side pressure drop

    surface_roughness = e_ss./D_i_tube;

    f_corrected = zeros(1,2);

    for j = 1:2

        f = 0;
        f_corrected(j) = 10;
        tick = 0;

        while abs(f_corrected(j)-f) > .01

            f = f+.001;
            f_corrected(j) = (1.14-2.*log10(surface_roughness+9.35./Re_i(j)./sqrt(f))).^-2;

        end
    end

    deltaP_tube = f_corrected .* L_tube.*num_tubes./D_i_tube .* 1/2.* rho_hot .* V_i .^2;

    % Return pressure loss (tube side)
    
    Gt = mdot_hot ./ (num_tubes .* D_i_tube.^2 / 4 * pi / tube_passes);
    deltaP_return = 4 * tube_passes .* Gt ./ 2 ./ rho_hot;

    deltaP = deltaP_return + deltaP_tube;
    head_tube= deltaP./rho_cold./g;

    % Shell side
    
    G_s = mdot_hot ./ cross_flow_area;
    N_b = round(L_final(i) ./ baffle_spacing);
    C_p = 0.86; % (triangular pitch)
    D_e = 4 * (C_p .* pitch - pi .* D_o_tube .^2 / 4) ./ ...
        (pi .* D_o_tube);
    Re_s = D_e .* G_s ./ mu_hot;

    f_corrected = zeros(1,2);

    for j = 1:2

        f = 0;
        f_corrected(j) = 10;
        tick = 0;

        while abs(f_corrected(j)-f) > .01

            f = f+.001;
            f_corrected(j) = (1.14-2.*log10(surface_roughness+9.35./Re_s(j)./sqrt(f))).^-2;

        end
    end

    deltaP_shell = 2 * f .* G_s .^2 .* D_i_shell .* (N_b + 1) ./ ...
    (rho_hot .* D_e .* .9^0.14);
    head_shell= deltaP_shell./rho_hot./g;

    % Total head

    head(i) = sum(head_tube) + sum(head_shell);

    head_cost(i) = head(i) .* 3.281 .* 9.29 + 833;
    
    % Total cost

    cost(i) = head_cost(i) + material_cost(i);

end

[cost,I] = min(cost)
L_final = L_final(I,:)
head = head(I)
head_cost = head_cost(I)
material_cost = material_cost(I)
mdot_hot = mdot_hot_out(I,2)
D_i_tube = D_i_tube_out(I)
D_i_shell = D_i_shell_out(I)

num_hex = num_hex_out(I)
num_tubes = 1;
  for j = 1:num_hex
        num_tubes = num_tubes + j* 6;
  end
 num_tubes