clear
clc



%% Problem Formulation
tStart = tic;

%% Parameters
Q = 3000;
Ta_in = 80;
rho_a = 1.225;
cp_air = 1005;


Tw_in = 32;
rho_w = 1000;
cp_water = 4.184;
Vdot_w = 0.00004921035319468;


%% Properties
Tube_Depth = 8/1000; %m
number_of_tubes = 8;

price = 2.50;
k = 237;
Tube_Length = 200/1000; %mm



%% Monte Carlo
n = input("Number of Points:");

%% Assigning Bounds
%Thickness
Tube_Height_lb = 1/1000; %m
Tube_Height_ub = 10/1000; %m

%Height
Tube_Width_lb = 1/1000; %m
Tube_Width_ub = 40/1000; %m

%% Randomly Sample n number of x1 and x2 within the bounds

Tube_Height = Tube_Height_lb + (Tube_Height_ub-Tube_Height_lb).*rand(n,1);
Tube_Width = Tube_Width_lb + (Tube_Width_ub-Tube_Width_lb).*rand(n,1);


%% Inquality Constraint

%Air Side
mdot_air = .015;
rho = 1.225;
Pfin = 2/1000;
Gap = 1.5/1000;
tfin = 0.178/1000;
twall = 0.35/1000;
Dhyd = 2.5/1000;
Pitch = 1;
Afr = 0.005609;
Aff = Afr*.82;
N = 10;
L_tube = 79/1000;
mu = 0.00001837;
cp_air = 1.006*1000;
kair = .026;
Pr = mu*cp_air/kair;
h_tube = 1.7/1000;
rho_water = 1000;
Acr = ((16-0.7)/1000)*(0.001);
G_air= mdot_air/Aff;
Re_air = G_air*Dhyd/mu;
j = [0.045 .023 .015 .011];
h_air = j*cp_air.*G_air/(Pr^(2/3));
Cair = mdot_air*cp_air;





%Water Side
mdot_water = rho_water*Vdot_w;
Dh_w = 4*(1/1000*15.3/1000)/(2*(1/1000+15.3/1000));
Pr_water = 3.56;
Cp_water = 4184;
A_tube = (15.3/1000)*1/1000*4;
v_water = Vdot_w/A_tube;
mu_water=0.0005383;
mu_w = 8.9*10^-4;
P_tube = ((15.3/1000)+1/1000)*8;
k_water = 0.64;
Dh_water = 2*A_tube/(P_tube);
Re_water = rho_water*Dh_w*v_water/mu_water;
Nu_water = (0.023*Re_water^0.8)*Pr_water^0.33;
Atot = 2*(Tube_Length*Tube_Height+Tube_Length*Tube_Width)*8;
h_water = k_water*Nu_water/Dh_w;
C_water = Cp_water*Vdot_w*1000;


c = Cair/C_water;
U = [((1/h_air(1))+(1/h_water))^-1 ((1/h_air(2))+(1/h_water))^-1 ((1/h_air(3))+(1/h_water))^-1 ((1/h_air(4))+(1/h_water))^-1];
NTU = Atot*U./Cair;

efficiency = 1-exp(NTU.^.22/c.*(exp(-c*NTU.^.78)-1)) ;
Qmax = Cair*(Tw_in-Ta_in);
Q = Qmax.*efficiency;
Ta_out = Q./Cair + Ta_in;
Tw_out = Tw_in - Q./C_water;


volume = Tube_Height.*Tube_Width*Tube_Length;
density = 2710;


Cost = price*volume*density;



%% Initiating Monte Carlo Sampling
Count = 0;
for i=1:n
    if Ta_out <= 80 
        Count = Count+1;
        
        %Create array of x1 and x2 that satisfy the inequality
        Height(Count) = Tube_Height(Count);
        Width(Count) = Tube_Width(Count);


        %Create array of resulting objective value
        P(Count) = Cost(i);

    end
end

[mincost, ind] = min(P);
BestHeight = Height(ind);
BestWidth = Width (ind);
Ta_out = Q./Cair + Ta_in;

fprintf("Optimal Heigh is %.2f. mm \n", BestHeight*1000)
fprintf("Optimal Thickness is %.2f mm \n", BestWidth*1000)
fprintf("Optimal Temp is %.2f\n", Ta_out(ind))
fprintf("Optimal Cost is $%.f.4 \n",mincost)
