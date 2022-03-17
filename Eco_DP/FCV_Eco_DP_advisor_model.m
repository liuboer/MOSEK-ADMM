function [X, C, I, out] = FCV_Eco_DP_advisor_model(inp,par)

%% define varibles
delta_s = par.delta_s;
veh_spd = inp.X{1};
veh_acc = inp.U{1};

% Cycle dynamics 
X{1} = sqrt(max(inp.X{1}.^2 + 2.*veh_acc.*delta_s,0));    % current velocity, m/s
veh_spd_mean = ((inp.X{1}+X{1})/2);
delta_t = delta_s ./ veh_spd_mean;
X{2} = inp.X{2} + delta_t;              % current time, s

traf_t = 0;
inf_traff = zeros(size(X{1}));
if inp.W{1} == 0
    inf_traff = (X{1} > inp.W{2});    % speed limit
elseif inp.W{1} == 1
    traf_t = (X{2} + inp.W{4} - floor((X{2}+inp.W{4})/inp.W{2})*inp.W{2}) .* inp.W{1};
    inf_traff = (traf_t <= inp.W{3}+1); % no running a red light
elseif inp.W{1} == 2
    inf_traff = X{1} > 1;             % speed limit
end
%% parametas
% Vehicle parametas
veh_whl_radius = 0.282;
veh_mass = 1380;
veh_rrc  = 0.009;
veh_air_density = 1.2;
veh_FA = 2.0;
veh_CD = 0.335;
veh_gravity = 9.81;
veh_fd_ratio = 6.67;

% Motor parameters
% efficiency map indexed vertically by mc_map_spd and horizontally by mc_map_trq
mc_map_spd = (0:1000:10000) * (2 * pi) / 60; % motor speed list (rad/s)
mc_map_trq = (-200:20:200) * 4.448 / 3.281; % motor torque list (Nm)
mc_eff_map=[...
0.7 	0.7     0.7     0.7     0.7 	0.7     0.7 	0.7     0.7     0.7     0.7	0.7     0.7     0.7 	0.7 	0.7     0.7 	0.7 	0.7     0.7 	0.7
0.78	0.78	0.79	0.8 	0.81	0.82	0.82	0.82	0.81	0.77	0.7	0.77	0.81	0.82	0.82	0.82	0.81	0.8 	0.79	0.78	0.78
0.85	0.86	0.86	0.86	0.87	0.88	0.87	0.86	0.85	0.82	0.7	0.82	0.85	0.86	0.87	0.88	0.87	0.86	0.86	0.86	0.85
0.86	0.87	0.88	0.89	0.9     0.9 	0.9     0.9     0.89	0.87	0.7	0.87	0.89	0.9 	0.9     0.9     0.9     0.89	0.88	0.87	0.86
0.81	0.82	0.85	0.87	0.88	0.9 	0.91	0.91	0.91	0.88	0.7	0.88	0.91	0.91	0.91	0.9     0.88	0.87	0.85	0.82	0.81
0.82	0.82	0.82	0.82	0.85	0.87	0.9     0.91	0.91	0.89	0.7	0.89	0.91	0.91	0.9     0.87	0.85	0.82	0.82	0.82	0.82
0.79	0.79	0.79	0.78	0.79	0.82	0.86	0.9     0.91	0.9     0.7	0.9     0.91	0.9     0.86	0.82	0.79	0.78	0.79	0.79	0.79
0.78	0.78	0.78	0.78	0.78	0.78	0.8     0.88	0.91	0.91	0.7	0.91	0.91	0.88	0.8     0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.8     0.9     0.92	0.7	0.92	0.9     0.8     0.78	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.88	0.92	0.7	0.92	0.88	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78
0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.8 	0.92	0.7	0.92	0.8     0.78	0.78	0.78	0.78	0.78	0.78	0.78	0.78];
% max torque curve of the motor indexed by mc_map_spd
mc_max_trq = [200, 200, 200, 175.2, 131.4, 105.1, 87.6, 75.1, 65.7, 58.4, 52.5] * 4.448/3.281; % (N*m)
mc_max_gen_trq = -1 * [200, 200, 200, 175.2, 131.4, 105.1, 87.6, 75.1, 65.7, 58.4, 52.5] * 4.448/3.281; % (N*m), estimate

%% Update
% Wheel speed (rad/s) & torque (Nm)
w_whl = veh_spd ./ veh_whl_radius;
F_roll = (veh_spd > 0) .* veh_mass .* veh_gravity .* veh_rrc;
F_drag = 0.5 .* veh_air_density .* veh_FA .* veh_CD .* (veh_spd.^2);
F_acc = veh_mass .* veh_acc;
T_whl = veh_whl_radius .* (F_acc + F_roll + F_drag);
P_dem_m = w_whl .* T_whl;

% Motor
mc_spd = w_whl .* veh_fd_ratio;
mc_trq = T_whl ./ veh_fd_ratio;
mc_eff = (mc_spd == 0) + (mc_spd ~= 0) .* interp2(mc_map_trq, mc_map_spd, mc_eff_map, mc_trq, mc_spd);
inf_mc = (isnan(mc_eff)) + (mc_trq < 0) .* (mc_trq < interp1(mc_map_spd, mc_max_gen_trq, mc_spd, 'linear', 'extrap')) + ...
    (mc_trq >= 0) .* (mc_trq > interp1(mc_map_spd, mc_max_trq, mc_spd, 'linear', 'extrap'));
mc_eff(isnan(mc_eff)) = 1;
mc_outpwr = mc_trq .* mc_spd;
mc_inpwr =  mc_outpwr .* (mc_eff.^(-sign((mc_outpwr))));

% Infeasiable
I = (inf_mc + inf_traff~= 0);

% Cost
cost_pwr = abs(mc_inpwr)/10^5;
cost_acc = abs(veh_acc)/3;
% C{1} = (F_acc + F_roll + F_drag./veh_acc*par.v_mean).* (delta_s); %
% C{1} = (cost_pwr+cost_acc) .* delta_t; % good
C{1} = abs(mc_inpwr);

% Output
out.P_dem_m = P_dem_m;
out.P_dem_e = mc_inpwr;
out.Mot_spd = mc_spd;
out.Mot_trq = mc_trq;
out.Mot_pwr = mc_outpwr;
out.Mot_eta = mc_eff;
out.veh_acc = veh_acc;
out.Inf_tot = I;
out.cost_pwr = cost_pwr;
out.cost_acc = cost_acc;
end
