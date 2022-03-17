% ADMM-based EMS of an FCV.
close all;clear all;clc

load NEDC.mat

speed_vector = NEDC / 3.6;
veh_spd = reshape(speed_vector,[],1);
veh_acc = [0;diff(veh_spd)];
N = size(veh_spd,1);
% Vehicle parametas
veh_whl_radius = 0.282;
veh_mass = 1380;
veh_rrc  = 0.009;
veh_air_density = 1.2;
veh_FA = 2.0;
veh_CD = 0.335;
veh_gravity = 9.81;
veh_fd_ratio = 6.67;
% Wheel speed (rad/s) & torque (Nm)
w_whl = veh_spd ./ veh_whl_radius;
F_roll = (veh_spd > 0) .* veh_mass .* veh_gravity .* veh_rrc;
F_drag = 0.5 .* veh_air_density .* veh_FA .* veh_CD .* (veh_spd.^2);
F_acc = veh_mass .* veh_acc;
T_whl = veh_whl_radius .* (F_acc + F_roll + F_drag);
P_dem_m = w_whl .* T_whl;

% motor
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
mc_spd = w_whl .* veh_fd_ratio;
mc_trq = T_whl ./ veh_fd_ratio;
mc_eff = (mc_spd == 0) + (mc_spd ~= 0) .* interp2(mc_map_trq, mc_map_spd, mc_eff_map, mc_trq, mc_spd);
mc_eff(isnan(mc_eff)) = 1;
mc_outpwr = mc_trq .* mc_spd;
mc_inpwr =  mc_outpwr .* (mc_eff.^(-sign((mc_outpwr))));% .* (1-0.1*(mc_outpwr<0));

% fuel cell system
fc_pwr_min = 0;
fc_pwr_max = 50*1000;
fc_pwr_map = [0, 2, 5, 7.5, 10, 20, 30, 40, 50] * 1000; % % kW (net) including parasitic losses
fc_fuel_map = [0.012, 0.05, 0.085, 0.117, 0.149, 0.280, 0.423, 0.594, 0.821]; % fuel use map (g/s)
% another war to calculate fc_fuel_map
fc_eff_map = [10, 33, 49.2, 53.3, 55.9, 59.6, 59.1, 56.2, 50.8] / 100; % % efficiency indexed by fc_pwr
fc_fuel_lhv = 120.0*1000; % (J/g), lower heating value of the fuel
fc_fuel_map2 = fc_pwr_map .* (1./fc_eff_map) / fc_fuel_lhv; % fuel consumption map (g/s)

[Eta_fcs,gof] = fit(fc_pwr_map',fc_fuel_map','poly2');
coeffs = [Eta_fcs.p1,Eta_fcs.p2,Eta_fcs.p3];

% Battery parameters
Num_cell = 25;
Q_bat = 26 * 3600; % coulombs, battery package capacity
% resistance and OCV list
ess_soc_map = 0:0.1:1;
% module's resistance to being discharged, indexed by ess_soc
ess_r_dis_map = [40.7, 37.0, 33.8, 26.9, 19.3, 15.1, 13.1, 12.3, 11.7, 11.8, 12.2] / 1000 * Num_cell; % (ohm)
% module's resistance to being charged, indexed by ess_soc
ess_r_chg_map = [31.6, 29.8, 29.5, 28.7, 28.0, 26.9, 23.1, 25.0, 26.1, 28.8, 47.2] / 1000 * Num_cell; % (ohm)
% module's open-circuit (a.k.a. no-load) voltage, indexed by ess_soc
ess_voc_map = [11.70, 11.85, 11.96, 12.11, 12.26, 12.37, 12.48, 12.59, 12.67, 12.78, 12.89] * Num_cell; % (V)
% Battery limitations
ess_min_volts = 9.5 * Num_cell; % 237.5 V
ess_max_volts = 16.5 * Num_cell; % 412.5 V
% V_OC and R_0
V_OC = interp1(ess_soc_map, ess_voc_map, 0.6, 'linear', 'extrap');
R_0 = (interp1(ess_soc_map, ess_r_dis_map, 0.6, 'linear', 'extrap')...
    + interp1(ess_soc_map, ess_r_chg_map, 0.6, 'linear', 'extrap')) / 2;

P_OC_min = (V_OC-ess_max_volts)*V_OC / R_0 .* ones(N,1);
P_OC_max = (V_OC-ess_min_volts)*V_OC / R_0 .* ones(N,1);

idx_pos = find(mc_inpwr>0);
idx_neg = find(mc_inpwr<=0);

P_OC_min = max(P_OC_min, f_g_inv(mc_inpwr-fc_pwr_max,V_OC,R_0));
P_OC_max = min(min(P_OC_max, f_g_inv(mc_inpwr-fc_pwr_min,V_OC,R_0)),V_OC.^2/(2*R_0));

SOC_0 = 0.6;
SOC_min = 0.3;
SOC_max = 0.9;

misc.epsilon = 1E3;
misc.maxIterations = 1000;
rho1_min = 1E-14;
rho1_max = 1E-10;
Flag = 1;
eps = 5e-5;

rho_1_list = 10.^(-14:0.2:-10);
SOC_f_list = zeros(size(rho_1_list));
for i = 1:length(rho_1_list)
rho1 = rho_1_list(i);
[SOC, P_OC, time, iters, iter2] = ADMM(coeffs,mc_inpwr,SOC_0,P_OC_min,P_OC_max,SOC_min,SOC_max,R_0,V_OC,Q_bat,misc,idx_pos,idx_neg,rho1);
SOC_f_list(i) = SOC(end);
h = waitbar(i/length(rho_1_list));
end
P_fcs = mc_inpwr - f_g(P_OC,V_OC,R_0);

delete(h)

SOC_f_rho1 = [SOC_f_list;rho_1_list];
save('SOC_f_rho1','SOC_f_rho1')

% plot
figure
semilogx(rho_1_list,SOC_f_list);grid on
xlabel('\rho_1 [-]')
ylabel('SOC_f [-]')

saveas(gcf,'SOC_f_rho1_NEDC.png')

function g_inv = f_g_inv(P_bat,V_OC,R_0)
g_inv = V_OC.^2./(2*R_0) .* (1-sqrt(1-4*R_0.*P_bat./V_OC.^2));
end
function g = f_g(P_OC,V_OC,R_0)
g = P_OC - R_0./V_OC^2 .* P_OC.^2;
end
