% DP-based Eco-driving of an FCV at signalized intersections.
% 2022-01-17, Liu Bo

clear;clc;close all

tic
model = 'FCV_EcoEMS_DP_advisor_model';
par.is_convex = true;

Results_all = cell(8,1);

for ii = [5,8,6]
par.delta_s = 1;
Scenario = get_scenarios_by_index(ii,par.delta_s);
s_ref = get_s_ref(Scenario);

spd_max = 60 / 3.6;
spd_min = 0;
acc_max = 1.4;
acc_min = -2;
SoC_0 = 0.60;
SoC_T = 0.60;
deta  = 0.001;
soc_min = 0.591;
soc_max = 0.601;
fc_pwr_min = 0;
fc_pwr_max = 50*1000;

% State variables
grd.Nx{1}    = 40 + 1;    % speed, m/s
grd.X0{1}    = 0;
grd.Xn{1}.hi = spd_max;
grd.Xn{1}.lo = spd_min;
grd.XN{1}.hi = spd_min+1;
grd.XN{1}.lo = spd_min;

grd.Nx{2}    = 100 + 1;    % distance, s
grd.X0{2}    = 0;
grd.Xn{2}.hi = [s_ref.s_upper_ref(1);s_ref.s_upper_ref];
grd.Xn{2}.lo = [0;s_ref.s_lower_ref];
grd.XN{2}.hi = s_ref.s_upper_ref(end);
grd.XN{2}.lo = s_ref.s_lower_ref(end);

grd.Nx{3}    = 100 + 1;    % SOC, -
grd.X0{3}    = SoC_0;
grd.Xn{3}.hi = soc_max; 
grd.Xn{3}.lo = soc_min;
grd.XN{3}.hi = SoC_T+deta;
grd.XN{3}.lo = SoC_T;

% Control variables
grd.Nu{1} = 20 + 1;   % acceleration, m/s^2
grd.Un{1}.hi = acc_max;
grd.Un{1}.lo = acc_min;

grd.Nu{2} = 10 + 1;   % P_fcs, W
grd.Un{2}.hi = fc_pwr_max;
grd.Un{2}.lo = fc_pwr_min;

% Problem
prb.Ts = par.delta_s;
prb.N  = s_ref.N;
prb.N0 = 1;

options = dpm();
options.MyInf = 1e10;
options.BoundaryMethod = 'none'; % also possible: 'none' or 'LevelSet';

[res,dyn] = dpm(model,par,grd,prb,options);

TOC = toc;
SOC = res.X{3};
H2 = sum(res.FC_fuel);
fprintf('SOC_0: %.4f , SOC_T: %.4f , H2 Consumption: %.4f g.\n', SoC_0, SOC(end), H2);

res.results = [TOC,SOC(end),H2];
Results_all{ii,1} = res;

save('FCV_EcoEMS_DP_advisor_results/Results_all','Results_all')

disp([sum(abs(res.P_dem_e)/10^5*1.5),sum(res.FC_fuel),sum(abs(res.veh_acc)/5)]);
%% plot Eco
LineWidth = 1;

figure
subplot(211) % distance-time
N = s_ref.N;
tt = 1:N;
vv = res.X{1}(2:end);
ss = res.X{2}(2:end);
soc = res.X{3}(2:end);
plot(tt,ss,'b-','linewidth',LineWidth); hold on
s_upper_ref = s_ref.s_upper_ref;
s_lower_ref = s_ref.s_lower_ref;
plot(tt,s_upper_ref,'r-','linewidth',LineWidth); hold on
plot(tt,s_lower_ref,'r-','linewidth',LineWidth); hold on
plot_scenario(Scenario,N)
axis([0 N 0 ss(end)]); grid on; hold off
ylabel('Distance (m)'); xlabel('Time (s)')

subplot(212) % speed-time
plot(tt,vv,'b-','linewidth',LineWidth)
axis([0 N 0 spd_max]); grid on
xlabel('Time (s)'); ylabel('Speed (m/s)')

saveas(gcf,['FCV_EcoEMS_DP_advisor_results/Scenario' num2str(ii) 'Eco.png'])
%% plot EMS
figure
subplot(211)
yyaxis left
plot(tt,vv,'b-','linewidth',LineWidth)
ylabel('Speed (m/s)')
yyaxis right
plot(tt,soc,'r-','linewidth',LineWidth)
ylabel('SOC (-)')
axis([0 N 0 10])
axis 'auto y'

subplot(212)
plot(tt,res.FCS_pwr/1000,'r-','linewidth',LineWidth)
hold on
plot(tt,res.Bat_pwr/1000,'b-','linewidth',LineWidth)
hold on
plot(tt,res.P_dem_m/1000, 'k-');
ylabel('power <kW>')
legend('FCS_{pwr}','Bat_{pwr}','P_{dem}','location','southwest');
axis([0 N 0 10])
axis 'auto y'
grid on

saveas(gcf,['FCV_EcoEMS_DP_advisor_results/Scenario' num2str(ii) 'EMS.png'])

end