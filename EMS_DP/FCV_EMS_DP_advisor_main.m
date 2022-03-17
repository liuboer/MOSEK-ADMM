% DP-based EMS of an FCV.
clear all;close all;clc

model = 'FCV_EMS_DP_advisor_model';
par.is_convex = false;

load ..\Eco_MOSEK\Eco_v_MOSEK.mat
load ..\Eco_DP\Eco_v_DP.mat

EMS_res = zeros(8,3);% time,SOC_f,HC
EMS_ux = cell(8,4); % SOC,P_fcs,P_mot_e,HC

% Eco_v = 'DP';
Eco_v = 'MOSEK';

for ii = [5,8,6]

if strcmp(Eco_v,'DP')
veh_spd = Eco_v_DP{ii};
else
veh_spd = Eco_v_MOSEK{ii};
end
veh_spd = reshape(veh_spd,[],1);
veh_acc = [0;diff(veh_spd)];
N = length(veh_spd);

SOC_0 = 0.60;
SOC_T = 0.60;
deta  = 0.001;
soc_min = 0.55;
soc_max = 0.65;
fc_pwr_min = 0;
fc_pwr_max = 50*1000;

% mesh grid
clear grd

grd.Nx{1}    = 1000 + 1;    % state variable, SOC
grd.X0{1}    = SOC_0;
grd.Xn{1}.hi = soc_max; 
grd.Xn{1}.lo = soc_min;
grd.XN{1}.hi = SOC_T+deta;
grd.XN{1}.lo = SOC_T;

grd.Nu{1}    = 100 + 1;   % control variable, P_fcs
grd.Un{1}.hi = fc_pwr_max;
grd.Un{1}.lo = fc_pwr_min;

% define problem
clear prb
prb.W{1} = veh_spd; 
prb.W{2} = veh_acc; 
prb.Ts = 1;
prb.N  = (N-1)*1/prb.Ts + 1;

tic
% set options
options = dpm();
options.MyInf = 1e4;
options.BoundaryMethod = 'none'; % also possible: 'none' or 'LevelSet';

[res, dyn] = dpm(model,par,grd,prb,options);

%% save resluts
TOC = toc;
SOC = res.X{1,1};
H2_tot = sum(res.FC_fuel);
fprintf('SOC_0: %.4f , SOC_T: %.4f , H2 Consumption: %.4f g.\n', SOC_0, SOC(end), H2_tot);

EMS_res(ii,:) = [TOC,SOC(end),H2_tot];
EMS_ux{ii,1} = SOC;
EMS_ux{ii,2} = res.FCS_pwr;
EMS_ux{ii,3} = res.P_dem_e;
EMS_ux{ii,4} = res.FC_fuel;

if ~par.is_convex
save(['Eco_',Eco_v,'_EMS_DP_res'],"EMS_res");
save(['Eco_',Eco_v,'_EMS_DP_ux'],'EMS_ux');
else
save(['Eco_',Eco_v,'_EMS_DP_res_convex'],"EMS_res");
save(['Eco_',Eco_v,'_EMS_DP_ux_convex'],'EMS_ux');
end
%% plot
figure
subplot(211)
yyaxis left
plot(veh_spd);
ylabel('Speed [m/s]')
yyaxis right
plot(SOC(2:end));grid on
xlabel('Time [-]')
ylabel('SOC [-]')
xlim([0 N])

subplot(212)
plot(res.FCS_pwr/1000, 'r-');hold on
plot(res.P_dem_e/1000, 'k-');hold off;grid on;
xlabel('Time [s]')
ylabel('power [kW]')
legend('FCS_{pwr}','P_{dem}');
xlim([0 N])

if ~par.is_convex
saveas(gcf,['Eco_Scenario' num2str(ii) '_' Eco_v '_EMS_DP.png'])
else
saveas(gcf,['Eco_Scenario' num2str(ii) '_' Eco_v '_EMS_DP_convex.png'])
end

end
