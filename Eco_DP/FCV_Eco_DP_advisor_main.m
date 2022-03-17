% DP-based Eco-driving of an FCV at signalized intersections.
% 2021-0520, Liu Bo
clear;clc;close all

model = 'FCV_Eco_DP_advisor_model';

Eco_v_DP = cell(8,1);
ComputTime = zeros(8,1);

for ii = [5,8,6]

tic

par.delta_s = 5;
Dtnc_list = [2200,nan,nan,nan,2200,3000,5000,2600];
Dtnc = Dtnc_list(ii) / par.delta_s;

spd_max = 60 / 3.6;
acc_max = 1.4;
acc_min = -2;

% Problem
prb.Ts = par.delta_s;  % x meters (distance) per step
prb.N  = Dtnc;
prb.N0 = 1;

Scenario = get_scenarios_by_index(ii,par.delta_s);
Tf_lt = Scenario.s_list / par.delta_s;
s_ref = get_s_ref(Scenario);
par.v_mean = Scenario.s_max / s_ref.N;

% reference IDM
% t_min_list = [0;s_ref.t_lower_ref(1:par.delta_s:end)]-10;
% t_max_list = t_min_list + 50;

% reference v_mean
t_ref = [0:par.delta_s:Scenario.s_max] ./ par.v_mean;
t_min_list = t_ref-30;
t_max_list = t_ref+30;

prb.W{1} = zeros(1,Dtnc);
prb.W{2} = spd_max*ones(1,Dtnc);
prb.W{3} = ones(1,Dtnc);
prb.W{4} = ones(1,Dtnc);

% prb.W1: 0-speed limit, 1-traffic light, 2-stop sign, 3-ped crossing
% prb.W2: 0-limit value, 1-light cycle (red to green)
% prb.W3:                1-red light duration
% prb.W4:                1-offfset

prb.W{1}(1,Tf_lt) = ones(1,length(Tf_lt));
prb.W{2}(1,Tf_lt) = Scenario.T;
prb.W{3}(1,Tf_lt) = Scenario.T_r;
prb.W{4}(1,Tf_lt) = Scenario.T_0;

% State variables
grd.Nx{1} = 40 + 1;   % velocity, m/s
grd.X0{1} = 0;
grd.Xn{1}.hi = spd_max;
grd.Xn{1}.lo = 0;
grd.XN{1}.hi = 1;
grd.XN{1}.lo = 0;

grd.Nx{2} = 300 + 1;   % time, s
grd.X0{2} = 0;
grd.Xn{2}.hi = t_max_list;
grd.Xn{2}.lo = t_min_list;
grd.XN{2}.hi = s_ref.N+0.5;
grd.XN{2}.lo = s_ref.N-0.5;

% Control variables
grd.Nu{1} = 10 + 1;   % acc, m/s^2
grd.Un{1}.hi = acc_max;
grd.Un{1}.lo = acc_min;

options = dpm();
options.MyInf = 1e10;
options.BoundaryMethod = 'none';
[res,dyn] = dpm(model,par,grd,prb,options);

%% save results
ComputTime(ii,1) = toc;
Eco_v_DP{ii,1} = interp1(res.X{1,2},res.X{1,1},0:round(res.X{1,2}(end)),"linear","extrap");

save('ComputTime','ComputTime');
save('Eco_v_DP',"Eco_v_DP");

% a = [sum(res.cost_acc),sum(res.cost_pwr)];
%% plot
LineWidth = 1;

figure
subplot(211) % distance-time
vv_s = res.X{1};
tt_s = res.X{2};
tt = 1:ceil(tt_s(end));
ss = interp1(tt_s,1:length(tt_s),tt,'linear','extrap') * par.delta_s;
plot(tt,ss,'b-','linewidth',LineWidth); hold on
plot_scenario(Scenario,s_ref.N)
axis([0 tt(end) 0 Scenario.s_max]); grid on; hold off
ylabel('Distance (m)'); xlabel('Time (s)')

subplot(212) % speed-time
plot(tt_s,vv_s,'b-','linewidth',LineWidth)
axis([0 tt(end) 0 spd_max+2]); grid on
xlabel('Time (s)'); ylabel('Speed (m/s)')

saveas(gcf,['Eco_Scenario' num2str(ii) '_DP.png'])
end
