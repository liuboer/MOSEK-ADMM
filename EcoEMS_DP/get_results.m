% get results for paper
clear;clc;close all

load FCV_EcoEMS_DP_advisor_results\Results_all.mat

EcoEMS_v_DP = cell(8,1);
EcoEMS_res = zeros(8,3); % time,SOC_f,HC
EcoEMS_ux = cell(8,4); % SOC,P_fcs,P_mot_e,HC

for ii = [5,8,6]
EcoEMS_v_DP{ii,1} = Results_all{ii,1}.X{1,1};
EcoEMS_res(ii,:) = Results_all{ii,1}.results;
EcoEMS_ux{ii,1} = [0.6,Results_all{ii,1}.X{1,3}];
EcoEMS_ux{ii,2} = [0,Results_all{ii,1}.FCS_pwr];
EcoEMS_ux{ii,3} = [0,Results_all{ii,1}.P_dem_e];
EcoEMS_ux{ii,4} = [0,Results_all{ii,1}.FC_fuel];
end

save('EcoEMS_v_DP_convex',"EcoEMS_v_DP");
save('EcoEMS_DP_res_convex',"EcoEMS_res");
save('EcoEMS_DP_ux_convex',"EcoEMS_ux");
