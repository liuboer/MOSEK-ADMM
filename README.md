# MOSEK-ADMM

Code for the paper "Bi-level Convex Optimization of Eco-driving for Connected Fuel Cell Hybrid Electric Vehicles Through Signalized Intersections".

## Dependencies

- Matlab
- CVX

## Code structure

```
|-- AMOSEK-ADMM
    |-- README.md
    |-- set_path.m
    |-- EcoEMS_DP
    |   |-- EcoEMS_DP_res_convex.mat
    |   |-- EcoEMS_DP_ux_convex.mat
    |   |-- EcoEMS_v_DP_convex.mat
    |   |-- FCV_EcoEMS_DP_advisor_main.m
    |   |-- FCV_EcoEMS_DP_advisor_model.m
    |   |-- get_results.m
    |   |-- FCV_EcoEMS_DP_advisor_results
    |       |-- Results_all.mat
    |       |-- Scenario5Eco.png
    |       |-- Scenario5EMS.png
    |       |-- Scenario6Eco.png
    |       |-- Scenario6EMS.png
    |       |-- Scenario8Eco.png
    |       |-- Scenario8EMS.png
    |-- Eco_DP
    |   |-- ComputTime.mat
    |   |-- Eco_Scenario5_DP.png
    |   |-- Eco_Scenario6_DP.png
    |   |-- Eco_Scenario8_DP.png
    |   |-- Eco_v_DP.mat
    |   |-- FCV_Eco_DP_advisor_main.m
    |   |-- FCV_Eco_DP_advisor_model.m
    |-- Eco_MOSEK
    |   |-- ComputTime.mat
    |   |-- Eco_Scenario1_MOSEK.png
    |   |-- Eco_Scenario5_MOSEK.png
    |   |-- Eco_Scenario6_MOSEK.png
    |   |-- Eco_Scenario8_MOSEK.png
    |   |-- Eco_v_MOSEK.mat
    |   |-- MOSEK_main.m
    |   |-- s_ref.mat
    |-- EMS_ADMM
    |   |-- ADMM.m
    |   |-- ADMM_main.m
    |   |-- ADMM_main_soc_rho1_for_plot.m
    |   |-- Eco_DP_EMS_ADMM_res.mat
    |   |-- Eco_DP_EMS_ADMM_ux.mat
    |   |-- Eco_MOSEK_EMS_ADMM_res.mat
    |   |-- Eco_MOSEK_EMS_ADMM_ux.mat
    |   |-- Eco_Scenario5_DP_EMS_ADMM.png
    |   |-- Eco_Scenario5_MOSEK_EMS_ADMM.png
    |   |-- Eco_Scenario6_DP_EMS_ADMM.png
    |   |-- Eco_Scenario6_MOSEK_EMS_ADMM.png
    |   |-- Eco_Scenario8_DP_EMS_ADMM.png
    |   |-- Eco_Scenario8_MOSEK_EMS_ADMM.png
    |   |-- f_BacktrackingNewtonVector.m
    |   |-- SOC_f_rho1.mat
    |   |-- SOC_f_rho1_NEDC.png
    |-- EMS_DP
    |   |-- Eco_DP_EMS_DP_res.mat
    |   |-- Eco_DP_EMS_DP_res_convex.mat
    |   |-- Eco_DP_EMS_DP_ux.mat
    |   |-- Eco_DP_EMS_DP_ux_convex.mat
    |   |-- Eco_MOSEK_EMS_DP_res.mat
    |   |-- Eco_MOSEK_EMS_DP_res_convex.mat
    |   |-- Eco_MOSEK_EMS_DP_ux.mat
    |   |-- Eco_MOSEK_EMS_DP_ux_convex.mat
    |   |-- Eco_Scenario5_DP_EMS_DP.png
    |   |-- Eco_Scenario5_DP_EMS_DP_convex.png
    |   |-- Eco_Scenario5_MOSEK_EMS_DP.png
    |   |-- Eco_Scenario5_MOSEK_EMS_DP_convex.png
    |   |-- Eco_Scenario6_DP_EMS_DP.png
    |   |-- Eco_Scenario6_DP_EMS_DP_convex.png
    |   |-- Eco_Scenario6_MOSEK_EMS_DP.png
    |   |-- Eco_Scenario6_MOSEK_EMS_DP_convex.png
    |   |-- Eco_Scenario8_DP_EMS_DP.png
    |   |-- Eco_Scenario8_DP_EMS_DP_convex.png
    |   |-- Eco_Scenario8_MOSEK_EMS_DP.png
    |   |-- Eco_Scenario8_MOSEK_EMS_DP_convex.png
    |   |-- FCV_EMS_DP_advisor_main.m
    |   |-- FCV_EMS_DP_advisor_model.m
    |-- Functions
        |-- get_scenarios_by_index.m
        |-- get_s_ref.m
        |-- IDM.m
        |-- NEDC.mat
        |-- plot_scenario.m
```