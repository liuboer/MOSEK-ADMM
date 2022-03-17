function Scenario = get_scenarios_by_index(S_idx,delta_s,varargin)
p = inputParser;
addParameter(p,'s_loss',0);
parse(p,varargin{:});
switch S_idx
    case 1
        s_unit = 200;
        num_light = 6;
        s_max = s_unit * (2*num_light-1) + p.Results.s_loss;
        Tf_lt = [s_unit:2*s_unit:s_max] / delta_s; % traffic light position
        Scenario.s_list = Tf_lt * delta_s;
        Scenario.s_max = Scenario.s_list(end);
        Scenario.T_0 = [10, 50, 10, 20, 40, 30];
        Scenario.T = 60 * ones(1,num_light);
        Scenario.T_r = Scenario.T - 30 * ones(1,num_light);
    case 5
        s_max = 2200 + p.Results.s_loss;
        Tf_lt = [250, 900, 1300, 1650, s_max] / delta_s; % traffic light position
        Scenario.s_list = Tf_lt * delta_s;
        Scenario.s_max = Scenario.s_list(end);
        Scenario.T_0 = [20, 30, 40, 10, 0];
        Scenario.T = [40, 60, 50, 55, 65];
        Scenario.T_r = Scenario.T - [25, 40, 20, 30, 30];
    case 6
        s_max = 3000 + p.Results.s_loss;
        Tf_lt = [300, 700, 1000, 1700, 2100, 2500, s_max] / delta_s; % traffic light position
        Scenario.s_list = Tf_lt * delta_s;
        Scenario.s_max = Scenario.s_list(end);
        Scenario.T_0 = [40, 10, 20, 50, 25, 15, 50];
        Scenario.T = [60, 50, 55, 70, 45, 60, 60];
        Scenario.T_r = Scenario.T - [35, 20, 30, 40, 25, 30, 20];
    case 7
        s_max = 5000 + p.Results.s_loss;
        Tf_lt = [300, 750, 1000, 1500, 1800, 2500, 3000, 3300, 4200, 4500, s_max] / delta_s; % traffic light position
        Scenario.s_list = Tf_lt * delta_s;
        Scenario.s_max = Scenario.s_list(end);
        Scenario.T_0 = [10, 50, 10, 20, 40, 30, 10, 50, 10, 20, 40];
        Scenario.T = [65, 60, 55, 60, 60, 60, 70, 60, 50, 60, 65];
        Scenario.T_r = Scenario.T - [40, 30, 35, 35, 30, 20, 35, 25, 30, 30, 40];
    case 8
        s_max = 2600 + p.Results.s_loss;
        Tf_lt = [300, 600, 1000, 1300, 2300, s_max] / delta_s; % traffic light position
        Scenario.s_list = Tf_lt * delta_s;
        Scenario.s_max = Scenario.s_list(end);
        Scenario.T_0 = [50, 10, 50, 10, 20, 40];
        Scenario.T = [60, 50, 55, 70, 45, 60];
        Scenario.T_r = Scenario.T - [35, 20, 30, 40, 25, 30];
end

end