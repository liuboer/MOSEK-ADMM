function plot_scenario(Scenario,N)
s_list = Scenario.s_list;
T = Scenario.T;
T_0 = Scenario.T_0;
T_r = Scenario.T_r;
for jj = 1:length(s_list)
    s = s_list(jj);
    for ii = 1:ceil(N/T(1,jj))+1
        plot([-T_0(1,jj)+T(1,jj)*ii T_r(1,jj)-T_0(1,jj)+T(1,jj)*ii],[s s],'r-','linewidth',1.5);
    end
end

end