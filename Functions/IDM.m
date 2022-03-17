function [s,v,a] = IDM(Scenario)
s_max = Scenario.s_max;
s_list = Scenario.s_list;
T = Scenario.T;
T_r = Scenario.T_r;
T_0 = Scenario.T_0;

v_max = 60/3.6; % m/s
a_max = 1.4; % m/s^2
b_com = 2.0; % m/s^2

s_0 = 0; % spacing, m
Ts = 0; % time spacing, s
v(1) = 0;
s(1) = 0;
ds_des(1) = s_0; 
ds(1) = s_0; % 
v_ref = v_max;
b_max = b_com;
t = 2;

mm = 0;
statu_T = [];

while t < s_max % /m
    id1 = sum(s_list - s(t-1) <= 0) + 1;
    [~,id2] = find(0 < (s_list - s(t-1)) & (s_list - s(t-1)) < 100);
    [~,id3] = find(0 < (s_list - s(t-1)) & (s_list - s(t-1)) <= 30);

    ds(t-1) = s_list(id1) - s(t-1);
    
    if t>2
        [~,id4] = find(30 <= (s_list - s(t-2)) & (s_list - s(t-2)) < 100);
        if ~isempty(id3) && ~isempty(id4)
            mm = mm + 1;
            statu_T(mm) = t-1 + T_0(id1) - floor( (t-1+T_0(id1)) / T(id1) ) * T(id1) > T_r(id1); % green 0-30m临界时刻是绿灯
        end
    end

    statu = t-1 + T_0(id1) - floor( (t-1+T_0(id1)) / T(id1) ) * T(id1) <= T_r(id1); % red
    
    if ~isempty(id2)
        v_ref = 0;
        b_max = 2 * b_com;
        if statu
            a(t-1) = max(-(v(t-1)^2)/(2*ds(t-1)));
        else
            a(t-1) = a_max * (1- (v(t-1)/v_max)^4);
        end
        
        if isempty(id3) && statu
            a(t-1) = a_max * (1- (v(t-1)/v_max)^4 - (ds_des(t-1) / ds(t-1))^2);
        end
        
%         if isempty(id3)
%             if statu
%                 v_ref = 0;
%             else
%                 v_ref = v_max / 3;
%             end
%             a(t-1) = a_max * (1- (v(t-1)/v_max)^4 - (ds_des(t-1) / ds(t-1))^2);
%         else
%             if statu
%                 a(t-1) = max(-(v(t-1)^2)/(2*ds(t-1)));
%             else
%                 a(t-1) = a_max * (1- (v(t-1)/v_max)^4);
%             end
%         end        
        
        if ~isempty(id3) && statu_T(mm)
            a(t-1) = a_max * (1- (v(t-1)/v_max)^4);
        end
        
    else
        v_ref = v_max;
        b_max = b_com; 
        a(t-1) = a_max * (1- (v(t-1)/v_max)^4 - (ds_des(t-1) / ds(t-1))^2);
    end

    
    v(t) = min(max(v(t-1)+a(t-1) , 0),v_max); % >=0
    s(t) = s(t-1) + v(t);
    deta_v(t) = v_ref - v(t); % 
    ds_des(t) = max(v(t)*Ts - v(t)*deta_v(t) / (2*(a_max*b_max)^0.5) , 0) + s_0;
    t = t+1;
    
    if (-1e-10 < (s(t-1) - s_max) && (s(t-1) - s_max) <= 0) || (s(t-1) - s_max) > 0
        t = 1e10;
    end
end

end