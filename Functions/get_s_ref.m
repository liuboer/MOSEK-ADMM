function par = get_s_ref(Scenario)
% get references for vehicle positon.

% Output:
% par.N: driving time
% par.s_upper_ref: indexd by time, length=par.N
% par.s_lower_ref: indexd by time, length=par.N
% par.t_lower_ref: indexd by position, length=par.s_max

 % Use IDM to get s_ref_upper
[s_upper_ref,~,~] = IDM(Scenario); % indexd by time

t_small = interp1(s_upper_ref,1:length(s_upper_ref),1:Scenario.s_max, 'linear', 'extrap'); % indexd by position
t_pass_light_small = interp1(1:Scenario.s_max,t_small,Scenario.s_list, 'linear', 'extrap'); % passing time
t_pass_light_big = zeros(size(t_pass_light_small)); % end time of the current red light phase
for i = 1:length(t_pass_light_small)
    t_pass_light_big(i) = ceil((Scenario.T_0(i)+t_pass_light_small(i))/Scenario.T(i)).*Scenario.T(i)-Scenario.T_0(i) - 1;  
end

par.N = round((length(s_upper_ref)-1+t_pass_light_big(end))/2); % determain final time
t_pass_light_big(end) = par.N;

for idx_i = length(t_pass_light_big):-1:2
    if (t_pass_light_big(idx_i) - t_pass_light_big(idx_i-1)) <= 0
        t_pass_light_big(idx_i-1) = t_pass_light_big(idx_i) - 1;
    end
end

% get s_ref_lower
s_lower_ref = interp1([0,1,5,10,t_pass_light_big],[0,0.5,1,10,Scenario.s_list],1:par.N, 'linear', 'extrap');

% final references
par.s_upper_ref = reshape([s_upper_ref(2:end),s_upper_ref(end)*ones(1,par.N-length(s_upper_ref(2:end)))],[],1);
par.s_lower_ref = reshape(s_lower_ref,[],1);
par.t_lower_ref = reshape(t_small,[],1);
end