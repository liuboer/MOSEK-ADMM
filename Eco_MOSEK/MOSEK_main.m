% QP-based Eco-driving of an FCV at signalized intersections.
close all;clear all;clc

% vehicle parameters
par.spd_max = 60/3.6;
par.spd_min = 0;
par.acc_min = -2;
par.acc_max = 1.4;
par.veh_mass = 1380;
par.veh_air_density = 1.2;
par.veh_CD = 0.335;
par.veh_FA = 2.0;
par.veh_gravity = 9.81;
par.veh_rrc = 0.009;

par.delta_s = 1;

% main results
Eco_v_MOSEK = cell(8,1);
ComputTime = zeros(8,1);

for ii = [1,5,8,6]

tic

Scenario = get_scenarios_by_index(ii,par.delta_s);

% get references
s_ref = get_s_ref(Scenario);
v_mean = Scenario.s_max / s_ref.N;

% Generate matrices
Phi = ones(s_ref.N,1);
Psi = tril(ones(s_ref.N));
I = eye(s_ref.N);

%% solve the QP problem
% Solvers => No.1: cvx; No.2: OSQP; No.3: matlab-quadprog; No.4: opti;
solver = 1;
switch solver
    case 1
        % cvx_solver —— Gurobi, Mosek, SDPT3, SeDuMi ---- faster
        cvx_solver Mosek
        cvx_begin
            cvx_precision low
            variable a(s_ref.N)
            minimize( cost_function(a,par,Phi,Psi,v_mean) )
            subject to
            a <= par.acc_max
            a >= par.acc_min
            Psi*a <= [par.spd_max*ones(s_ref.N-1,1); 1] % final speed < 1m/s
            Psi*a >= par.spd_min
            Psi*(Psi*a) <= s_ref.s_upper_ref
            Psi*(Psi*a) >= s_ref.s_lower_ref
        cvx_end
    case 2
        % OSQP solver ---- fastest
        H = sparse(2*Psi + par.veh_air_density*par.veh_CD.*par.veh_FA*v_mean/par.veh_mass * (Psi'*Psi));
        H = (H'+H)/2;
        f = par.veh_gravity*par.veh_rrc*Psi'*Phi;
        A = sparse([I;Psi;Psi*Psi]);
        l = [Phi*par.acc_min;Phi*par.spd_min;s_ref.s_lower_ref];
        u = [Phi*par.acc_max;Phi*par.spd_max;s_ref.s_upper_ref;];
        prob = osqp;
        prob.setup(H, f, A, l, u, 'alpha', 1);
        res = prob.solve();
        a = res.x;
    case 3
        % QP solver in matlab ---- slower
        H = 2*Psi + par.veh_air_density*par.veh_CD.*par.veh_FA*v_mean/par.veh_mass * (Psi'*Psi);
        H = (H'+H)/2;
        f = par.veh_gravity*par.veh_rrc*Psi'*Phi;
        A = [Psi;-Psi;Psi*Psi;-Psi*Psi];
        b = [Phi*par.spd_max;-Phi*par.spd_min;s_ref.s_upper_ref;-s_ref.s_lower_ref];
        lb = Phi*par.acc_min;
        ub = Phi*par.acc_max;
        Aeq = [];
        beq = [];
        a = quadprog(H,f,A,b,Aeq,beq,lb,ub);
    case 4
        % opti tools   optiSolver('QP') OOQP, CLP, IPOPT, MATLAB ---- faster
        optiSolver OOQP
        H = 2*Psi + par.veh_air_density*par.veh_CD.*par.veh_FA*v_mean/par.veh_mass * (Psi'*Psi);
        H = (H'+H)/2;
        f = par.veh_gravity*par.veh_rrc*Psi'*Phi;
        A = [Psi;-Psi;Psi*Psi;-Psi*Psi];
        b = [Phi*par.spd_max;-Phi*par.spd_min;s_ref.s_upper_ref;-s_ref.s_lower_ref];
        lb = Phi*par.acc_min;
        ub = Phi*par.acc_max;
        Opt = opti('H',H,'f',f,'ineq',A,b,'lb',lb,'ub',ub);
        [a,fval,exitflag,info] = solve(Opt);
end

v = [0;Psi*a];
s = [0;Psi*(Psi*a)];

%% save results
ComputTime(ii) = toc;
Eco_v_MOSEK{ii} = v;

save('ComputTime','ComputTime');
save('Eco_v_MOSEK',"Eco_v_MOSEK");
%% plot
figure
subplot(211) % distance-time
plot(s_ref.s_upper_ref); hold on
plot(s_ref.s_lower_ref); hold on
plot(s); hold on
plot_scenario(Scenario,s_ref.N)
legend('Upper','Lower','Optimal','Location','best')
axis([0 s_ref.N 0 Scenario.s_max]); grid on; hold off
ylabel('Distance (m)'); xlabel('Time (s)')

subplot(212) % speed-time
plot(v)
axis([0 s_ref.N 0 par.spd_max+2]); grid on
xlabel('Time (s)'); ylabel('Speed (m/s)')

saveas(gcf,['Eco_Scenario' num2str(ii) '_MOSEK.png'])
end

function J = cost_function(a,par,Phi,Psi,v_mean)
cost1 = par.veh_mass*a'*(Psi*a);
cost2 = par.veh_mass*par.veh_gravity*par.veh_rrc*Phi'*(Psi*a);
cost3 = par.veh_air_density*par.veh_CD.*par.veh_FA*v_mean/2*(Psi*a)'*(Psi*a);
J = cost1 + cost2 + cost3;
end