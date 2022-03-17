function [ x, u, time, iterations, iter2 ] = ADMM(coeffs,P_dem_e,SOC_0,P_OC_min,P_OC_max,SOC_min,SOC_max,R_0,V_OC,Q_bat,misc,idx_pos,idx_neg,rho1)

N = length(P_dem_e);

% rho1 = 1E-9;
rho2 = 3E-1;

%% Feasibility check

%% Initialization (these operations can be performed offline and are not timed)

I = eye(N);
Psi = tril(ones(N,N));
M1 = (rho1 * I + rho2/(Q_bat*V_OC)^2 * (Psi' * Psi));

u = P_OC_max;
% u(idx_neg) = V_OC.^2./(2*R_0) .* (1-sqrt(1-4*R_0.*P_dem_e(idx_neg)./V_OC.^2));
zeta = -u;
x = min(SOC_max,max(SOC_min,SOC_0 + cumsum(zeta)/(Q_bat*V_OC)));
lambda1 = zeros(N,1);
lambda2 = SOC_0 + cumsum(zeta)/(Q_bat*V_OC) - x;

%% Algorithm
tic

iterations = 0;
flag = 1;

while flag
    
    % The u update for k in P is the solution of a convex optimization
    % problem. This is solved using a newton method that is included in the
    % function f_BacktrackingNewtonVector
    [u,iter2] = f_BacktrackingNewtonVector(coeffs,P_OC_min,P_OC_max,R_0,V_OC,P_dem_e, rho1, zeta, lambda1);
%     [u(idx_pos),iter2] = f_BacktrackingNewtonVector(coeffs,P_OC_min(idx_pos),P_OC_max(idx_pos),R_0,V_OC,P_dem_e(idx_pos), rho1, zeta(idx_pos), lambda1(idx_pos));
    
    %hold the current value of zeta for residual calculations
    zeta_old = zeta;
    
    % The zeta update is solved
    vec = - rho1 .* (u + lambda1) - rho2 * cumsum(SOC_0 - x + lambda2)/(Q_bat*V_OC);
    zeta = (vec \ M1)';

    % The x update is trivial
    x_old = x;
    x = SOC_0 + cumsum(zeta)/(Q_bat*V_OC) + lambda2;
    x = min(SOC_max,max(SOC_min,x));
    
    lambda1 = lambda1 + (u + zeta);
    lambda2 = lambda2 + (SOC_0 + cumsum(zeta)/(Q_bat*V_OC) - x);
    
    % The residual and Lagrange multiplier updates are also trivial
%     r = [u + zeta; SOC_0 + cumsum(zeta)/(Q_bat*V_OC) - x];
%     s = [rho1 * (zeta_old - zeta); rho2*cumsum(zeta_old - zeta)/(Q_bat*V_OC) - rho2*(x_old-x)];
    r = [u + zeta; cumsum(zeta) - (x-SOC_0)*(Q_bat*V_OC)];
    s = [rho1 .* (zeta_old - zeta); rho2*cumsum(zeta_old - zeta) - rho2*(x_old-x)*(Q_bat*V_OC)];
    
    iterations = iterations + 1;
    
    % Termination criteria
    if iterations > misc.maxIterations || max(norm(r), norm(s)) < misc.epsilon
        flag = 0;
    end
    
end

time = toc;

x = SOC_0 - cumsum(u)/(Q_bat*V_OC);

return
