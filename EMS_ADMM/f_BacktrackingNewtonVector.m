function [u, iter] = f_BacktrackingNewtonVector(coeffs,P_OC_min,P_OC_max,R_0,V_OC,P_dem_e, rho1, zeta, lambda1)

alpha2 = coeffs(1);
alpha1 = coeffs(2);

u = P_OC_max;
iter = 1;
delta_t = 0.5;

flag = 1;

while flag
    g = u - R_0./V_OC^2 .* u.^2;
    g_d1 = 1 - 2.*R_0./V_OC^2 .* u;
    g_d2 = - 2.*R_0./V_OC^2;

    w = P_dem_e - g;
    
    %f = alpha2.*w.^2 + alpha1.*w + alpha0;
    f_d1 = 2.*alpha2.*w + alpha1;
    f_d2 = 2.*alpha2;

    %y = f + rho1/2*(u + zeta + lambda_1).^2;
    y_d1 = -g_d1.*f_d1 + rho1.*(u + zeta + lambda1);
    y_d2 = (g_d1).^2.*f_d2 - g_d2.*f_d1 + rho1;
    
    %Calculate Newton step
    deltax = -y_d1./y_d2;
    
    %Take Newton Step
    u_old = u;
    u = u + delta_t.*deltax;
    u = min(P_OC_max, max(P_OC_min, u));
    
    iter = iter + 1;
    
    if norm(u - u_old) < 1E-5 || iter > 100
        flag = 0;
    end
    
end

end

