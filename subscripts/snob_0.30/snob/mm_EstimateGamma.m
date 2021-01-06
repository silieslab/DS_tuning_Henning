function [phi, mu, M] = mm_EstimateGamma(S, L, n)

%% Initial guess at 'phi'
mu = S/n;
if (n > 4)
    s_phi = log(S/n) - L/n;
    phi = exp(log(3 - s_phi + sqrt((s_phi - 3)^2 + 24*s_phi)) - log(12) - log(s_phi));
else
    phi = 0.5 / (log(S/n) - L/n);
end

%% Iterative generalized Newton loop
c = L + n;
err = 1e-8; M = 0;

while (1)
    phi_old = phi;
    M = M+1;
    
    % update mu
    C = [1 + n*phi, -phi*S, -1 + n*phi, -phi*S];    
    [x1,x2,x3] = cubicroots(C(2)/C(1), C(3)/C(1), C(4)/C(1));
    v = [x1,x2,x3];
    mu=v(v>0);

    % update phi
    a = -phi*n - phi^2/(1 + phi)^2 + n*phi^2*psi(1,phi);
    b = c - 1/(1 + phi) - S/mu - n*log(mu/phi) - n*psi(0,phi) - a/phi;
    
    phi = -a/b; 
    
    % Done?
    if (abs(phi_old / phi - 1) < err)
        break;
    end
end

end