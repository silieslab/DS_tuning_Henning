%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function mm = mm_EstimateTheta(mm, data, wClasses)
   
%% Estimate parameters
% For each class
for k = wClasses

    r = mm.r(:,k);  % posterior probabiity of item belonging to the class    
    
    % For each model (across columns), estimate parameters
    for i = 1:mm.nModelTypes
        
        model = mm.class{k}.model{i};           % model                
        y = data(:, mm.class{k}.model{i}.Ivar); % data 
        ix = ~isnan(y);                         % ignore missing data
        
        switch model.type
            
            %% Univariate Weibull
            case 'weibull'

            theta = fminunc(@(X) weibull_msglen(y(ix), r(ix), X(1), X(2)), [log(mean(y(ix))), 0], mm.opts.SearchOptions);
            model.theta = exp(theta(:));
            
            %% Univariate exponential
            case 'exp'
                
            % Sufficient statistics
            S  = sum(r(ix) .* y(ix));           % sum y_i

            % Estimate parameters
            Nk  = sum(r(ix));   
            Coeff = [-(Nk+1), S, 1-Nk, S];
            Coeff = Coeff ./ Coeff(1);            
            
            lambda = cubicroots(Coeff(2), Coeff(3), Coeff(4));
            
            model.theta = lambda;
            
            %% Negative binomial distribution
            case 'negb'
            
            theta = fminunc(@(X) negb_msglen(y(ix), r(ix), X(1), X(2)), [0, 0], mm.opts.SearchOptions);
            model.theta = exp(theta(:));            
            
            %% Laplace distribution
            case 'Laplace'
            
            Nk = sum(r(ix));
            [mu,fval] = fminunc(@(MU) sum(r(ix) .* abs(y(ix) - MU)), 0, mm.opts.SearchOptions);
            %mu = medianw(y(ix),r(ix));
            %fval = sum(r(ix) .* abs(y(ix) - mu));
            b  = max(fval / (Nk - 1) , 1e-3 );
            
            model.theta = [mu; b];
            
            %% Gamma distribution
            case 'gamma'
                
            % Sufficient statistics
            S  = sum(r(ix) .* y(ix));
            L  = sum(r(ix) .* log(y(ix)));
            Nk = sum(r(ix)); 
            
%             % solve for phi
%             s_phi = log(S/Nk) - L/Nk ;
%             logphi_init = log(3 - s_phi + sqrt((s_phi - 3)^2 + 24*s_phi)) - log(12) - log(s_phi);
%             phi = exp( fminunc(@(X) gamma_msglen(Nk,S,L,X), logphi_init, mm.opts.SearchOptions) );
%             
%             % solve for mu
%             Coeff = [(1 + Nk*phi), -S*phi, (-1 + Nk*phi), -S*phi];
%             Coeff = Coeff ./ Coeff(1);                       
%             [x1,x2,x3] = cubicroots(Coeff(2), Coeff(3), Coeff(4));            
%             v = [x1,x2,x3];
%             mu=v(v>0);
%                 
%             model.theta = [mu; phi];

            [phi, mu, ~] = mm_EstimateGamma(S, L, Nk);
            model.theta = [mu; phi];
            
            %% Univariate k-nomial model
            case 'multi'
                
            % Hyperparameters
            alpha = mm.ModelTypes{i}.alpha;
            M = mm.ModelTypes{i}.nStates;
            A = mm.ModelTypes{i}.A;
            
            % Sufficient statistics
            Cnts = zeros(M, 1);
            R = r(ix);
            for j = 1:M
                Cnts(j) = sum(R(y(ix) == j));
            end
            Nk = sum(Cnts);

            % Estimate parameters
            model.theta = (Cnts + alpha - 1/2) ./ (Nk + A - M/2);
            
            %% Univariate Gaussian model
            case 'Gaussian'
        
            % Sufficient statistics
            s  = sum(r(ix) .* y(ix));           % sum y_i
            s2 = sum(r(ix) .* (y(ix).^2));      % sum y_i^2

            % Estimate parameters
            Nk  = sum(r(ix));                   % n
            mu  = s/Nk;
            tau = max( (s2 - 2*mu*s + Nk*mu^2) / (Nk-1), 1e-3);

            model.theta = [mu; tau];
            
            %% Single factor analysis
            case 'sfa'

            % Parameters
            ix  = ~any(isnan(y),2);        
            Nk  = sum(r(ix));                   % n            
            mu  = sum(bsxfun(@times, y(ix,:), r(ix))) / Nk;    
            xmu = bsxfun(@minus, y(ix,:), mu);
            rxmu = bsxfun(@times, xmu, sqrt(r(ix)));
            [a_sfa, sigma, v, collapsed] = mmlsfa(rxmu, Nk);
            
            model.theta = [mu(:); sigma(:); a_sfa(:)];   
            model.v     = v;
            model.collapsed = collapsed;

            
            %% Multivariate Gaussian model
            case 'mvg'

            % Estimate parameters
            d = mm.ModelTypes{i}.nDim;
            ix  = ~any(isnan(y),2);                                 
            Nk  = sum(r(ix));                   % n
            mu  = sum(bsxfun(@times, y(ix,:), r(ix))) / Nk;
            xmu = bsxfun(@minus, y(ix,:), mu);
            Sigma = (bsxfun(@times, xmu, r(ix))'*xmu + eye(d))/ (Nk+d);
            
            model.theta = [mu(:); Sigma(:)];                
                            
            %% Poisson model
            case 'Poisson'
                
            % Sufficient statistics
            s  = sum(r(ix) .* y(ix));           % sum y_i
            
            % Estimate parameters
            Nk  = sum(r(ix));                   % n            
            lambda = (sqrt(s^2 + 2*s*(Nk-1) + (Nk+1)^2) + s - Nk - 1) / 2 / Nk;
            
            model.theta = lambda;
            
            %% Geometric distribution
            case 'geometric'
                
            s  = sum(r(ix) .* y(ix));           % sum y_i
            Nk  = sum(r(ix));                   % n            
            r = roots( [Nk+s, 1-4*Nk-3*s, 1+6*Nk+3*s, -4-5*Nk-2*s, 2+2*Nk] );
            
            model.theta = r(imag(r) == 0 & r > 0 & r < 1);
            
            %% Inverse Gaussian
            case 'invGaussian'
                
            % Sufficient statstics
            S1 = sum(r(ix) .* y(ix));
            S2 = sum((r(ix) ./ y(ix)));
            
            % Estimate parameters
            Nk = sum(r(ix));
            mu = S1 / Nk;
            lambda = max(1e-5, (S1*S2 - Nk^2) / (Nk - 1) / S1);
            
            model.theta = [mu; lambda];
            
            %% Linear regression
            case 'linreg'
                
            CovIx = mm.ModelTypes{i}.CovIx; 
            P   = length(CovIx);
            ix  = ix & ~any(isnan(data(:,CovIx)),2);         
            X   = data(ix, CovIx);
            y   = data(ix, mm.class{k}.model{i}.Ivar);

            % Estimate regression coefficients
            yr = y .* sqrt(r(ix));
            Xr = bsxfun(@times, X, sqrt(r(ix)));
            b = wridge(X, y, r(ix), 1);
            %b = lscov([ones(length(y),1), X], y, r(ix));

            % Estimate tau
            Nk  = sum(r(ix));   
            mu_z = Xr * b(2:end);
            Kconst = mu_z'*mu_z;
            e2 = sum( (yr - mu_z - sqrt(r(ix))*b(1)).^2 );
%            e2 = sum( (yr - mu_z - b(1)).^2 );
            log_kappa_P = -P*log(2) + log(P) + (1-P)*log(pi) + 2*psi(1)-P;
            logKmult = ( log_kappa_P + P * log(pi * Kconst) - 2*gammaln(P/2 + 1) );    
            
            logtau_init = log( max(1e-5, e2 / (Nk - P - 1)) );
            % uncomment below if actual MML estimate is required
            %logtau = fminunc(@(LOGTAU) linreg_msglentau(Nk, P, e2, logKmult, LOGTAU), logtau_init, opts.SearchOptions);
            logtau = logtau_init;
            tau = exp(max([logtau, -10]));
            
            model.theta = [tau; b];
            
            %% Logistic regression
            case 'logreg'
                
            CovIx = mm.ModelTypes{i}.CovIx; 
            P   = length(CovIx);
            ix  = ix & ~any(isnan(data(:,CovIx)),2);         
            X   = data(ix, CovIx);
            y   = data(ix, mm.class{k}.model{i}.Ivar);   
            
            % Estimate regression coefficients
            Xr = bsxfun(@times, X, sqrt(r(ix)));            
            yr = y .* sqrt(r(ix));
            
            b = fminunc(@(Z) logreg_msglen(r(ix), Xr, yr, X, y, Z(1), Z(2:end)), zeros(P+1,1), mm.opts.SearchOptions);            
            model.theta = b;
            
            
        end
        
        mm.class{k}.model{i} = model;
        
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = linreg_msglentau(n, p, e2, logKmult, logtau)

tau = exp(logtau);
f = n/2*logtau + e2/2/tau - logtau/2 + log1p(exp(logKmult - p*logtau))/2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = weibull_msglen(x, r, log_lambda, log_k)

k = exp(log_k);
lambda = exp(log_lambda);

F = -log_lambda;
h = log(1+lambda^2) + log(1+k^2);
L = -log_k + log_lambda + (k-1)*log_lambda - (k-1)*log(x) + (x./lambda).^k;

f = F + h + r'*L;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = gamma_msglen(n, S, P, logphi)

phi = exp(logphi);
Coeff = [(1 + n*phi), -S*phi, (-1 + n*phi), -S*phi];
Coeff = Coeff ./ Coeff(1);                       
[x1,x2,x3] = cubicroots(Coeff(2), Coeff(3), Coeff(4));
v = [x1,x2,x3];
mu=v(v>0);

f = n*gammaln(phi) + n*phi*log(mu/phi) - (phi-1)*P + phi/mu*S;  % neg ll
f = f + (1/2)*log( phi*psi(1,phi) - 1 ) + log(1 + phi) + log(phi)/2;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, sigma, v, all_zeros] = mmlsfa(w, N)

%% setup
K = size(w,2);
V = w'*w;

% Initial guess for beta
[v,d] = eig(V/N);
beta = v(:,end);
b2 = d(end);
beta = beta * sqrt(b2);

%% Estimate beta
done = false;
all_zeros = false;
while(~done)
    b2 = beta'*beta;    
    
    if((N-1)*b2 < K)
       beta = zeros(K,1);
       sigma = std(w)';
       v = zeros(size(w,1),1);       
       all_zeros = true;
       done = true;
    else
       beta_old = beta;
       
       sigma = sqrt( sum(w.^2, 1)' / (N-1) ./ (1+beta.^2) );
       Y = V ./ (sigma * sigma');
       beta = Y*beta * (1 - K/(N-1)/b2) / (N-1) / (1+b2);      
       done = norm(beta - beta_old) < 1e-3;
   end
end

%% Other parameters
b2 = beta'*beta;
a = sigma .* beta;
if(~all_zeros)
    y = bsxfun(@rdivide, w, sigma');    
    v = (y * beta) * (1 - K/(N-1)/b2) / (1 + b2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = negb_msglen(x, r, logmu, logphi)

n = sum(r);
mu = exp(logmu);
phi = exp(logphi);

L = -gammaln(x+phi) + gammaln(phi) + gammaln(x+1) - x .* log(mu/(mu+phi)) - phi*log(phi/(mu+phi));
h = -2*log(2) + 2*log(pi) + log(1+phi^2) + log(1+mu^2);

numerator = mu + phi*(mu + phi) * psi(1,phi)*((phi/(mu+phi))^phi - 1);
J = log(n) - log(mu)/2 - log(mu+phi) + log(-numerator)/2;

f = r'*L + h + J;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msglen = logreg_msglen(r, Xr, yr, X, y, b0, b)

% number of parameters
P = length(b);
K = P + 1;

% negative log-likelihood
[L,mu] = logregnll(r,X,y,b0,b);

% assertion
v = mu.*(1-mu);
Xv = bsxfun(@times,Xr,sqrt(v));
J = Xv'*Xv;
yhat = Xr*b;
Kb = yhat'*yhat;
logh = gammaln(P/2+1) + logdet(Xr'*Xr)/2 - (P/2)*log(pi) - (P/2)*log(Kb);
log_kappa_K = -K*log(2) + log(K) + (1-K)*log(pi) + 2*psi(1) - K;
logratio = log_kappa_K + logdet(J) - 2*logh;

assertion = log1p(exp(logratio))/2;

% msglen
msglen = assertion + L + K/2;

end

function [f,mu,yhat] = logregnll(r,X,y,b0,b)

lowerBnd = log(eps); 
upperBnd = -lowerBnd;
muLims = [eps, 1-eps];

%% negative log-likelihood
yhat=constrain(b0 + X*b, lowerBnd, upperBnd);
mu=1./(1 + exp(-yhat));

if any(mu < muLims(1) | muLims(2) < mu)
    mu = max(min(mu,muLims(2)),muLims(1));
end

f = -(y.*log(mu) + (1.0-y).*log(1.0-mu));
f = r'*f;


end