function mm = mm_MsgLen(mm, data)

% Model and data structure
% ------------------------
n = size(data, 1);                  % Data set size
K = mm.nClasses;                    % Number of mixtures
a = mm.a;                           % Mixing proportions

% Assertion length for K (uniform prior)
% --------------------------------------
Ak = K*log(2);    % p(K) = 2^(-K)

% Assertion length for the mixing proportions
% -------------------------------------------
Aa = (K-1)*log(n)/2 - sum(log(a))/2 - gammaln(K);

% Detail length and assertion length for kn
% -----------------------------------------
p = exp(-mm_Likelihood(mm, data, 1:mm.nModelTypes));
p(p==0) = realmin; 
r = bsxfun(@rdivide, p, sum(p,2));
Nk = sum(r,1);      % how many things in each class

An_L = -sum(log(sum(p,2)));
mm.L = An_L;

% Assertion length for the hyperparameters
% ----------------------------------------
Atheta = 0;
for i = 1:mm.nModelTypes
    switch mm.ModelTypes{i}.type
        %% Gaussian hyperparameters mu \in [mu0,mu1], tau \in [exp(-a),exp(+a)]
        case 'Gaussian'
            tau = zeros(K,1);
            for k = 1:K    
                tau(k) = mm.class{k}.model{i}.theta(2);
            end
            a_tau = FindPriorRange(tau);
            % mu hyperparameters; each coded as log(n)/2
            % tau hyperparameter coded as logstar(a)
            Atheta = Atheta + sum(log(Nk)) + K*log(2*a_tau) + logstar(a_tau);   
            
        %% Laplace hyperparameters mu \in [mu0,mu1], tau \in [exp(-a),exp(+a)]
        case 'Laplace'
            tau = zeros(K,1);
            for k = 1:K    
                tau(k) = mm.class{k}.model{i}.theta(2);
            end
            a_tau = FindPriorRange(tau);
            % mu hyperparameters; each coded as log(n)/2
            % tau hyperparameter coded as logstar(a)
            Atheta = Atheta + sum(log(Nk)) + K*log(2*a_tau) + logstar(a_tau);             
            
        case 'mvg'
            
            d = mm.ModelTypes{i}.nDim;
            Atheta = Atheta + sum(log(Nk)) * d;
            
        case 'sfa'
                        
            d = mm.ModelTypes{i}.nDim;
            sigma = zeros(K,d);
            for k =1:K
                sigma(k,:) = mm.class{k}.model{i}.theta(d+1:2*d);
            end
            a_sigma = FindPriorRange(sigma(:));            
            
            Atheta = Atheta + sum(log(Nk))*d + K*d*log(a_sigma) + logstar(a_sigma);   
            
            
        %% Inverse Gaussian hyperparameters lambda \in [exp(-a), exp(+a)]
        case 'invGaussian'
            lambda = zeros(K,1);
            for k =1:K
                lambda(k) = mm.class{k}.model{i}.theta(2);
            end
            a_lambda = FindPriorRange(lambda);
            % mu hyperparameter; coded as log(n)/2
            % lambda hyperparameter coded as logstar(a)
            Atheta = Atheta + sum(log(Nk))/2 + K*log(2*a_lambda) + logstar(a_lambda);   
            
        %% Gaussian linear regression
        case 'linreg'     
            tau = zeros(K,1);
            for k = 1:K    
                tau(k) = mm.class{k}.model{i}.theta(1);
            end
            a_tau = FindPriorRange(tau);
            % two mu hyperparameters; each coded as log(n)/2
            % K region parameters coded as 1/2*log(n) each
            % tau hyperparameter coded as logstar(a)
            Atheta = Atheta + 3/2*sum(log(Nk)) + K*log(2*a_tau) + logstar(a_tau);               
            
        %% Logistic regression
        case 'logreg'
            
            Atheta = Atheta + 3/2*sum(log(Nk));
    end
end

% Assertion length for theta
% --------------------------
D = (K-1);      % number of mixing proportions
totalParams = (K-1);
nParams = 0;
for k = 1:K
    for i = 1:mm.nModelTypes

        model = mm.class{k}.model{i};           % model                
        switch model.type

            %% Weibull distribution
            case 'weibull'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                lambda = model.theta(1);
                k_wbl = model.theta(2);
                h_theta = -log(2) + log(pi) + log(1+lambda^2) -log(2) + log(pi) + log(1+k_wbl^2);
                F_theta = log(pi) - log(6)/2 - log(lambda) + log(Nk(k));
                AssLen = h_theta + F_theta;   
                
            %% Laplace distribution
            case 'Laplace'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                b_lap = model.theta(2);   
                               
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;
                
                h_theta = log(Rmu) + log(b_lap);
                F_theta = log(Nk(k)) - 2*log(b_lap);
                AssLen = h_theta + F_theta;                   
            
            %% Univariate exponential
            case 'exp'
                nParams = 1;
                totalParams = totalParams + nParams;
                
                lambda = model.theta;
                h_theta = -log(2) + log(pi) + log(1+lambda^2);
                F_theta = log(Nk(k))/2 - log(lambda);
                AssLen = h_theta + F_theta;
                
            %% Univariate gamma
            case 'gamma'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                mu = model.theta(1);
                phi = model.theta(2);
                
                h_theta = -log(2) + 2*log(pi) + log(1+mu^2) + log(phi)/2 + log(1 + phi);
                F_theta =  log(Nk(k))- log(mu) + 1/2*log(phi*psi(1,phi) - 1);
                AssLen = h_theta + F_theta;
            
            %% Univariate k-nomial model
            case 'multi'
                
                % Hyperparameters
                alpha = mm.ModelTypes{i}.alpha;
                M = mm.ModelTypes{i}.nStates;
                A = mm.ModelTypes{i}.A;
                
                nParams = M-1;
                totalParams = totalParams + nParams;                
                theta = model.theta;
                
                h_theta = -gammaln(A) + sum(gammaln(alpha)) - sum( (alpha-1) .* log(theta) );
                F_theta = (M-1)/2*log(Nk(k)) - sum(log(theta))/2;
                
                AssLen = h_theta + F_theta;
            
            %% Single factor analysis model
            case 'sfa'
                % Parameters
                d = mm.ModelTypes{i}.nDim;                
                theta = model.theta;
                sigma = theta(d+1:2*d);
                if(~model.collapsed)    % correlated normal
                    nParams = 3*d + Nk(k);                    
                    a_sfa = theta(2*d+1:end);
                    beta = a_sfa ./ sigma;
                    b2 = beta'*beta;
                    v = model.v;
                               
                    Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;

                    Bk = zeros(d,1);
                    Bk(1:2)=[pi/2, pi];
                    for T=3:d
                        Bk(T)=2*pi * Bk(T-2)/(d-1);
                    end
                    Bk=Bk(d); 

                    F_theta = d/2*log(2*Nk(k)) + d/2*abs(log(Nk(k)*sum(v.^2)) - sum(v)^2) + (Nk(k)-2)/2*log(1+b2) - 3*sum(log(sigma));
                    h_theta = sum(log(Rmu)) + sum(log(sigma)) + (d+1)/2*log(1+b2) + log(Bk) + (Nk(k)/2)*log(2*pi) + v'*v/2;
                
                else    % uncorrelated model
                    nParams = 2*d;                    
                    
                    Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;

                    h_theta = sum(log(Rmu)) + sum(log(sigma));
                    F_theta = d*log(Nk(k)) - 2*sum(log(sigma)) + d*log(2)/2;                                       
                end
                
                totalParams = totalParams + nParams;                
                AssLen = h_theta + F_theta;                
                
            %% Multivariate Gaussian model
            case 'mvg'
                
                d = mm.ModelTypes{i}.nDim;
                nParams = d + d*(d+1)/2;
                totalParams = totalParams + nParams;
                
                theta = model.theta;
                Sigma = reshape(theta(d+1:end),d,d);                
                
                R = cholcov(Sigma,0);              
                logDetSigma = 2*sum(log(diag(R)));
                Rinv = R\eye(d);
                
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;
                
                h_theta = sum(log(Rmu)) + (d+1)*logDetSigma + trace(Rinv*Rinv')/2 + log(2)*d*(d+1)/2 + logmvgamma(d,(d+1)/2);
                F_theta = d*(d+3)/4*log(Nk(k)) - d/2*log(2) -(d+2)/2*logDetSigma;
                AssLen = h_theta + F_theta;                
                
                
            %% Univariate Gaussian model
            case 'Gaussian'
                nParams = 2;
                totalParams = totalParams + nParams;
                        
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;
                
                tau = model.theta(2);
                h_theta = log(Rmu) + log(tau);
                F_theta = log(Nk(k)) - 3*log(tau)/2 - log(2)/2;
                AssLen = h_theta + F_theta;
                
            %% Poisson model
            case 'Poisson'
                nParams = 1;
                totalParams = totalParams + nParams;
                
                lambda = model.theta;
                h_theta = log(pi) + log(lambda)/2 + log(1+lambda);
                F_theta = log(Nk(k))/2 - log(lambda)/2;
                AssLen = h_theta + F_theta;
                
            %% Negative binomial
            case 'negb'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                mu = model.theta(1);
                phi = model.theta(2);
                
                h_theta = -2*log(2) + 2*log(pi) + log(1+phi^2) + log(1+mu^2);
                numerator = mu + phi*(mu + phi) * psi(1,phi)*((phi/(mu+phi))^phi - 1);
                F_theta = log(Nk(k)) - log(mu)/2 - log(mu+phi) + log(-numerator)/2;                
                AssLen = h_theta + F_theta;                
                
            %% geometric model
            case 'geometric'
                nParams = 1;
                totalParams = totalParams + nParams;
                
                theta = model.theta;
                h_theta = -log(2-theta) + log(pi) + log(1-theta)/2 + log(theta^2 - theta + 1);
                F_theta = log(Nk(k))/2 - log(theta) - log(1-theta)/2;
                AssLen = h_theta + F_theta;
                
            %% Inverse Gaussian
            case 'invGaussian'
                nParams = 2;
                totalParams = totalParams + nParams;
                
                % Parameters
                mu = model.theta(1);
                lambda = model.theta(2);                         
                
                % Hyperparameters
                mu0 = mm.ModelTypes{i}.mu0;        
                
                h_theta = -log(mu0)/2 + log(2) + 3*log(mu)/2 + log(lambda);
                F_theta = log(Nk(k)) - log(2)/2 - 3*log(mu)/2 - 3*log(lambda)/2;                
                AssLen = h_theta + F_theta;
                
                
            %% Gaussian linear regression
            case 'linreg'
                nParams = 2; % beta0 and tau; the others are handled within                
                CovIx = mm.ModelTypes{i}.CovIx; 
                Ivar  = mm.class{k}.model{i}.Ivar;
                totalParams = totalParams + nParams + length(CovIx);                
                
                
                P   = length(CovIx);
                ix  = ~isnan(data(:, Ivar)) & ~any(isnan(data(:,CovIx)),2);         
                X   = data(ix,mm.ModelTypes{i}.CovIx);                
                
                beta = model.theta(3:end);               
                logtau = log(model.theta(1));
                
                Xr = bsxfun(@times, X, sqrt(r(ix,k)));
                Kconst = beta'*(Xr'*Xr)*beta;
    
                log_kappa_P = -P*log(2) + log(P) + (1-P)*log(pi) + 2*psi(1)-P;
                logKmult = ( log_kappa_P + P * log(pi * Kconst) - 2*gammaln(P/2 + 1) );
                
                Rmu = mm.ModelTypes{i}.mu1 - mm.ModelTypes{i}.mu0;               
                AssLen = log1p(exp(logKmult - P*logtau))/2 + log(Rmu) + logtau + log(Nk(k)) - log(2)/2 - 3*logtau/2;
                
            %% Logistic regression
            case 'logreg'
                nParams = 0; % all parameters (b0,b) coded within
                CovIx = mm.ModelTypes{i}.CovIx; 
                Ivar  = mm.class{k}.model{i}.Ivar;
                totalParams = totalParams + nParams + length(CovIx);        
                
                P   = length(CovIx);
                ix  = ~isnan(data(:, Ivar)) & ~any(isnan(data(:,CovIx)),2);         
                X   = data(ix,mm.ModelTypes{i}.CovIx);                
                Xr = bsxfun(@times, X, sqrt(r(ix,k)));
                
                beta0 = model.theta(1);
                beta  = model.theta(2:end);      
                
                % get mu
                lowerBnd = log(eps); 
                upperBnd = -lowerBnd;
                muLims = [eps, 1-eps];

                %% negative log-likelihood
                yhat = constrain(beta0 + X*beta, lowerBnd, upperBnd);
                mu = 1./(1 + exp(-yhat));

                if any(mu < muLims(1) | muLims(2) < mu)
                    mu = max(min(mu,muLims(2)),muLims(1));
                end                
                
                % assertion
                v = mu.*(1-mu);
                Xv = bsxfun(@times,Xr,sqrt(v));
                J = Xv'*Xv;
                yhat = Xr*beta;
                Kb = yhat'*yhat;

                logh = gammaln(P/2+1) + logdet(X'*X)/2 - (P/2)*log(pi) - (P/2)*log(Kb);
                log_kappa_P = -(P+1)*log(2) + log(P+1) + (1-(P+1))*log(pi) + 2*psi(1)-(P+1);
                logratio = log_kappa_P + logdet(J) - 2*logh;
                
                AssLen = log1p(exp(logratio))/2;                    

        end

        % Total assertion length and number of parameters
        Atheta = Atheta + AssLen;
        D = D + nParams;
    end
end

% Constant terms
% --------------
gamma_d = -psi(1);
if D > 0
    c_d = -D/2*log(2*pi) + log(D*pi)/2 - gamma_d; 
else
    c_d = 0;
end
constant = c_d  - gammaln(K+1);

% Total Msglen for a mixture model
% --------------------------------
mm.Ak = Ak;
mm.Aa = Aa;
mm.Atheta = Atheta;
mm.constant = constant;
mm.nParams = totalParams;
mm.msglen = Ak + Aa + Atheta + An_L + constant;

% Compute BIC and AIC
% -------------------
mm.AIC = mm.L + (mm.nParams/2);
mm.BIC = mm.L + (mm.nParams/2)*log(n);

end