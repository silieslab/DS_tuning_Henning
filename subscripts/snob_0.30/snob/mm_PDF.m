%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [p,L] = mm_PDF(mm, Y, wModel, wClass)

%% If no model specified, use all 
if (~exist('wModel','var'))
    wModel = 1:mm.nModelTypes;
end

%% If no class specified, use all
if (~exist('wClass','var'))
    wClass = 1:mm.nClasses;
end

%%
n = size(Y,1);
K = mm.nClasses;   % number of mixtures
a = mm.a;          % mixing proportions

L = zeros(n,K);

%% Get likelihoods
% For each class
for k = wClass
    % For each model (across columns), get likelihoods
    for i = wModel
        
        subL = zeros(n, 1);
                    
        m = mm.class{k}.model{i};            
        I = ~any(isnan( Y(:,m.Ivar) ), 2);
        
        switch m.type            
            
            %% Weibull
            case 'weibull'
                lambda = m.theta(1);
                k_wbl = m.theta(2);
                
                subL(I) = -log(k_wbl/lambda) + (k_wbl-1)*log(lambda) - (k_wbl-1)*log(Y(I, m.Ivar)) + (Y(I, m.Ivar)./lambda).^k_wbl;
            
            %% exponential
            case 'exp'
                
                lambda = m.theta;
                subL(I) = log(lambda) +  Y(I, m.Ivar) ./ lambda;
                
            case 'Laplace'
                mu = m.theta(1);
                b_lap = m.theta(2);
                
                subL(I) = log(2*b_lap) + abs(Y(I, m.Ivar) - mu)/b_lap;            
                
            %% gamma
            case 'gamma'
                mu = m.theta(1);
                phi = m.theta(2);
                
                subL(I) = gammaln(phi) + phi*log(mu/phi) - (phi-1)*log(Y(I, m.Ivar)) + phi/mu*Y(I, m.Ivar);
            
            %% Univaraite k-nomial
            case 'multi'
                
                % Parameters
                theta = m.theta;
                subL(I) = -log( theta(Y(I, m.Ivar)) );
                
            %% Single factor analysis model
            case 'sfa'
                
                % Parameters
                d = mm.ModelTypes{i}.nDim;
                theta = m.theta;
                mu = theta(1:d);
                sigma = theta(d+1:2*d);
                a_sfa = theta(2*d+1:end);
                
                x0 = bsxfun(@minus, Y(I, m.Ivar), mu');
                vrep = repmat(m.v, [1,d]);
                w = x0 - bsxfun(@times,vrep,a_sfa');
                e2 = sum(bsxfun(@rdivide, w.^2, sigma.^2'),2);
                subL(I) = (d/2)*log(2*pi) + sum(log(sigma)) + e2/2;
                
                
            %% Multivariate Gaussian model
            case 'mvg'
                
                % Parameters
                d = mm.ModelTypes{i}.nDim;
                theta = m.theta;
                mu = theta(1:d);
                Sigma = reshape(theta(d+1:end),d,d);
                
                X0 = bsxfun(@minus, Y(I, m.Ivar), mu');
                R = cholcov(Sigma,0);              
                logSqrtDetSigma = sum(log(diag(R)));
                xRinv = X0 / R;
                quadform = sum(xRinv.^2, 2);
                
                % Negative log-likelihood
                subL(I) = 0.5*quadform + logSqrtDetSigma + d*log(2*pi)/2;
                
            %% Univariate Gaussian model
            case 'Gaussian'

                % Parameters
                mu  = m.theta(1);
                tau = m.theta(2);

                % Negative log-likelihood
                subL(I) = (1/2)*log(2*pi*tau) + (Y(I, m.Ivar) - mu).^2 / 2 / tau;                
                
            %% Poisson distribution
            case 'Poisson'
                % Parameters
                lambda = m.theta(1);
                
                % Negative log-likelihood
                subL(I) = lambda - Y(I, m.Ivar)*log(lambda) + gammaln(Y(I,m.Ivar)+1);
                
                
            %% Inverse Gaussian
            case 'invGaussian'
                % Parameters
                mu = m.theta(1);
                lambda = m.theta(2);
                
                % Negative log-likelihood
                subL(I) = (1/2)*log(2*pi*lambda) + (1/2)*log(Y(I, m.Ivar).^3) - 1/lambda/mu + Y(I, m.Ivar)/2/mu^2/lambda + (1./Y(I, m.Ivar))/2/lambda;
                
            %% Gaussian linear regression
            case 'linreg'
                
                % Parameters
                tau = m.theta(1);
                b0  = m.theta(2);
                b   = m.theta(3:end);
                
                X   = Y(:,mm.ModelTypes{i}.CovIx);
                I   = I & ~any(isnan(X),2);         % missing covs
                mu  = b0 + X(I,:)*b;
                
                % Negative log-likelihood
                subL(I) = (1/2)*log(2*pi*tau) + (Y(I, m.Ivar) - mu).^2 / 2 / tau;                
                
            otherwise
                error('Model not available');
        end
        
        L(:,k) = L(:,k) + subL;
    end
    
    % Finally, weight by mixing proportions
    L(:,k) = L(:,k) - log(a(k));
end


%% Finally, combine them to get our probability density function
p = sum(exp(-L(:,wClass)),2);
L = -log(p);

end