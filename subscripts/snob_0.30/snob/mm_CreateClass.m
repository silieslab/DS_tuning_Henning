function class = mm_CreateClass(ModelTypes)

%% Process all models for eeach class
for i = 1:length(ModelTypes)
    
    %% Type of model
    switch ModelTypes{i}.type

        case 'sfa'
            d = ModelTypes{i}.nDim;
            mu = randn(d,1);
            sigma = rand(d,1);
            a = randn(d,1);
            
            class.model{i}.type  = 'sfa';                % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = [mu; sigma(:); a(:)]; % [mu,sigma,a,v]       
            class.model{i}.v     = 0;
            class.model{i}.collapsed = true;
        
        case 'mvg'
            d = ModelTypes{i}.nDim;
            mu = randn(d,1);
            Sigma = eye(d);
            
            class.model{i}.type  = 'mvg';                % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = [mu; Sigma(:)];      % [mu,Sigma]       
            
        case 'weibull'
            class.model{i}.type  = 'weibull';            % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(2,1);            % [lambda,k]            

        case 'Laplace'
            class.model{i}.type  = 'Laplace';            % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(2,1);            % [mu,b]               
        
        case 'negb'
            class.model{i}.type  = 'negb';               % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(2,1);            % [mu,phi]    
            
        case 'exp'
            class.model{i}.type  = 'exp';                % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(1);              % [mean]            
        
        %% Create a k-nomial model
        case 'multi'
            class.model{i}.type  = 'multi';              % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            
            M = ModelTypes{i}.nStates;
            theta = rand(M,1);
            class.model{i}.theta = theta ./ sum(theta);  % [p1,p2,...,pM]
        
        %% Create a gamma model
        case 'gamma'
            
            class.model{i}.type  = 'gamma';           % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(2,1);            % [mean, dispersion]
            
        %% Create a Gaussian model
        case 'Gaussian'
            
            class.model{i}.type  = 'Gaussian';           % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(2,1);            % [mean, variance]            
            
        %% Create a Poisson model
        case 'Poisson'
            class.model{i}.type  = 'Poisson';            % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(1);              % [mean]
            
        %% Create a geometric model
        case 'geometric'
            class.model{i}.type  = 'geometric';          % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(1);              % [probability]            
        
        
        %% Create an inverse Gaussian model
        case 'invGaussian'
            class.model{i}.type  = 'invGaussian';        % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?
            class.model{i}.theta = rand(2,1);            % [mu, lambda]        
        
        %% Create a Gaussian linear regression model
        case 'linreg'
            
            class.model{i}.type  = 'linreg';             % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?    
            
            P = length(ModelTypes{i}.CovIx) + 2;
            class.model{i}.theta = rand(P,1);            % [tau, b0, b]
            
        %% Create a logistic regression model
        case 'logreg'
            
            class.model{i}.type  = 'logreg';             % type
            class.model{i}.Ivar  = ModelTypes{i}.Ivar;   % which variable in the data?    
            
            P = length(ModelTypes{i}.CovIx) + 1;
            class.model{i}.theta = zeros(P,1);           % [b0, b]              
            
    end

end

end