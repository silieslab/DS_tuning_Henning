%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [i,Ix,ModelTypes,VarsUsed] = mm_ProcessModelTypes(model_list, i, Ix, ModelTypes, data, VarsUsed)

k = length(ModelTypes) + 1;
cols = model_list{i+1};     
if(any(cols<1) || any(cols>size(data,2)))
    error('Model specification not possible')
end
if(length(cols) < 1)
    error('Must specify at least one data column for each model type')
end
if(length(unique(cols)) ~= length(cols))
    error('Duplicate model entries found')
end

%% Check that model type i was specified correctly
switch lower(model_list{i})
    
    %% Univariate exponential distribution
    case {'exp','exponential'}
                        
        %% Basic data  
        for j = 1:length(cols)
            ModelTypes{k}.type = 'exp';
            ModelTypes{k}.Ivar = cols(j);            
            ModelTypes{k}.MinMembers = 2;        
            
            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end               
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));     
            y = data(ix,ModelTypes{k}.Ivar);
            if(min(y) < 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end 
            if(std(y) == 0)
                error(['Data column ', int2str(cols(j)), ': zero variance']);
            end  
            
            k = k + 1;
        end      
        
    %% Univariate gamma distribution        
    case {'gamma'}
        %% Basic data  
        for j = 1:length(cols)
            ModelTypes{k}.type = 'gamma';
            ModelTypes{k}.Ivar = cols(j);            
            ModelTypes{k}.MinMembers = 4;        
            
            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end               
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));     
            y = data(ix,ModelTypes{k}.Ivar);
            if(min(y) < 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end 
            if(std(y) == 0)
                error(['Data column ', int2str(cols(j)), ': zero variance']);
            end  
            
            k = k + 1;
        end            

    %% Univariate inverse Gaussian
    case {'igauss','invgaussian','igaussian','invg'}
        
        %% Basic data
        for j = 1:length(cols)
            ModelTypes{k}.type = 'invGaussian';
            ModelTypes{k}.Ivar = cols(j);    
            ModelTypes{k}.MinMembers = 4;

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end               
                        
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));
            y = data(ix,ModelTypes{k}.Ivar);            
            if(min(y) <= 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end     
            if(std(y) == 0)
                error(['Data column ', int2str(cols(j)), ': zero variance'])
            end
                            
            ModelTypes{k}.mu0 = min(y);
            
            k = k + 1;
        end

    %% Univariate multinomial distribution
    case {'multinomial','multi'}
        
        %% Basic data
        for j = 1:length(cols)                
            ModelTypes{k}.type = 'multi';
            ModelTypes{k}.Ivar = cols(j);         
            ModelTypes{k}.MinMembers = 1;
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));
            y = data(ix,ModelTypes{k}.Ivar);                     
            uy = sort( unique( y ) )';            
            ModelTypes{k}.nStates = length(uy);

            %% Error checking     
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end               
            
            if(ModelTypes{k}.nStates < 2  || ModelTypes{k}.nStates > 15)
                error(['Data column ', int2str(cols(j)), ': number of states for multinomial distribution is out of range']);
            end                   
            if( length(uy) ~= ModelTypes{k}.nStates || any(uy ~= 1:ModelTypes{k}.nStates) )
                error(['Data column ', int2str(cols(j)), ': a K-nomial variable must be coded as {1,2,...,K}']);
            end
                            
            ModelTypes{k}.alpha = ones(ModelTypes{k}.nStates,1);
            ModelTypes{k}.A     = sum(ModelTypes{k}.alpha);

            k = k + 1;   
        end

    %% Multivariate Gaussian distribution
    case {'mvg','mvn'}
        ModelTypes{k}.type = 'mvg';
        ModelTypes{k}.Ivar = cols;
        
        %% Error checking
        CovIx = ModelTypes{k}.Ivar;
        if( any(CovIx<1) || any(CovIx > size(data,2)) )
            error('All covariates must be included in data matrix');
        end
        if(length(CovIx) < 2)
            error('Multivariate Gaussian must use at least two data columns');
        end
        if(any(VarsUsed(CovIx)))
            error('Cannot use multiple models for the same data column');
        end
        
        d = length(CovIx);
        ModelTypes{k}.nDim = d;
        ModelTypes{k}.MinMembers = d + d*(d+1)/2;
        
        ix = ~any(isnan(data(:,CovIx)),2);
        y = data(ix, CovIx);
        ModelTypes{k}.mu0 = min(y);
        ModelTypes{k}.mu1 = max(y);             
        
    %% Negative binomial distribution
    case {'negb','nbin'}

         %% Basic data
        for j = 1:length(cols)        
            ModelTypes{k}.type = 'negb';
            ModelTypes{k}.Ivar = cols(j);
            ModelTypes{k}.MinMembers = 4;

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end     
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));            
            y = data(ix,ModelTypes{k}.Ivar);            
            if(std(y) == 0)
                error(['Data column ', int2str(cols(j)), ': zero variance']);
            end            
            if(min(y) < 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end            
            if(any(mod(y,1) ~= 0))
                error(['Data column ', int2str(cols(j)), ': data must consist of only positive integers and zero']);
            end            
            
            k = k + 1;
        end

        
    %% Univariate normal distribution
    case {'gaussian','normal','norm','gauss'}
                
        %% Basic data
        for j = 1:length(cols)        
            ModelTypes{k}.type = 'Gaussian';
            ModelTypes{k}.Ivar = cols(j);
            ModelTypes{k}.MinMembers = 4;

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end     
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));            
            y = data(ix,ModelTypes{k}.Ivar);            
            if(std(y) == 0)
                error(['Data column ', int2str(cols(j)), ': zero variance']);
            end            
            
            %% Prior hyperparameters
            % Range for uniform distribution
            ModelTypes{k}.mu0 = min(y);
            ModelTypes{k}.mu1 = max(y);
            
            k = k + 1;
        end
        
    %% Univariate Laplace distribution
    case {'laplace'}
                
        %% Basic data
        for j = 1:length(cols)        
            ModelTypes{k}.type = 'Laplace';
            ModelTypes{k}.Ivar = cols(j);
            ModelTypes{k}.MinMembers = 4;

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end     
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));            
            y = data(ix,ModelTypes{k}.Ivar);            
            if(std(y) == 0)
                error(['Data column ', int2str(cols(j)), ': zero variance']);
            end            
            
            %% Prior hyperparameters
            % Range for uniform distribution
            ModelTypes{k}.mu0 = min(y);
            ModelTypes{k}.mu1 = max(y);
            
            k = k + 1;
        end        

    %% Univariate Poisson distribution
    case {'poisson'}
        
        %% Basic data
        for j = 1:length(cols)
            ModelTypes{k}.type = 'Poisson';
            ModelTypes{k}.Ivar = cols(j);
            ModelTypes{k}.MinMembers = 2; 

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end                
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));   
            y = data(ix,ModelTypes{k}.Ivar);            
            if(min(y) < 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end        
            if(any(mod(y,1) ~= 0))
                error(['Data column ', int2str(cols(j)), ': data must consist of only positive integers and zero']);
            end
            
            ModelTypes{k}.index = ix;   
            
            k = k + 1;
        end
        
    %% Univariate geometric distribution
    case {'geometric'}
        
        %% Basic data
        for j = 1:length(cols)
            ModelTypes{k}.type = 'geometric';
            ModelTypes{k}.Ivar = cols(j);
            ModelTypes{k}.MinMembers = 2; 

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined']);
            end                
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));   
            y = data(ix,ModelTypes{k}.Ivar);            
            if(min(y) < 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end        
            if(any(mod(y,1) ~= 0))
                error(['Data column ', int2str(cols(j)), ': data must consist of only positive integers and zero']);
            end
            
            ModelTypes{k}.index = ix;   
            
            k = k + 1;
        end        

    %% Single-factor analysis
    case {'sfa'}
        ModelTypes{k}.type = 'sfa';
        ModelTypes{k}.Ivar = cols;        
        
        %% Error checking
        CovIx = ModelTypes{k}.Ivar;
        if( any(CovIx<1) || any(CovIx > size(data,2)) )
            error('All covariates must be included in data matrix');
        end
        if(length(CovIx) < 2)
            error('Multivariate Gaussian must use at least two data columns');
        end
        if(any(VarsUsed(CovIx)))
            error('Cannot use multiple models for the same data column');
        end
        
        d = length(CovIx);
        ModelTypes{k}.nDim = d;
        ModelTypes{k}.MinMembers = 3*d;
        
        ix = ~any(isnan(data(:,CovIx)),2);        
        y = data(ix,CovIx);
        ModelTypes{k}.mu0 = min(y);
        ModelTypes{k}.mu1 = max(y);        
        
                    
    %% Univariate Weibull distribution
    case {'weibull','wbl'}
        
        %% Basic data
        for j = 1:length(cols)        
            ModelTypes{k}.type = 'weibull';
            ModelTypes{k}.Ivar = cols(j);       
            ModelTypes{k}.MinMembers = 4;        

            %% Error checking
            if(VarsUsed(cols(j)))
                error(['Data column ', int2str(cols(j)), ': multiple models defined'])
            end               
            
            ix = ~isnan(data(:,ModelTypes{k}.Ivar));
            y = data(ix,ModelTypes{k}.Ivar);            
            if(min(y) < 0)
                error(['Data column ', int2str(cols(j)), ': data cannot be negative']);
            end
            
            ModelTypes{k}.index = ix;

            k = k + 1;
        end
             
    %% Gaussian linear regression
    case {'linreg'}
        
        if(length(cols) < 2)
            error('Must have at least one covariate specified for Gaussian linear regression');
        end
        
        %% Basic data
        TargetIx = cols(1);
        CovIx = cols(2:end);
        
        ModelTypes{k}.type = 'linreg';
        ModelTypes{k}.Ivar = TargetIx;
        ModelTypes{k}.CovIx = CovIx;
        ModelTypes{k}.MinMembers = length(CovIx) + 3;        
        
        %% Error checking
        if(VarsUsed(TargetIx))
            error(['Data column ', int2str(TargetIx), ': multiple models defined'])
        end                                

        %% Prior hyperparameters
        % Range for uniform distribution
        ix = ~any(isnan(data(:,cols)),2);
        y = data(ix,ModelTypes{k}.Ivar);
        ModelTypes{k}.mu0 = min(y);
        ModelTypes{k}.mu1 = max(y);
        
        
        %% check that the covariates are valid
        if( any(ModelTypes{k}.Ivar == CovIx) )
            error('Target cannot also be a covariate');
        end
        if( any(CovIx<1) || any(CovIx > size(data,2)) )
            error('All covariates must exist in the data');
        end
        
        %% Ensure only target column is flagged as 'used'
        cols(2:end) = [];     
        
    %% Logistic regression
    case {'logreg'}
        
        if(length(cols) < 2)
            error('Must have at least one covariate specified for Gaussian linear regression');
        end
        
        %% Basic data
        TargetIx = cols(1);
        CovIx = cols(2:end);
        
        ModelTypes{k}.type = 'logreg';
        ModelTypes{k}.Ivar = TargetIx;
        ModelTypes{k}.CovIx = CovIx;
        ModelTypes{k}.MinMembers = length(CovIx) + 3;        
        
        %% Error checking
        if(VarsUsed(TargetIx))
            error(['Data column ', int2str(TargetIx), ': multiple models defined'])
        end                                
                
        ix = ~any(isnan(data(:,cols)),2);
        uy = unique(data(ix,ModelTypes{k}.Ivar));               
        if( length(uy) ~= 2 || any(uy ~= [0:1]') )
            error(['Data column ', int2str(TargetIx), ': logistic regression target variable must be coded as {0,1}']);
        end
        
        
        %% check that the covariates are valid
        if( any(ModelTypes{k}.Ivar == CovIx) )
            error('Target cannot also be a covariate');
        end
        if( any(CovIx<1) || any(CovIx > size(data,2)) )
            error('All covariates must exist in the data');
        end
        
        %% Ensure only target column is flagged as 'used'
        cols(2:end) = [];           
           
        
    %% Otherwise
    otherwise
        error('Model specification unknown.');
        
end

%% update list
VarsUsed(cols) = true;                            
i = i + 2;        

end