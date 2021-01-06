function y = mm_ImputeModel(mm, y, iModel, allData)

% where are the missing values?
j = find( sum(isnan(y),2) );
probs = zeros(length(j), length(mm.class{1}.model{iModel}.theta) );   % only used for k-nomials

% set all missing values to 0 to start with
ind = isnan(y);
y(ind) = 0;

%% impute all missing entries for this model
for k = 1:mm.nClasses        

    switch mm.class{k}.model{iModel}.type

        % Exponential distribution
        case 'exp'
            impVal = mm.class{k}.model{iModel}.theta(1);
            y(j) = y(j) + mm.r(j,k) * impVal;                     

        % Gamma distribution
        case 'gamma'
            impVal = mm.class{k}.model{iModel}.theta(1);
            y(j) = y(j) + mm.r(j,k) * impVal;     

        % Gaussian distribution
        case 'Gaussian'
            impVal = mm.class{k}.model{iModel}.theta(1);
            y(j) = y(j) + mm.r(j,k) * impVal;                     

        % Inverse Gaussian distribution
        case 'invGaussian'
            impVal = mm.class{k}.model{iModel}.theta(1);
            y(j) = y(j) + mm.r(j,k) * impVal;                     

        % Laplace distribution
        case 'Laplace'
            impVal = mm.class{k}.model{iModel}.theta(1);
            y(j) = y(j) + mm.r(j,k) * impVal;     

        % Linear regression
        case 'linreg'
            b0  = mm.class{k}.model{iModel}.theta(2);
            b   = mm.class{k}.model{iModel}.theta(3:end);

            x   = allData(j, mm.ModelTypes{iModel}.CovIx);
            impVal  = b0 + x*b;                
            y(j) = y(j) + mm.r(j,k) * impVal;     

        % Multinomial distribution
        case 'multi'
            probs = probs + bsxfun( @times, repmat(mm.class{k}.model{iModel}.theta',length(j),1), mm.r(j,k) );

        % Multivariate Gaussian distribution
        case 'mvg'
            d = mm.ModelTypes{iModel}.nDim;
            mu = mm.class{k}.model{iModel}.theta(1:d)';
            for sample = 1:length(j)
                ix = ind(j(sample),:);
                y(j(sample),ix) = y(j(sample),ix) + mm.r(j(sample),k) * mu(ix);
            end

        % Poisson distribution
        case 'Poisson'
            impVal = mm.class{k}.model{iModel}.theta(1);
            y(j) = y(j) + mm.r(j,k) * impVal;     

        % Single factor analysis
        case 'sfa'
            d = mm.ModelTypes{iModel}.nDim;
            mu = mm.class{k}.model{iModel}.theta(1:d)';
            for sample = 1:length(j)
                ix = ind(j(sample),:);
                y(j(sample),ix) = y(j(sample),ix) + mm.r(j(sample),k) * mu(ix);
            end
            
        % Weibull distribution
        case 'weibull'
            lwbl = mm.class{k}.model{iModel}.theta(1);
            kwbl = mm.class{k}.model{iModel}.theta(2);
            impVal = lwbl * gamma(1 + 1/kwbl);
            y(j) = y(j) + mm.r(j,k) * impVal;     

        otherwise
            error('Model not available');
    end        

end    

%% Special imputation procedure for k-nomials
if( strcmp(mm.class{k}.model{iModel}.type,'multi') )
    [~,y(j)] = max(probs,[],2);
end


end