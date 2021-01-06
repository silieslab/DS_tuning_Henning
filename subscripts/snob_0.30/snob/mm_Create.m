function mm = mm_Create(data, ModelTypes, opts)

%% Structure to holde all mixture model information
mm.name = 'Mixture Model';
mm.nClasses = opts.nClasses;
mm.ModelTypes = ModelTypes;
mm.nModelTypes = length(ModelTypes);
mm.N = size(data,1);
mm.opts = opts;

%% Get minimum members required for each class
minmembers = zeros(mm.nModelTypes, 1);
for i = 1:mm.nModelTypes
    minmembers(i) = mm.ModelTypes{i}.MinMembers;
end
mm.MinMembers = max(minmembers);    % minimum items for each class

K = opts.nClasses;

%% Initialise all parameters using specified algorithm
switch opts.Initialisation
    
    %% Random assignment of data
    case 'random'
        
    % Mixing proportions
    mm.a = rand(K,1); 
    mm.a = mm.a ./ sum(mm.a);

    % Create parameters for each model type and class
    for k = 1:K
        mm.class{k} = mm_CreateClass(ModelTypes);
    end
    
    % Randomly assign data equally to all K classes
    mm.r  = rand(size(data,1), K);
    mm.r  = bsxfun(@rdivide, mm.r, sum(mm.r,2));
    mm.Nk = sum(mm.r,1)';
    
    %% Assignment based on the k-means++ algorithm
    case 'kmeans++'
        
        % run kmeans first
        warning('off', 'stats:kmeans:FailedToConverge');    % supress kmeans++ warnings
        warning('off', 'stats:kmeans:MissingDataRemoved');
        [~,~,~,Dist] = kmeans(data, K);
        warning('on', 'stats:kmeans:FailedToConverge');    
        warning('on', 'stats:kmeans:MissingDataRemoved');
        
        R  = bsxfun(@rdivide, Dist, sum(Dist,2));
        ix = sum(isnan(Dist),2) > 0;
        if(any(ix))                
            t       = rand(sum(ix), K);
            R(ix,:) = bsxfun(@rdivide, t, sum(t,2));
        end
        
        mm.r  = R;
        mm.Nk = sum(mm.r,1)';
        mm.a  = mm.Nk / sum(mm.Nk);
        
        for k = 1:K
            mm.class{k} = mm_CreateClass(ModelTypes);
        end        
    
end

%% negative log-likeliood and message length of mixture model
mm.L      = inf;
mm.msglen = inf;        

%% Run estimation functions for all the classes to seed models with parameters (using random assignments)
if(min(mm.Nk) >= mm.MinMembers) % do this only if we have enough items in each class
    mm = mm_EstimateTheta(mm, data, 1:K);
end

end