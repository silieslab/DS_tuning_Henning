%MM_PLOTMODEL1D    Creates a 1D plot of a model.
%  mm_PlotModel1d(.) creats a plot of one of the models within the
%  structure mm.
%  
%  The input arguments are:
%   mm      - structure respresenting the complete mixture model
%   x       - [n x p] data set that was used to train the model (n samples; p variables)
%   wModel  - the model to be displayed
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
%% Plot a 1-d model from the given mixture model
function mm_PlotModel1d(mm, data, wModel)

%% error checking
if(wModel < 1 || wModel > mm.nModelTypes)
    error(['The model to plot must be an integer between 1 and ', num2str(mm.nModelTypes), '.']);
end
if(length(wModel) ~= 1)
    error(['The model to plot must be an integer between 1 and ', num2str(mm.nModelTypes), '.']);
end

%% asthetics
n = 1e3; % how many points to use for graph
style = {'r--','g--','k--','m--','y--','c--','r:','g:','k:','m:','y:','c:'}; % line styles

%% Determine min/max plot range as [mean - 3*std, mean + 3*std]
K = mm.nClasses;
mu = zeros(K,1);
tau = zeros(K,1);
switch mm.ModelTypes{wModel}.type    
    case 'mvg'
        if(mm.ModelTypes{wModel}.nDim > 2)
            error('Can only plot 2D Gaussians');
        end
        min_mvg = [inf, inf];
        max_mvg = [-inf, -inf];
        for k = 1:K
            d = mm.ModelTypes{wModel}.nDim;
            theta = mm.class{k}.model{wModel}.theta;
            mu_mvg = theta(1:d);
            Sigma_mvg = reshape(theta(d+1:end),d,d);   
            
            min_mvg = min(min_mvg, mu_mvg' - 3*sqrt(Sigma_mvg(1:3:4)));
            max_mvg = max(max_mvg, mu_mvg' + 3*sqrt(Sigma_mvg(1:3:4)));
        end        
        
    case 'exp'
        for k = 1:K
            lam = mm.class{k}.model{wModel}.theta(1);
            mu(k) = lam;
            tau(k) = lam^2;
        end 
        
    case 'weibull'
        for k = 1:K
            lam = mm.class{k}.model{wModel}.theta(1);
            kwbl = mm.class{k}.model{wModel}.theta(2) ;
            mu(k) = lam*(log(2))^(1/kwbl);
            tau(k) = lam^2*(gamma(1+2/kwbl) - gamma(1+1/kwbl)^2);            
        end    
        
    case 'Gaussian'
        for k = 1:K
            mu(k) = mm.class{k}.model{wModel}.theta(1);
            tau(k) = mm.class{k}.model{wModel}.theta(2);
        end
        
    case 'Laplace'
        for k = 1:K
            mu(k) = mm.class{k}.model{wModel}.theta(1);
            tau(k) = 2 * mm.class{k}.model{wModel}.theta(2)^2;
        end        
                
      case 'invGaussian'
        for k = 1:K
            mu(k) = mm.class{k}.model{wModel}.theta(1);
            lam = mm.class{k}.model{wModel}.theta(2);
            tau(k) = mu(k)^3 / lam;
        end
        
    case 'Poisson'
        for k = 1:K
            lam = mm.class{k}.model{wModel}.theta(1);
            mu(k) = lam;
            tau(k) = lam;
        end
        
    case 'gamma'
        for k = 1:K
            mu(k) = mm.class{k}.model{wModel}.theta(1);
            tau(k) = mu(k)^2 / mm.class{k}.model{wModel}.theta(2);
        end
        
    case 'multi'
        probs = zeros(mm.ModelTypes{wModel}.nStates, K);
        for k = 1:K
            probs(:,k) = mm.class{k}.model{wModel}.theta;
        end        

    otherwise
        error('This type of model cannot be graphed');
end

min_val = mu - 3*sqrt(tau);
max_val = mu + 3*sqrt(tau);
min_val = min(min_val);
max_val = max(max_val);

switch mm.ModelTypes{wModel}.type
    case 'mvg'
        
    case 'weibull'
        min_val = max(min_val, 1e-3);
    case 'exp'
        min_val = max(min_val, 1e-3);        
    case 'invGaussian'
        min_val = max(min_val, 1e-3);  
        max_val = max(data(:,wModel))*1.5;
    case 'Poisson'
        min_val = max(min_val, 1e-3);             
end

%% Plot data 
if(~strcmp(mm.ModelTypes{wModel}.type,'linreg') && ~strcmp(mm.ModelTypes{wModel}.type,'multi') && ~strcmp(mm.ModelTypes{wModel}.type,'mvg'))
    clf; 
    hold on; 
    grid on;
    
    %% histogram
    histogram(data(:,wModel), 'normalization', 'pdf', 'facealpha', 0.05, 'linewidth', 0.4);
    
    %% Plot each class, then overall
    y = ones(n, mm.nModelTypes) * NaN;
    y(:,wModel) = linspace(min_val,max_val,n)';
    for k = 1:K
        plot(y, real(mm_PDF(mm, y, wModel, k)), style{mod(k,length(style)-1)+1},'LineWidth',2);
    end

    % Then, plot complete mixture PDF
    mm_pdf = real(mm_PDF(mm, y));
    plot(y, mm_pdf);

    %% Finally, add the legend
    s = 'legend(''Data'',';
    for k = 1:mm.nClasses
        s = strcat(s,sprintf('%cClass %d%c,',39,k,39));
    end
    s = strcat(s,sprintf('%cPDF%c);',39,39));
    eval(s);

    axis([min_val,max_val,0,max(mm_pdf)*1.1]);
    set(gca, 'FontSize', 16)
    set(gcf,'color','w');
    ylabel('Density','Fontsize',18);
    xlabel('Data','Fontsize',18);
    box;
    
elseif(strcmp(mm.ModelTypes{wModel}.type,'mvg'))
       
    g1 = linspace(min_mvg(1), max_mvg(1), 1e2);
    g2 = linspace(min_mvg(2), max_mvg(2), 1e2);
    [X1,X2] = meshgrid(g1,g2);
    F = zeros(1e4,1);
    for k = 1:K
        d = mm.ModelTypes{wModel}.nDim;
        theta = mm.class{k}.model{wModel}.theta;
        mu = theta(1:d);
        Sigma = reshape(theta(d+1:end),d,d);   
        F = F + mm.a(k)*mvnpdf([X1(:) X2(:)], mu', Sigma);
    end
    F = reshape(F,1e2,1e2);
    surf(g1,g2,F);
    caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
    set(gcf,'color','w'); 
    set(gca, 'FontSize', 16)
    xlabel('x1','Fontsize',18); 
    ylabel('x2','Fontsize',18); 
    zlabel('Probability Density','Fontsize',18);    

elseif(strcmp(mm.ModelTypes{wModel}.type,'multi'))
    if(mm.nClasses > 2)
        bar(probs', 0.7, 'grouped', 'FaceColor', [0.65 0.65 0.65]);    
    else
        bar([probs'; nan(1,mm.ModelTypes{wModel}.nStates)], 0.7, 'grouped', 'FaceColor', [0.65 0.65 0.65]);    
        a = axis; axis([a(1) a(2)-0.8 a(3:4)]);
    end
    grid;
    set(gca, 'FontSize', 16)
    set(gcf,'color','w'); 
    ylabel('Probability','Fontsize',18);
end

end