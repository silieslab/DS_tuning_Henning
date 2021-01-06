%% Example - 1D Function smoothing
clear;

% Create some data and plot the underlying model
n = 500;
x = linspace(-pi,pi,n)';
y = 1./(3*x.^2 - 2*x + 7) - cos(x.^2) + randn(n,1)*0.2;
plot(x,y,'.');
grid;
hold;

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate Gaussian distribution: {'norm',1}
%     (2) column 2 is a Gaussian linear regression; covariate is column 1: 'linreg',[2,1]    
%     (3) initial number of classes is 5: 'k',5
%     (3) Snob will automatically attempt to discover the optimal number of mixtures (subpopulations)
mm = snob([x, y], {'norm',1,'linreg',[2,1]}, 'k',5);

% Compute predictions and plot them
ypred=zeros(n,mm.nClasses); 
for i = 1:mm.nClasses
    b0 = mm.class{i}.model{2}.theta(2);
    b  = mm.class{i}.model{2}.theta(3:end);
    ypred(:,i) = b0 + b*x; 
end
r = mm_EstimateR(mm, [x,nan(n,1)]);
ypred = sum( ypred .* r, 2);
plot(x,ypred,'r.')
