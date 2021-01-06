%% Example - Acidity data (real data)
clear;

% Load the data. 
% The data consists of a single numerical continuous variable.
%
% Reference: 
% Richardson, S. and Green, P. J. (1997). 
% On Bayesian analysis of mixtures with unknown number of components (with discussion). Journal of the Royal Statistical Society, Series B, 59, 731--792
load data/acidity;  

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate Gaussian distribution: {'norm',1}
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
mm = snob(acidity, {'norm',1}, 'k',5);

% Print a summary of all the components (parameters and structure) of the
% mixture model we have discovered. Snob discovered two classes; one of the
% classes [N(mu = 4.3, sigma = 0.4)] has a sample size n~92, the other
% class [N(mu = 6.3, sigma = 0.5)] has a smaller sample size n~63. 
%
% The total message length of this model is ~209 nits.
mm_Summary(mm);

% Plot the mixture distribution for the acidity data
mm_PlotModel1d(mm, acidity, 1);

% We now force Snob to use 3 classes and turn off subpopulation discovery. 
mm2 = snob(acidity, {'norm',1}, 'k',3,'fixedstructure',true);

% Print a summary of the new model. 
% The total message length of the 3-class model is ~215 nits; 6 nits longer than
% the two class model. 
mm_Summary(mm2);

% The two class model is preferred as it has a smalled message length.
% Specifically, the two class model is approximately
exp( -(mm.msglen - mm2.msglen) )
% times more likely a posteriori than the three class model.