%% Example - Simulated data (mixtures of exponential distributions)
clear;

% Seed the random number generator
rng(1);

% Generate some data from a mixture of two exponential distributions 
% x ~ 0.6 Exp(5) + 0.4 Exp(1)
x = [exprnd(5, 60, 1); exprnd(1, 40, 1)];

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate exponential distribution: {'exp',1}
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
mm = snob(x, {'exp',1});

% Print a summary of the discovered mixture model.
% The total message length of this model is 239.64 nits. 
% The model discovered is x ~ 0.43 Exp(0.9) + 0.57 Exp(6.5)
mm_Summary(mm);

% What if we fit three sub-populations to this data?
% Run Snob with the following options: 
%     (1) the data is modelled using a univariate exponential distribution: {'exp',1}
%     (2) initial number of classes is 3: 'k',3
%     (2) Snob will NOT search for the best model structure: 'fixedstructure',true
mm2 = snob(x, {'exp',1},'k',3,'fixedstructure',true);

% The total message length of the three class model is 243.19 nits.
mm_Summary(mm2);

% The two class model is preferred as it has a smalled message length.
% Specifically, the two class model is approximately
exp( -(mm.msglen - mm2.msglen) )
% times more likely a posteriori than the three class model.

% Next, we use a mixture of gamma distributions instead of exponentials.
mm3 = snob(x, {'gamma',1});

% The two class exponential model is preferred as it has a smalled message length.
% Specifically, the two class exponential model is approximately
exp( -(mm.msglen - mm3.msglen) )
% times more likely a posteriori than the one-class gamma model.
