%% Example - Thyroid data (url: http://archive.ics.uci.edu/ml/datasets/thyroid+disease)
clear;

% Load the data
% The data consists of 5 continuous variables.
load data/thyroid;

% In this example we know that there exist three sub-populations in the
% sample. The labels of the three populations are stored in the variable
% "label" (Class attribute (1 = normal, 2 = hyper, 3 = hypo)). We will use
% snob to attempt to discover these subpopulations and compare the results to
% the actual labels.

% We will try and fit a few different distributions to this data.
% multivariate Gaussian
mm_mvg = snob(x, {'mvg',1:5},'k',1);
% Single factor analysis
mm_sfa = snob(x, {'sfa',1:5},'k',1);
% Univariate normal; i.e., multivariate Gaussian with a diagonal covariance matrix
mm_norm = snob(x, {'norm',1:5},'k',1);
% Univariate Laplace distribution
mm_lapl = snob(x, {'laplace',1:5},'k',1);

% What is the best fitting model according to MML?
msglen = [mm_mvg.msglen, mm_sfa.msglen, mm_norm.msglen, mm_lapl.msglen];
[val,I] = min(msglen);

% The model with the smallest message length is the single factor analysis model
mm_Summary(mm_sfa);
