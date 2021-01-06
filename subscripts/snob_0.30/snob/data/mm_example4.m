%% Example - Simulated data, single factor analysis (SFA) model
clear;

% Seed the random number generator
rng(1);

% Generate data from a mixture of two single factor analysis models.
% A single factor analysis model is defined as
%               x_nk             = mu_k + v_n a_k + sigma_k*r_nk,
%                                   {v_n, {r_nk, k=1,...,K},n=1,...,N} ~ N(0,1)
%               where a_k are the factor loadings, v_n are the factor
%               scores.
a1     = randn(5,1);
a1     = 2*(a1 / norm(a1));
sigma1 = ones(5,1);
mu1    = 2*ones(5,1);
Sigma1 = a1*a1' + diag(sigma1.^2);
x1     = mvnrnd(mu1, Sigma1, 1e2);

a2     = randn(5,1);
a1     = 2*(a2 / norm(a2));
sigma2 = ones(5,1);
Sigma2 = a2*a2' + diag(sigma2.^2);
mu2    = -3*ones(5,1);
x2     = mvnrnd(mu2, Sigma2, 3e2);

x = [x1; x2];

% Fit a multivariate Gaussian distribution to the data
mm_mvg = snob(x, {'mvg',1:5},'k',1,'display',false);
mm_Summary(mm_mvg);

% Fit a single factor analysis model to the data
mm_sfa = snob(x, {'sfa',1:5},'k',1,'display',false);
mm_Summary(mm_sfa);

% The SFA model is strongly preferred according to the message length. 
% Specifically, the SFA model is
exp( -(mm_sfa.msglen - mm_mvg.msglen) )
% times more likely a posteriori than the one-class gamma model.