function [a, r, Nk, p] = mm_EstimateMixingWeights(mm, data)

K = mm.nClasses;         % number of mixtures
n = size(data, 1);

% Estimate 
[r, p] = mm_EstimateR(mm, data);
Nk = sum(r,1)';

% Estimate the mixing weights using the MML87 multinomial estimate
a = (Nk+1/2) ./ (n+K/2);

end