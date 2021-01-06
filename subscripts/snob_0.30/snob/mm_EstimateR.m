function [r, p] = mm_EstimateR(mm, data)

% Normalize to get the posterior probability of assignment to classes
p = exp(-mm_Likelihood(mm, data, 1:mm.nModelTypes));
p(p==0) = realmin; 
r = bsxfun(@rdivide, p, sum(p,2));

end