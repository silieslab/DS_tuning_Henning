%LOGMVGAMMA    Log of the multivariate gamma function.
%   LOGMVGAMMA(p,a) computes log( Gamma_p( x ) ).
%   See, e.g., https://en.wikipedia.org/wiki/Multivariate_gamma_function 
%
%  The input arguments are:
%   p    - [1 x 1] 
%   x    - [1 x 1] 
%
%  Returns:
%   f    - [1 x 1] log( Gamma_p( x ) )
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function f = logmvgamma(p, x)

j = 1:p;
f = p*(p-1)/4*log(pi) + sum( gammaln(x + (1-j)./2) );

end