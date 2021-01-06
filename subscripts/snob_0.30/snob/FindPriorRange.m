%FINDPRIORRANGE    Return a, such that x \in [exp(-a), exp(+a)]
%
%   FINDPRIORRANGE(x) returns the smallest number a such that all numbers in vector x are
%   between exp(-a) and exp(+a).
%
%  The input arguments are:
%   x    - [p x 1] 
%
%  Returns:
%   a    - [1 x 1] 
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function a = FindPriorRange(x)

a = max(abs([floor(log(min(x))), ceil(log(max(x)))]));

end