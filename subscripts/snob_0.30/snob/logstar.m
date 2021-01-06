%LOGSTAR    Logstar code for an integer
%   LOGSTAR(x) computes Rissanen's logstar codelength for integer x.
%
%  The input arguments are:
%   x    - [1 x 1] 
%
%  Returns:
%   f    - [1 x 1] log* codelength of x
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function v = logstar(x)

%%
v = 0;
while(1)
    if (log(x) > 0)
        v = v + log(x);
        x = log(x);
    else
        break;
    end
end

end