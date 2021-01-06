%% Solve a cubic analytically
function [x1,x2,x3] = cubicroots(a, b, c)

% Determine if it is the three real roots case or one real, two complex
Q  = (a^2 - 3*b)/9;
R  = (2*a^3 - 9*a*b + 27*c)/54; 

Q3 = Q^3;
R2 = R^2;

% If three real roots, use trigonometric solutions
if (R2 < Q3)
    theta = acos(R/sqrt(Q3));
    SQ    = sqrt(Q);
    x1  = -2*SQ*cos(theta/3) - a/3;
    x2  = -2*SQ*cos((theta + 2*pi)/3) - a/3;
    x3  = -2*SQ*cos((theta - 2*pi)/3) - a/3;
    
    return
end

% Otherwise we have one real, two complex
A = -(abs(R) + sqrt(R2-Q3))^(1/3);
if (R < 0)
    A = -A; 
end
if (A==0)
    B = 0;
else
    B = Q/A;
end

AB   = A + B;
x1   = AB - a/3;
comp = 1i*sqrt(3)*(A-B)/2;
x2   = -0.5*AB - a/3 + comp;
x3   = -0.5*AB - a/3 - comp;

end