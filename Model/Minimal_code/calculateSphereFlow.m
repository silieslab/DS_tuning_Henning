%% function [phi, theta, u, v] = CALCULATESPHEREFLOW(phi, theta, r, vT, vR)
% Takes phi, and theta, as azimuth and elevation coordinates in radians
% and radius to calculate the flow caused by a translation vT and rotation vR
% vT and vR are defines in Cartesian coordinates.
% phi is in the range [-pi, pi]
% theta is in the range [-pi/2, pi/2]
% r is in the range [0, inf]
% Example:
% vT = [0, 0, 0]; % m/s
% vR = [1, 0, 0]; % rad/s
% phi = linspace(-pi, pi, 20);
% theta = linspace(-pi/2, pi/2, 10);
% [phi, theta] = meshgrid(phi, theta);
% r = ones(size(phi));

%%
%
%  Copyright 2020 Giordano Ramos-Traslosheros
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [phi, theta, u, v] = calculateSphereFlow(phi, theta, r, vT, vR)

[X, Y, Z] = sph2cart(phi, theta, r);
phi   = rad2deg(phi);
theta = rad2deg(theta);

p = zeros(numel(X), 3);
pSpherical = zeros(size(p));
for iP = 1: numel(X)
    p(iP, :) = getFlowVector([X(iP) Y(iP) Z(iP)], vT, vR);
    pSpherical(iP, :) = cart2sphvec(p(iP, :)', phi(iP), theta(iP))';
end

phi = phi(:);
theta = theta(:);
u = pSpherical(:, 1);
v = pSpherical(:, 2);

end

function p = getFlowVector(d, vT, vR)
    r = norm(d);
    p = -1 / r * (vT - dot(vT, d) * d) - cross(vR, d);
end
