% fitFlowField_CV
% Function to fit the optic flow fields to the T4T5 population tuning data from
% Henning, Ramos-Traslosheros et al.
% Called from fitT4T5FlowField, within the script T4T5_fitFlowField_expanded_CV.
% Inputs:
% - Vector components of visual landmark position in a sphere given
%   in spherical coordinates (T4T5 receptive field position):
%   - phiData, thetaData, rData
% - Flow vector components (T4T5 direction selectivity index), defined as
%   tangetial to the spherical surface at visual landmark point given by
%   (phiData, thetaData, rData):
%   - uData, vData
% - flowType: Integer scalar defining the type of flow field to be fitted
%   - 0: translational plus rotational flow 'T+R'
%   - 1: pure translational flow 'T'
%   - 2: pure roational flow 'R'
%   - 3: no optic flow field, just a uniform vector field 'Uni'


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


function [vT, vR, vUni, loss] = fitFlowField_CV(phiData, thetaData, rData, uData, vData, flowType)

% TODO include normalization options and interpolated and raw data fits.
lossFunction = @(x) -getLinearProjection(deg2rad(phiData.train), deg2rad(thetaData.train), ...
                                         rData.train, uData.train, vData.train, x(1:3), x(4:6));
lossFunctionUni = @(x) -getLinearProjectionUniformField(deg2rad(phiData.train), ...
                                                        deg2rad(thetaData.train), ...
                                                        rData.train, uData.train, vData.train, x);

lossFunctionTest = @(x) -getLinearProjection(deg2rad(phiData.test), deg2rad(thetaData.test), ...
                                             rData.test, uData.test, vData.test, x(1:3), x(4:6));
lossFunctionUniTest = @(x) -getLinearProjectionUniformField(deg2rad(phiData.test), ...
                                                            deg2rad(thetaData.test), ...
                                                            rData.test, uData.test, vData.test, x);
% Model params initial conditions.
meanU = nanmean(uData.train(:)) * (1 + 0.5 * rand);
meanV = nanmean(vData.train(:)) * (1 + 0.5 * rand);

vT0 = [meanU, meanV, 0]; % m/s
vT0 = 2 * rand(1, 3) - 1;
vT0 = vT0 / norm(vT0);
vR0 = [meanU, meanV, 0]; % rad/s
vR0 = 2 * rand(1, 3) - 1;
vR0 = vR0 / norm(vR0);
vUni0 = [meanU, meanV];
vUni0 = 2 * rand(1, 2) - 1;
vUni0 = vUni0 / norm(vUni0);

%%
vR = nan(1, 3);
vT = nan(1, 3);
vUni = nan(1, 2);

switch flowType
    case 0 % General flow
        lBound = -ones(1, 6);
        uBound =  ones(1, 6);
        flowVector = fmincon(lossFunction, [vT0 vR0], [], [], [], [], ...
                             lBound, uBound, @isUnitVector2);
        vT = flowVector(1: 3);
        vR = flowVector(4: 6);
        loss = lossFunctionTest(flowVector);
    case 1 % Pure translation
        lossFunctionT = @(x) lossFunction([x(1:3), 0 0 0]);
        lossFunctionTTest = @(x) lossFunction([x(1:3), 0 0 0]);
        lBound = -ones(1, 3);
        uBound =  ones(1, 3);
        flowVector = fmincon(lossFunctionT, vT0, [], [], [], [], ...
                             lBound, uBound, @isUnitVector);
        vT = flowVector(1: 3);
        vR = [0 0 0];
        loss = lossFunctionTTest(flowVector);
    case 2 % Pure rotation
        lossFunctionR = @(x) lossFunction([0 0 0, x(1:3)]);
        lossFunctionRTest = @(x) lossFunction([0 0 0, x(1:3)]);
        lBound = -ones(1, 3);
        uBound =  ones(1, 3);
        flowVector = fmincon(lossFunctionR, vR0, [], [], [], [], ...
                             lBound, uBound, @isUnitVector);
        vT = [0 0 0];
        vR = flowVector(1: 3);
        loss = lossFunctionRTest(flowVector);
    case 3 % Uniform vector
        lBound = -ones(1, 2);
        uBound =  ones(1, 2);
        flowVector = fmincon(lossFunctionUni, vUni0, [], [], [], [], ...
                             lBound, uBound, @isUnitVector);
        vUni = flowVector(1: 2);
        loss = lossFunctionUniTest(flowVector);
    otherwise % General flow
        lBound = -ones(1, 6);
        uBound =  ones(1, 6);
        flowVector = fmincon(lossFunction, [vT0 vR0], [], [], [], [], ...
                             lBound, uBound, @isUnitVector2);
        vT = flowVector(1: 3);
        vR = flowVector(4: 6);
        loss = lossFunctionTest(flowVector);
end

end

function [c, ceq] = isUnitVector(x)
         ceq = 1 - norm(x);
         c = [];
end
function [c, ceq] = isUnitVector2(x)
         ceq(1) = 1 - norm(x(1:3));
         ceq(2) = 1 - norm(x(4:6));
         c = [];
end



function linearProjection = getLinearProjectionUniformField(...
    phiData, thetaData, rData, uData, vData, vUni)
    u = vUni(1) * ones(size(phiData));
    v = vUni(2) * ones(size(thetaData));
    validPoints = ~isnan(uData) & ~isnan(vData);
    u = reshape(u(validPoints), [], 1);
    uData = reshape(uData(validPoints), [], 1);
    v = reshape(v(validPoints), [], 1);
    vData = reshape(vData(validPoints), [], 1);
    linearProjection = sum(u .* uData) + sum(v .* vData);
    % Normalize by the total norm of all vectors
    totalNorm = sum(sqrt((u.^2 + v.^2) .* (uData.^2 + vData.^2)));
    linearProjection = linearProjection / totalNorm;
end


function linearProjection = getLinearProjection(phiData, thetaData, rData, ...
                                                uData, vData, vT, vR)
    [~, ~, u, v] = calculateSphereFlow(phiData, thetaData, rData, vT, vR);
    validPoints = ~isnan(uData) & ~isnan(vData);
    u = reshape(u(validPoints), [], 1);
    uData = reshape(uData(validPoints), [], 1);
    v = reshape(v(validPoints), [], 1);
    vData = reshape(vData(validPoints), [], 1);
    linearProjection = sum(u .* uData) + sum(v .* vData);
    % Normalize by the total norm of all vectors
    totalNorm = sum(sqrt((u.^2 + v.^2) .* (uData.^2 + vData.^2)));
    linearProjection = linearProjection / totalNorm;
end
