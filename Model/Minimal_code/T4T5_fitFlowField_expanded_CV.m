% T4T5_fitFlowField_expanded_CV
% Script to fit the optic flow fields to the T4T5 population tuning data from
% Henning, Ramos-Traslosheros et al.
% Requires experimental data in the file:
% '../Input_data/t4t5_flowFields.mat'
% Saves data to 'Model_output' directory, output data is used by the scripts:
% - AnalyzeFlowFits_expanded
% - MakeFigureBestFlowFieldFits

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


close all;
clearvars;
clc;

%% move to parent directory
[parentDir, ~] = fileparts(pwd);
cd(parentDir);
% addpath(genpath(parentDir));
%% Load data
load(['Input_data' filesep 't4t5_flowFields.mat']);
%% set output path
outputDir = 'Model_output';
if ~exist(['.' filesep outputDir], 'dir'); mkdir(outputDir); end

%% Set figure settings
set(groot, 'DefaultFigureVisible', 'on');
set(groot, 'DefaultFigureWindowStyle', 'normal');

%%
DataCell = {[Coord_T4A_I Coord_T5A_I], [Quiver_T4A_I Quiver_T5A_I];
            [Coord_T4A_II Coord_T5A_II], [Quiver_T4A_II Quiver_T5A_II];
			[Coord_T4B_I Coord_T5B_I], [Quiver_T4B_I Quiver_T5B_I];
			[Coord_T4B_II Coord_T5B_II], [Quiver_T4B_II Quiver_T5B_II];
			[Coord_T4C Coord_T5C], [Quiver_T4C Quiver_T5C];
			[Coord_T4D Coord_T5D], [Quiver_T4D Quiver_T5D]};
FilenameCell = {'LayerA1', 'LayerA2', 'LayerB1', 'LayerB2', 'LayerC', 'LayerD'};
%%
nCellTypes = size(DataCell, 1);
nModels = 4;
nKCrossValidation = 10;
normalizeData = 0;
interpData = 0;

vTCell = cell(nCellTypes, nModels, nKCrossValidation);
vRCell = cell(nCellTypes, nModels, nKCrossValidation);
vUniCell = cell(nCellTypes, nModels, nKCrossValidation);
lossArray = nan(nCellTypes, nModels, nKCrossValidation);

nRepeats = 10;
%%
for iLayer = 1: nCellTypes
    XAll = DataCell{iLayer, 1};
    VAll = DataCell{iLayer, 2};
    crossValInds = crossvalind('Kfold', VAll, nKCrossValidation);
    crossValIndsCell{iLayer} = crossValInds;
    for jFlowType = 0: nModels - 1
        for kRun = 1: nKCrossValidation
            testInds = (crossValInds == kRun);
            X.train = XAll(:, ~testInds);
            X.test = XAll(:, testInds);
            V.train = VAll(~testInds);
            V.test = VAll(testInds);
            flowType = jFlowType;
            fileName = [FilenameCell{iLayer} ...
                        '-Interp' num2str(interpData, '%d') ...
                        '-Norm' num2str(normalizeData, '%d') ...
                        '-FlowType' num2str(flowType, '%d')];
            tmpLoss = 100;
            loss = 100;
            for lRepeat = 1: nRepeats
                fprintf('Running... %s, fold (%d/%d), run (%d/%d), loss %0.4f->%0.4f \n', ...
                        fileName, kRun, nKCrossValidation, lRepeat, nRepeats, tmpLoss, loss)
                [vT, vR, vUni, loss] = fitT4T5FlowField(X, V, interpData, ...
                                                        normalizeData, ...
                                                        flowType, fileName);

                if loss < tmpLoss
                    display(loss)
                    display(tmpLoss)
                    tmpLoss = loss;

                    vTCell{iLayer, jFlowType + 1, kRun} = vT;
                    vRCell{iLayer, jFlowType + 1, kRun} = vR;
                    vUniCell{iLayer, jFlowType + 1, kRun} = vUni;
                    lossArray(iLayer, jFlowType + 1, kRun) = loss;
                end
                clc
            end
        end
        close all;
    end
end
%% Store the motion vectors, fit quality and data partition used in each CV.
outputFileName = ['flowFitsCV' num2str(nKCrossValidation, '%d')  'FoldRand.mat'];
outputFilePath = [outputDir filesep outputFileName];
save(outputFilePath, 'vTCell', 'vRCell', 'vUniCell', 'lossArray', 'crossValIndsCell', '-v7.3')
%%

function [vT, vR, vUni, loss] = fitT4T5FlowField(X, V, interpData, normalizeData, flowType, fileName)
    %% Fit flow
    [phiData.train, thetaData.train, rData.train, uData.train, vData.train] = ...
    	getT4T5Data(X.train, V.train, interpData, normalizeData);
    [phiData.test, thetaData.test, rData.test, uData.test, vData.test] = ...
    	getT4T5Data(X.test, V.test, interpData, normalizeData);
    [vT, vR, vUni, loss] = fitFlowField_CV(phiData, thetaData, rData, uData, vData, flowType);
end

function [phiData, thetaData, rData, uData, vData] = getT4T5Data(X, V, interpData, normalizeData)
    if interpData
        [Xq, Yq, Vq] = interpolateGridQuiver(X, V);
    else
        validInds = ~any(isnan(X), 1);
        Xq = X(2, validInds) - 34;
        Yq = X(1, validInds) + 36;
        Vq = V(validInds);
    end

    phiData = Xq;
    thetaData = Yq;
    uData = real(Vq);
    vData = imag(Vq);
    rData = ones(size(phiData));

    if normalizeData
        normData = sqrt(uData .^ 2 + vData .^ 2);
        uData = uData ./ normData;
        vData = vData ./ normData;
    end
end

function [Xq, Yq, Vq] = interpolateGridQuiver(X, V)
[Xq, Yq] = meshgrid(-35: 3: 45);
validInds = ~any(isnan(X), 1);
Vq = griddata(X(2, validInds) - 34, X(1, validInds) + 36, ...
              V(validInds), Xq, Yq, 'cubic');

end
