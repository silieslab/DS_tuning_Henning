% AnalyzeFlowFits_expanded
% Script to generate figure panels related to the optic flow model in
% Henning, Ramos-Traslosheros et al.
% Requires data generated by the script:
% - T4T5_fitFlowField_expanded_CV

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
codeDir = pwd;
[parentDir, ~] = fileparts(pwd);
cd(parentDir);
addpath(genpath(codeDir));
%% Load fit data
fileName = 'flowFitsCV10FoldRand';
load(['Model_output' filesep fileName '.mat']);
%% set output path
outputDir = 'Model_figures';
if ~exist(['.' filesep outputDir], 'dir'); mkdir(outputDir); end
figDir = [outputDir filesep 'flowAnalysis'];
if ~exist(figDir, 'dir'); mkdir(figDir); end
set(groot, 'DefaultFigureVisible', 'on');
set(groot, 'DefaultFigureWindowStyle', 'normal');

FilenameCell = {'LayerA1', 'LayerA2', 'LayerB1', 'LayerB2', 'LayerC', 'LayerD'};

%% Compare losses for cross validated fits
% Used in extended data figure 5
figure('Position', [0 0 500 1000]);
legendStr = cell(1, 4);
for iLayer = 1: size(vRCell, 1)                
    hSubAx(iLayer) = subplot(ceil(size(vRCell, 1) / 2), 2, iLayer);
    iBar = 1;
    for jFlowType = 0: 3
        loss = -squeeze(lossArray(iLayer, jFlowType + 1, :));
        cmap = colormap(tab20c);
        col = cmap(jFlowType * 4 + 1, :);
        hold on;
        line(iBar + 0.3 * [-1 1], mean(loss) * [1 1], 'Color', col, 'LineWidth', 3);
        scatter(iBar*ones(size(loss)) + 0.25 * (rand(size(loss)) - 0.5), loss, [], col);
        iBar = iBar + 1;
    end   
    uniLoss = squeeze(lossArray(iLayer, 4, :));
    lossDiff = (uniLoss' - squeeze(lossArray(iLayer, 1: 3, :)));
    meanDiff = mean(lossDiff, 2)
    uniVs = arrayfun(@(x) signrank(uniLoss, squeeze(lossArray(iLayer, x, :)), 'tail', 'right'), 1:3);
    for iModel = 1: 3
        scatter(iModel, -mean(lossArray(iLayer, iModel, :), 3), 50, ...
                'MarkerFaceColor', [1 1 1] * (uniVs(iModel) < 0.05 / 3));
        hold on
    end
    title(FilenameCell{iLayer});
    yLims = get(gca, 'YLim');
    ylabel('Linear Projection');
    xticks(1: 1: 4);
    xticklabels({'T+R', 'T', 'R', 'Uni'});
    xlim([0 5]);
    set(gca, 'TickDir', 'out');
end

%%
setEqualAxes = 0;
if setEqualAxes
    linkaxes(hSubAx, 'y');
    print(gcf, [figDir filesep 'fitQualityEqAx-' fileName '.pdf'], '-dpdf');
else
    print(gcf, [figDir filesep 'fitQuality-' fileName '.pdf'], '-dpdf');
end
%% Compare losses for cross validated fits, by differences
figure('Position', [0 0 500 1000]);
legendStr = cell(1, 4);
for iLayer = 1: size(vRCell, 1)             
    hSubAx(iLayer) = subplot(ceil(size(vRCell, 1) / 2), 2, iLayer);
    iBar = 1;
    uniLoss = -squeeze(lossArray(iLayer, 4, :));
    for jFlowType = 0: 2
        loss = -squeeze(lossArray(iLayer, jFlowType + 1, :));
        lossDiff = loss - uniLoss;
        cmap = colormap(tab20c);
        col = cmap(jFlowType * 4 + 1, :);
        hold on;
        line(iBar + 0.3 * [-1 1], mean(lossDiff) * [1 1], 'Color', col, 'LineWidth', 3);
        scatter(iBar*ones(size(loss)) + 0.25 * (rand(size(loss)) - 0.5), lossDiff, [], col);
        
        uniVs = signrank(loss, uniLoss, 'tail', 'right');
        
        scatter(iBar, mean(lossDiff), 50, 'square',...
                'MarkerFaceColor', [1 1 1] * (uniVs < 0.05 / 3));
        iBar = iBar + 1;
    end
    refline(0, 0);
    title(FilenameCell{iLayer});
    ylabel('\Delta (Linear Projection)');
    xticks(1: 1: 3);
    xticklabels({'T+R', 'T', 'R', 'Uni'});
    xlim([0 4]);
end
%%
print(gcf, [figDir filesep 'fitQualityStatisticsVsUni-' fileName '.pdf'], '-dpdf');
%% Pooled statistics across layers and cross-validation runs
% minus because algorithm was minimizing, return to positive linear
% projection value.
pooledLossesAcrossLayers = -reshape(permute(lossArray, [2 1 3]), 4, []);

%% Used in figure 4
uniVs = arrayfun(@(x) signrank(pooledLossesAcrossLayers(x, :), ...
                               pooledLossesAcrossLayers(4, :), ...
                               'tail', 'right'), 1: 3);

pooledRelativeLosses = pooledLossesAcrossLayers(1: 3, :) - ...
                       pooledLossesAcrossLayers(4, :);

figure('Position', [0 0 900 400]);
hAx(1) = subplot(1, 3, 1);
x = repmat(1: 3, 1, size(pooledRelativeLosses, 2)) + 0.4 * (rand(1, numel(pooledRelativeLosses)) - 0.5);
scatter(x, pooledRelativeLosses(:), 'filled', 'MarkerFaceAlpha', 0.3);
% refline(0,0)
for iBar = 1: 3
    line(iBar + 0.4 * [-1 1], mean(pooledRelativeLosses(iBar, :)) * [1 1], 'Color', 'k', 'LineWidth', 3);
end

ylabel('\Delta (Linear Projection)')
xticks(1: 1: 3);
xticklabels({'T+R', 'T', 'R', 'Uni'});


hAx(2) = subplot(1, 3, 2);
hBar = bar(mean(pooledRelativeLosses, 2), 'EdgeColor', 'none');
hold on
errorbar(mean(pooledRelativeLosses, 2), ...
         std(pooledRelativeLosses, [], 2) ./ ...
         sqrt(size(pooledRelativeLosses, 2)), ...
         'LineStyle', 'none', 'LineWidth', 2, 'Color', [1 1 1] * 0.25, 'CapSize', 0)
xticks(1: 1: 3);
xticklabels({'T+R', 'T', 'R', 'Uni'});
arrayfun(@(x) text(2, (1-0.07*x)*hAx(2).YLim(2), ['p = ' num2str(uniVs(x), '%0.1e')]), 1: 3);
% legend(hBar, arrayfun(@(x) ['p = ' num2str(uniVs(x), '%0.2f')], 1:3, 'uni', 0))

hAx(3) = subplot(1, 3, 3);
bar(mean(pooledLossesAcrossLayers, 2))

ylabel('Linear Projection')
xticks(1: 1: 4);
xticklabels({'T+R', 'T', 'R', 'Uni'});


x = repmat(1: 4, 1, size(pooledLossesAcrossLayers, 2)) + 0.4 * (rand(1, numel(pooledLossesAcrossLayers)) - 0.5);
scatter(x, pooledLossesAcrossLayers(:), 'filled', 'MarkerFaceAlpha', 0.3);
% refline(0,0)
for iBar = 1: 4
    line(iBar + 0.4 * [-1 1], mean(pooledLossesAcrossLayers(iBar, :)) * [1 1], 'Color', 'k', 'LineWidth', 3);
end

ylabel('Linear Projection')
xticks(1: 1: 4);
xticklabels({'T+R', 'T', 'R', 'Uni'});

set(hAx, ...
    'Box', 'off', ...
    'tickdir', 'out', ...
    'yminortick', 'off');
%%
print(gcf, [figDir filesep 'fitQualityPooled-' fileName '.pdf'], '-dpdf');