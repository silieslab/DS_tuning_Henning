
%%%%%%%%%%%%%%
% This skript plots all figures shown in Figure1 and Extended Data Figure1
%%%%%%%%%%%%%% 

addpath(genpath('subscripts'))

% Load preprocessed Data matrix
Conditions.Control=load('Data/Data_Edges/processed_Data_SIMA_CS5_sh.mat');
% Load text file with imaging conditins (e.g. z-depth and orientation to the screen)
[fname,turn,Zdepth]=textread('Data/Data_Edges/Turn_info.txt','%s %f %f','headerlines',0,'delimiter','\t');

for NF=1:size(Conditions.Control.T4T5_mb,2)
    Conditions.Control.T4T5_mb(NF).turn=turn(NF);
    Conditions.Control.T4T5_mb(NF).z_depth=Zdepth(NF);
end


%% Circular Rayleigh test: 
% extract all tuning vectors from all flies: 
Z = averageDirectionVectors(Conditions.Control.T4T5_mb);
% Test T4 and T5 separatly:
[pval, z]=circ_rtest([Z.T4A.ALL,Z.T4B.ALL,Z.T4C.ALL,Z.T4D.ALL])
[pval, z]=circ_rtest([Z.T5A.ALL,Z.T5B.ALL,Z.T5C.ALL,Z.T5D.ALL])
% Test T4 and T5 together:
[pval, z]=circ_rtest([Z.T4A.ALL,Z.T4B.ALL,Z.T4C.ALL,Z.T4D.ALL,Z.T5A.ALL,Z.T5B.ALL,Z.T5C.ALL,Z.T5D.ALL])

% ---> In all cases P is highly significant, which means that the
% distribution is not uniform!

%% Figure 1a: Compass Plots with z-depth Information

CONi='Control'; 
ConditionI=eval(['Conditions.',CONi]);
    
Average=false; 

[F1]=CompassPlot_Zdepth_8dirsEdge(ConditionI,CONi,Average);



%% Figure 1b,c and Extended Data Figure 1d: Circular Histogram, 

    
[F3,F4,F5]=CircHist_8dirsEdge_sh(ConditionI);




%%%%% NOTE %%%
% to plot results (Gaussians) from SNOB analysis as shown in Extended Data Figure 1c,
% open Find_underl_distr_SNOB.m




