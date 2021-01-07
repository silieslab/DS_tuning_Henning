%%%%%%%%%%%%%%
% This skript automatically finds ROIs of single T4 and T5 axon terminals
% - It uses the selective responses of ON and OFF polarity to distinguish T4 and T5
% - and the direction selective response to distinguish between subtypes,
% - to cluster pixels belonging to one cell it utilizes the response timingand Pixel coordinates
%
%%%%%%%%%%%%%%


% This skript loads the imaging data or p-data if available
% Then it compares the stimulus that was used for recording and if it was the desired stimulus
% (Drifting edges into four or eight directions) it calculates clusters using this function:  FindCanROIs_v9
% if it was another stimulus used it calls clusters from a previous imaging series and merges
% them on the data


addpath(genpath('subscripts'))

folder='Data/RawDataExample/200729_Fly4'; %specify pdata folder
% NOTE: pdata files contain the almost raw image data, which were only
% preprocessed using SIMA python package for Motion alignment. (These raw
% data were too big for the upload on github, only one example is given.
% For more raw data files please contact: mhenning@uni-minz.de or msilies@uni-mainz.de



Crop=[];

SelectLayers=1;    %Put 1 if you want to select Layers manually, 0 if not!
% CalcClusters=1;  %Put 1 if you want to use cluster algorithm and 2 if you want to use existing clusters
Merge=0; %Put a zero if you do not want to use the function that merges ROIs of the same cell terminal

ROISize=2.5;  %expected cluster size based on approximate size of cell terminal
%%
for image=[4] % If more T-series were recorded for one imaging plane it will loop through all 'Images' and calculate cluster/ROIs if the correct stimulus was shown, or assign clusters from a previous recording..
    
    
    close all
    file=[folder, '/Image',num2str(image)];
    
    
    disp(file);
    curDir = pwd;
    cd(file);
    
    m = dir('data_file_SIMA_only_m*');
    p = dir('*pData_SIMA_only_m.mat');
    
    if ( ~ (isempty(p))) % If pData are available load them and just add new information
        load(p.name);
        out=strct;
    elseif (~ (isempty(m))) && ((isempty(p))) % If no pData but only data_file, then load this.
        load('data_file_SIMA_only_m');
    end
    
    
    DesiredStim='DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs_optimized';
    DesiredStim2='Search_DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs_optimized'; %Sometimes I accidentally pick the search stimulus, which does not matter, analysis stays the same
    DesiredStim3='DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_8Dirs_optimized';
    
    if strcmp(DesiredStim, out.stim_type)
        CalcClusters=1;
    elseif strcmp(DesiredStim2, out.stim_type)
        CalcClusters=1;
    elseif strcmp(DesiredStim3, out.stim_type)
        CalcClusters=1;
    else
        CalcClusters=2;
    end
    
    % % If use previos clusters, then copy CA_information of previous Image into
    % the actual one
    
    if CalcClusters==2
        if SelectLayers==0
            copyfile([folder, '/Image',num2str(image-1),'/CA_information.mat'],...
                [folder, '/Image',num2str(image),'/CA_information.mat'])
        elseif SelectLayers==1
            copyfile([folder, '/Image',num2str(image-1),'/CA_information_ManuallySelect.mat'],...
                [folder, '/Image',num2str(image),'/CA_information_ManuallySelect.mat'])
        end
    end
    
    % % %
    
    
    
    [out,Crop] = FindCanROIs_v9(out,CalcClusters,SelectLayers,Crop,ROISize); %2nd input=1 if clusters should be calculated
    %=2 if masks are already existing
    %3rd input=1 (only for V5 upwards) if you want to select Layers manually
    %input=0 if not
    %4th input= Size of expected ROI in uM,so if you cell type of interest is app. 5 uM big in x and y then put 5
    
    
    
    
    if SelectLayers==1
        Masks= out.ClusterInfo_ManuallySelect.masks_CA;
    else
        Masks=out.ClusterInfo.masks_CA;
    end
    
    Clusters= figure;
    Try1=nan(size(Masks{1,1}));
    for i=1:length(Masks)
        ClusterType=i;
        x=find(Masks{1,i});
        Try1(x)=ClusterType;
    end
    Try1(isnan(Try1))=0;
    imagesc(Try1)
    
    colormap('colorcube');
end



%% Plot results: 

[h2,h3,h4,out] = FFFlash_res_display_1ch_mh_forV6(out,SelectLayers);  %Plots the results of ROI selection

% Bright traces belong to T4 cells, Dark traces to T5 cells

if isfield(out,'ch1a_old')
    out=rmfield(out, 'ch1a_old');
end

save_processed_data_1ch_eni_SIMA(out, file);%


% To check which stimtype:
d = dir('_stimulus_*');
fid = fopen(d.name,'r');
currline = fgetl(fid);
ind = strfind(currline,'\');
fprintf('stimulus: %s\n',currline(ind(end)+1:end));
fclose(fid);

