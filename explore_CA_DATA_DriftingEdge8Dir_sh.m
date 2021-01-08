%%%%%%%%%%%%%%
% This skript collects information from all p-data for all flies that match
% the criteria, e.g. correct genotype and stimulus shown (in this case for
% all T4\T5 flies that were shown drifting edges moving into 8 directions
% (big data set shown in the paper Henning et al., 2021)

%--> calculates the tuning vector Z for each cell
%--> plots average responses across flies
%--> saves Data on server in subfolders
%%%%%%%%%%%%%%


Foldertosave='Data/Data_Edges' ;
Homepath='/Volumes/SILIESLAB/MiriH/Github_Reps/DS_tuning_Henning'; 
addpath(genpath(Homepath))


% CONSTANTS

Condition='i_T4T5Rec_Control.*i_mb_movEdgesONOFF_optimized_8Dir.*i_responsivefly.*~moving.*i_PDATA_SIMA';

% Define path where the pdata are saved (pdata for the whole data set can
% be given upon request to msilies@uni-mainz.de or mhenning@uni-mainz.de

pdatapath='/Volumes/Seagate/Documents/2p-imaging/Miriam_pData';
cd(pdatapath);

database_select_samples_mh_DStuning; % get the right stuff...
f_mb_T4T5 = find(eval(Condition));

cd(pdatapath);
T4T5_mb = load_neuron_data10Hz_ks_CA_average_Edges_8Dir_mh(f_mb_T4T5,pdatapath,'ClusterInfo_ManuallySelect');
%Data structure with the following fields:
%T4T5_mb(#Recording).Flyname= Identity of the recording;   
%T4T5_mb(#Recording).NROIs= Number of ROIs for each layer identity
%T4T5_mb(#Recording).Z= Tuning vector for each ROI
%T4T5_mb(#Recording).Masks= Masks for each ROI;
%T4T5_mb(#Recording).MAXdeg=Maximal response for each neuron and the direction it responded most to;
%T4T5_mb(#Recording).AV= Average intensity projection of the LP_cropped version
%T4T5_mb(#Recording).PixelSize= Size of the Pixel im microns;

%%
% SAVE Data structure

save([Homepath,'/',Foldertosave, '/processed_Data_SIMA_CS5_sh'],'T4T5_mb');


