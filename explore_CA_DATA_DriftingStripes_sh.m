
%%%%%%%%%%%%%%
% This skript collects information from all p-data for all flies that match
% the criteria, e.g. correct genotype and stimulus shown (in this case for
% all T4\T5 flies that were shown stripes moving into 8 directions
% (data shown in Fig.1, Fig.4, and supplementary Fig.4 shown in the manuscript Henning et al., 2021)

%--> calculates the tuning vector Z for each cell

%--> saves Data on server in subfolders


%%%%%%%%%%%%%%

clear all;
close all;

Homepath='/Volumes/SILIESLAB/MiriH/Github_Reps/DS_tuning_Henning'; 
addpath(genpath(Homepath))


pdatapath='/Volumes/Seagate/Documents/2p-imaging/Miriam_pData';
cd(pdatapath);
    
  
database_select_samples_mh_DStuning; % get the right stuff...

%%
% Cellect data for dark (OFF) bars
Foldertosave='Data/Data_Stripes/Dark_Stripes' ;

% This selects data shown in Fig.1 and 4 
Condition='i_T4T5Rec_Kir_Control.*i_mb_lumdec_20degs_5deg_8dir_FullContrast.*i_responsivefly.*~moving.*i_PDATA_SIMA';

f_mb_T4T5 = find(eval(Condition));
    
%     fT4T5_mb = create_neuron_structure_all_CA_data_mh(f_mb_T4T5,'ClusterInfo_ManuallySelect'); % or just 'ClusterInfo'
   
T4T5_mb = load_neuron_data10Hz_ks_CA_average_Stripes_mh(f_mb_T4T5,pdatapath,'ClusterInfo_ManuallySelect');
        
    
save([Homepath,'/',Foldertosave, '/processed_Data_SIMA_CS5_sh'],'T4T5_mb');





% Cellect data for bright (ON) bars
Foldertosave='Data/Data_Stripes/Bright_Stripes' ;

% This selects data shown in Fig.1 and 4 
Condition='i_T4T5Rec_Kir_Control.*i_mb_luminc_20degs_5deg_8dir_FullContrast.*i_responsivefly.*~moving.*i_PDATA_SIMA';

f_mb_T4T5 = find(eval(Condition));
    
%     fT4T5_mb = create_neuron_structure_all_CA_data_mh(f_mb_T4T5,'ClusterInfo_ManuallySelect'); % or just 'ClusterInfo'
   
T4T5_mb = load_neuron_data10Hz_ks_CA_average_Stripes_mh(f_mb_T4T5,pdatapath,'ClusterInfo_ManuallySelect');
        
    
save([Homepath,'/',Foldertosave, '/processed_Data_SIMA_CS5_sh'],'T4T5_mb');




%% Now again for the data set, where each fly was recorded at different z-depth 



% Collect data for dark (OFF) bars
Foldertosave='Data/Data_Stripes/Dark_Stripes' ;

% This selects data shown in Fig.1 and 4 
Condition='i_T4T5Rec_LayerControl.*i_mb_lumdec_20degs_5deg_8dir_FullContrast.*i_responsivefly.*~moving.*i_PDATA_SIMA';

f_mb_T4T5 = find(eval(Condition));
    
%     fT4T5_mb = create_neuron_structure_all_CA_data_mh(f_mb_T4T5,'ClusterInfo_ManuallySelect'); % or just 'ClusterInfo'
   
T4T5_mb = load_neuron_data10Hz_ks_CA_average_Stripes_mh(f_mb_T4T5,pdatapath,'ClusterInfo_ManuallySelect');
        
    
save([Homepath,'/',Foldertosave, '/processed_Data_SIMA_CS5_Layer'],'T4T5_mb');





%% Cellect data for bright (ON) bars
Foldertosave='Data/Data_Stripes/Bright_Stripes' ;

% This selects data shown in Fig.1 and 4 
Condition='i_T4T5Rec_LayerControl.*i_mb_luminc_20degs_5deg_8dir_FullContrast.*i_responsivefly.*~moving.*i_PDATA_SIMA';

f_mb_T4T5 = find(eval(Condition));
    
%     fT4T5_mb = create_neuron_structure_all_CA_data_mh(f_mb_T4T5,'ClusterInfo_ManuallySelect'); % or just 'ClusterInfo'
   
T4T5_mb = load_neuron_data10Hz_ks_CA_average_Stripes_mh(f_mb_T4T5,pdatapath,'ClusterInfo_ManuallySelect');
        
    
save([Homepath,'/',Foldertosave, '/processed_Data_SIMA_CS5_Layer'],'T4T5_mb');



