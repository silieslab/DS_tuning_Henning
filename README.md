# DS_tuning_Henning
This repository contains data and code shown in the manuscript Henning et al., 2021


## Code

##NOTE: to use this code download the circular statistics toolbox from Philipp Berens 
P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009 
As well as: SNOB which is a Matlab implementation of finite mixture models: https://de.mathworks.com/matlabcentral/fileexchange/72310-flexible-mixture-models-for-automatic-clustering


### Plot_DS_tuning_Analysis_CompassPlots:
Plots Tuning vectors of all recorded neurons in compass plots and histograms, also splits data from layer A and B into subtypes based on results from SNOB analysis (**Fig.1**)


### Find_underlying_distr_SNOB.m: 
Finds underlying Gaussian distributions in the tuning data of neurons from each layer. It uses the SNOB analysis, a finite Gaussian mixture model: (Wallace et al., 2005) *downloaded from: https://de.mathworks.com/matlabcentral/fileexchange/72310-flexible-mixture-models-for-automatic-clustering* (**Fig.1**)


### Plot_ClassIdent_back_on_brain.m: 
Plots the ROIs back on to the average calcium image of the brain and color codes them by subtype identity based on SNOB Analysis (**Fig.3, fig S3**)
Plots the examples shown in the Manuscript, but changing the IMAGES to a range of [1:114] it will plot all examples.


### Plot_GlobalVSLocal_Tuning_plusMedial.m:
Plots the compass plots and circular histograms shown in (**Fig.3b,c**). 
Additionally plots the tuning distribution (std of angular preference) for each subtype as shown in (**fig.S4c**)


### Plot_DS_Tuning_per_ROI_on_brain.m
This skript plots all figures shown in (**Fig.4**)


### Plot_DS_Tuning_per_ROI_on_brain_LayerControl.m
This skript plots all figures shown in (**fig.S4**)


### Plot_DS_tuning_on_screen_stitched.m
Plots the tuning maps for each T4\T5 subtype as shown in (**Fig.5**, **fig.S5c**), as well as single fly data shown in (**Fig.2**) and the mean tuning vectors across 10deg wide bins from tuning maps of Fig 5 as shown in (**fig.S5a,b**). In addition to that it plots compass plots for different z-depth shown in (**Fig.3d**).


### Plot_RFCenter_Methods.m
This scripts plots the RF location on screen of one recording (**fig.S2c**)


### today_single_1ch_eni_CA_Kir_SIMA.m:
script for automatic ROI selection of T4/T5 axon terminals based on their unique response properties to ON and OFF moving edges
	
	Loads: 
	* data_file_SIMA_only_m.mat

	Saves: 
	* **CA_information_ManuallySelect** (contains information about the extracted ROIs that can be used for 	subsequent recordings to use ROI location 

	* **..pData_SIMA_only_m.mat** (contains a structure with all relevant imaging information and the 	extracted ROIs
	* * ch1: original imaging frames 
	* * ch1a_crop: cropped image 
	* * CLusterInfo_ManuallySelect: contains Clustering information (for each ROI the **Layer** identity, **T4_T5** identity, the average Response of each ROI to each direction 	of the stimulus, averaged across epochs (**avSignal1_CA**), background subtracted response for each stimulus direction (**dSignal1_CA**), **DSI**- and **CSI**-thresholds, the 	**masks** for each ROI, and the the Background from Otsu thresholding (**Real_Background**) 

	Uses:*FindCanROIs_v9.m, 


### explore_CA_DATA_DriftingEdge8Dir_sh.m:
Collects data from PData structure for all recorded flies in one structure named: **processed_Data_SIMA_CS5_sh.mat** , for more information see Data
Calculates the tuning vectors *Z* for each cell


### Spatial_RFs_8DirEdge.m: 
This scripts calculates the RF field centers of T4T5 cells based on responses to an 8 dir moving stimulus. 


### delay_analysis_T45_RF.py:
Finds the delay between moving and static stimulus. Loads the pickle files which contain responses for T4T5 ROIs to static or moving stimuli. Uses the henning_helper.py 


### all other functions not listed here are used by any of the above listed main scripts.



________________________________________________
## Data 


#### --> processed_Data_SIMA_CS5_sh_Edges.mat: 
contains the structure T4T5_mb: This structure contains all calcium imaging data from all cells recorded in a total of 14 flies (114 recordings) with the following fields: 	
* **Flyname**: date of recording, Fly ID, Image ID 
* **NROIs**: number of ROIs for T4 and T5 cells for each lobula plate layer 
* **Masks**: masks for each cell/ROI (T4/T5 axon terminal) 
* **Z**: tuning vector calculated after Mazurek 2014 
* **MAXdeg**: Preferred direction, based on maximal response amplitude
* **AV**: Average calcium image of the recoding (averaged across frames)
* **PixelSize**: Size of the Pixel im microns

#### --> processed_Data_ROI_rf_Edges.mat:
contains the structure T4T5_mb_new: This structure contains the receptive field locations for all cells (114 recordings) with the following fields: 
* **RFCenter**: X and Y coordinates for the center of the neuron's receptive field 			based on back transformation of responses to 8 directions to the 			screen(see methods in the paper) 
* **CellID**

#### --> Turn_info_Edges.txt:
Flyname, orientation to the screen (0, 45 or 90 deg), z-depth location



#### --> Snob_Cluster_Info_Edges.mat: 
contains the structure ClusterR:
* TA_T4: For T4 cells from Layer A: Class(Subtype) Identity from SNOB analysis, 		assigned based on highest probability of falling into one of the 			underlying distributions(see methods)
*TA_T5: same for T5 cells from layer A, and so on 
* mm_A: Model statistics for Layer A, and so on 


#### --> processed_Data_SIMA_CS5_Layer_Stripe_Bright.mat: processed data for responses to Bright Stripes in different depth
#### --> processed_Data_SIMA_CS5_Layer_Stripe_Dark.mat: processed data for responses to Dark Stripes in different depth
#### --> processed_Data_SIMA_CS5_Stripe_Bright.mat: processed data for responses to Bright Stripes
#### --> processed_Data_SIMA_CS5_Stripe_Dark.mat: processed data for responses to Dark 



#### --> data_file_SIMA_only_m.mat
Contains the almost raw imaging data from one example recording, which were only preprocessed using SIMA python package for Motion alignment. 
Data are too big to upload them for each recording, thus only one example is given. 
For more raw data files please contact: mhenning@uni-mainz.de or msilies@uni-mainz.de

#### --> pickle files
Contain response data of T4T5 cells to static and moving stimuli used for the calculation of the delay in: 'delay_analysis_T45_RF.py'





__________________________________________________

## Model of optic flow: 
Fits optic flow fields to the T4T5 data 


### --> AnalyzeFlowFits_expanded.m
### --> calculateSphereFlow.m
### --> fitFlowField_CV.m
### --> MakeFigureBestFlowFieldFits.m
### --> T4T5_fitFlowField_expanded_CV.m



## Input data: 

# t4t5_flowFields.mat
# Flow_field_Data_new_new.mat

