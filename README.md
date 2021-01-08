# DS_tuning_Henning
This repository contains data and code shown in the manuscript Henning et al., 2021


## Code

### Plot_DS_tuning_Analysis_CompassPlots:
Plots Tuning vectors of all recorded neurons in compass plots and histograms, also splits data from layer A and B into subtypes based on results from SNOB analysis (**Fig.1, Extended Data Fig.1**)


### Find_underlying_distr_SNOB.m: 
Finds underlying Gaussian distributions in the tuning data of neurons from each layer. It uses the SNOB analysis, a finite Gaussian mixture model: (Wallace et al., 2005) *downloaded from: https://de.mathworks.com/matlabcentral/fileexchange/72310-flexible-mixture-models-for-automatic-clustering* (**Extended Data Fig. 1**)


### Plot_ClassIdent_back_on_brain.m: 
Plots the ROIs back on to the average calcium image of the brain and color codes them by subtype identity based on SNOB Analysis (**Figure 2, Extended Data Figure 2**)
Plots the examples shown in the Manuscript, but changing the IMAGES to a range of [1:114] it will plot all examples.


### Plot_GlobalVSLocal_Tuning.m:
Plots the compass plots as shown for one example in (**Fig. 2c,d,e**) for all flies. 
Additionally plots the tuning distribution (std of angular preference) for each subtype as shown in (**Extended Data Fig.3**)


### today_single_1ch_eni_CA_Kir_SIMA.m:
script for automatic ROI selection of T4/T5 axon terminals based on their unique response properties to ON and OFF moving edges

Saves: 
* **CA_information_ManuallySelect** (contains information about the extracted ROIs that can be used for subsequent recordings to use ROI location 

* ** *pData_SIMA_only_m.mat (contains a structure with all relevant imaging information and the extracted ROIs
* * ch1: original imaging frames 
* * ch1a_crop: cropped image 
* * CLusterInfo_ManuallySelect: contains Clustering information (for each ROI the **Layer** identity, **T4_T5** identity, the average Response of each ROI to each direction of the stimulus, averaged across epochs (**avSignal1_CA**), background subtracted response for each stimulus direction (**dSignal1_CA**), **DSI**- and **CSI**-thresholds, the **masks** for each ROI, and the the Background from Otsu thresholding (**Real_Background**) 


### explore_CA_DATA_DriftingEdge8Dir_sh.m:
Collects data from PData structure for all recorded flies in one structure named: **processed_Data_SIMA_CS5**
Calculates the tuning vectors *Z* for each cell
%Data structure with the following fields:
%T4T5_mb(#Recording).Flyname= Identity of the recording;   
%T4T5_mb(#Recording).NROIs= Number of ROIs for each layer identity
%T4T5_mb(#Recording).Z= Tuning vector for each ROI
%T4T5_mb(#Recording).Masks= Masks for each ROI;
%T4T5_mb(#Recording).MAXdeg=Maximal response for each neuron and the direction it responded most to;
%T4T5_mb(#Recording).AV= Average intensity projection of the LP_cropped version
%T4T5_mb(#Recording).PixelSize= Size of the Pixel im microns;





________________________________________________
## Data 

### Data/Data_Edges 

#### --> processed_Data_SIMA_Cs5_sh.mat: 
contains the structure T4T5_mb: This structure contains all calcium imaging data from all cells recorded in a total of 14 flies (114 recordings) with the following fields: 	
* **Flyname**: date of recording, Fly ID, Image ID 
* **NROIs**: number of ROIs for T4 and T5 cells for each lobula plate layer 
* **Masks**: masks for each cell/ROI (T4/T5 axon terminal) 
* **Z**: tuning vector calculated after Mazurek 2014 
* **MAXdeg**: Preferred direction, based on maximal response amplitude
* **AV**: Average calcium image of the recoding (averaged across frames)

#### --> processed_Data_ROI_rf.mat:
contains the structure T4T5_mb_new: This structure contains the receptive field locations for all cells (114 recordings) with the following fields: 
* **RFCenter**: X and Y coordinates for the center of the neuron's receptive field 			based on back transformation of responses to 8 directions to the 			screen(see methods in the paper) 
* **CellID**

#### --> Turn_info.txt:
Flyname, orientation to the screen (0, 45 or 90 deg), z-depth location



#### --> Snob_Cluster_Info.mat: 
contains the structure ClusterR:
* TA_T4: For T4 cells from Layer A: Class(Subtype) Identity from SNOB analysis, 		assigned based on highest probability of falling into one of the 			underlying distributions(see methods)
*TA_T5: same for T5 cells from layer A, and so on 
* mm_A: Model statistics for Layer A, and so on 








