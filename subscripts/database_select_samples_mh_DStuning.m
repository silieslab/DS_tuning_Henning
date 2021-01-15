% script to read in database and make various categories...

[fname,datetime,frames,zdepth,tottime,stimbouts,locnum,...
    activecells,inverted,stimcode,quality,driver,...
    moving,layer,wavelength,flyID, responsivefly,Layercheck,SingleFly,PTX_5_100,Washout5mM,Washout100mM]...
    =textread('MASTER_foldersummary_Miriam_sh.txt',...
    '%s %s %f %f %f %f %f %s %f %s %f %s %f %s %f %f %f %f %f %f %f %f\r',...
    'headerlines',1,'delimiter','\t');%MS: %s reads a white space or delimiter separated string, %f a floating-point value
% reads now string for layers


for ll=1:length(fname)
    S=strfind(fname{ll,1},'SIMA');
    if ~isempty(S)
    i_PDATA_SIMA(ll)=1; 
    else 
    i_PDATA_SIMA(ll)=0;  
    end 
end 
i_PDATA_SIMA=i_PDATA_SIMA'; 


% Stimulus
i_mb_luminc_20degs_5deg_8dir = strcmpi('DriftingStripe_LumInc_1period_20degPerSec_10deg_BlankRand_8Dir',(stimcode));
i_mb_lumdec_20degs_5deg_8dir = strcmpi('DriftingStripe_LumDec_1period_20degPerSec_10deg_BlankRand_8Dir',(stimcode));
i_mb_movEdgesONOFF = strcmpi('DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs',(stimcode));


i_mb_luminc_20degs_5deg_8dir_FullContrast = strcmpi('DriftingStripe_LumInc_1period_20degPerSec_10deg_BlankRand_8Dir_fullContrast',(stimcode));
i_mb_lumdec_20degs_5deg_8dir_FullContrast = strcmpi('DriftingStripe_LumDec_1period_20degPerSec_10deg_BlankRand_8Dir_fullContrast',(stimcode));
i_mb_movEdgesONOFF_optimized = strcmpi('DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs_optimized',(stimcode));
i_mb_movEdgesONOFF_optimized_8Dir = strcmpi('DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_8Dirs_optimized',(stimcode));
i_mb_movEdgesONOFF_optimized_16Dir = strcmpi('DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_16Dirs_optimized',(stimcode));

i_mb_DriftingSquare_static=strcmpi('DriftingSquare_3secMoving_1secStanding_30degPerSec_30degPerCycle_1cyclePerSec_BlankRand_4Dir',(stimcode));
i_fff5s  = strcmpi('LocalCircle_5sec_220deg_0degAz_0degEl_Sequential_LumDec_LumInc',(stimcode)); 

% Quality of recording
i_moving = (moving == 0);
i_responsivefly = (responsivefly == 1);  

% Driver Line / Genotype
i_T4T5Rec_Control =strcmpi('w+_R59E08lexA-homo>>GCaMP6f',driver);
i_T4T5Rec_Kir_Control = strcmpi('R59E08lexA>>GCaMP6f_UASKirControl',driver);
i_T4T5Rec_LayerControl =strcmpi('w+_R59E08lexA-homo>>GCaMP6f_LayerControl',driver);

