% script to read in database and make various categories...

[fname,datetime,frames,zdepth,tottime,stimbouts,locnum,...
    activecells,inverted,stimcode,quality,driver,...
    moving,layer,wavelength,flyID, responsivefly,Layercheck,SingleFly,PTX_5_100,Washout5mM,Washout100mM]...
    =textread('MASTER_foldersummary_Miriam.txt',...
    '%s %s %f %f %f %f %f %s %f %s %f %s %f %s %f %f %f %f %f %f %f %f\r',...
    'headerlines',1,'delimiter','\t');%MS: %s reads a white space or delimiter separated string, %f a floating-point value
% reads now string for layers

%i_fff2s  = strcmpi('LocalCircle_2sec_220deg_0degAz_0degEl_Sequential_LumDec_LumInc',(stimcode)); 
%i_fff5s  = strcmpi('LocalCircle_5sec_220deg_0degAz_0degEl_Sequential_LumDec_LumInc',(stimcode)); 


for ll=1:length(fname)
    S=strfind(fname{ll,1},'SIMA');
    if ~isempty(S)
    i_PDATA_SIMA(ll)=1; 
    else 
    i_PDATA_SIMA(ll)=0;  
    end 
end 
i_PDATA_SIMA=i_PDATA_SIMA'; 

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


%i_15bouts = stimbouts>=15;
%i_inverted = (inverted == 1);
i_moving = (moving == 0);
%i_active = (quality ~= 0);
i_responsivefly = (responsivefly == 1);  %ms addition

%changed to read strings
% i_L4M2 = strcmpi('L4M2',layer);
% i_L4M4 = strcmpi('L4M4',layer);
% i_lob1 = strcmpi('lob1',layer);


% i_482 = (wavelength == 482);
iT4C= strcmpi('T4C>>GCaMP6f',driver);
iT5= strcmpi('T5>>GCaMP6f',driver);

i_UASGCaMP6m_T4T5GMR79D04Gal4 = strcmpi('UASGCaMP6m_T4T5GMR79D04Gal4',driver);
i_UASGCaMP6m_T4T5GMR59E08lexA = strcmpi('UASGCaMP6m_T4T5GMR59E08lexA',driver);

i_T4T5Rec_shi_Control_1hincub = strcmpi('R59E08lexA>>GCaMP6f_UASshiControl_1h_37degWB',driver);
i_T4T5Rec_shi_C2_1hincub = strcmpi('R59E08lexA>>GCaMP6f_UASshiC2_1h_37degWB',driver);
i_T4T5Rec_shi_C3_1hincub = strcmpi('R59E08lexA>>GCaMP6f_UASshiC3_1h_37degWB',driver);
i_T4T5Rec_shi_Control_noincub = strcmpi('R59E08lexA>>GCaMP6f_UASshiControl_noincub',driver);
i_T4T5Rec_shi_C2_noincub = strcmpi('R59E08lexA>>GCaMP6f_UASshiC2_noincub',driver);
i_T4T5Rec_shi_C3_noincub = strcmpi('R59E08lexA>>GCaMP6f_UASshiC3_noincub',driver);

i_T4T5Rec_Kir_Control = strcmpi('R59E08lexA>>GCaMP6f_UASKirControl',driver);

i_T4T5Rec_Kir_C2 = strcmpi('R59E08lexA>>GCaMP6f_UASKirC2',driver);
i_T4T5Rec_Kir_C3 = strcmpi('R59E08lexA>>GCaMP6f_UASKirC3',driver);
i_T4T5Rec_Kir_C2C3 = strcmpi('R59E08lexA>>GCaMP6f_UASKirC2C3',driver);
% i_T4T5Rec_Light =strcmpi('R59E08lexA-homo>>GCaMP6f_light',driver);
% i_T4T5Rec_Dark =strcmpi('R59E08lexA-homo>>GCaMP6f_dark',driver);
 i_T4T5Rec_Dark =strcmpi('w+_R59E08lexA-homo>>GCaMP6f_dark',driver);
i_T4T5Rec_Adapt =strcmpi('R59E08lexA-homo>>GCaMP6f_light_adaptation',driver);
i_T4T5Rec_LayerControl =strcmpi('w+_R59E08lexA-homo>>GCaMP6f_LayerControl',driver);
i_T4T5Rec_Control =strcmpi('w+_R59E08lexA-homo>>GCaMP6f',driver);

i_T4T5Rec_PosFront =strcmpi('w+_R59E08lexA-homo>>GCaMP6f_Position_front',driver);
i_T4T5Rec_PosSide =strcmpi('w+_R59E08lexA-homo>>GCaMP6f_Position_side',driver);
i_T4T5Rec_PosFlat =strcmpi('w+_R59E08lexA-homo>>GCaMP6f_Position_flat',driver);


i_T4T5Rec_Kir_Control_RF = strcmpi('R59E08lexA>>GCaMP6f_UASKirControl_forRF',driver);
i_T4T5Rec_Kir_C2_RF = strcmpi('R59E08lexA>>GCaMP6f_UASKirC2_forRF',driver);
i_T4T5Rec_Kir_C3_RF = strcmpi('R59E08lexA>>GCaMP6f_UASKirC3_forRF',driver);


i_T4T5Reiser_shi_Control_1hincub = strcmpi('ReiserT4/T5_UASshiControl_1h_37degWB',driver);
% i_T4T5Reiser_shi_C2_1hincub = strcmpi('UASGCaMP6m_T4T5GMR59E08lexA',driver);
% i_T4T5Reiser_shi_C3_1hincub = strcmpi('UASGCaMP6m_T4T5GMR59E08lexA',driver);
i_T4T5Reiser_shi_Control_noincub = strcmpi('ReiserT4/T5_UASshiControl_noincub',driver);
% i_T4T5Reiser_shi_C2_noincub = strcmpi('UASGCaMP6m_T4T5GMR59E08lexA',driver);
% i_T4T5Reiser_shi_C3_noincub = strcmpi('UASGCaMP6m_T4T5GMR59E08lexA',driver);



i_T4T5_Rdl_Control_25_before = strcmpi('R59E08lexA>>GCaMP6f_Rdl_Control_2-5uMPTX_before',driver);
i_T4T5_Rdl_Control_25_3min = strcmpi('R59E08lexA>>GCaMP6f_Rdl_Control_2-5uMPTX_3min',driver);
i_T4T5_Rdl_Control_25_10min = strcmpi('R59E08lexA>>GCaMP6f_Rdl_Control_2-5uMPTX_10min',driver);

i_T4T5_Rdl_Control_noPTX = strcmpi('R59E08lexA>>GCaMP6f_Rdl_Control_noPTX',driver);
i_T4T5_Rdl_Control_1uMPTX = strcmpi('R59E08lexA>>GCaMP6f_Rdl_Control_1uMPTX',driver);

i_T4T5_RdlMDRR_noPTX = strcmpi('R59E08lexA>>GCaMP6f_RdlMDRR_noPTX',driver);
i_T4T5_RdlMDRR_1uMPTX = strcmpi('R59E08lexA>>GCaMP6f_RdlMDRR_1uMPTX',driver);

i_T4T5_GluCl_Control_noPTX = strcmpi('R59E08lexA>>GCaMP6f_GluCl_Control_noPTX',driver);
i_T4T5_GluCl_Control_1uMPTX = strcmpi('R59E08lexA>>GCaMP6f_GluCl_Control_1uMPTX',driver);

i_T4T5_GluCl_S278T_noPTX = strcmpi('R59E08lexA>>GCaMP6f_GluCl_S278T_noPTX',driver);
i_T4T5_GluCl_S278T_1uMPTX = strcmpi('R59E08lexA>>GCaMP6f_GluCl_S278T_1uMPTX',driver);



i_nframes = (frames == 3600);