% explore Ca --> Drifting Stripes 

%--> plots single Fly Plots 
%--> plots average responses across flies 
%--> saves Data on server in subfolders 

% clc;
clear all;
close all;
Foldertosave='/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning' ;  % Results_T4T5_Imaging_DS_tuning
addpath(genpath('/Users/mhenning/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis'));

% CONSTANTS
savein= 'Control';
Condition='i_T4T5Rec_Control.*i_mb_movEdgesONOFF_optimized_8Dir.*i_responsivefly.*~moving.*i_PDATA_SIMA';


pdatapath='/Volumes/Seagate/Documents/2p-imaging/Miriam_pData';
cd(pdatapath);

database_select_samples_mh; % get the right stuff...
f_mb_T4T5 = find(eval(Condition));


T4T5_mb = load_neuron_data10Hz_ks_CA_average_Edges_8Dir_mh(f_mb_T4T5,pdatapath,'ClusterInfo_ManuallySelect');

    
 %%  
    
    %***********************************************
    %  PLOTS - by FLY
    %***********************************************    

            COLOR_OF_PLOT_A = [0 .5 0]; % GREEN
            COLOR_CLOUD_A =[.5 .7 .5];  % GREEN cloud

            COLOR_OF_PLOT_B = [0 0 1];  %BLUE;
            COLOR_CLOUD_B = [.5 .5 1];  %BLUE cloud

            COLOR_OF_PLOT_C = [1 0 0];  %RED
            COLOR_CLOUD_C = [1 .5 .5];  %RED cloud
    
            COLOR_OF_PLOT_D = [.7 .7 0];%YELLOW
            COLOR_CLOUD_D = [.9 .9 .5]; %YELLOW cloud
            
    TRACE_OFFSET_PLOT = 4; % spacing in dF/F between the traces for the 4 directions
    t=[0:0.1:4.9];
    
    ALL_Flies.T5A=[];
    ALL_Flies.T5B=[];
    ALL_Flies.T5C=[];
    ALL_Flies.T5D=[];
    
    ALL_Flies.T4A=[];
    ALL_Flies.T4B=[];
    ALL_Flies.T4C=[];
    ALL_Flies.T4D=[];

   %% \\ T5 
 F1=  figure(1); hold on;

   subplot(2,2,1) %Layer 1/A    
   hold on  
   for iFly=1:length(T4T5_mb) 
       try 
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5A(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5A(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_A); 
       end 
             
       ALL_Flies.T5A=cat(3,ALL_Flies.T5A,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5A);
       catch 
       end 
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5A)==0;
          ALL_Flies.T5A=cat(3,ALL_Flies.T5A,nan(16,80));
       end 
       
   end 
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
   title(['T5  NFly=', num2str(length(T4T5_mb)), '\newline LayerA '])
   
   
   subplot(2,2,2) %Layer 2/B
   hold on
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5B(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5B(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_B); 
       end 
       ALL_Flies.T5B=cat(3,ALL_Flies.T5B,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5B);

       catch 
       end 
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5B)==0;
          ALL_Flies.T5B=cat(3,ALL_Flies.T5B,nan(16,80));
       end
   end 
   
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
   title(['LayerB '])

   
  subplot(2,2,3) %Layer 3/C
  hold on  
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5C(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5C(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_C); 
       end   
       ALL_Flies.T5C=cat(3,ALL_Flies.T5C,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5C);

       catch 
       end 
        if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5C)==0;
          ALL_Flies.T5C=cat(3,ALL_Flies.T5C,nan(16,80));
       end
   end 
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
   title(['LayerC '])

   
   subplot(2,2,4) %Layer 4/D
   hold on  
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5D(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5D(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_D); 
       end 
       ALL_Flies.T5D=cat(3,ALL_Flies.T5D,T4T5_mb(iFly).iAV_ROI_Resp.iM_T5D);

       catch 
       end 
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T5D)==0;
          ALL_Flies.T5D=cat(3,ALL_Flies.T5D,nan(16,80));
       end
   end 
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
      title(['LayerD '])

   
   
 %% \\ T4 
 F2=   figure(2); hold on;

   subplot(2,2,1) %Layer 1/A
   hold on  
   
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4A(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4A(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_A); 
       end 
       
       ALL_Flies.T4A=cat(3,ALL_Flies.T4A,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4A);

       catch 
       end 
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4A)==0;
          ALL_Flies.T4A=cat(3,ALL_Flies.T4A,nan(16,80));
       end
   end 
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
   title(['T4  NFly=', num2str(length(T4T5_mb)), '\newline LayerA'])
   
   
   subplot(2,2,2) %Layer 1/A
   hold on 
  
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4B(1,:))-1)/10]; 
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4B(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_B); 
       end          
       ALL_Flies.T4B=cat(3,ALL_Flies.T4B,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4B);
       catch
       end
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4B)==0;
          ALL_Flies.T4B=cat(3,ALL_Flies.T4B,nan(16,80));
       end
   end 
   
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
   title(['LayerB '])

   
   subplot(2,2,3) %Layer 3/C
   hold on  
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4C(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4C(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_C); 
       end 
       ALL_Flies.T4C=cat(3,ALL_Flies.T4C,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4C);
       catch 
       end 
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4C)==0;
          ALL_Flies.T4C=cat(3,ALL_Flies.T4C,nan(16,80));
       end
   end 
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
   title(['LayerC '])

   subplot(2,2,4) %Layer 4/D
     hold on  
   for iFly=1:length(T4T5_mb) 
       try
       t=[0:0.1:(length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4D(1,:))-1)/10];
       for NDir=1:16;
            plot(t,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4D(NDir,:)+(NDir-1)*TRACE_OFFSET_PLOT,'Color',COLOR_OF_PLOT_D); 
       end 
       ALL_Flies.T4D=cat(3,ALL_Flies.T4D,T4T5_mb(iFly).iAV_ROI_Resp.iM_T4D);
       catch 
       end 
       if length(T4T5_mb(iFly).iAV_ROI_Resp.iM_T4D)==0;
          ALL_Flies.T4D=cat(3,ALL_Flies.T4D,nan(16,80));
       end
   end 
   set(gca, 'YLim', [-0.5 65])
   set(gca, 'XLim', [0 8])
      title(['LayerD '])

   
      
 %%     
  
    
    %***********************************************
    %  PLOTS - average across FLY
    %***********************************************        
    
    % \\\\\ T5 

    F3= figure(3);
     subplot(2,2,1)
   
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T5A(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T5A(i,:,:))')'),COLOR_OF_PLOT_A,COLOR_CLOUD_A);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['T5  LayerA - NFlies = ', num2str(size(ALL_Flies.T5A(i,:,:),3))])
     
            
    subplot(2,2,2)
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T5B(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T5B(i,:,:))')'),COLOR_OF_PLOT_B,COLOR_CLOUD_B);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['LayerB - NFlies = ', num2str(size(ALL_Flies.T5B(i,:,:),3))])
     
   
    subplot(2,2,3)
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T5C(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T5C(i,:,:))')'),COLOR_OF_PLOT_C,COLOR_CLOUD_C);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['LayerC - NFlies = ', num2str(size(ALL_Flies.T5C(i,:,:),3))])


       subplot(2,2,4)
   for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T5D(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T5D(i,:,:))')'),COLOR_OF_PLOT_D,COLOR_CLOUD_D);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['LayerD - NFlies = ', num2str(size(ALL_Flies.T5D(i,:,:),3))])

   
   
   
   
       % \\\\\ T4 

   F4=  figure(4);
   subplot(2,2,1)
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T4A(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T4A(i,:,:))')'),COLOR_OF_PLOT_A,COLOR_CLOUD_A);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['T4  LayerA - NFlies = ', num2str(size(ALL_Flies.T4A(i,:,:),3))])
     
            
    subplot(2,2,2)
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T4B(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T4B(i,:,:))')'),COLOR_OF_PLOT_B,COLOR_CLOUD_B);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['LayerB - NFlies = ', num2str(size(ALL_Flies.T4B(i,:,:),3))])
     
   
    subplot(2,2,3)
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T4C(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T4C(i,:,:))')'),COLOR_OF_PLOT_C,COLOR_CLOUD_C);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['LayerC - NFlies = ', num2str(size(ALL_Flies.T4C(i,:,:),3))])


       subplot(2,2,4)
    for i=1:16;
       plot_err_patch_v2(t, TRACE_OFFSET_PLOT*(i-1)+(nanmean(ALL_Flies.T4D(i,:,:),3)'),...
          ( nanstd(squeeze(ALL_Flies.T4D(i,:,:))')'),COLOR_OF_PLOT_D,COLOR_CLOUD_D);
       line([t(1) t(end)],[TRACE_OFFSET_PLOT*(i-1)  TRACE_OFFSET_PLOT*(i-1)],...
       'color',[0 0 0],'LineStyle','--');
%        text(0.5+mod_t(range(1)),0.25 + TRACE_OFFSET_PLOT*spacer_counter,num2str(dirs(idirs)));
    end         
   set(gca, 'YLim', [-0.5 63])
   set(gca, 'XLim', [0 8])
   title(['LayerD - NFlies = ', num2str(size(ALL_Flies.T4D(i,:,:),3))])

   
% /size(ALL_Flies.T4A(i,:,:),3)
Folders=dir([Foldertosave, '/Responses_to_Edges16Dir/', savein]);
if length(Folders)<1
    mkdir([Foldertosave, '/Responses_to_Edges16Dir/',savein]);
end 


saveas(F1, [Foldertosave, '/Responses_to_Edges16Dir/',savein, '/T5 responses to Edges16Dir_SIMA_CS5.pdf'])    
saveas(F2, [Foldertosave, '/Responses_to_Edges16Dir/',savein,  '/T4 responses to Edges16Dir_SIMA_CS5.pdf'])    
saveas(F3, [Foldertosave, '/Responses_to_Edges16Dir/',savein,  '/Average T5 responses to Edges16Dir_SIMA_CS5.pdf'])    
saveas(F4, [Foldertosave, '/Responses_to_Edges16Dir/',savein,  '/Average T4 responses to Edges16Dir_SIMA_CS5.pdf'])    
    

save([Foldertosave, '/Responses_to_Edges16Dir/',savein,  '/processed_Data_SIMA_CS5'],'ALL_Flies','T4T5_mb'); 
