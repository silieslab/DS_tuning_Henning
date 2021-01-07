% explore_CA_DATA_DriftingStripes_tuning width 


close all
clear all
clc
% 

addpath('/Users/mhennin2/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/ClusterAnalysisData')

Control_Dark=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/LayerControl/Dark_Stripes/processed_Data_SIMA_CS5.mat');
Control_Bright=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/LayerControl/Bright_Stripes/processed_Data_SIMA_CS5.mat');


CondName= 'LayerControl';


load('/Users/mhennin2/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/MyColormap.mat')


%Find same Flies, because for some flies i didnt get Off clusters, so I
%have less Flies in the Control_Dark than in the Control_Bright matrix,
%but plot the clusters for the same fly together I first need to
%identify which Data belong to the same Fly!! 

    
    FlynameBright={};
    for ii=1:size(Control_Bright.T4T5_mb,2)
    FlynameBright{ii}=Control_Bright.T4T5_mb(ii).Flyname(1:13);
    end 
    
    FlynameDark={};
    for ii=1:size(Control_Dark.T4T5_mb,2)
    FlynameDark{ii}=Control_Dark.T4T5_mb(ii).Flyname(1:13);
    end 

    HIT=[]; %HIT contains the indice of Dark(first row) and the belonging indice of Bright(2nd row)
    count=1;
    for ii=1:length(FlynameDark)   
     COMP=strcmp(FlynameDark{ii},FlynameBright);
        if sum(COMP)==1
            HIT(1,count)=ii;
            HIT(2,count)=find(COMP);
        count=count+1;
        end
    end
    
        
    
%Plot Flies for Control 

for i=1:length(HIT)
    IFly=Control_Bright.T4T5_mb(i);
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    
    F1=figure
    subplot(1,2,1)
    
%     imagesc(IFly.AV)
    
    
    ColROI_T4=nan(size(IFly.Masks.T4B{1,1}));
    for ii=1:length(IFly.Masks.T4A)
        ColROI_T4(find(IFly.Masks.T4A{1,ii}))= angle(IFly.Z.T4A(ii)); 
    end
    
    for ii=1:length(IFly.Masks.T4B)
        ColROI_T4(find(IFly.Masks.T4B{1,ii}))= angle(IFly.Z.T4B(ii)); 
    end
    
    for ii=1:length(IFly.Masks.T4C)
        ColROI_T4(find(IFly.Masks.T4C{1,ii}))= angle(IFly.Z.T4C(ii)); 
    end
    
    for ii=1:length(IFly.Masks.T4D)
        ColROI_T4(find(IFly.Masks.T4D{1,ii}))= angle(IFly.Z.T4D(ii)); 
    end
    ng=find(ColROI_T4<=0);
    ColROI_T4(ng)=(ColROI_T4(ng))+2*pi;
    ColROI_T4=ColROI_T4*180/pi;
%     ColROI=cat(3,ones(size(ColROI)),ones(size(ColROI)),ColROI);
    
h = imagesc(ColROI_T4);
    
load('/Users/mhennin2/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/MyColormap.mat')
colormap(cmap)
title(IFly.Flyname(1:11))


IFly2=Control_Dark.T4T5_mb(i);
    
  
    subplot(1,2,2)
    
%     imagesc(IFly.AV)
    
    
    ColROI_T5=nan(size(IFly2.Masks.T4B{1,1}));
    for ii=1:length(IFly2.Masks.T5A)
        ColROI_T5(find(IFly2.Masks.T5A{1,ii}))= angle(IFly2.Z.T5A(ii)); 
    end
    
    for ii=1:length(IFly2.Masks.T5B)
        ColROI_T5(find(IFly2.Masks.T5B{1,ii}))= angle(IFly2.Z.T5B(ii)); 
    end
    
    for ii=1:length(IFly2.Masks.T5C)
        ColROI_T5(find(IFly2.Masks.T5C{1,ii}))= angle(IFly2.Z.T5C(ii)); 
    end
    
    for ii=1:length(IFly2.Masks.T5D)
        ColROI_T5(find(IFly2.Masks.T5D{1,ii}))= angle(IFly2.Z.T5D(ii)); 
    end
    ng=find(ColROI_T5<=0);
    ColROI_T5(ng)=(ColROI_T5(ng))+2*pi;
    ColROI_T5=ColROI_T5*180/pi;
%     ColROI=cat(3,ones(size(ColROI)),ones(size(ColROI)),ColROI);
    
 h = imagesc(ColROI_T5);
    
colormap(cmap)
colorbar

% saveas(F1, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorangle_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'_sep.pdf'])    

 
    
    
    try
        ColROI=nan(size(IFly.Masks.T4A{1,1}));
        IntROI=nan(size(IFly.Masks.T4A{1,1}));
    catch
        ColROI=nan(size(IFly.Masks.T4B{1,1}));
        IntROI=nan(size(IFly.Masks.T4B{1,1}));
    end 
    
    
    
    for ii=1:length(IFly.Masks.T4A)
        ColROI(find(IFly.Masks.T4A{1,ii}))= angle(IFly.Z.T4A(ii)); 
        IntROI(find(IFly.Masks.T4A{1,ii}))= abs(IFly.Z.T4A(ii)); 
    end

    for ii=1:length(IFly.Masks.T4B)
        ColROI(find(IFly.Masks.T4B{1,ii}))= angle(IFly.Z.T4B(ii)); 
        IntROI(find(IFly.Masks.T4B{1,ii}))= abs(IFly.Z.T4B(ii)); 
    end
    
    for ii=1:length(IFly.Masks.T4C)
        ColROI(find(IFly.Masks.T4C{1,ii}))= angle(IFly.Z.T4C(ii)); 
        IntROI(find(IFly.Masks.T4C{1,ii}))= abs(IFly.Z.T4C(ii)); 
    end
    
    for ii=1:length(IFly.Masks.T4D)
        ColROI(find(IFly.Masks.T4D{1,ii}))= angle(IFly.Z.T4D(ii)); 
        IntROI(find(IFly.Masks.T4D{1,ii}))= abs(IFly.Z.T4D(ii));
    end

 
    
    for ii=1:length(IFly2.Masks.T5A)
        ColROI(find(IFly2.Masks.T5A{1,ii}))= angle(IFly2.Z.T5A(ii)); 
        IntROI(find(IFly2.Masks.T5A{1,ii}))= abs(IFly2.Z.T5A(ii)); 
    end
    
    for ii=1:length(IFly2.Masks.T5B)
        ColROI(find(IFly2.Masks.T5B{1,ii}))= angle(IFly2.Z.T5B(ii)); 
        IntROI(find(IFly2.Masks.T5B{1,ii}))= abs(IFly2.Z.T5B(ii)); 
    end
    
    for ii=1:length(IFly2.Masks.T5C)
        ColROI(find(IFly2.Masks.T5C{1,ii}))= angle(IFly2.Z.T5C(ii)); 
        IntROI(find(IFly2.Masks.T5C{1,ii}))= abs(IFly2.Z.T5C(ii));
    end
    
    for ii=1:length(IFly2.Masks.T5D)
        ColROI(find(IFly2.Masks.T5D{1,ii}))= angle(IFly2.Z.T5D(ii)); 
        IntROI(find(IFly2.Masks.T5D{1,ii}))= abs(IFly2.Z.T5D(ii)); 
    end
    ng=find(ColROI<=0);
    ColROI(ng)=(ColROI(ng))+2*pi;
    ColROI=ColROI*180/pi;
    
    
F1=figure    
h = imagesc(ColROI);
colorbar
colormap(cmap)
title(['Vector direction: ', IFly.Flyname(1:11)])
c=colorbar;
% c.Limits=[0 360];
caxis([0 360])

% F2=figure    
% h = imagesc(IntROI);
% 
% colormap(gray)
% title(['Vector length: ', IFly.Flyname(1:11)])
% c=colorbar;
% c.Limits=[0 1];
% 
%  % to plot this in Hue stile I need to add a third dimension
% % Colors need to be in a range between 0 and 1 
%  
% ColROI_Hue= cat(3,ColROI/360, ones(size(ColROI))*0.7 , IntROI);
% 
% ColROI_Hue=hsv2rgb(ColROI_Hue);
% 
% F3=figure;
% h = imagesc(ColROI_Hue);
% 
% CondName='Control';


% saveas(F2, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorlength_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% saveas(F3, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Hue_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% % 
 
 
end 

%% __________________________________________   
    
% Plot T4T5 distribution

for i=1:length(HIT)
%     IFly=Control_Bright.T4T5_mb(i);
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    
load('/Users/mhennin2/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/MyColormap.mat')

    
    try
        T4T5ROI=zeros(size(IFly.Masks.T4A{1,1}));
    catch
        T4T5ROI=zeros(size(IFly.Masks.T4B{1,1}));
    end 
    
    
    
    for ii=1:length(IFly.Masks.T4A)
        T4T5ROI(find(IFly.Masks.T4A{1,ii}))= 2; 
    end

    for ii=1:length(IFly.Masks.T4B)
        T4T5ROI(find(IFly.Masks.T4B{1,ii}))= 2; 
    end
    
    for ii=1:length(IFly.Masks.T4C)
        T4T5ROI(find(IFly.Masks.T4C{1,ii}))= 2; 
    end
    
    for ii=1:length(IFly.Masks.T4D)
        T4T5ROI(find(IFly.Masks.T4D{1,ii}))= 2; 
    end

%     IFly2=Control_Dark.T4T5_mb(i);   
    
    for ii=1:length(IFly2.Masks.T5A)
        T4T5ROI(find(IFly2.Masks.T5A{1,ii}))= 1; 
    end
    
    for ii=1:length(IFly2.Masks.T5B)
        T4T5ROI(find(IFly2.Masks.T5B{1,ii}))= 1; 
    end
    
    for ii=1:length(IFly2.Masks.T5C)
        T4T5ROI(find(IFly2.Masks.T5C{1,ii}))= 1; 
    end
    
    for ii=1:length(IFly2.Masks.T5D)
        T4T5ROI(find(IFly2.Masks.T5D{1,ii}))= 1; 
    end
    
    
F4=figure    
h = imagesc(T4T5ROI);
colorbar
colormap(copper)
title(['T4 T5 distribution: ', IFly.Flyname(1:11)])
c=colorbar;


 
%  saveas(F4, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/', CondName, '/Single_FlyPlots/T4T5_Distribution_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
%   saveas(F4, ['/Users/mhennin2/Desktop/Henning et al/Figures_Tuning Paper/Responses_to_Stripes/', CondName, '/Single_FlyPlots/T4T5_Distribution_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    

 
end 








