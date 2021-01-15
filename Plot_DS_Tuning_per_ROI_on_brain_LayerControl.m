% explore_CA_DATA_DriftingStripes_tuning width 

%%%%%%%%%%%%%%
% This skript plots all figures shown in Supplementary Figure4
%%%%%%%%%%%%%%

addpath(genpath('subscripts'))
addpath(genpath('Data/Data_Edges'))

% addpath('/Users/mhenning/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/ClusterAnalysisData')

% Load preprocessed Data matrix

Control_Dark=load('Data/Data_Stripes/Dark_Stripes/processed_Data_SIMA_CS5_Layer.mat');
Control_Bright=load('Data/Data_Stripes/Bright_Stripes/processed_Data_SIMA_CS5_Layer.mat');

CondName= 'Control';


load('MyColormap.mat')
%Fly1
Plot=[3:7]; % to plot the examples shown in Supplementary Fig. 4
Z_depth=[15,30,45,60,75]; 

% %Fly2
% Plot=[24:27]; % to plot the examples shown in Supplementary Fig. 4
% Z_depth=[30,45,60,75];

% %Fly3
% Plot=[13:16]; % to plot the examples shown in Supplementary Fig. 4
% Z_depth=[15,30,45,60];



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
countPlot=0;
for i=1:length(HIT)
    IFly=Control_Bright.T4T5_mb(i);
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
  
    
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
    
    
    if sum(i==Plot) 
        countPlot=countPlot+1; 
        F1=figure    
        h = imagesc(ColROI);
        colorbar
        colormap(cmap)
        title(['Vector direction: ', IFly.Flyname(1:11), '  Z-depth: ', num2str(Z_depth(countPlot))])
        c=colorbar;
        % c.Limits=[0 360];
        caxis([0 360])
    end 
 
 
end 







