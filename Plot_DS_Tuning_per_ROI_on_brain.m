
%%%%%%%%%%%%%%
% This skript plots all figures shown in Figure3 
%%%%%%%%%%%%%% 

addpath(genpath('subscripts'))
addpath(genpath('Data/Data_Edges'))

% addpath('/Users/mhenning/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/ClusterAnalysisData')

% Load preprocessed Data matrix

Control_Dark=load('Data/Data_Stripes/Dark_Stripes/processed_Data_SIMA_CS5_sh.mat');
Control_Bright=load('Data/Data_Stripes/Bright_Stripes/processed_Data_SIMA_CS5_sh.mat');

CondName= 'Control';


load('MyColormap.mat')

      

    %Find same Flies, because for some flies i didnt get Off clusters, so I
    %have less Flies in the Control_Dark than in the Control_Bright matrix,
    %but plot the clusters for the same fly together I first need to
    %identify which Data belong to the same Fly!! 
    
%     FlynameBright={};
%     for ii=1:size(Control_Bright.T4T5_mb,2)
%     FlynameBright{ii}=Control_Bright.T4T5_mb(ii).Flyname(1:11);
%     end 
%     
%     FlynameDark={};
%     for ii=1:size(Control_Dark.T4T5_mb,2)
%     FlynameDark{ii}=Control_Dark.T4T5_mb(ii).Flyname(1:11);
%     end 
    
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
    for ii=1:length(FlynameDark);    
     COMP=strcmp(FlynameDark{ii},FlynameBright);
        if sum(COMP)==1
            HIT(1,count)=ii;
            HIT(2,count)=find(COMP);
        count=count+1;
        end
    end
    
    
    
    
%Plot Flies for Control 

for i=4%1:length(HIT)
    IFly=Control_Bright.T4T5_mb(i);
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    
    
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



IFly2=Control_Dark.T4T5_mb(i);
    
    
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

load('/Users/mhenning/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/MyColormap.mat')




F0=figure

subplot(1,2,1)
% imagesc(IFly.AV)
h = imagesc(ColROI_T4);
colormap(cmap)
caxis([0 360])
title(IFly.Flyname(1:11))

subplot(1,2,2)
% imagesc(IFly.AV)
h = imagesc(ColROI_T5);
colormap(cmap)
% colormap(M)
colorbar
caxis([0 360])


% saveas(F0, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorangle_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'_sep.pdf'])    

   
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

F2=figure    
h = imagesc(IntROI);
colormap(gray)
title(['Vector length: ', IFly.Flyname(1:11)])
c=colorbar;
c.Limits=[0 1];

 % to plot this in Hue stile I need to add a third dimension
% Colors need to be in a range between 0 and 1 
 
ColROI_Hue= cat(3,ColROI/360, ones(size(ColROI))*0.7 , IntROI);

ColROI_Hue=hsv2rgb(ColROI_Hue);

F3=figure;
h = imagesc(ColROI_Hue);

CondName='Control';

%  saveas(F1, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorangle_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
%  saveas(F2, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorlength_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
%  saveas(F3, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Hue_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% saveas(F1, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorangle_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% saveas(F1, ['/Users/mhennin2/Desktop/Henning et al/Figures_Tuning Paper/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorangle_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% saveas(F2, ['/Users/mhennin2/Desktop/Henning et al/Figures_Tuning Paper/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorlength_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% saveas(F3, ['/Users/mhennin2/Desktop/Henning et al/Figures_Tuning Paper/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Hue_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    


% saveas(F2, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Vectorlength_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% saveas(F3, ['/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/', CondName, '/Single_FlyPlots/DS_Distribution_Hue_',IFly.Flyname(1:6), '_',IFly.Flyname(8:13),'.pdf'])    
% % 
 
 
end 

%% __________________________________________
% Control_Dark=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/Control/Dark_Stripes/processed_Data_SIMA_CS5.mat');
% Control_Dark=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/LayerControl/Dark_Stripes/processed_Data_SIMA_CS5.mat');
% Control_Dark=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/C2Silence/Dark_Stripes/processed_Data_SIMA_CS5.mat');
% Control_Dark=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/C3Silence/Dark_Stripes/processed_Data_SIMA_CS5.mat');

% Control_Bright=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/Control/Bright_Stripes/processed_Data_SIMA_CS5.mat');
% Control_Bright=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/LayerControl/Bright_Stripes/processed_Data_SIMA_CS5.mat');
% Control_Bright=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/C2Silence/Bright_Stripes/processed_Data_SIMA_CS5');
% Control_Bright=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4_T5_Imaging_C2_C3Kir/Responses_to_Stripes/C3Silence/Bright_Stripes/processed_Data_SIMA_CS5');


% Control_Bright=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/PositionControl/Pos_Flat/Bright_Stripes/processed_Data_SIMA_CS5.mat');
% Control_Dark=load('/Volumes/ukme06/mhennin2/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Stripes/PositionControl/Pos_Flat/Dark_Stripes/processed_Data_SIMA_CS5.mat');
% 
% CondName= 'PositionControl/Pos_Flat';

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
    for ii=1:length(FlynameDark);    
     COMP=strcmp(FlynameDark{ii},FlynameBright);
        if sum(COMP)==1
            HIT(1,count)=ii;
            HIT(2,count)=find(COMP);
        count=count+1;
        end
    end
    
    
    
    
% Plot Flies for Control 

for i=1:length(HIT)
%     IFly=Control_Bright.T4T5_mb(i);
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    
load('/Users/mhenning/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/MyColormap.mat')

    
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




%% Plott tuning with respect to quantifyied position- proximal vs. distal

for i=4 %1:length(HIT) do this just for one example fly
    
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    Tuning_T4_A=[];
    X_axis_pos_T4_A=[];
    for ii=1:length(IFly.Masks.T4A)
            Tuning_T4_A=[Tuning_T4_A,angle(IFly.Z.T4A(ii))]; 
            [~,X]=find(IFly.Masks.T4A{1,ii});
            X_axis_pos_T4_A=[X_axis_pos_T4_A,mean(X)];
    end
    
 
   Tuning_T4_B=[];
    X_axis_pos_T4_B=[];
    for ii=1:length(IFly.Masks.T4B)
            Tuning_T4_B=[Tuning_T4_B,angle(IFly.Z.T4B(ii))]; 
            [~,X]=find(IFly.Masks.T4B{1,ii});
            X_axis_pos_T4_B=[X_axis_pos_T4_B,mean(X)];
    end
    
    
    Tuning_T4_C=[];
    X_axis_pos_T4_C=[];
    for ii=1:length(IFly.Masks.T4C)
            Tuning_T4_C=[Tuning_T4_C,angle(IFly.Z.T4C(ii))]; 
            [~,X]=find(IFly.Masks.T4C{1,ii});
            X_axis_pos_T4_C=[X_axis_pos_T4_C,mean(X)];
    end
    
 
   Tuning_T4_D=[];
    X_axis_pos_T4_D=[];
    for ii=1:length(IFly.Masks.T4D)
            Tuning_T4_D=[Tuning_T4_D,angle(IFly.Z.T4D(ii))]; 
            [~,X]=find(IFly.Masks.T4D{1,ii});
            X_axis_pos_T4_D=[X_axis_pos_T4_D,mean(X)];
    end
    
%%    
    Tuning_T5_A=[];
    X_axis_pos_T5_A=[];
    for ii=1:length(IFly2.Masks.T5A)
            Tuning_T5_A=[Tuning_T5_A,angle(IFly2.Z.T5A(ii))]; 
            [~,X]=find(IFly2.Masks.T5A{1,ii});
            X_axis_pos_T5_A=[X_axis_pos_T5_A,mean(X)];
    end
    
 
   Tuning_T5_B=[];
    X_axis_pos_T5_B=[];
    for ii=1:length(IFly2.Masks.T5B)
            Tuning_T5_B=[Tuning_T5_B,angle(IFly2.Z.T5B(ii))]; 
            [~,X]=find(IFly2.Masks.T5B{1,ii});
            X_axis_pos_T5_B=[X_axis_pos_T5_B,mean(X)];
    end
    
    
    Tuning_T5_C=[];
    X_axis_pos_T5_C=[];
    for ii=1:length(IFly2.Masks.T5C)
            Tuning_T5_C=[Tuning_T5_C,angle(IFly2.Z.T5C(ii))]; 
            [~,X]=find(IFly2.Masks.T5C{1,ii});
            X_axis_pos_T5_C=[X_axis_pos_T5_C,mean(X)];
    end
    
 
   Tuning_T5_D=[];
    X_axis_pos_T5_D=[];
    for ii=1:length(IFly2.Masks.T5D)
            Tuning_T5_D=[Tuning_T5_D,angle(IFly2.Z.T5D(ii))]; 
            [~,X]=find(IFly2.Masks.T5D{1,ii});
            X_axis_pos_T5_D=[X_axis_pos_T5_D,mean(X)];
    end

    
    Tuning_T4_A=convert_angle(Tuning_T4_A);
    Tuning_T4_B=convert_angle(Tuning_T4_B);
    Tuning_T4_C=convert_angle(Tuning_T4_C);
    Tuning_T4_D=convert_angle(Tuning_T4_D);
    
      
    Tuning_T5_A=convert_angle(Tuning_T5_A);
    Tuning_T5_B=convert_angle(Tuning_T5_B);
    Tuning_T5_C=convert_angle(Tuning_T5_C);
    Tuning_T5_D=convert_angle(Tuning_T5_D);


    
end 


figure 
subplot(2,2,1)
scatter(X_axis_pos_T4_A,Tuning_T4_A)
hold on 
scatter(X_axis_pos_T5_A,Tuning_T5_A)
set(gca, 'YLim', [0 60])
set(gca, 'XLim', [70 180])

title('LayerA')

subplot(2,2,2)
scatter(X_axis_pos_T4_B,Tuning_T4_B)
hold on 
scatter(X_axis_pos_T5_B,Tuning_T5_B)
set(gca, 'YLim', [180 240])
set(gca, 'XLim', [70 180])
title('LayerB')

subplot(2,2,3)
scatter(X_axis_pos_T4_C,Tuning_T4_C)
hold on 
scatter(X_axis_pos_T5_C,Tuning_T5_C)
set(gca, 'YLim', [80 140])
set(gca, 'XLim', [70 180])
title('LayerC')

subplot(2,2,4)
scatter(X_axis_pos_T4_D,Tuning_T4_D)
hold on 
scatter(X_axis_pos_T5_D,Tuning_T5_D)
set(gca, 'YLim', [250 310])
set(gca, 'XLim', [70 180])
title('LayerD')


figure 
scatter(X_axis_pos_T4_A,Tuning_T4_A,'g')
hold on 
scatter(X_axis_pos_T5_A,Tuning_T5_A,'g')

scatter(X_axis_pos_T4_B,Tuning_T4_B,'b')
hold on 
scatter(X_axis_pos_T4_B,Tuning_T4_B,'b')

scatter(X_axis_pos_T4_C,Tuning_T4_C,'r')
hold on 
scatter(X_axis_pos_T4_C,Tuning_T4_C,'r')

scatter(X_axis_pos_T4_D,Tuning_T4_D,'y')
hold on 
scatter(X_axis_pos_T4_D,Tuning_T4_D,'y')


% calculate change of deg 

%Layer A 
X_A=[X_axis_pos_T4_A,X_axis_pos_T5_A];
Tuning_A=[Tuning_T4_A,Tuning_T5_A];
[Min,min_pos]=min(X_A);
Prox_tun=Tuning_A(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE

Tuning_T4_A_lin=Tuning_T4_A;
Tuning_T4_A_lin(Tuning_T4_A<100)=Tuning_T4_A_lin(Tuning_T4_A<100)+360;
Change_Tun_T4_A=Tuning_T4_A_lin-Prox_tun;
Tuning_T5_A_lin=Tuning_T5_A;
Tuning_T5_A_lin(Tuning_T5_A<100)=Tuning_T5_A_lin(Tuning_T5_A<100)+360;
Change_Tun_T5_A=Tuning_T5_A_lin-Prox_tun;

%Layer B 
X_B=[X_axis_pos_T4_B,X_axis_pos_T5_B];
Tuning_B=[Tuning_T4_B,Tuning_T5_B];
[Min,min_pos]=min(X_B);
Prox_tun=Tuning_B(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_T4_B=Tuning_T4_B-Prox_tun;
Change_Tun_T5_B=Tuning_T5_B-Prox_tun;


%Layer C 
X_C=[X_axis_pos_T4_C,X_axis_pos_T5_C];
Tuning_C=[Tuning_T4_C,Tuning_T5_C];
[Min,min_pos]=min(X_C);
Prox_tun=Tuning_C(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_T4_C=Tuning_T4_C-Prox_tun;
Change_Tun_T5_C=Tuning_T5_C-Prox_tun;


%Layer D 
X_D=[X_axis_pos_T4_D,X_axis_pos_T5_D];
Tuning_D=[Tuning_T4_D,Tuning_T5_D];
[Min,min_pos]=min(X_D);
Prox_tun=Tuning_D(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_T4_D=Tuning_T4_D-Prox_tun;
Change_Tun_T5_D=Tuning_T5_D-Prox_tun;


figure 
scatter(X_axis_pos_T4_A,Change_Tun_T4_A,[],[0, 153/255, 0])
hold on 
scatter(X_axis_pos_T5_A,Change_Tun_T5_A,[],[0, 77/255, 0])

scatter(X_axis_pos_T4_B,Change_Tun_T4_B,[],[0, 0, 1])
scatter(X_axis_pos_T5_B,Change_Tun_T5_B,[],[0, 0, 179/255])

scatter(X_axis_pos_T4_C,Change_Tun_T4_C,[],[1, 0, 0])
scatter(X_axis_pos_T5_C,Change_Tun_T5_C,[],[179/255, 0, 0])

scatter(X_axis_pos_T4_D,Change_Tun_T4_D,[],[204/255, 204/255, 0])
scatter(X_axis_pos_T5_D,Change_Tun_T5_D,[],[153/255, 153/255, 0])



%Now do this again for all brains, 

Max_change_A=[];
Max_change_B=[];
Max_change_C=[];
Max_change_D=[];


for i=1:length(HIT) 
    
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    Tuning_T4_A=[];
    X_axis_pos_T4_A=[];
    for ii=1:length(IFly.Masks.T4A)
            Tuning_T4_A=[Tuning_T4_A,angle(IFly.Z.T4A(ii))]; 
            [~,X]=find(IFly.Masks.T4A{1,ii});
            X_axis_pos_T4_A=[X_axis_pos_T4_A,mean(X)];
    end
    
 
   Tuning_T4_B=[];
    X_axis_pos_T4_B=[];
    for ii=1:length(IFly.Masks.T4B)
            Tuning_T4_B=[Tuning_T4_B,angle(IFly.Z.T4B(ii))]; 
            [~,X]=find(IFly.Masks.T4B{1,ii});
            X_axis_pos_T4_B=[X_axis_pos_T4_B,mean(X)];
    end
    
    
    Tuning_T4_C=[];
    X_axis_pos_T4_C=[];
    for ii=1:length(IFly.Masks.T4C)
            Tuning_T4_C=[Tuning_T4_C,angle(IFly.Z.T4C(ii))]; 
            [~,X]=find(IFly.Masks.T4C{1,ii});
            X_axis_pos_T4_C=[X_axis_pos_T4_C,mean(X)];
    end
    
 
   Tuning_T4_D=[];
    X_axis_pos_T4_D=[];
    for ii=1:length(IFly.Masks.T4D)
            Tuning_T4_D=[Tuning_T4_D,angle(IFly.Z.T4D(ii))]; 
            [~,X]=find(IFly.Masks.T4D{1,ii});
            X_axis_pos_T4_D=[X_axis_pos_T4_D,mean(X)];
    end
    
%%    
    Tuning_T5_A=[];
    X_axis_pos_T5_A=[];
    for ii=1:length(IFly2.Masks.T5A)
            Tuning_T5_A=[Tuning_T5_A,angle(IFly2.Z.T5A(ii))]; 
            [~,X]=find(IFly2.Masks.T5A{1,ii});
            X_axis_pos_T5_A=[X_axis_pos_T5_A,mean(X)];
    end
    
 
   Tuning_T5_B=[];
    X_axis_pos_T5_B=[];
    for ii=1:length(IFly2.Masks.T5B)
            Tuning_T5_B=[Tuning_T5_B,angle(IFly2.Z.T5B(ii))]; 
            [~,X]=find(IFly2.Masks.T5B{1,ii});
            X_axis_pos_T5_B=[X_axis_pos_T5_B,mean(X)];
    end
    
    
    Tuning_T5_C=[];
    X_axis_pos_T5_C=[];
    for ii=1:length(IFly2.Masks.T5C)
            Tuning_T5_C=[Tuning_T5_C,angle(IFly2.Z.T5C(ii))]; 
            [~,X]=find(IFly2.Masks.T5C{1,ii});
            X_axis_pos_T5_C=[X_axis_pos_T5_C,mean(X)];
    end
    
 
   Tuning_T5_D=[];
    X_axis_pos_T5_D=[];
    for ii=1:length(IFly2.Masks.T5D)
            Tuning_T5_D=[Tuning_T5_D,angle(IFly2.Z.T5D(ii))]; 
            [~,X]=find(IFly2.Masks.T5D{1,ii});
            X_axis_pos_T5_D=[X_axis_pos_T5_D,mean(X)];
    end

    
    Tuning_T4_A=convert_angle(Tuning_T4_A);
    Tuning_T4_B=convert_angle(Tuning_T4_B);
    Tuning_T4_C=convert_angle(Tuning_T4_C);
    Tuning_T4_D=convert_angle(Tuning_T4_D);
    
      
    Tuning_T5_A=convert_angle(Tuning_T5_A);
    Tuning_T5_B=convert_angle(Tuning_T5_B);
    Tuning_T5_C=convert_angle(Tuning_T5_C);
    Tuning_T5_D=convert_angle(Tuning_T5_D);


    
    
    
    % calculate change of deg 

%Layer A 
X_A=[X_axis_pos_T4_A,X_axis_pos_T5_A];
Tuning_A=[Tuning_T4_A,Tuning_T5_A];
[Min,min_pos]=min(X_A);
[Max,max_pos]=max(X_A);

Tuning_A_lin=Tuning_A;
Tuning_A_lin(Tuning_A<100)=Tuning_A_lin(Tuning_A<100)+360;

Prox_tun=Tuning_A_lin(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_A=Tuning_A_lin-Prox_tun;

Max_change_A=[Max_change_A,abs(Change_Tun_A(max_pos))];   % TUNING ANGLE AT MOST DISTAL SIDE



%Layer B 
X_B=[X_axis_pos_T4_B,X_axis_pos_T5_B];
Tuning_B=[Tuning_T4_B,Tuning_T5_B];
[Min,min_pos]=min(X_B);
[Max,max_pos]=max(X_B);
Prox_tun=Tuning_B(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_B=Tuning_B-Prox_tun;
Max_change_B=[Max_change_B,abs(Change_Tun_B(max_pos))];   % TUNING ANGLE AT MOST DISTAL SIDE



%Layer C 
X_C=[X_axis_pos_T4_C,X_axis_pos_T5_C];
Tuning_C=[Tuning_T4_C,Tuning_T5_C];
[Min,min_pos]=min(X_C);
[Max,max_pos]=max(X_C);
Prox_tun=Tuning_C(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_C=Tuning_C-Prox_tun;
Max_change_C=[Max_change_C,abs(Change_Tun_C(max_pos))];   % TUNING ANGLE AT MOST DISTAL SIDE


%Layer D 
X_D=[X_axis_pos_T4_D,X_axis_pos_T5_D];
Tuning_D=[Tuning_T4_D,Tuning_T5_D];
[Min,min_pos]=min(X_D);
[Max,max_pos]=max(X_D);
Prox_tun=Tuning_D(min_pos);   % TUNING ANGLE AT MOST PROXIMAL SIDE
Change_Tun_D=Tuning_D-Prox_tun;
Max_change_D=[Max_change_D,abs(Change_Tun_D(max_pos))];   % TUNING ANGLE AT MOST DISTAL SIDE

       
    
end 




figure
Tuning_change=[Max_change_A;Max_change_B;Max_change_C;Max_change_D];
boxplot(Tuning_change', 'Notch', 'on') 
title('Tuning change from proximal to distal')





% Plot the average projection of the brains 
for i=[2,4,11,5]
AV=Control_Bright.T4T5_mb(i).AV;

figure 
imagesc(mean(AV,3))
colormap('gray')

end 




%% Plot ColorCoded compass plots 

 Deg_cat=1:360/64:360; 

for i=4%1:length(HIT)
% 
    IFly=Control_Bright.T4T5_mb(HIT(2,i));
    IFly2=Control_Dark.T4T5_mb(HIT(1,i)); 
    
    figure
    subplot(1,2,1)
    P=compass(1);
    set(P, 'Visible', 'off')
    hold on

    T4A_ANGLE=convert_angle(angle(IFly.Z.T4A)); 
    Comp=compass(IFly.Z.T4A);
    for nn=1:length(Comp)
        T4A_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T4A_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    hold on 

    T4B_ANGLE=convert_angle(angle(IFly.Z.T4B)); 
     Comp=compass(IFly.Z.T4B);
    for nn=1:length(Comp)
        T4B_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T4B_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
      T4C_ANGLE=convert_angle(angle(IFly.Z.T4C)); 
     Comp=compass(IFly.Z.T4C);
    for nn=1:length(Comp)
        T4C_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T4C_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
      T4D_ANGLE=convert_angle(angle(IFly.Z.T4D)); 
     Comp=compass(IFly.Z.T4D);
    for nn=1:length(Comp)
        T4D_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T4D_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    
    
    
    subplot(1,2,2)
 
    P=compass(1);
    set(P, 'Visible', 'off')
    hold on

    T5A_ANGLE=convert_angle(angle(IFly2.Z.T5A)); 
    Comp=compass(IFly2.Z.T5A);
    for nn=1:length(Comp)
        T5A_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T5A_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    hold on   
    T5B_ANGLE=convert_angle(angle(IFly2.Z.T5B)); 
     Comp=compass(IFly2.Z.T5B);
    for nn=1:length(Comp)
        T5B_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T5B_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
      T5C_ANGLE=convert_angle(angle(IFly2.Z.T5C)); 
     Comp=compass(IFly2.Z.T5C);
    for nn=1:length(Comp)
        T5C_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T5C_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
      T5D_ANGLE=convert_angle(angle(IFly2.Z.T5D)); 
     Comp=compass(IFly2.Z.T5D);
    for nn=1:length(Comp)
        T5D_ANGLE(nn)
        [~, Pos]=min(abs(Deg_cat-T5D_ANGLE(nn))); 
        set(Comp(nn),'color',cmap(Pos,:));
    end
   
   
end 

