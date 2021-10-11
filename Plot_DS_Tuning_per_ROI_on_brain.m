
%%%%%%%%%%%%%%
% This skript plots all figures shown in Figure4
%%%%%%%%%%%%%%

addpath(genpath('subscripts'))
addpath(genpath('Data/Data_Edges'))

% addpath('/Users/mhenning/Documents/MATLAB/TP_code_Miriam/T4_T5 Analysis/ClusterAnalysisData')

% Load preprocessed Data matrix

Control_Dark=load('Data/Data_Stripes/Dark_Stripes/processed_Data_SIMA_CS5_sh.mat');
Control_Bright=load('Data/Data_Stripes/Bright_Stripes/processed_Data_SIMA_CS5_sh.mat');

CondName= 'Control';


load('MyColormap.mat')


HIT=4; %HIT contains number of all flies to plot
%HIT=1:11; for all flies


for i=1:length(HIT)
    
    IFly=Control_Bright.T4T5_mb(HIT(i));
    IFly2=Control_Dark.T4T5_mb(HIT(i));
    
    
    
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
    
    
    
    F0=figure('Position', [200 200 1000 400]);
    
    subplot(1,2,1)
    h = imagesc(ColROI_T4);
    colormap(cmap)
    caxis([0 360])
    title(IFly.Flyname(1:11))
    text(5,5,'T4','Color', [1 1 1], 'FontSize', 16)
    
    subplot(1,2,2)
    h = imagesc(ColROI_T5);
    colormap(cmap)
    colorbar
    caxis([0 360])
    text(5,5,'T5','Color', [1 1 1], 'FontSize', 16)
    
    
    % Repeat for combination of T4 and T5
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
    
    
    F1=figure;
    h = imagesc(ColROI);
    colorbar
    colormap(cmap)
    title(['Tuning direction: ', IFly.Flyname(1:11)])
    c=colorbar;
    % c.Limits=[0 360];
    caxis([0 360])
    text(5,5,'T4/T5','Color', [1 1 1], 'FontSize', 16)
    
    
    F2=figure;
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
    title(['Hue (tuning selectivity + tuning direction): ', IFly.Flyname(1:11)])
    
    
    
end




%% Plot ColorCoded compass plots

Deg_cat=1:360/64:360;

for i=1:length(HIT)
    %
    IFly=Control_Bright.T4T5_mb(HIT(i));
    IFly2=Control_Dark.T4T5_mb(HIT(i));
    
    figure
    subplot(1,3,1)
    P=compass(1);
    set(P, 'Visible', 'off')
    hold on
    
    T4A_ANGLE=convert_angle(angle(IFly.Z.T4A));
    Comp=compass(IFly.Z.T4A);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4A_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    hold on
    
    T4B_ANGLE=convert_angle(angle(IFly.Z.T4B));
    Comp=compass(IFly.Z.T4B);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4B_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    T4C_ANGLE=convert_angle(angle(IFly.Z.T4C));
    Comp=compass(IFly.Z.T4C);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4C_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    T4D_ANGLE=convert_angle(angle(IFly.Z.T4D));
    Comp=compass(IFly.Z.T4D);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4D_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    title('Tuning T4')
    
    
    
    subplot(1,3,2)
    
    P=compass(1);
    set(P, 'Visible', 'off')
    hold on
    
    T5A_ANGLE=convert_angle(angle(IFly2.Z.T5A));
    Comp=compass(IFly2.Z.T5A);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5A_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    hold on
    T5B_ANGLE=convert_angle(angle(IFly2.Z.T5B));
    Comp=compass(IFly2.Z.T5B);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5B_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    T5C_ANGLE=convert_angle(angle(IFly2.Z.T5C));
    Comp=compass(IFly2.Z.T5C);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5C_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    T5D_ANGLE=convert_angle(angle(IFly2.Z.T5D));
    Comp=compass(IFly2.Z.T5D);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5D_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    title('Tuning T5')
    
    
    
    
    subplot(1,3,3)
    P=compass(1);
    set(P, 'Visible', 'off')
    hold on
    
    T4A_ANGLE=convert_angle(angle(IFly.Z.T4A));
    Comp=compass(IFly.Z.T4A);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4A_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    hold on
    T5A_ANGLE=convert_angle(angle(IFly2.Z.T5A));
    Comp=compass(IFly2.Z.T5A);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5A_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    
    T4B_ANGLE=convert_angle(angle(IFly.Z.T4B));
    Comp=compass(IFly.Z.T4B);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4B_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    T5B_ANGLE=convert_angle(angle(IFly2.Z.T5B));
    Comp=compass(IFly2.Z.T5B);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5B_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    T4C_ANGLE=convert_angle(angle(IFly.Z.T4C));
    Comp=compass(IFly.Z.T4C);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4C_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    T5C_ANGLE=convert_angle(angle(IFly2.Z.T5C));
    Comp=compass(IFly2.Z.T5C);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5C_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    
    T4D_ANGLE=convert_angle(angle(IFly.Z.T4D));
    Comp=compass(IFly.Z.T4D);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T4D_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    T5D_ANGLE=convert_angle(angle(IFly2.Z.T5D));
    Comp=compass(IFly2.Z.T5D);
    for nn=1:length(Comp)
        [~, Pos]=min(abs(Deg_cat-T5D_ANGLE(nn)));
        set(Comp(nn),'color',cmap(Pos,:));
    end
    title('Tuning T4 and T5')
    
    
end


%% Plot tuning with respect to quantifyied position- proximal vs. distal

for i=1:length(HIT) %do this just for one example fly
    
    IFly=Control_Bright.T4T5_mb(HIT(i));
    IFly2=Control_Dark.T4T5_mb(HIT(i)); 
    
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




%Load P-data for PixelResolution
load(['/Volumes/SILIESLAB/MiriH/2p-imaging/Raw_Data/Miriam_pData/',IFly.Flyname])
PixelSI=str2num(strct.xml.micronsPerPixel.XAxis);
Deg_cat=1:360/64:360;

%now with tuning color
figure 
subplot(2,2,1)
X_axis_pos_T4_AmM=(X_axis_pos_T4_A-min([X_axis_pos_T4_A,X_axis_pos_T5_A]))*PixelSI;
X_axis_pos_T5_AmM=(X_axis_pos_T5_A-min([X_axis_pos_T4_A,X_axis_pos_T5_A]))*PixelSI;

for NC=1:length(X_axis_pos_T4_AmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T4_A(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T4_AmM(NC),Tuning_T4_A(NC), [], Colori)
    hold on
end 
hold on 

for NC=1:length(X_axis_pos_T5_AmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T5_A(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T5_AmM(NC),Tuning_T5_A(NC), [], Colori)
    hold on
end 
title('LayerA')
set(gca, 'Ylim', [10 70])



X_axis_pos_T4_BmM=(X_axis_pos_T4_B-min([X_axis_pos_T4_B,X_axis_pos_T5_B]))*PixelSI;
X_axis_pos_T5_BmM=(X_axis_pos_T5_B-min([X_axis_pos_T4_B,X_axis_pos_T5_B]))*PixelSI;


subplot(2,2,2)
for NC=1:length(X_axis_pos_T4_BmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T4_B(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T4_BmM(NC),Tuning_T4_B(NC), [], Colori)
    hold on
end 

for NC=1:length(X_axis_pos_T5_BmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T5_B(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T5_BmM(NC),Tuning_T5_B(NC), [], Colori)
    hold on
end
title('LayerB')
set(gca, 'Ylim', [190 250])


subplot(2,2,3)

X_axis_pos_T4_CmM=(X_axis_pos_T4_C-min([X_axis_pos_T4_C,X_axis_pos_T5_C]))*PixelSI;
X_axis_pos_T5_CmM=(X_axis_pos_T5_C-min([X_axis_pos_T4_C,X_axis_pos_T5_C]))*PixelSI;

for NC=1:length(X_axis_pos_T4_CmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T4_C(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T4_CmM(NC),Tuning_T4_C(NC), [], Colori)
    hold on
end

for NC=1:length(X_axis_pos_T5_CmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T5_C(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T5_CmM(NC),Tuning_T5_C(NC), [], Colori)
    hold on
end
title('LayerC')
set(gca, 'Ylim', [80 140])


subplot(2,2,4)

X_axis_pos_T4_DmM=(X_axis_pos_T4_D-min([X_axis_pos_T4_D,X_axis_pos_T5_D]))*PixelSI;
X_axis_pos_T5_DmM=(X_axis_pos_T5_D-min([X_axis_pos_T4_D,X_axis_pos_T5_D]))*PixelSI;

for NC=1:length(X_axis_pos_T4_DmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T4_D(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T4_DmM(NC),Tuning_T4_D(NC), [], Colori)
    hold on
end

for NC=1:length(X_axis_pos_T5_DmM)
    [~, Pos]=min(abs(Deg_cat-Tuning_T5_D(NC)));
    Colori=cmap(Pos,:);
    scatter(X_axis_pos_T5_DmM(NC),Tuning_T5_D(NC), [], Colori)
    hold on
end
title('LayerD')
set(gca, 'Ylim', [250 310])



