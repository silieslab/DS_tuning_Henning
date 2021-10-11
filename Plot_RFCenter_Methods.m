% This scripts plots the RF location on screen of one recording 
% fig. S2C

addpath(genpath('subscripts'))

% Load preprocessed Data matrix

Control=load('Data/Data_Edges/processed_Data_SIMA_CS5_sh.mat');
Control_add=load('Data/Data_Edges/processed_Data_ROI_rf.mat');
load('Data/Data_Edges/Snob_Cluster_Info.mat') % load subtype identity from SNOB analysis

% Load text file with imaging conditins (e.g. z-depth and orientation to the screen)
[fname,turn,Zdepth]=textread('Data/Data_Edges/Turn_info.txt','%s %f %f','headerlines',0,'delimiter','\t');

for NF=1:size(Control.T4T5_mb,2)
%     Control.T4T5_mb(NF).RF_loc=Control_add.T4T5_mb_new(NF).RFloc;
    Control.T4T5_mb(NF).RFCenter=Control_add.T4T5_mb_new(NF).RFCenter;
    Control.T4T5_mb(NF).CellID=Control_add.T4T5_mb_new(NF).CellID;
    Control.T4T5_mb(NF).turn=turn(NF);
    Control.T4T5_mb(NF).Zdepth=Zdepth(NF);
end

%% 

for NFlyi=97:98
    h1=figure;
    Ifly=Control.T4T5_mb(NFlyi);
    
    try
        CMask = zeros(size(Ifly.Masks.T4A{1,1},1), size(Ifly.Masks.T4A{1,1},2), 3);
        cm=colormap;
        
        Turn=Ifly.turn;
        
        CoordX=-Ifly.RFCenter.T4A(2,:);
        CoordY=Ifly.RFCenter.T4A(1,:)+Turn;
        
        Quiver=[Ifly.Z.T4A]*10;
        
        Center=[];
        
        for NROi=1:size(Ifly.Masks.T4A,2)
            Masks=Ifly.Masks.T4A{NROi};
            [Y,X]=find(Masks);
            Center(NROi)=round(mean(X));
        end
        
        [CenterI,Seq]=sort(Center);
        CJ=floor(64/size(Ifly.Masks.T4A,2));
        
        for NROi=1:size(Ifly.Masks.T4A,2)
            NRO=Seq(NROi);
            
            if ~isnan(CoordX(NRO))
                Masks=Ifly.Masks.T4A;
                curColor = cm(NROi*CJ,:);
                curMask = cat(3,curColor(1).*Masks{NRO},curColor(2).*Masks{NRO},curColor(3).*Masks{NRO});
                CMask = CMask + curMask;
            end
            
        end
        
        
%         h1=figure;
        subplot(3,1,1)
        image=imshow(Ifly.AV,[],'InitialMagnification',600)
        hold on
        hcllus=imshow(CMask,'InitialMagnification',600,'Parent', h1.Children);
        set(hcllus,'AlphaData',0.5);
        h1.Children.Title.String=[Ifly.Flyname(1:12),'  Zdepth:',num2str(Ifly.Zdepth), '  turn:', num2str(Ifly.turn)];
        
        
        subplot(3,1,2)
        for NROi=1:size(Ifly.Masks.T4A,2)
            NRO=Seq(NROi);
            curColor = cm(NROi*CJ,:);
            
            plot(CoordY(NRO)-34,CoordX(NRO)+36, 'o', 'Color', curColor)
            hold on
            
        end
        axis('equal')
        set(gca,'XLim', [-34,44+2*45])
        set(gca,'YLim', [-17,36])
        title(['RF Loc'])
        
        
        subplot(3,1,3)
        for NROi=1:size(Ifly.Masks.T4A,2)
            NRO=Seq(NROi);
            curColor = cm(NROi*CJ,:);
            
            quiver(CoordY(NRO)-34,CoordX(NRO)+36, real(Quiver(NRO)), imag(Quiver(NRO)), 'Color',  curColor,  'AutoScale','off')
            hold on
            
        end
        
        axis('equal')
        set(gca,'XLim', [-34,44+2*45])
        set(gca,'YLim', [-17,36])
    end
    
end

