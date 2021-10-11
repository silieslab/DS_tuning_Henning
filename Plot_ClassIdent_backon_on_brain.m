
%%%%%%%%%%%%%%
% This skript plots some figures shown in Figure3 and Extended Data Figure3
%%%%%%%%%%%%%%

addpath(genpath('subscripts'))
addpath(genpath('Data'))
%% Plot Subgroups of Layer A and B back to the brain

% to plot only the examples shown in the paper
IMAGES=[16,17,18,20,21,22,58,59,60];
% to plot all 114 example recordings
%IMAGES=[1:114];

%% load Data
% Load preprocessed Data matrix
Conditions.Control=load('Data/Data_Edges/processed_Data_SIMA_CS5_sh.mat');
% Load text file with imaging conditins (e.g. z-depth and orientation to the screen)
[fname,turn,Zdepth]=textread('Data/Data_Edges/Turn_info.txt','%s %f %f','headerlines',0,'delimiter','\t');

for NF=1:size(Conditions.Control.T4T5_mb,2)
    Conditions.Control.T4T5_mb(NF).turn=turn(NF);
    Conditions.Control.T4T5_mb(NF).z_depth=Zdepth(NF);
end

load('Data/Data_Edges/Snob_Cluster_Info.mat') % load subtype identity from SNOB analysis

%% Initialize Data Matrices

NCellsT4A=[];
NCellsT4B=[];
NCellsT5A=[];
NCellsT5B=[];
NCellsT4C=[];
NCellsT5C=[];
NCellsT4D=[];
NCellsT5D=[];

count_T4A=1;
count_T5A=1;
count_T4B=1;
count_T5B=1;
count_T4C=1;
count_T5C=1;
count_T4D=1;
count_T5D=1;

Layer={'A','B', 'C', 'D'};
T4T5={'T4', 'T5'};

cm_A=[151,217,0; 170, 170, 170; 23,106,0]/355;
cm_A=[166,217,63; 170, 170, 170; 0,127,0]/355;

% cm_B=[ 0,0,154; 117,0,154; 170, 170, 170]/355; purple color
cm_B=[ 0,0,154; 72,144,203; 170, 170, 170]/355;
cm_B=[ 0,0,222; 123,191,217; 170, 170, 170]/355;


cm_C=[ 237,32,39;  170, 170, 170; 237,32,39;]/355;
cm_D=[ 214,214,37; 170, 170, 170]/355;


ClusterIdent_A=[];
Z_depthI_A=[];
Roiloc_A=[];
FlyID_A=[];

ClusterIdent_B=[];
Z_depthI_B=[];
Roiloc_B=[];
FlyID_B=[];

ClusterIdent_C=[];
Z_depthI_C=[];
Roiloc_C=[];
FlyID_C=[];

ClusterIdent_D=[];
Z_depthI_D=[];
Roiloc_D=[];
FlyID_D=[];

FlyIDs=nan(1,114);
FlyIDs(1:9)=1; 
FlyIDs(10:18)=2;
FlyIDs(19:26)=3;
FlyIDs(27:38)=4;
FlyIDs(39:48)=5;
FlyIDs(49:56)=6;
FlyIDs(57:64)=7;
FlyIDs(65:70)=8;
FlyIDs(71:79)=9;
FlyIDs(80:84)=10;
FlyIDs(85:88)=11;
FlyIDs(89:93)=12;
FlyIDs(94:103)=13;
FlyIDs(104:114)=14;


%% Main
% loop through each Recording and plot ROIs back onto the brain
for II=1:size(Conditions.Control.T4T5_mb,2)
    
    IFly=Conditions.Control.T4T5_mb(II);
    ZDepth=IFly.z_depth;
    % Initialize mask, which will later contain the ROI locations
    CMask=[];
    try
        CMask = zeros(size(IFly.Masks.T4B{1,1},1), size(IFly.Masks.T4B{1,1},2), 3);% Miri 09.11.2018
    catch
        try
            CMask = zeros(size(IFly.Masks.T4A{1,1},1), size(IFly.Masks.T4A{1,1},2), 3);% Miri 09.11.2018
        catch
            try
                CMask = zeros(size(IFly.Masks.T5A{1,1},1), size(IFly.Masks.T5A{1,1},2), 3);% Miri 09.11.2018
            catch
                try
                    CMask = zeros(size(IFly.Masks.T5B{1,1},1), size(IFly.Masks.T5B{1,1},2), 3);% Miri 09.11.2018
                catch
                    CMask=[];
                end
            end
        end
    end
    
    %now loop through all cells of the current recording
    if ~isempty(CMask)
        for La= 1:4 % for each Layer
            Layeri= Layer{La};
            
            for T=1:2 % for T4 and T5 cells
                T4T5i=T4T5{T};
                
                
                Masks=eval(['IFly.Masks.', T4T5i,Layeri]);
                Clusteri=eval(['ClusterR.T',Layeri,'_', T4T5i]);
                count=eval(['count_',T4T5i,Layeri]); %For each 'Celltype' I count the number already assigned, so that I keep track of the SubtypeIdenty
                nMasks=size(Masks,2);
                Clusterii=Clusteri(count:count+nMasks-1);
                
                cm2=eval(['cm_', Layeri]); % use correct color code based on SubtypeIdentity
                
                
                if nMasks>0
                    for i = 1:nMasks
                        curColor = cm2(Clusterii(i),:);
                        curMask = cat(3,curColor(1).*Masks{i},curColor(2).*Masks{i},curColor(3).*Masks{i});
                        
                        CMask = CMask + curMask;
                        [~,col]=find(Masks{i});
                        
                        if La==1
                            ClusterIdent_A=[ClusterIdent_A,Clusterii(i)];
                            Z_depthI_A=[Z_depthI_A,ZDepth];
                            Roiloc_A=[Roiloc_A,mean(col)];
                            FlyID_A=[FlyID_A,FlyIDs(II)];
                            
                        elseif La==2
                            ClusterIdent_B=[ClusterIdent_B,Clusterii(i)];
                            Z_depthI_B=[Z_depthI_B,ZDepth];
                            Roiloc_B=[Roiloc_B,mean(col)];
                            FlyID_B=[FlyID_B,FlyIDs(II)];
                            
                        elseif La==3
                            ClusterIdent_C=[ClusterIdent_C,Clusterii(i)];
                            Z_depthI_C=[Z_depthI_C,ZDepth];
                            Roiloc_C=[Roiloc_C,mean(col)];
                            FlyID_C=[FlyID_C,FlyIDs(II)];
                            
                        elseif La==4
                            ClusterIdent_D=[ClusterIdent_D,Clusterii(i)];
                            Z_depthI_D=[Z_depthI_D,ZDepth];
                            Roiloc_D=[Roiloc_D,mean(col)];
                            FlyID_D=[FlyID_D,FlyIDs(II)];
                            
                        end
                        
                    end
                end
                count=count+nMasks;
                
                if strcmp(T4T5i, 'T4') && strcmp(Layeri, 'A')
                    count_T4A=count;
                elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'A')
                    count_T5A=count;
                elseif strcmp(T4T5i, 'T4') && strcmp(Layeri, 'B')
                    count_T4B=count;
                elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'B')
                    count_T5B=count;
                elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'C')
                    count_T5C=count;
                elseif strcmp(T4T5i, 'T4') && strcmp(Layeri, 'C')
                    count_T4C=count;
                elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'D')
                    count_T5D=count;
                elseif strcmp(T4T5i, 'T4') && strcmp(Layeri, 'D')
                    count_T4D=count;
                    
                end
                
            end
        end
    end
    
    % Now plot the ROIs on the brain, if it is in the range of recordings
    % that are supposed to be plotted: Define above
    if sum(II == IMAGES)
        
        AV = IFly.AV;
        pause(5) 
        h3=figure;
%         image=imshow(AV,[],'InitialMagnification',600);
%         h3.Children.Title.String=[IFly.Flyname(1:12),'  Zdepth:',num2str(IFly.z_depth), '  turn:', num2str(IFly.turn)];
        
        hold on
        CZ=CMask==0;
        for row=1:size(CMask,1)
            for col=1:size(CMask,2)
                
                if sum(CZ(row,col,:)==0)==0  % only if all three color indices are zero for an individual pixel, I put a one, such that the background would be white instead of black!! 
                    CMask(row,col,:)=1; 
                end 
                
            end 
        end 
        
        hcllus=imshow(CMask,'InitialMagnification',600);%, 'Parent', h3.Children);
%          set(hcllus,'AlphaData',0.5);
        h3.Children.Title.String=[IFly.Flyname(1:12),'  Zdepth:',num2str(IFly.z_depth), '  turn:', num2str(IFly.turn)];
        
        %add scale bar
        PixelSize=str2num(IFly.PixelSize);
        NumPixel= round(10/PixelSize);
        [YS,XS]=size(AV);
        XSt=XS-NumPixel-5; %to have the scale bar within the image
        YSt=YS-5;
        line([XSt XSt+NumPixel], [YSt YSt], 'Color', [1 1 1], 'LineWidth', 4 )
        
        
    end
end

%% Plot Histograms Counts of subtypes along dorso-ventral axis
% (old Figure 2  before revision)

% Zhist_A=nan(2,max(length(Z_depthI_A(ClusterIdent_A==1)),length(Z_depthI_A(ClusterIdent_A==3))));
% Zhist_A(1,1:length(Z_depthI_A(ClusterIdent_A==3)))=Z_depthI_A(ClusterIdent_A==3);
% Zhist_A(2,1:length(Z_depthI_A(ClusterIdent_A==1)))=Z_depthI_A(ClusterIdent_A==1);
% Zhist_B=nan(2,max(length(Z_depthI_B(ClusterIdent_B==1)),length(Z_depthI_B(ClusterIdent_B==2))));
% Zhist_B(1,1:length(Z_depthI_B(ClusterIdent_B==1)))=Z_depthI_B(ClusterIdent_B==1);
% Zhist_B(2,1:length(Z_depthI_B(ClusterIdent_B==2)))=Z_depthI_B(ClusterIdent_B==2);
% 
% figure
% subplot(2,1,1)
% hist(Zhist_A',4)
% h = findobj(gca,'Type','patch');
% h(1).FaceColor=cm_A(1,:);
% h(2).FaceColor=cm_A(3,:);
% title('LayerA')
% ylabel('# Count')
% legend({'subtype A.I', 'subtype A.II'})
% 
% subplot(2,1,2)
% hist(Zhist_B',4)
% h = findobj(gca,'Type','patch');
% h(1).FaceColor=cm_B(2,:);
% h(2).FaceColor=cm_B(1,:);
% legend({'subtype B.I', 'subtype B.II'})
% title('LayerB')
% xlabel('Z-depth')
% ylabel('# Count')
% 
% 
% 
% Zhist_C=nan(2,length(Z_depthI_C(ClusterIdent_C==1))+length(Z_depthI_C(ClusterIdent_C==3)));
% Zhist_C(1,:)=[Z_depthI_C(ClusterIdent_C==1),Z_depthI_C(ClusterIdent_C==3)];
% Zhist_D=nan(2,length(Z_depthI_D(ClusterIdent_D==1)));
% Zhist_D(1,:)=Z_depthI_D(ClusterIdent_D==1);
% 
% figure
% subplot(2,1,1)
% hist(Zhist_C',4)
% h = findobj(gca,'Type','patch');
% h(2).FaceColor=cm_C(1,:);
% legend({'subtype C'})
% title('LayerC')
% ylabel('# Count')
% 
% 
% subplot(2,1,2)
% hist(Zhist_D',4)
% h = findobj(gca,'Type','patch');
% h(2).FaceColor=cm_D(1,:);
% legend({'subtype D'})
% title('LayerD')
% xlabel('Z-depth')
% ylabel('# Count')


% %% Plot Histograms Counts of subtypes along dorso-ventral axis average per Fly and single Fly plots
% Zdepths=[30,35;45,50;60,65;75,80];
% CellCounts_AI=nan(14,4);
% CellCounts_AII=nan(14,4);
% CellCounts_BI=nan(14,4);
% CellCounts_BII=nan(14,4);
% 
% CellCounts_C=nan(14,4);
% CellCounts_D=nan(14,4);
% 
%  
% for FLy=1:14
%     
%     for Zdepth=1:size(Zdepths,1)
%         
%      CellCounts_AI(FLy,Zdepth)=length(find(Z_depthI_A(find((ClusterIdent_A==3).*(FlyID_A==FLy)))==Zdepths(Zdepth,1)))+...
%         length(find(Z_depthI_A(find((ClusterIdent_A==3).*(FlyID_A==FLy)))==Zdepths(Zdepth,2)));
%     
%      CellCounts_AII(FLy,Zdepth)=length(find(Z_depthI_A(find((ClusterIdent_A==1).*(FlyID_A==FLy)))==Zdepths(Zdepth,1)))+...
%         length(find(Z_depthI_A(find((ClusterIdent_A==1).*(FlyID_A==FLy)))==Zdepths(Zdepth,2)));
%     
%      CellCounts_BI(FLy,Zdepth)=length(find(Z_depthI_B(find((ClusterIdent_B==1).*(FlyID_B==FLy)))==Zdepths(Zdepth,1)))+...
%         length(find(Z_depthI_B(find((ClusterIdent_B==1).*(FlyID_B==FLy)))==Zdepths(Zdepth,2)));
%     
%      CellCounts_BII(FLy,Zdepth)=length(find(Z_depthI_B(find((ClusterIdent_B==2).*(FlyID_B==FLy)))==Zdepths(Zdepth,1)))+...
%         length(find(Z_depthI_B(find((ClusterIdent_B==2).*(FlyID_B==FLy)))==Zdepths(Zdepth,2)));
%     
%     
%     
%      CellCounts_C(FLy,Zdepth)=length(find(Z_depthI_C(find((((ClusterIdent_C==1)+ (ClusterIdent_C==3)).*(FlyID_C==FLy))))==Zdepths(Zdepth,1)))+...
%         length(find(Z_depthI_C(find((((ClusterIdent_C==1)+ (ClusterIdent_C==3)).*(FlyID_C==FLy))))==Zdepths(Zdepth,2)));
%     
% 
%      CellCounts_D(FLy,Zdepth)=length(find(Z_depthI_D(find((ClusterIdent_D==1).*(FlyID_D==FLy)))==Zdepths(Zdepth,1)))+...
%         length(find(Z_depthI_D(find((ClusterIdent_D==1).*(FlyID_D==FLy)))==Zdepths(Zdepth,2)));
% 
%     
%     
%     end
% end 
% 
% %%
% 
% 
% figure('Position', [200 200 1000 500])
% subplot(2,2,1)
% bar(Zdepths(:,2), mean(CellCounts_AI,1),0.3,'FaceColor',cm_A(3,:))
% hold on
% errorbar(Zdepths(:,2), mean(CellCounts_AI,1), std(CellCounts_AI,1)/sqrt(14), 'ko')
% 
% bar(Zdepths(:,1), mean(CellCounts_AII,1),0.3,'FaceColor',cm_A(1,:))
% errorbar(Zdepths(:,1), mean(CellCounts_AII,1), std(CellCounts_AII,1)/sqrt(14), 'ko')
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 30])
% 
% 
% 
% subplot(2,2,3)
% bar(Zdepths(:,2), mean(CellCounts_BI,1),0.3,'FaceColor',cm_B(1,:))
% hold on 
% errorbar(Zdepths(:,2), mean(CellCounts_BI,1), std(CellCounts_BI,1)/sqrt(14), 'ko')
% 
% bar(Zdepths(:,1), mean(CellCounts_BII,1),0.3,'FaceColor',cm_B(2,:))
% errorbar(Zdepths(:,1), mean(CellCounts_BII,1), std(CellCounts_BII,1)/sqrt(14), 'ko')
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% xlabel('Z-depth')
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 30])
% 
% 
% 
% 
% subplot(2,2,2)
% 
% bar(Zdepths(:,2)-2.5, mean(CellCounts_C,1),0.3,'FaceColor',cm_C(1,:))
% hold on 
% errorbar(Zdepths(:,2)-2.5, mean(CellCounts_C,1), std(CellCounts_C,1)/sqrt(14), 'ko')
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 30])
% 
% 
% 
% subplot(2,2,4)
% 
% bar(Zdepths(:,2)-2.5, mean(CellCounts_D,1),0.3,'FaceColor',cm_D(1,:))
% hold on 
% errorbar(Zdepths(:,2)-2.5, mean(CellCounts_D,1), std(CellCounts_D,1)/sqrt(14), 'ko')
% set(gca, 'XLim', [0 100])
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% xlabel('Z-depth')
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 30])
% 
% 
% 
% 
% %Now for single flies
% 
% 
% figure('Position', [200 200 1000 500])
% 
% subplot(2,2,1)
% hold on 
% plot(Zdepths(:,2), CellCounts_AI, 'o-', 'Color',cm_A(3,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
% hold on
% plot(Zdepths(:,1), CellCounts_AII, 'o-', 'Color',cm_A(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 50])
% 
% 
% 
% subplot(2,2,3)
% hold on
% plot(Zdepths(:,2), CellCounts_BI, 'o-', 'Color',cm_B(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
% hold on 
% plot(Zdepths(:,1), CellCounts_BII, 'o-', 'Color',cm_B(2,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% xlabel('Z-depth')
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 50])
% 
% 
% 
% subplot(2,2,2)
% 
% plot(Zdepths(:,2), CellCounts_C, 'o-', 'Color',cm_C(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 50])
% 
% 
% subplot(2,2,4)
% 
% plot(Zdepths(:,2), CellCounts_D, 'o-', 'Color',cm_D(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
% 
% set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
% set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
% xlabel('Z-depth')
% ylabel('Average Cell Count +/- ste')
% set(gca, 'Xlim', [20 90])
% set(gca, 'Ylim', [0 50])

 %% Plot percentage of cells
 
 Zdepths=[30,35;45,50;60,65;75,80];
CellCounts_AI=nan(14,4);
CellCounts_AII=nan(14,4);
CellCounts_BI=nan(14,4);
CellCounts_BII=nan(14,4);

CellCounts_C=nan(14,4);
CellCounts_D=nan(14,4);

 
for FLy=1:14
    
    for Zdepth=1:size(Zdepths,1)
        
     CellCounts_AI(FLy,Zdepth)=((length(find(Z_depthI_A(find((ClusterIdent_A==3).*(FlyID_A==FLy)))==Zdepths(Zdepth,1)))+...
        length(find(Z_depthI_A(find((ClusterIdent_A==3).*(FlyID_A==FLy)))==Zdepths(Zdepth,2))))/length(find(FlyID_A==FLy)))*100;
    
     CellCounts_AII(FLy,Zdepth)=((length(find(Z_depthI_A(find((ClusterIdent_A==1).*(FlyID_A==FLy)))==Zdepths(Zdepth,1)))+...
        length(find(Z_depthI_A(find((ClusterIdent_A==1).*(FlyID_A==FLy)))==Zdepths(Zdepth,2))))/length(find(FlyID_A==FLy)))*100;
    
     CellCounts_BI(FLy,Zdepth)=((length(find(Z_depthI_B(find((ClusterIdent_B==1).*(FlyID_B==FLy)))==Zdepths(Zdepth,1)))+...
        length(find(Z_depthI_B(find((ClusterIdent_B==1).*(FlyID_B==FLy)))==Zdepths(Zdepth,2))))/length(find(FlyID_B==FLy)))*100;
    
     CellCounts_BII(FLy,Zdepth)=((length(find(Z_depthI_B(find((ClusterIdent_B==2).*(FlyID_B==FLy)))==Zdepths(Zdepth,1)))+...
        length(find(Z_depthI_B(find((ClusterIdent_B==2).*(FlyID_B==FLy)))==Zdepths(Zdepth,2))))/length(find(FlyID_B==FLy)))*100;
    
    
    
     CellCounts_C(FLy,Zdepth)=((length(find(Z_depthI_C(find((((ClusterIdent_C==1)+(ClusterIdent_C==3)).*(FlyID_C==FLy))))==Zdepths(Zdepth,1)))+...
        length(find(Z_depthI_C(find((((ClusterIdent_C==1)+(ClusterIdent_C==3)).*(FlyID_C==FLy))))==Zdepths(Zdepth,2))))/length(find(FlyID_C==FLy)))*100;
    

     CellCounts_D(FLy,Zdepth)=((length(find(Z_depthI_D(find((ClusterIdent_D==1).*(FlyID_D==FLy)))==Zdepths(Zdepth,1)))+...
        length(find(Z_depthI_D(find((ClusterIdent_D==1).*(FlyID_D==FLy)))==Zdepths(Zdepth,2))))/length(find(FlyID_D==FLy)))*100;
    

    
    
    end
end 




figure('Position', [200 200 1000 500])
subplot(2,2,1)
bar(Zdepths(:,2), nanmean(CellCounts_AI,1),0.3,'FaceColor',cm_A(3,:))
hold on
errorbar(Zdepths(:,2), nanmean(CellCounts_AI,1), nanstd(CellCounts_AI,1)/sqrt(14), 'ko')

bar(Zdepths(:,1), nanmean(CellCounts_AII,1),0.3,'FaceColor',cm_A(1,:))
errorbar(Zdepths(:,1), nanmean(CellCounts_AII,1), nanstd(CellCounts_AII,1)/sqrt(14), 'ko')

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 40])



subplot(2,2,3)
bar(Zdepths(:,2), nanmean(CellCounts_BI,1),0.3,'FaceColor',cm_B(1,:))
hold on 
errorbar(Zdepths(:,2), nanmean(CellCounts_BI,1), nanstd(CellCounts_BI,1)/sqrt(14), 'ko')

bar(Zdepths(:,1), nanmean(CellCounts_BII,1),0.3,'FaceColor',cm_B(2,:))
errorbar(Zdepths(:,1), nanmean(CellCounts_BII,1), nanstd(CellCounts_BII,1)/sqrt(14), 'ko')

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
xlabel('Z-depth')
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 40])


subplot(2,2,2)

bar(Zdepths(:,2)-2.5, nanmean(CellCounts_C,1),0.3,'FaceColor',cm_C(1,:))
hold on 
errorbar(Zdepths(:,2)-2.5, nanmean(CellCounts_C,1), nanstd(CellCounts_C,1)/sqrt(14), 'ko')

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 40])



subplot(2,2,4)

bar(Zdepths(:,2)-2.5, nanmean(CellCounts_D,1),0.3,'FaceColor',cm_D(1,:))
hold on 
errorbar(Zdepths(:,2)-2.5, nanmean(CellCounts_D,1), nanstd(CellCounts_D,1)/sqrt(14), 'ko')

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
xlabel('Z-depth')
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 40])




%Now for single flies


% figure('Position', [200 200 1000 500])

subplot(2,2,1)
hold on 
plot(Zdepths(:,2), CellCounts_AI, 'o', 'Color',cm_A(3,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
hold on
plot(Zdepths(:,1), CellCounts_AII, 'o', 'Color',cm_A(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 70])



subplot(2,2,3)
hold on
plot(Zdepths(:,2), CellCounts_BI, 'o', 'Color',cm_B(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))
hold on 
plot(Zdepths(:,1), CellCounts_BII, 'o', 'Color',cm_B(2,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
xlabel('Z-depth')
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 70])



subplot(2,2,2)

plot(Zdepths(:,2)-2.5, CellCounts_C, 'o', 'Color',cm_C(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 70])


subplot(2,2,4)

plot(Zdepths(:,2)-2.5, CellCounts_D, 'o', 'Color',cm_D(1,:))% (FlyID,:),0.3,'FaceColor',cm_A(3,:))

set(gca, 'XTick', [32.5; 47.5; 62.5; 77.5])
set(gca, 'XTicklabel', [{'30-35'}, {'45-50'},{'60-65'},{'75-80'}])
xlabel('Z-depth')
ylabel('Percentage of cells +/- ste')
set(gca, 'Xlim', [20 90])
set(gca, 'Ylim', [0 70])





