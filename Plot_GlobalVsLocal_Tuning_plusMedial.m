

%%%%%%%%%%%%%%
% This skript plots compass plots for local vs. global tuning shown in
% Figure 3c,d,e and Extended Data Figure 3c
% %% NOTE!: This skript loops through all 14 flies and plots Z30 and Z60 (or whatever specified in 'z-depth'),
%           but not all Flies were recorded at that z-depth planes, therefore some local plots are empty!   
%%%%%%%%%%%%%%

addpath(genpath('subscripts'))
addpath(genpath('Data/Data_Edges'))

% Load preprocessed Data matrix
Conditions.Control=load('Data/Data_Edges/processed_Data_SIMA_CS5_sh.mat');
% Load text file with imaging conditins (e.g. z-depth and orientation to the screen)
[fname,turn,Zdepth]=textread('Data/Data_Edges/Turn_info.txt','%s %f %f','headerlines',0,'delimiter','\t');

for NF=1:size(Conditions.Control.T4T5_mb,2)
    Conditions.Control.T4T5_mb(NF).turn=turn(NF);
    Conditions.Control.T4T5_mb(NF).z_depth=Zdepth(NF);
end

load('Data/Data_Edges/Snob_Cluster_Info.mat') % load subtype identity from SNOB analysis

Foldertosave='/Volumes/SILIESLAB/MiriH/2p-imaging/Results_T4T5_Imaging_DS_tuning/Responses_to_Edges8Dir/'
%% Define

% Z_depth=[45,60,75]; % choose z-depth planes to compare: here Z-depth of 30(most dorsal plane) vs 60(ventral plane)
Z_depth=[30,45,60]; % choose z-depth planes to compare: here Z-depth of 30(most dorsal plane) vs 60(ventral plane)

Turni=0;
%specifiy Colors for plotting
ColorAI=[23,106,0]/255;
ColorAII=[151,217,0]/255;
ColorBI=[0,0,154]/255;
ColorBII=[117,0,154]/255;
ColorC=[240,0,0]/255;
ColorD=[214,214,0]/255;

Layer={'A','B', 'C', 'D'};
T4T5={'T4', 'T5'};

%% Initializing
count_T4A=1;
count_T5A=1;
count_T4B=1;
count_T5B=1;
count_T4C=1;
count_T5C=1;
count_T4D=1;
count_T5D=1;


DSI_byFly_T4_AII=nan(14,1000);
DSI_byFly_T4_AI=nan(14,1000);
DSI_byFly_T5_AII=nan(14,1000);
DSI_byFly_T5_AI=nan(14,1000);
DSI_byFly_T4_BII=nan(14,1000);
DSI_byFly_T4_BI=nan(14,1000);
DSI_byFly_T5_BII=nan(14,1000);
DSI_byFly_T5_BI=nan(14,1000);
DSI_byFly_T4_C=nan(14,1000);
DSI_byFly_T5_C=nan(14,1000);
DSI_byFly_T4_D=nan(14,1000);
DSI_byFly_T5_D=nan(14,1000);

%% Main

FlytoPlot='200729_Fly4';
count_Flyn=1;
Flyname_before=['noname'];
for II=1:size(Conditions.Control.T4T5_mb,2)
    
    IFly=Conditions.Control.T4T5_mb(II);
    ZDepth=IFly.z_depth;
    TURN=IFly.turn;
    Flyname=IFly.Flyname;
    Ind=strfind(Flyname, '_');
    Flynamei=Flyname(1:Ind(2)-1);
    
    
    if ~strcmp(Flynamei, Flyname_before)
        
        if II>1
            DSI_byFly_T4_AI(count_Flyn,1:length(DSI_T4_AI))=DSI_T4_AI;
            DSI_byFly_T4_AII(count_Flyn,1:length(DSI_T4_AII))=DSI_T4_AII;
            DSI_byFly_T5_AI(count_Flyn,1:length(DSI_T5_AI))=DSI_T5_AI;
            DSI_byFly_T5_AII(count_Flyn,1:length(DSI_T5_AII))=DSI_T5_AII;
            DSI_byFly_T4_BI(count_Flyn,1:length(DSI_T4_BI))=DSI_T4_BI;
            DSI_byFly_T4_BII(count_Flyn,1:length(DSI_T4_BII))=DSI_T4_BII;
            DSI_byFly_T5_BI(count_Flyn,1:length(DSI_T5_BI))=DSI_T5_BI;
            DSI_byFly_T5_BII(count_Flyn,1:length(DSI_T5_BII))=DSI_T5_BII;
            DSI_byFly_T4_C(count_Flyn,1:length(DSI_T4_C))=DSI_T4_C;
            DSI_byFly_T5_C(count_Flyn,1:length(DSI_T5_C))=DSI_T5_C;
            DSI_byFly_T4_D(count_Flyn,1:length(DSI_T4_D))=DSI_T4_D;
            DSI_byFly_T5_D(count_Flyn,1:length(DSI_T5_D))=DSI_T5_D;
            
            count_Flyn=count_Flyn+1;
            
            %Save previous figure
            if strcmp(FlytoPlot, Flyname_before)
            
            F1=figure('Position',[200 200 1180 700]);
            subplot(2,4,1)
            P=compass(1);
            set(P, 'Visible', 'off')
            hold on
            title(['Local Tuning Z-depth:', num2str(Z_depth(1)), ' _',Flyname_before])
            Comp=compass(DSI_AI_Z1);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAI);
            end
            Comp=compass(DSI_AII_Z1);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAII);
            end
            Comp=compass(DSI_BI_Z1);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBI);
            end
            Comp=compass(DSI_BII_Z1);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBII);
            end
            Comp=compass(DSI_C_Z1);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorC);
            end
            Comp=compass(DSI_D_Z1);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorD);
            end
            
            
            subplot(2,4,2)
            P=compass(1);
            set(P, 'Visible', 'off')
            hold on
            title(['Local Tuning Z-depth:', num2str(Z_depth(2))])
            Comp=compass(DSI_AI_Z2);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAI);
            end
            Comp=compass(DSI_AII_Z2);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAII);
            end
            Comp=compass(DSI_BI_Z2);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBI);
            end
            Comp=compass(DSI_BII_Z2);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBII);
            end
            Comp=compass(DSI_C_Z2);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorC);
            end
            Comp=compass(DSI_D_Z2);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorD);
            end
            
            
            subplot(2,4,3)
            P=compass(1);
            set(P, 'Visible', 'off')
            hold on
            title(['Local Tuning Z-depth:', num2str(Z_depth(3))])
            Comp=compass(DSI_AI_Z3);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAI);
            end
            Comp=compass(DSI_AII_Z3);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAII);
            end
            Comp=compass(DSI_BI_Z3);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBI);
            end
            Comp=compass(DSI_BII_Z3);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBII);
            end
            Comp=compass(DSI_C_Z3);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorC);
            end
            Comp=compass(DSI_D_Z3);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorD);
            end
            
            
            subplot(2,4,4)
            P=compass(1);
            set(P, 'Visible', 'off')
            hold on
            title(['Global Tuning'])
            Comp=compass([DSI_T4_AI,DSI_T5_AI]);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAI);
            end
            Comp=compass([DSI_T4_AII,DSI_T5_AII]);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorAII);
            end
            Comp=compass([DSI_T4_BI,DSI_T5_BI]);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBI);
            end
            Comp=compass([DSI_T4_BII,DSI_T5_BII]);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorBII);
            end
            Comp=compass([DSI_T4_C,DSI_T5_C]);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorC);
            end
            Comp=compass([DSI_T4_D,DSI_T5_D]);
            for ic=1:length(Comp)
                set(Comp(ic),'color',ColorD);
            end
            
            
            nbins= 72;
            subAx1 = subplot(2, 4, 5, polaraxes);
            
            if ~isempty([DSI_C_Z1])
                obj5 = Plot_circHist([DSI_C_Z1], ColorC, subAx1, [-0.5, 4],nbins);
                delete(obj5.scaleBar)
            end
            
            if ~isempty([DSI_D_Z1])
                obj6 = Plot_circHist([DSI_D_Z1], ColorD, subAx1, [-0.5, 4],nbins);
                delete(obj6.scaleBar)
            end
            
            if ~isempty([DSI_AI_Z1])
                obj1 = Plot_circHist([DSI_AI_Z1], ColorAI, subAx1, [-0.5, 4],nbins);
                delete(obj1.scaleBar)
            end
            
            if ~isempty([DSI_AII_Z1])
                obj2 = Plot_circHist([DSI_AII_Z1], ColorAII, subAx1, [-0.5, 4],nbins);
            end
            
            if ~isempty([DSI_BI_Z1])
                obj3 = Plot_circHist([DSI_BI_Z1], ColorBI, subAx1, [-0.5, 4],nbins);
                delete(obj3.scaleBar)
            end
            
            if ~isempty([DSI_BII_Z1])
                obj4 = Plot_circHist([DSI_BII_Z1], ColorBII, subAx1, [-0.5, 4],nbins);
                delete(obj4.scaleBar)
            end
            
            
            
            
            subAx2 = subplot(2, 4, 6, polaraxes);
            
            if ~isempty([DSI_C_Z2])
                obj5 = Plot_circHist([DSI_C_Z2], ColorC, subAx2, [-0.5,4],nbins);
                delete(obj5.scaleBar)
            end
            
            if ~isempty([DSI_D_Z2])
                obj6 = Plot_circHist([DSI_D_Z2], ColorD, subAx2, [-0.5, 4],nbins);
                delete(obj6.scaleBar)
            end
            
            if ~isempty([DSI_AI_Z2])
                obj1 = Plot_circHist([DSI_AI_Z2], ColorAI, subAx2, [-0.5, 4],nbins);
            end
            
            if ~isempty([DSI_AII_Z2])
                obj2 = Plot_circHist([DSI_AII_Z2], ColorAII, subAx2, [-0.5, 4],nbins);
                delete(obj2.scaleBar)
            end
            
            if ~isempty([DSI_BI_Z2])
                obj3 = Plot_circHist([DSI_BI_Z2], ColorBI, subAx2, [-0.5, 4],nbins);
                delete(obj3.scaleBar)
            end
            
            if ~isempty([DSI_BII_Z2])
                obj4 = Plot_circHist([DSI_BII_Z2], ColorBII, subAx2, [-0.5, 4],nbins);
                delete(obj4.scaleBar)
            end
            
            
             subAx3 = subplot(2, 4, 7, polaraxes);
            
            if ~isempty([DSI_C_Z3])
                obj5 = Plot_circHist([DSI_C_Z3], ColorC, subAx3, [-0.5,6],nbins);
                delete(obj5.scaleBar)
            end
            
            if ~isempty([DSI_D_Z3])
                obj6 = Plot_circHist([DSI_D_Z3], ColorD, subAx3, [-0.5, 6],nbins);
                delete(obj6.scaleBar)
            end
            
            if ~isempty([DSI_AI_Z3])
                obj1 = Plot_circHist([DSI_AI_Z3], ColorAI, subAx3, [-0.5,6],nbins);
            end
            
            if ~isempty([DSI_AII_Z3])
                obj2 = Plot_circHist([DSI_AII_Z3], ColorAII, subAx3, [-0.5, 6],nbins);
                delete(obj2.scaleBar)
            end
            
            if ~isempty([DSI_BI_Z3])
                obj3 = Plot_circHist([DSI_BI_Z3], ColorBI, subAx3, [-0.5,6],nbins);
                delete(obj3.scaleBar)
            end
            
            if ~isempty([DSI_BII_Z3])
                obj4 = Plot_circHist([DSI_BII_Z3], ColorBII, subAx3, [-0.5, 6],nbins);
                delete(obj4.scaleBar)
            end
            
            
            subAx4 = subplot(2, 4, 8, polaraxes);
            if ~isempty([DSI_T4_C,DSI_T5_C])
                obj5 = Plot_circHist([DSI_T4_C,DSI_T5_C], ColorC, subAx4, [-0.5, 12],nbins);
                delete(obj5.scaleBar)
                
            end
            
            if ~isempty([DSI_T4_D,DSI_T5_D])
                obj6 = Plot_circHist([DSI_T4_D,DSI_T5_D], ColorD, subAx4, [-0.5, 12],nbins);
                delete(obj6.scaleBar)
            end
            
            if ~isempty([DSI_T4_AI,DSI_T5_AI])
                obj1 = Plot_circHist([DSI_T4_AI,DSI_T5_AI], ColorAI, subAx4, [-0.5, 12],nbins);
            end
            
            if ~isempty([DSI_T4_AII,DSI_T5_AII])
                obj2 = Plot_circHist([DSI_T4_AII,DSI_T5_AII], ColorAII, subAx4, [-0.5, 12],nbins);
                delete(obj2.scaleBar)
            end
            
            if ~isempty([DSI_T4_BI,DSI_T5_BI])
                obj3 = Plot_circHist([DSI_T4_BI,DSI_T5_BI], ColorBI, subAx4, [-0.5, 12],nbins);
                delete(obj3.scaleBar)
            end
            
            if ~isempty([DSI_T4_BII,DSI_T5_BII])
                obj4 = Plot_circHist([DSI_T4_BII,DSI_T5_BII], ColorBII, subAx4, [-0.5, 12],nbins);
                delete(obj4.scaleBar)
                
            end
            
            
            
            set(F1,'PaperSize', [50 40])
            F1.Renderer='Painters';
%                   saveas(F1, [Foldertosave, 'Global_vs_LocalTuning_',Flyname_before,'_plusMedial.pdf'])
%             
%                           close all
           end 
        end
        
        DSI_T4_AI=[];
        DSI_T5_AI=[];
        DSI_T4_AII=[];
        DSI_T5_AII=[];
        DSI_T4_BI=[];
        DSI_T5_BI=[];
        DSI_T4_BII=[];
        DSI_T5_BII=[];
        DSI_T4_C=[];
        DSI_T5_C=[];
        DSI_T4_D=[];
        DSI_T5_D=[];
        
        
        DSI_AI_Z1=[];
        DSI_AI_Z2=[];
        DSI_AI_Z3=[];
        DSI_AII_Z1=[];
        DSI_AII_Z2=[];
        DSI_AII_Z3=[];
        DSI_BI_Z1=[];
        DSI_BI_Z2=[];
        DSI_BI_Z3=[];
        DSI_BII_Z1=[];
        DSI_BII_Z2=[];
        DSI_BII_Z3=[];
        DSI_C_Z1=[];
        DSI_C_Z2=[];
        DSI_C_Z3=[];
        DSI_D_Z1=[];
        DSI_D_Z2=[];
        DSI_D_Z3=[];
        
    end
    
    
    
    for La= 1:4
        Layeri= Layer{La};
        
        for T=1:2
            T4T5i=T4T5{T};
            
            
            Z=eval(['IFly.Z.', T4T5i,Layeri]);
            
            nCells=size(Z,2);
            
            count=eval(['count_',T4T5i,Layeri]);
            Clusteri=eval(['ClusterR.T',Layeri,'_', T4T5i]);
            Clusterii=Clusteri(count:count+nCells-1);
            
            
            if nCells>0
                for i = 1:nCells
                    
                    if La==1 
                        
                        if Clusterii(i)==3
                            if strcmp(T4T5i, 'T4')
                                DSI_T4_AI=[DSI_T4_AI,Z(i)];
                            elseif strcmp(T4T5i, 'T5')
                                DSI_T5_AI=[DSI_T5_AI,Z(i)];
                            end
                            
                            if ZDepth==Z_depth(1) && TURN==Turni
                                DSI_AI_Z1=[DSI_AI_Z1,Z(i)]; %save again only for certain Z_depth
                            elseif ZDepth==Z_depth(2) && TURN==Turni
                                DSI_AI_Z2=[DSI_AI_Z2,Z(i)];
                            elseif ZDepth==Z_depth(3) && TURN==Turni
                                DSI_AI_Z3=[DSI_AI_Z3,Z(i)];
                            end
                            
                            
                            
                        elseif Clusterii(i)==1
                            if strcmp(T4T5i, 'T4')
                                DSI_T4_AII=[DSI_T4_AII,Z(i)];
                            elseif strcmp(T4T5i, 'T5')
                                DSI_T5_AII=[DSI_T5_AII,Z(i)];
                            end
                            
                            
                            if ZDepth==Z_depth(1) && TURN==Turni
                                DSI_AII_Z1=[DSI_AII_Z1,Z(i)];
                            elseif ZDepth==Z_depth(2)  && TURN==Turni
                                DSI_AII_Z2=[DSI_AII_Z2,Z(i)];
                            elseif ZDepth==Z_depth(3)  && TURN==Turni
                                DSI_AII_Z3=[DSI_AII_Z3,Z(i)];
                            end
                            
                        end
                        
                    elseif La==2 
                        
                        
                        if Clusterii(i)==1
                            if strcmp(T4T5i, 'T4')
                                DSI_T4_BI=[DSI_T4_BI,Z(i)];
                            elseif strcmp(T4T5i, 'T5')
                                DSI_T5_BI=[DSI_T5_BI,Z(i)];
                            end
                            
                            if ZDepth==Z_depth(1)
                                DSI_BI_Z1=[DSI_BI_Z1,Z(i)];
                            elseif ZDepth==Z_depth(2)
                                DSI_BI_Z2=[DSI_BI_Z2,Z(i)];
                            elseif ZDepth==Z_depth(3)
                                DSI_BI_Z3=[DSI_BI_Z3,Z(i)];
                            end
                            
                        elseif Clusterii(i)==2
                            if strcmp(T4T5i, 'T4')
                                DSI_T4_BII=[DSI_T4_BII,Z(i)];
                            elseif strcmp(T4T5i, 'T5')
                                DSI_T5_BII=[DSI_T5_BII,Z(i)];
                            end
                            
                            
                            if ZDepth==Z_depth(1)
                                DSI_BII_Z1=[DSI_BII_Z1,Z(i)];
                            elseif ZDepth==Z_depth(2)
                                DSI_BII_Z2=[DSI_BII_Z2,Z(i)];
                            elseif ZDepth==Z_depth(3)
                                DSI_BII_Z3=[DSI_BII_Z3,Z(i)];
                            end
                        end
                        
                        
                    elseif La==3 
                        
                        if Clusterii(i)==1 || Clusterii(i)==3
                            if strcmp(T4T5i, 'T4')
                                DSI_T4_C=[DSI_T4_C,Z(i)];
                            elseif strcmp(T4T5i, 'T5')
                                DSI_T5_C=[DSI_T5_C,Z(i)];
                            end
                            
                            if ZDepth==Z_depth(1) && TURN==Turni
                                DSI_C_Z1=[DSI_C_Z1,Z(i)];
                            elseif ZDepth==Z_depth(2) && TURN==Turni
                                DSI_C_Z2=[DSI_C_Z2,Z(i)];
                            elseif ZDepth==Z_depth(3) && TURN==Turni
                                DSI_C_Z3=[DSI_C_Z3,Z(i)];
                            end
                        end
                        
                        
                    elseif La==4 
                        if Clusterii(i)==1
                            if strcmp(T4T5i, 'T4')
                                DSI_T4_D=[DSI_T4_D,Z(i)];
                            elseif strcmp(T4T5i, 'T5')
                                DSI_T5_D=[DSI_T5_D,Z(i)];
                            end
                            
                            if ZDepth==Z_depth(1) && TURN==Turni
                                DSI_D_Z1=[DSI_D_Z1,Z(i)];
                            elseif ZDepth==Z_depth(2) && TURN==Turni
                                DSI_D_Z2=[DSI_D_Z2,Z(i)];
                            elseif ZDepth==Z_depth(3) && TURN==Turni
                                DSI_D_Z3=[DSI_D_Z3,Z(i)];
                            end
                            
                        end
                    end
                    
                end
            end
            
            
            count=count+nCells;
            
            
            
            if strcmp(T4T5i, 'T4') && strcmp(Layeri, 'A')
                count_T4A=count;
            elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'A')
                count_T5A=count;
            elseif strcmp(T4T5i, 'T4') && strcmp(Layeri, 'B')
                count_T4B=count;
            elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'B')
                count_T5B=count;
            elseif strcmp(T4T5i, 'T4') && strcmp(Layeri, 'C')
                count_T4C=count;
            elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'C')
                count_T5C=count;
            elseif strcmp(T4T5i, 'T4') && strcmp(Layeri, 'D')
                count_T4D=count;
            elseif strcmp(T4T5i, 'T5') && strcmp(Layeri, 'D')
                count_T5D=count;
            end
            
        end % for T
    end % for La
    %     end
    Flyname_before=Flynamei;
    
    if count_Flyn==14
        DSI_byFly_T4_AI(count_Flyn,1:length(DSI_T4_AI))=DSI_T4_AI;
        DSI_byFly_T4_AII(count_Flyn,1:length(DSI_T4_AII))=DSI_T4_AII;
        DSI_byFly_T5_AI(count_Flyn,1:length(DSI_T5_AI))=DSI_T5_AI;
        DSI_byFly_T5_AII(count_Flyn,1:length(DSI_T5_AII))=DSI_T5_AII;
        DSI_byFly_T4_BI(count_Flyn,1:length(DSI_T4_BI))=DSI_T4_BI;
        DSI_byFly_T4_BII(count_Flyn,1:length(DSI_T4_BII))=DSI_T4_BII;
        DSI_byFly_T5_BI(count_Flyn,1:length(DSI_T5_BI))=DSI_T5_BI;
        DSI_byFly_T5_BII(count_Flyn,1:length(DSI_T5_BII))=DSI_T5_BII;
        DSI_byFly_T4_C(count_Flyn,1:length(DSI_T4_C))=DSI_T4_C;
        DSI_byFly_T5_C(count_Flyn,1:length(DSI_T5_C))=DSI_T5_C;
        DSI_byFly_T4_D(count_Flyn,1:length(DSI_T4_D))=DSI_T4_D;
        DSI_byFly_T5_D(count_Flyn,1:length(DSI_T5_D))=DSI_T5_D;
    end
end
%%
Variance_all_T4=nan(6,14);
Variance_all_T5=nan(6,14);

for NFlies=1:14
    
    DSI=DSI_byFly_T4_AI(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T4(1,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T4_AII(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T4(2,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T4_BI(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T4(3,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T4_BII(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T4(4,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T4_C(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T4(5,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T4_D(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T4(6,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    
    
    
    DSI=DSI_byFly_T5_AI(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T5(1,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T5_AII(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T5(2,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T5_BI(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T5(3,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T5_BII(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T5(4,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T5_C(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T5(5,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
    DSI=DSI_byFly_T5_D(NFlies,:);
    DSI=DSI(~isnan(DSI));
    Variance_all_T5(6,NFlies)=circ_std(convert_angle(angle(DSI'), 'rad'));
    
end


%% (Extended Data Fig.4c)

% plot the std (variance of angular preference [deg]) 
colors=[ColorD;ColorC;ColorBII;ColorBI;ColorAII;ColorAI];

figure
subplot(2,1,1)
boxplot(rad2deg(Variance_all_T4'),'label', {'LayerA.I', 'LayerA.II', 'LayerB.I', 'LayerB.II', 'LayerC', 'LayerD'})
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',0.7);
end
ylabel('Std of angular preference')
title('T4')
set(gca, 'YLim', [0 45])
subplot(2,1,2)
boxplot(rad2deg(Variance_all_T5'),'label', {'LayerA.I', 'LayerA.II', 'LayerB.I', 'LayerB.II', 'LayerC', 'LayerD'})
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',0.7);
end
ylabel('Std of angular preference')
title('T5')
set(gca, 'YLim', [0 45])

%Statistics, this scripts tests for normal distributions of data 
% if data are normalized it does a ANOVA if not KKW test 
Multcomp_Stats(Variance_all_T4,  {'LayerA.I', 'LayerA.II', 'LayerB.I', 'LayerB.II', 'LayerC', 'LayerD'})
Multcomp_Stats(Variance_all_T5,  {'LayerA.I', 'LayerA.II', 'LayerB.I', 'LayerB.II', 'LayerC', 'LayerD'})



