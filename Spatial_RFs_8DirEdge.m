%%
% Created on Apr 29 2020
% @author: miriam henning

%%% NOTE: This scripts only works with the raw p-data files, available upon
%%% request to msilies@uni-mainz.de or mhenning@uni-mainz.de

%%% The script can be used to analyze responses to the Edge or Stripe stimulus moving into 8 different directions
% Initialize

close all
clear all
clc
%% %%%%% Variables to set:
SAVE=0;   %Put 1 if you want to save!!
LPF=false; % use low pass filter
Delay=true; % correct for Calcium delay
Z_Transform=true; % Z-transform Data to equalize PD responses
PLOT=false; %Plot single ROIs and traces
%% %%%%%

addpath(genpath('subscripts'))

Condition='i_T4T5Rec_Control.*i_mb_movEdgesONOFF_optimized_8Dir.*i_responsivefly.*~moving.*i_PDATA_SIMA';

Foldertosave='Data/Data_Edges' ;

pdatapath='/Volumes/Seagate/Documents/2p-imaging/Miriam_pData';


cd(pdatapath);

database_select_samples_mh_DStuning; % get the right stuff...
f_mb_T4T5 = find(eval(Condition));


[fname,datetime,frames,zdepth,tottime,stimbouts,locnum,...
    activecells,inverted,stimcode,quality,driver,...
    moving,layer,wavelength,flyID, responsivefly,Layercheck,SingleFly,PTX_5_100,Washout5mM,Washout100mM]...
    =textread('MASTER_foldersummary_Miriam_sh.txt',...
    '%s %s %f %f %f %f %f %s %f %s %f %s %f %s %f %f %f %f %f %f %f %f\r',...
    'headerlines',1,'delimiter','\t');


currd=pwd;
% cd(locdir);
ScreenDimension=[53,74,78,72,53,74,78,72]; % These values describe the size of the screen in degress of the flies visual field
%Measured by Burak: Horizontal: 78deg
%Vertical: 53deg
%Diagonal:72deg

Cut=[-17,36;-24,50;-34,44;-44,28;-36,17;-50,24;-44,34;-28,44];


for iFLY=1:length(f_mb_T4T5) % for each Fly
    
    FLYname=fname{f_mb_T4T5(iFLY)};
    
    x=load(FLYname);
    
    
    
    % Now take the responses that are not yet averaged across repetition, so I
    % can compute the std
    ROIS_resp=x.strct.ClusterInfo_ManuallySelect.dSignal1_CA;
    Stimulus=x.strct.stim_type;
    fps=x.strct.xml.framerate;
    
    if contains(Stimulus, 'LumDecLumInc')
        CellT=[4,5];
    elseif contains(Stimulus, 'LumInc')
        CellT=4;
    elseif contains(Stimulus, 'LumDec')
        CellT=5;
    else
        disp('Wrong Stimulus')
    end
    
    T4T5=x.strct.ClusterInfo_ManuallySelect.T4_T5;
    Layer=x.strct.ClusterInfo_ManuallySelect.Layer;
    Stim=x.strct.fstimval;
    fps=x.strct.xml.framerate;
    %     center=x.strct.ClusterInfo_ManuallySelect.centers;
    masks=x.strct.ClusterInfo_ManuallySelect.masks_CA;
    AV = squeeze(sum(x.strct.ch1a_crop,3))/size(x.strct.ch1a_crop,3);
    
    

    %%
    
    for i=1:4 %for each Layer
        if length(CellT)>1
            CellInd=find((Layer==i));
        else
            CellInd=find((T4T5==CellT).*(Layer==i));
        end
        
        All_resp_dir_ACells=[];
        CellID_T5=[];
        RFloc_T5={};
        CellID_T4=[];
        RFloc_T4={};
        RFcenter_T4=[];
        RFcenter_T5=[];
        
        counT5=1;
        counT4=1;
        
        if ~isempty(CellInd)
            NUM_EPOCHS=8;
            
            %Find Epochs:
            DF=diff(Stim);
            Startepoch=find(DF>0);
            Startepoch=Startepoch(1:end-1); %Discard last epoch
            Endepoch=find(DF<0);
            if length(Endepoch)>length(Startepoch)
                Endepoch=Endepoch(1:end-1);
            end
            LengthEpoch=round(mean(Endepoch-Startepoch));
            
            StimNumber=Stim(Startepoch+1);
            Epochs=nan(LengthEpoch,length(StimNumber));
            for NEpoch=1:length(StimNumber)
                Epochs(:,NEpoch)=[Startepoch(NEpoch)+1:Startepoch(NEpoch)+LengthEpoch]';
            end
            
            
            %Get the response of each ROI
            %%
            kc=1;
            
            for k=1:length(CellInd) %for each cell
                All_resp_dir=[];
                ID=CellInd(k);
                
                ROI_Respi=ROIS_resp(ID,:);
                Layeri=Layer(ID);
                T4T5i=T4T5(ID);
                maski=masks{1,ID}; %Mask of this cell
                %                 centeri=center(:,ID);
                
                BaseLine_raw=x.strct.BaseLine;
                BaseLine_raw=BaseLine_raw/max(max(max(BaseLine_raw)));
                BaseLine=mean(BaseLine_raw,3);
                Baselineresp=BaseLine(find(maski));
                avBaselineresp=mean(Baselineresp); % Average Baseline response of this Cluster to grey Image
                
                ROI_Resp_DF=(ROI_Respi-avBaselineresp)/avBaselineresp;
                
                if Delay
                    % add delay of 9.6 degrees, measured by Burak
                    % (delay between a 20dps edge and the center of the
                    % RF probed with standing stripes)
                    StimulusSpeed=20; %deg/s
                    
                    shift=round(9.6/StimulusSpeed*fps);
                    
                    ROI_Resp_DF_delayed=circshift(ROI_Resp_DF,-shift);
                else
                    
                    ROI_Resp_DF_delayed=ROI_Resp_DF;
                end
                
 
                
                %Average responses across epoch repetitions
                All_maps=[];
                for kk=1:8 % loop through 8 directions of movement, 1= upward, 2=to right, 3=downward, 4=to left
                    
                    
                    Double=find(StimNumber==kk);
                    if length(Double)>1
                        allepochs=[];
                        for pp=1:length(Double)
                            RespEpi=ROI_Resp_DF_delayed(Epochs(:,Double(pp)));
                            allepochs=[allepochs;RespEpi];
                        end
                        RespEpi=mean(allepochs,1); %average Roi response across repetitions of stimulus direction
                        
                    else
                        RespEpi=ROI_Resp_DF_delayed(Epochs(:,Double));
                        
                    end
                    
                    % interpolate to 20Hz, because the bar moves with
                    % 20deg/s, so we now have the axis in deg
                    
                    t=[0:size(RespEpi,2)-1]/fps;
                    if contains(Stimulus, 'DriftingEdge')
                        it=[0:0.05:t(length(t))];
                        StartPos=[-60,-75,-75, -85, -75, -85, -85, -75] + 40; %for edge stimulus (80 deg wide bar)
                        
                    else
                        it=[0:0.05:t(length(t))]; % should be correct for both stim types, check again!
                        %                             it=[0.5/fps:0.05:(0.5)/fps+4.95];
                        StartPos=[-50, -50, -50, -50, -50, -50, -50, -50] + 5; % The bar is 10 deg wide, the starting position is determined by the middle of the bar: (10/2)
                        
                    end
                    
                    iRespEpi=interp1(t,RespEpi,it,'linear','extrap');
                    
                    if LPF
                        % low pass filter
                        [b,a]=butter(3,0.2);
                        iRespEpi_lp=filter(b,a,iRespEpi);
                        iRespEpi_lp=iRespEpi_lp-min(iRespEpi_lp);
                    else
                        iRespEpi_lp=iRespEpi;%-min(iRespEpi);
                    end
                    
                   
                    iRespEpi_lp_delayed=iRespEpi_lp;
                    
                    All_resp_dir=[All_resp_dir;iRespEpi_lp_delayed];
                    
                    
                    
                    % Now correct for time where stim is seen on the
                    % screen --> cut the trace
                    if contains(Stimulus, 'DriftingEdge') % This drifting edge stimulus_optimized is designed in a way that the edge directly enters the screen
                        %thus i do not need to cut anything but I need
                        %to devide into the part of ON or OFF edge
                        
                        if T4T5i==4 %Take on edge (second half)
                            iRespEpi_lp_delayed_half=iRespEpi_lp_delayed(round(size(iRespEpi_lp_delayed,2)/2)+1:end);
                            
                        elseif T4T5i==5 %Take first half
                            iRespEpi_lp_delayed_half=iRespEpi_lp_delayed(1:round(size(iRespEpi_lp_delayed,2)/2));
                            
                        end
                        Cuti=Cut(kk,:);
                        CutBeg=abs(StartPos(kk)-Cuti(1))+1;
                        CutEnd=CutBeg+sum(abs(Cuti))-1; %frames traveled before visible on the screen
                        
                        if CutEnd>length(iRespEpi_lp_delayed_half)
                            iRespEpi_lp_delayed_CUT=iRespEpi_lp_delayed_half(CutBeg:end);
                            iRespEpi_lp_delayed_CUT=[iRespEpi_lp_delayed_CUT,nan(1,sum(abs(Cuti))-length(iRespEpi_lp_delayed_CUT))];
                            
                        else
                            iRespEpi_lp_delayed_CUT=iRespEpi_lp_delayed_half(CutBeg:CutEnd);
                        end
                        
                    else
                        
                        Cuti=Cut(kk,:);
                        CutBeg=abs(StartPos(kk)-Cuti(1))+1;
                        CutEnd=CutBeg+sum(abs(Cuti))-1; %frames traveled before visible on the screen
                        iRespEpi_lp_delayed_CUT=iRespEpi_lp_delayed(CutBeg:CutEnd);
                        
                        
                    end
                    
                    % standardize responses so that DS responses dominate less
                    if Z_Transform
                        if max(iRespEpi_lp_delayed_CUT)>3*std(ROI_Resp_DF_delayed)+mean(ROI_Resp_DF_delayed)
                            SD=sqrt(nansum((iRespEpi_lp_delayed_CUT(:)).^2)/length(iRespEpi_lp_delayed_CUT(:)-1));
                            iRespEpi_z=iRespEpi_lp_delayed_CUT./SD;
                            %                                 disp(['Z-Transf- ', num2str(kk)])
                        else
                            iRespEpi_z=iRespEpi_lp_delayed_CUT;
                        end
                    else
                        iRespEpi_z=iRespEpi_lp_delayed_CUT;
                    end
                    
                    %map to the screen
                    orthogonalDir= kk+2;
                    if orthogonalDir>8
                        orthogonalDir=orthogonalDir-8;
                    end
                    
                    
                    DIM1=ScreenDimension(kk);
                    DIM2=ScreenDimension(orthogonalDir);
                    rotDIM=kk*45-45;
                    
                    
                    ScMAP=repmat(iRespEpi_z',[1,DIM2]);
                    ScMAP=flipdim(ScMAP,1);
                    
                    ScMAP=ScMAP+1;
                    ScMAP_rot=imrotate(ScMAP,-rotDIM);%,'nearest', 'crop');
                    ScMAP_rot(ScMAP_rot==0)=nan;
                    ScMAP_rot=ScMAP_rot-1; %This is done, because when rotating to diagonal directions it creates 0 at the egdes,...
                    % which I want to convert to nans, but to be sure, that I do not convert 0 within the RF I add a 1 before and...
                    %subtract it afterwards again
                    
                    if size(ScMAP_rot,1)>53 % Rotating to diagonal directions creates bigger dimensions, so I have to cut the result ...
                        % to the screen dimensions
                        
                        ScreenDIM1=53;
                        ScreenDIM2=78;
                        
                        Diff=ceil((size(ScMAP_rot,1)-ScreenDIM1)/2);
                        Diff2=ceil((size(ScMAP_rot,2)-ScreenDIM2)/2);
                        
                        
                        ScMAP_rot=ScMAP_rot(Diff:Diff+ScreenDIM1-1,Diff2:Diff2+ScreenDIM2-1);
                        
                    end
                    
                    All_maps=cat(3,All_maps,ScMAP_rot);
                    
                end
                %%
                ROI=mean(All_maps,3);
                
                
                
                % normalize
                ROI=ROI-nanmin(nanmin(ROI)); %min(mean(All_maps(:,:,[1,3,5,7])))));
                ROI=ROI./nanmax(nanmax(ROI));
                [ROIfit,A]=fit_2D_Gauss_on_RF(ROI);
                
                MAXi=nanmax(nanmax(nanmax(All_maps)));
                MINi=nanmin(nanmin(nanmin(All_maps)));
                
                if PLOT
                   
                    COLOR_OF_PLOT_GREEN = [0 .5 0];% GREEN
                    COLOR_OF_PLOT_RED = [1 0 0]; %RED
                    COLOR_OF_PLOT_YELLOW = [.7 .7 0];% YELLOW
                    COLOR_OF_PLOT_BLUE = [0 0 1]; %BLUE;
                    Color=[COLOR_OF_PLOT_GREEN;COLOR_OF_PLOT_BLUE;COLOR_OF_PLOT_RED;COLOR_OF_PLOT_YELLOW];
                    
                    
                    F3=figure('Position', [300 300 1400 400]);
                    
                    h3=subplot(1,3,1);
                    CMask = zeros(size(masks{1,1},1), size(masks{1,1},2), 3);% Miri 09.11.2018
                    curColor = Color(i,:);
                    curMask = cat(3,curColor(1).*maski,curColor(2).*maski,curColor(3).*maski);
                    CMask = CMask + curMask;
                    imshow(AV,[],'InitialMagnification',600);
                    hold(h3,'on');
                    hcllus=imshow(CMask,'InitialMagnification',600, 'Parent', h3);
                    set(hcllus,'AlphaData',0.5);
                    hold(h3,'off')
                    hold off
                    
                    subplot(1,3,2)
                    imagesc(ROI)
                    set(gca, 'XTick',[])
                    set(gca, 'YTick',[])
                    
                    
                    subplot(1,3,3)
                    imagesc(ROIfit)
                    set(gca, 'XTick',[])
                    set(gca, 'YTick',[])
                    
                    NameStop=strfind(FLYname,'_')-1;
                    
                    if SAVE
                        Savename=['RF_Center_MAP-', FLYname(1:NameStop(3)), '_T', num2str(CellT),' Layer', num2str(i),' CellID',num2str(ID)];
                        set(F3,'PaperSize', [50 20])
                        
                        saveas(F3, [Foldertosave,'/', Savename,'.pdf'])
                        
                        
                    end
                end
                
                
                
                ROIind=ROIfit>0.95;
                [X,Y]=find(ROIind);
                XDIM=max(X)-min(X);
                YDIM=max(Y)-min(Y);
                
                
                if k==1
                    CurrColMap=colormap;
                    ROIs_ON_Sc=zeros(size(ROI,1),size(ROI,2),3);
                    
                    
                    if size(CurrColMap,1)>length(CellInd)
                        Col_steps=1:floor(size(CurrColMap,1)/length(CellInd)):size(CurrColMap,1);
                    end
                    CMask = zeros(size(masks{1,1},1), size(masks{1,1},2), 3);% Miri 09.11.2018
                    
                end
                
                if XDIM<15 && YDIM<15
                    curColor = CurrColMap(Col_steps(k),:);
                    curMask = cat(3,curColor(1).*maski,curColor(2).*maski,curColor(3).*maski);
                    
                    try
                        CMask = CMask + curMask;
                    catch
                        CMask = zeros(size(masks{1,1},1), size(masks{1,1},2), 3);% Miri 09.11.2018
                        CMask = CMask + curMask;
                        
                    end
                    
                    curScrInd=cat(3,curColor(1).*ROIind,curColor(2).*ROIind,curColor(3).*ROIind);
                    ROIs_ON_Sc=ROIs_ON_Sc+curScrInd;
                    
                    if T4T5i==4
                        RFloc_T4{1,counT4}=ROIind;
                        CellID_T4(counT4)=ID;
                        RFcenter_T4(:,counT4)=[A(2);A(3)];  %[muX; muY]
                        counT4=counT4+1;
                    elseif T4T5i==5
                        CellID_T5(counT5)=ID;
                        RFloc_T5{1,counT5}=ROIind;
                        RFcenter_T5(:,counT5)=[A(2);A(3)];  %[muX; muY]
                        counT5=counT5+1;
                    end
                else
                    
                    
                    if T4T5i==4
                        RFloc_T4{1,counT4}=[];
                        CellID_T4(counT4)=ID;
                        RFcenter_T4(:,counT4)=[nan,nan];  %[muX; muY]
                        counT4=counT4+1;
                    elseif T4T5i==5
                        RFloc_T5{1,counT5}=[];
                        CellID_T5(counT5)=ID;
                        RFcenter_T5(:,counT5)=[nan,nan];  %[muX; muY]
                        counT5=counT5+1;
                    end
                end
                
            end
            
            %%
            
            NameStop=strfind(FLYname,'_')-1;
            
            Title=['RF Center MAP ', FLYname(1:NameStop(3)), ' AllCells  Layer', num2str(i)];
            
            F4=figure('Position', [300 300 1000 400]);
            h3=subplot(1,2,1);
            imshow(AV,[],'InitialMagnification',600);
            hold(h3,'on');
            hcllus=imshow(CMask,'InitialMagnification',600, 'Parent', h3);
            set(hcllus,'AlphaData',0.5);
            hold(h3,'off')
            hold off
            title(Title)
            
            subplot(1,2,2)
            imagesc(ROIs_ON_Sc)
            
            if SAVE
                Savename=['All ROIs_ON_Sc-', FLYname(1:NameStop(3)), '_T', num2str(CellT),' Layer', num2str(i)];
                set(F4,'PaperSize', [30 20])
                saveas(F4, [Foldertosave, '/Responses_to_Edges8Dir/',savein,'/', Savename,'.pdf'])
            end
            
        end
        
        if i==1
            RFloc.T4A=RFloc_T4;
            RFloc.T5A=RFloc_T5;
            RFCenter.T4A=RFcenter_T4;
            RFCenter.T5A=RFcenter_T5;
            CellID.T4A=CellID_T4;
            CellID.T5A=CellID_T5;
        elseif i==2
            RFloc.T4B=RFloc_T4;
            RFloc.T5B=RFloc_T5;
            RFCenter.T4B=RFcenter_T4;
            RFCenter.T5B=RFcenter_T5;
            CellID.T4B=CellID_T4;
            CellID.T5B=CellID_T5;
        elseif i==3
            RFloc.T4C=RFloc_T4;
            RFloc.T5C=RFloc_T5;
            RFCenter.T4C=RFcenter_T4;
            RFCenter.T5C=RFcenter_T5;
            CellID.T4C=CellID_T4;
            CellID.T5C=CellID_T5;
        elseif i==4
            RFloc.T4D=RFloc_T4;
            RFloc.T5D=RFloc_T5;
            RFCenter.T4D=RFcenter_T4;
            RFCenter.T5D=RFcenter_T5;
            CellID.T4D=CellID_T4;
            CellID.T5D=CellID_T5;
        end
        
        
        
        %% plot the ROIs together On screen
    end
%     T4T5_mb_new(iFLY).RFloc=RFloc;
    T4T5_mb_new(iFLY).RFCenter=RFCenter;
    T4T5_mb_new(iFLY).CellID=CellID;
    
    close all
end

% save([Foldertosave, '/Responses_to_Edges8Dir/',savein,  '/processed_Data_ROI_rf'],'T4T5_mb_new');

