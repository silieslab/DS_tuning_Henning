function [out,crop] = FindCanROIs_v9(in,mode, ManuallySelect,crop,ROISize)


%% Automatically find candidate ROIs
% Miriam Henning (15.05.17)

%%%INPUT%%%
% in  =strct containing imaging information from Batch Alignment 
% mode=1 if clusters should be calculated
%     =2 if masks are already existing
% ManuallySelect=1 it will ask you to select the four Layers manually 
%               =0 it will assign clusters to Layer identity based on their preferred direction 
% crop=[]; if the recorded image has not been cropped before 
% ROISize= desired size of the clusters/ROIs, based on how big the axon terminals are supposed to be 
%            for T4 and T5 I choose a size of 2.5uM 



%New in v3: correct calculation of ON/OFF split
%           new threshold calculation, whole response trace instead of averaged across epoch rep.
% New in v5: "Manually select Layer" mode AND correct use of treshold including the mean!!!!
% New in v6: " Correction after new alignment 09.02.18 &  Data differently saved for Layers selected manually
%          ClusterInfo are saved in out in clusterInfo structure

%New in v8: Adapted to new stimulus optimized!!
% v7 is also adapted (by sebas) but the way how data are saved was not
% changed by then.

%New in v9: 
% --> It additionally asks to specify region of lobula plate, to
%     prevent loosing Pixels from otsu thresholding due to bright cell bodies !!
% --> And added crop variable, this is used from cropping previous Images of the
%     same fly, to prevent haveing different image sizes for the different
%     stimuli recorded (if f.e.  responses to bright and dark stripes are going
%     to be compared e.g. plot PD on the brain, I need the same image sizes,
%     however this will only work if the image didnt drift too much between sequences of recording)
% --> real Background subtraction (by two times thresholding on the
%     original image size (before cropping)
% --> accepts also 8 direction edges
%%
DSIthreshold=-2000; % I do not want to include a DSI threshold, thats why I make it really small
CSIthreshold=0.2; % for all data concerning DS for stripe stimuli and so on I used 0.3 and 0.5 for STRF-data
timesSTD=2; % threshold for determining a real response, so only if the maximum is higher then this value x std of the trace 
%%

DesiredStim='DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs_optimized';
DesiredStim2='Search_DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs_optimized'; %Sometimes I accidentally pick the search stimulus, which does not matter, analysis stays the same
DesiredStim3='DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_8Dirs_optimized'; %Sometimes I accidentally pick the search stimulus, which does not matter, analysis stays the same


nframes = in.xml.frames;
fps = in.xml.framerate;

try
    ch1a=in.ch1a;
catch
    ch1a=in.ch1;
end


if max(max(max(ch1a)))>1 
    %convert to values between 0 and 1, this is normaly done in the batch
    %aligning 
    ch1a=ch1a/max(max(max(ch1a))); 
end 
    


if mode==1 % clusters should be calculated
    
    % Check if the Stimulus for calculating Cluster/ROIs was shown in that
    % recording:
    
    GoON=0;
    if strcmp(DesiredStim,in.stim_type)==1
        GoON=1;
        StimInd=[1:4];
        
    elseif strcmp(DesiredStim2,in.stim_type)==1
        GoON=1;
        StimInd=[1:4];

        
    elseif strcmp(DesiredStim3,in.stim_type)==1
        GoON=1;
        StimInd=[1:2:8];

    else
        disp(['The data are not recorded with the correct stimulus type. You need: DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs! Your stimulus used is: ', in.stim_type])
        
    end
    
    if GoON==1 % Stimulus is correct, continue!
        
        
        
        %% thresholding using Otsu' method, to identify pixels eligible for further analysis
        %(non background pixels)
        %   The algorithm assumes that the image contains two classes of pixels following bi-modal
        %   histogram (foreground pixels and background pixels), it then calculates the optimum
        %   threshold separating the two classes so that their combined spread (intra-class variance)
        %   is minimal, or equivalently (because the sum of pairwise squared distances is constant),
        %   so that their inter-class variance is maximal.
        
        if ~isfield(in,'ch1a_crop')            
            AV = squeeze(sum(ch1a,3))/nframes;
          
            if isempty(crop) % if crop is empty because for this fly the image has not been cropped yet
                F1=figure; 
                h=imshow(AV,[],'InitialMagnification',600)
                title('Crop the image')
                [I2,crop]=imcrop(h);
                ch1a_crop=ch1a(crop(2):crop(2)+crop(4),crop(1):crop(1)+crop(3),:); 

            else
                I2=AV(crop(2):crop(2)+crop(4),crop(1):crop(1)+crop(3));
                ch1a_crop=ch1a(crop(2):crop(2)+crop(4),crop(1):crop(1)+crop(3),:); 
            end 
            
        else % if this Image had been cropped before, I want to use the cropped image instead of cropping it again
            AV = squeeze(sum(ch1a,3))/nframes;
            ch1a_crop=in.ch1a_crop; 
            crop=in.crop; 
            I2=AV(crop(2):crop(2)+crop(4),crop(1):crop(1)+crop(3));
        end 
        

        GF=imgaussfilt(I2,1.5);
        [level,~] = graythresh(GF); 
        BW = im2bw(GF,level);
        
        figure;  
        imshow(BW,'InitialMagnification',600)
        title('Result of Gaussfilter + Otsu thresholding')

        
        
        %For getting real background I do another thresholding But on the
        %origial Image not the cropped one,because it will contain more
        %true background pixel. This background information will later be used for Background
        %subtraction of responsesbefore for dF/F calculation
        GF=imgaussfilt(AV,1.5);
        [level, ~] = graythresh(GF);
        BW2 = im2bw(GF,level);
        Background=nan(size(BW2));
        Background(BW2==0)=AV(BW2==0);
             
        %Second time threshold 
        V_BG=reshape(Background,1,size(Background,1)*size(Background,2));
        level = prctile(V_BG,10);
              
        Real_Background=zeros(size(BW2));
        Real_Background((Background<level))=1;
        
        
        %% Separate recorded trace into traces of certain directions of moving stim
        if ManuallySelect ==1
            Succeed=0;
            
            try load([in.fileloc,'/Layer_masks.mat'])
                if sum(size(masks_L{1,1}) == size(I2))==2  % If the size of the Layer masks fits the cropped image size
                    Succeed=1;
                end 
            catch      
            end
            
            if Succeed==1 % If Layers have been selected before
                curMask = cat(3,0.4.*masks_L{1},0.8.*masks_L{1},0.*masks_L{1});
                curMask2 = cat(3,0.*masks_L{2},0.4.*masks_L{2},0.8.*masks_L{2});
                curMask3 = cat(3,0.8.*masks_L{3},0.*masks_L{3},0.*masks_L{3});
                curMask4 = cat(3,0.8.*masks_L{4},0.8.*masks_L{4},0.*masks_L{4});
                curM=curMask+curMask2+curMask3+curMask4;
                AV = squeeze(sum(ch1a_crop,3))/nframes;
                
                figure;
                imshow(AV,[],'InitialMagnification',600);
                hold on;h = imshow(curM,'InitialMagnification',600);
                set(h,'AlphaData',0.2);
                title('Layers have been selected before')
                
            elseif Succeed==0
                curM=[];
                AV = squeeze(sum(ch1a_crop,3))/nframes;
                
                figure;imshow(AV,[],'InitialMagnification',600)
                title('Select Layer 1');
                masks_L{1} = roipoly;
                
                curMask = cat(3,0.4.*masks_L{1},0.8.*masks_L{1},0.*masks_L{1});
                curM=curMask;
                
                
                hold on;h = imshow(curM,'InitialMagnification',600);
                set(h,'AlphaData',0.2);
                title('Select Layer 2');
                hold off
                
                
                masks_L{2} = roipoly;
                
                curMask2 = cat(3,0.*masks_L{2},0.4.*masks_L{2},0.8.*masks_L{2});
                curM=curMask2+curMask;
                
                imshow(AV,[],'InitialMagnification',600);
                hold on;h = imshow(curM,'InitialMagnification',600);
                set(h,'AlphaData',0.2);
                title('Select Layer 3');
                hold off
                
                
                masks_L{3} = roipoly;
                
                curMask3 = cat(3,0.8.*masks_L{3},0.*masks_L{3},0.*masks_L{3});
                curM=curM+curMask3;
                
                imshow(AV,[],'InitialMagnification',600);
                hold on;h = imshow(curM,'InitialMagnification',600);
                set(h,'AlphaData',0.2);
                title('Select Layer 4 ');
                hold off
                
                
                masks_L{4} = roipoly;
                
                curMask4 = cat(3,0.8.*masks_L{4},0.8.*masks_L{4},0.*masks_L{4});
                curM=curM+curMask4;
                
                imshow(AV,[],'InitialMagnification',600);
                hold on;h = imshow(curM,'InitialMagnification',600);
                set(h,'AlphaData',0.2);
                
                curDir = pwd;
                
                save('Layer_masks','masks_L');
            end
            
            
            % Now use the mask info of the selected Layers and create a Pdir
            % matrix based on to which Layer this Pixel has been selected
            %Layer 1 has a preferred direction of right (2)
            %Layer 2 has a preferred direction of left (4)
            %Layer 3 has a preferred direction of up (1)
            %Layer 4 has a preferred direction of down (3)
            Pdir_select=2.*masks_L{1}+4.*masks_L{2}+1.*masks_L{3}+3.*masks_L{4};
        end
        
        fstimval=in.fstimval; % average stimulus value per recorded frame
        
        %Calculate average duration of each stimulus epoch
        changes=diff(fstimval);
        
        start_epoch=find(changes>0)+1;
        end_epoch=find(changes<0);
        
        dur=diff([start_epoch(1:length(end_epoch)),end_epoch]');
        avdur=min(dur);
        ep_all=nan(size(end_epoch,1), avdur); %ep_all contains the range of frames belonging to each stimulus epoch
        
        
        for pp=1:length(end_epoch)
            epoch=start_epoch(pp):start_epoch(pp)+avdur-1;
            ep_all(pp,:)=epoch;
        end
        ep_all=ep_all';
        dur=avdur;
        
        
        %% Calculate DSI, Pdir, ON/OFF resp
        % Compare each pixel reponse strength to one of the 4 directions of movements
        
        DSI=nan((size(ch1a_crop(:,:,1)))); % will contain a Direction Selectivity Index for each Pixel
        CSI=nan((size(ch1a_crop(:,:,1)))); % will contain a Contrast Selectivity Index for each Pixel
        Pdir=nan((size(ch1a_crop(:,:,1)))); % will contain the preferred Direction of each Pixel
        TRespPD=nan(size(ch1a_crop(:,:,1))); % will contain the timing of the Peak response to PD of each Pixel (relative to Stimulus onset)
        MaxRespPD=nan(size(ch1a_crop(:,:,1))); % will contain the Maximum of the Peak response to PD of each Pixel
        ONorOFFpixel=nan(size(ch1a_crop(:,:,1))); % will contain 0 for OFF (T5) and 1 for ON (T4) pixel
        Pdir_select_new=nan((size(ch1a_crop(:,:,1)))); %will contain the selected Layer info for all Pixels that have been declared to be foreground pixels (Otsu) and that overcome the selected threshold of response
        nomovement=find(diff(in.fstimpos1((ep_all(:,1))))==0);
        
        RespDir1_perPixel=nan(size(ch1a_crop,1),size(ch1a_crop,2),dur);
        RespDir2_perPixel=nan(size(ch1a_crop,1),size(ch1a_crop,2),dur);
        RespDir3_perPixel=nan(size(ch1a_crop,1),size(ch1a_crop,2),dur);
        RespDir4_perPixel=nan(size(ch1a_crop,1),size(ch1a_crop,2),dur);
        RT=[]; %nan((size(ch1a_crop(:,:,1)))); % will contain a Response Timing for each Pixel
        
        
        for p=1:size(ch1a_crop,1)
            for pp=1:size(ch1a_crop,2)
                
                if BW(p,pp)==1  %Calculate a DSI for each Pixel
                    
                    
                    presp=squeeze(ch1a_crop(p,pp,:));
                    
                    % average traces belonging to the same direction of movement
                    Stim=fstimval(ep_all(2,:)); %Stimulus direction during all epochs
                    
                    presp_dir=nan(avdur,4); %initiate matrix
                    
                    for kk=1:4 % loop through 4 directions of moving stim
                        IND=StimInd(kk);
                        Double=find(Stim==IND);
                        if length(Double)>1 %if I have more than 1 one repetition of one Stimulus direction...
                            allepochs=[];   %then create a matrix with all repetitions....
                            for k=1:length(Double)
                                allepochs=[allepochs, presp(ep_all(:,Double(k)))];  %...
                            end
                            presp_dir(:,kk)=mean(allepochs,2);  % ... and take the mean response of all repetitions
                        elseif length(Double)==0   % if any stimtype of 1 to 4 did not appear we did not record long enough
                            disp(['Error: No Stiumulus Type ' num2str(IND),' found, Stimulus length was too short!!!']);
                        else
                            presp_dir(:,kk)= presp(ep_all(:,Double));  % if just one repetition just take the response to this
                            
                        end
                    end
                    
                    %threshold=std([presp_dir(:,1);presp_dir(:,2);presp_dir(:,3);presp_dir(:,4)]); %threshold that defines if maximum response of Pixel counts as a 'Response'
                    threshold=std(presp);
                    MEAN=mean(presp);
                    % Check for responsivity
                    
                    
                    responsivity = zeros(1,4);    % Do I find a response of my Pixel to any of the Stim directions?
                    for u=1:4
                        %if max(presp_dir(:,u))> 5*std(presp_dir(:,u))
                        if max(presp_dir(:,u))> MEAN+timesSTD*threshold  % !!!!!!!This I have to validate if 3*std is good here!!!!!!!!!!!!!!
                            responsivity(u)=1;
                            [Max,RT_i]=max(presp_dir(:,u));
                            RT=[RT,RT_i];
                        end
                    end
                    
                    if sum(responsivity)>0 % Only if a response to moving edge is at least visible to one direction
                        % direction
                        if ManuallySelect ==1
                            Pdir(p,pp)=Pdir_select(p,pp);
                            PD=Pdir_select(p,pp);
                            MAXdir=max(presp_dir);
                            [~,rPD]=max(MAXdir);  % here PD=Pixels real PD (not based on Layer assignment)
                        elseif ManuallySelect ==0
                            MAXdir=max(presp_dir);
                            [MAX,PD]=max(MAXdir);
                            Pdir(p,pp)= PD;
                            rPD=PD;
                        end
                        
                        %DSI(p,pp)=(max(MAXdir)-min(MAXdir))/(max(MAXdir));
                        if     PD==1
                            ND=3;
                        elseif PD==2
                            ND=4;
                        elseif PD==3
                            ND=1;
                        elseif PD==4
                            ND=2;
                        end
                        
                        if PD>0
                            DSI(p,pp)=(MAXdir(PD)-MAXdir(ND))/MAXdir(PD);
                            [M,Timing]=max(presp_dir(:,PD));
                            TRespPD(p,pp)=Timing;
                            MaxRespPD(p,pp)=M;
                            
                            
                            stimuluspos=in.fstimpos1;
                            
                            Presp_nomov_new=presp_dir;
                            
                            RespDir1_perPixel(p,pp,:)=Presp_nomov_new(:,1);
                            RespDir2_perPixel(p,pp,:)=Presp_nomov_new(:,2);
                            RespDir3_perPixel(p,pp,:)=Presp_nomov_new(:,3);
                            RespDir4_perPixel(p,pp,:)=Presp_nomov_new(:,4);
                            
                            OFFedge=Presp_nomov_new(1:round(length(Presp_nomov_new)/2),:);
                            ONedge= Presp_nomov_new(round(length(Presp_nomov_new)/2+1):length(Presp_nomov_new),:);
                            
                            %\\\\\\\\\\\\\\\ Now separate pixels based on ON or OFF responses
                            % Write a matrix containing a 0 for OFF and 1 for ON response
                            
                            % Compare response to PD between ON and OFF
                            PDresp_OFF=OFFedge(:,rPD); % take out preferred direction response
                            PDresp_ON=ONedge(:,rPD);
                            
                            if max(PDresp_OFF)>max(PDresp_ON) %OFF pixel
                                ONorOFFpixel(p,pp)=0;
                                CSI(p,pp)=(max(PDresp_OFF)-max(PDresp_ON))/max(PDresp_OFF);
                                
                            elseif max(PDresp_OFF)< max(PDresp_ON)
                                ONorOFFpixel(p,pp)=1;
                                CSI(p,pp)=(max(PDresp_ON)-max(PDresp_OFF))/max(PDresp_ON);
                                
                            end
                            
                            %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                            
                        end
                        
                    end
                end
            end
        end
        
        
        
        
        
        
        
        
        %% Now we use the Components that decribe each Pixel Response (calculated above)
        % DSI = Direction Selectivity Index
        % CSI = Contrast Selectivity Index
        % Pdir = Preferred Direction
        % TRespPD = Timing of Response to PD
        % MaxRespPD = Maximum of Response to PD
        % ONorOFFpixel = T4 ot T5 cell
        % X &
        % Y position
        
        % to define a valuable size of ROIs fitting to the visual field of one T4/T5 cell
        % we need to know the pixelsize\resolution of our images:
        
        micronsX=round(str2num(in.xml.micronsPerPixel.XAxis),4);
        micronsY=round(str2num(in.xml.micronsPerPixel.YAxis),4);
        
        % Define ROI size
        if micronsX==micronsY
            N_Pi_perROI= round(ROISize/micronsX)* round(ROISize/micronsY);
            %looking at confocal images of terminals it looks like they are
            %between 5 and 10um big! 
            %I chose 2.5 because I want a minimal area of approximately 6um 
            
        else
            disp('Pixel resolution in X and Y is not equal ')
            
        end
        
        % calculate cluster analysis separated for Layers and ON/OFF pixel
        
        
        A=DSI>DSIthreshold;
        C=CSI>CSIthreshold;
        B_L1=(Pdir==2); %right
        B_L2=(Pdir==4); %left
        B_L3=(Pdir==1); %Up
        B_L4=(Pdir==3); %Down
        D_ON=(ONorOFFpixel==1);
        D_OFF=(ONorOFFpixel==0);
        
        
        ALL_L1_ON=A.*C.*B_L1.*D_ON;
        [relevantPixel_X_L1_ON,relevantPixel_Y_L1_ON]=find(ALL_L1_ON);
        
        ALL_L2_ON=A.*C.*B_L2.*D_ON;
        [relevantPixel_X_L2_ON,relevantPixel_Y_L2_ON]=find(ALL_L2_ON);
        
        ALL_L3_ON=A.*C.*B_L3.*D_ON;
        [relevantPixel_X_L3_ON,relevantPixel_Y_L3_ON]=find(ALL_L3_ON);
        
        ALL_L4_ON=A.*C.*B_L4.*D_ON;
        [relevantPixel_X_L4_ON,relevantPixel_Y_L4_ON]=find(ALL_L4_ON);
        
        
        ALL_L1_OFF=A.*C.*B_L1.*D_OFF;
        [relevantPixel_X_L1_OFF,relevantPixel_Y_L1_OFF]=find(ALL_L1_OFF);
        
        ALL_L2_OFF=A.*C.*B_L2.*D_OFF;
        [relevantPixel_X_L2_OFF,relevantPixel_Y_L2_OFF]=find(ALL_L2_OFF);
        
        ALL_L3_OFF=A.*C.*B_L3.*D_OFF;
        [relevantPixel_X_L3_OFF,relevantPixel_Y_L3_OFF]=find(ALL_L3_OFF);
        
        ALL_L4_OFF=A.*C.*B_L4.*D_OFF;
        [relevantPixel_X_L4_OFF,relevantPixel_Y_L4_OFF]=find(ALL_L4_OFF);
        
        % -------------------------------------------------------------------------
        
        Components_L1_OFF=nan(length(relevantPixel_X_L1_OFF),3);
        Components_L1_OFF(:,1)= relevantPixel_X_L1_OFF;
        Components_L1_OFF(:,2)= relevantPixel_Y_L1_OFF;
        for ll=1:length(relevantPixel_X_L1_OFF)
            Components_L1_OFF(ll,3)=TRespPD(relevantPixel_X_L1_OFF(ll),relevantPixel_Y_L1_OFF(ll));
        end
        
        Components_L2_OFF=nan(length(relevantPixel_X_L2_OFF),3);
        Components_L2_OFF(:,1)= relevantPixel_X_L2_OFF;
        Components_L2_OFF(:,2)= relevantPixel_Y_L2_OFF;
        for ll=1:length(relevantPixel_X_L2_OFF)
            Components_L2_OFF(ll,3)=TRespPD(relevantPixel_X_L2_OFF(ll),relevantPixel_Y_L2_OFF(ll));
        end
        
        Components_L3_OFF=nan(length(relevantPixel_X_L3_OFF),3);
        Components_L3_OFF(:,1)= relevantPixel_X_L3_OFF;
        Components_L3_OFF(:,2)= relevantPixel_Y_L3_OFF;
        for ll=1:length(relevantPixel_X_L3_OFF)
            Components_L3_OFF(ll,3)=TRespPD(relevantPixel_X_L3_OFF(ll),relevantPixel_Y_L3_OFF(ll));
        end
        
        Components_L4_OFF=nan(length(relevantPixel_X_L4_OFF),3);
        Components_L4_OFF(:,1)= relevantPixel_X_L4_OFF;
        Components_L4_OFF(:,2)= relevantPixel_Y_L4_OFF;
        for ll=1:length(relevantPixel_X_L4_OFF)
            Components_L4_OFF(ll,3)=TRespPD(relevantPixel_X_L4_OFF(ll),relevantPixel_Y_L4_OFF(ll));
        end
        
        Components_L1_ON=nan(length(relevantPixel_X_L1_ON),3);
        Components_L1_ON(:,1)= relevantPixel_X_L1_ON;
        Components_L1_ON(:,2)= relevantPixel_Y_L1_ON;
        for ll=1:length(relevantPixel_X_L1_ON)
            Components_L1_ON(ll,3)=TRespPD(relevantPixel_X_L1_ON(ll),relevantPixel_Y_L1_ON(ll));
        end
        
        Components_L2_ON=nan(length(relevantPixel_X_L2_ON),3);
        Components_L2_ON(:,1)= relevantPixel_X_L2_ON;
        Components_L2_ON(:,2)= relevantPixel_Y_L2_ON;
        for ll=1:length(relevantPixel_X_L2_ON)
            Components_L2_ON(ll,3)=TRespPD(relevantPixel_X_L2_ON(ll),relevantPixel_Y_L2_ON(ll));
        end
        
        Components_L3_ON=nan(length(relevantPixel_X_L3_ON),3);
        Components_L3_ON(:,1)= relevantPixel_X_L3_ON;
        Components_L3_ON(:,2)= relevantPixel_Y_L3_ON;
        for ll=1:length(relevantPixel_X_L3_ON)
            Components_L3_ON(ll,3)=TRespPD(relevantPixel_X_L3_ON(ll),relevantPixel_Y_L3_ON(ll));
        end
        
        Components_L4_ON=nan(length(relevantPixel_X_L4_ON),3);
        Components_L4_ON(:,1)= relevantPixel_X_L4_ON;
        Components_L4_ON(:,2)= relevantPixel_Y_L4_ON;
        for ll=1:length(relevantPixel_X_L4_ON)
            Components_L4_ON(ll,3)=TRespPD(relevantPixel_X_L4_ON(ll),relevantPixel_Y_L4_ON(ll));
        end
        
        %
        % ClusterAnalysis is a function that calculates the best suitable 'maxclust'variable
        % based on the defined Size of the ROIs (size of axon Terminals)
        if size(Components_L1_OFF,1)>1
            ValidClusters_L1_OFF= ClusterAnalysis_3(Components_L1_OFF,N_Pi_perROI,in);
        else
            ValidClusters_L1_OFF=[];
            disp('Not enough L1_OFF Pixels for cluster Analysis')
        end
        
        if size(Components_L2_OFF,1)>1
            ValidClusters_L2_OFF= ClusterAnalysis_3(Components_L2_OFF,N_Pi_perROI,in);
        else
            ValidClusters_L2_OFF=[];
            disp('Not enough L2_OFF Pixels for cluster Analysis')
        end
        
        if size(Components_L3_OFF,1)>1
            ValidClusters_L3_OFF= ClusterAnalysis_3(Components_L3_OFF,N_Pi_perROI,in);
        else
            ValidClusters_L3_OFF=[];
            disp('Not enough L3_OFF Pixels for cluster Analysis')
        end
        
        if size(Components_L4_OFF,1)>1
            ValidClusters_L4_OFF= ClusterAnalysis_3(Components_L4_OFF,N_Pi_perROI,in);
        else
            ValidClusters_L4_OFF=[];
            disp('Not enough L1_OFF Pixels for cluster Analysis')
        end
        
        if size(Components_L1_ON,1)>1
            ValidClusters_L1_ON= ClusterAnalysis_3(Components_L1_ON,N_Pi_perROI,in);
        else
            ValidClusters_L1_ON=[];
            disp('Not enough L1_ON Pixels for cluster Analysis')
        end
        
        if size(Components_L2_ON,1)>1
            ValidClusters_L2_ON= ClusterAnalysis_3(Components_L2_ON,N_Pi_perROI,in);
        else
            ValidClusters_L2_ON=[];
            disp('Not enough L2_ON Pixels for cluster Analysis')
        end
        
        if size(Components_L3_ON,1)>1
            ValidClusters_L3_ON= ClusterAnalysis_3(Components_L3_ON,N_Pi_perROI,in);
        else
            ValidClusters_L3_ON=[];
            disp('Not enough L3_ON Pixels for cluster Analysis')
        end
        
        if size(Components_L4_ON,1)>1
            ValidClusters_L4_ON= ClusterAnalysis_3(Components_L4_ON,N_Pi_perROI,in);
        else
            ValidClusters_L4_ON=[];
            disp('Not enough L4_ON Pixels for cluster Analysis')
        end
        
        
        
        %% Plots for T5 (OFF Pixels)

%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L2_OFF)
%             
%             for i=1:length(ValidClusters_L2_OFF(:,1))
%                 ClusterType=ValidClusters_L2_OFF(i,3);
%                 Try1(ValidClusters_L2_OFF(i,1),ValidClusters_L2_OFF(i,2))=ClusterType;
%             end
%             masks_L2_OFF=cell(1,max(ValidClusters_L2_OFF(:,3)));
%             
%             for ii=1:max(ValidClusters_L2_OFF(:,3))
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L2_OFF{1,ii}=mask;
%             end
%             
%             
%             figure('Color', [1 1 1])
%             subplot(2,2,2)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 2 - T5 ROIs')
%         else
%             
%             figure('Color', [1 1 1])
%             subplot(2,2,2)
%             text(0.2,0.5,'No valid clusters found')
%             masks_L2_OFF=[];
%             
%         end
%         %-----------------------------------------
%         
%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L4_OFF);
%             for i=1:length(ValidClusters_L4_OFF(:,1));
%                 ClusterType=ValidClusters_L4_OFF(i,3);
%                 Try1(ValidClusters_L4_OFF(i,1),ValidClusters_L4_OFF(i,2))=ClusterType;
%             end
%             masks_L4_OFF=cell(1,max(ValidClusters_L4_OFF(:,3)));
%             
%             for ii=1:max(ValidClusters_L4_OFF(:,3))
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L4_OFF{1,ii}=mask;
%             end
%             
%             subplot(2,2,4)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 4 - T5 ROIs')
%             
%         else
%             subplot(2,2,4)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L4_OFF=[];
%             
%             
%         end
%         
%         %-----------------------------------------
%         
%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L1_OFF);
%             
%             for i=1:length(ValidClusters_L1_OFF(:,1));
%                 ClusterType=ValidClusters_L1_OFF(i,3);
%                 Try1(ValidClusters_L1_OFF(i,1),ValidClusters_L1_OFF(i,2))=ClusterType;
%             end
%             masks_L1_OFF=cell(1,max(ValidClusters_L1_OFF(:,3)));
%             
%             for ii=1:max(ValidClusters_L1_OFF(:,3));
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L1_OFF{1,ii}=mask;
%             end
%             
%             subplot(2,2,1)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 1 - T5 ROIs')
%         else
%             subplot(2,2,1)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L1_OFF=[];
%             
%         end
%         %-----------------------------------------
%         
%         if ~isempty(ValidClusters_L3_OFF);
%             
%             Try1=nan(size(DSI));
%             for i=1:length(ValidClusters_L3_OFF(:,1));
%                 ClusterType=ValidClusters_L3_OFF(i,3);
%                 Try1(ValidClusters_L3_OFF(i,1),ValidClusters_L3_OFF(i,2))=ClusterType;
%             end
%             masks_L3_OFF=cell(1,max(ValidClusters_L3_OFF(:,3)));
%             
%             for ii=1:max(ValidClusters_L3_OFF(:,3));
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L3_OFF{1,ii}=mask;
%             end
%             
%             
%             subplot(2,2,3)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 3 - T5 ROIs')
%         else
%             subplot(2,2,4)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L3_OFF=[];
%             
%             
%         end
%         %-----------------------------------------
%         
%         Try1=ones(size(DSI));
%         if ~isempty(ValidClusters_L2_OFF);
%             for i=1:length(ValidClusters_L2_OFF(:,1));
%                 ClusterType=ValidClusters_L2_OFF(i,3);
%                 Try1(ValidClusters_L2_OFF(i,1),ValidClusters_L2_OFF(i,2))=2;
%             end
%         end
%         
%         if ~isempty(ValidClusters_L4_OFF);
%             for i=1:length(ValidClusters_L4_OFF(:,1));
%                 ClusterType=ValidClusters_L4_OFF(i,3);
%                 Try1(ValidClusters_L4_OFF(i,1),ValidClusters_L4_OFF(i,2))=3;
%             end
%         end
%         
%         if ~isempty(ValidClusters_L1_OFF);
%             for i=1:length(ValidClusters_L1_OFF(:,1));
%                 ClusterType=ValidClusters_L1_OFF(i,3);
%                 Try1(ValidClusters_L1_OFF(i,1),ValidClusters_L1_OFF(i,2))=4;
%             end
%         end
%         
%         if ~isempty(ValidClusters_L3_OFF);
%             for i=1:length(ValidClusters_L3_OFF(:,1));
%                 ClusterType=ValidClusters_L3_OFF(i,3);
%                 Try1(ValidClusters_L3_OFF(i,1),ValidClusters_L3_OFF(i,2))=5;
%             end
%         end
%         
%         figure('Color', [1 1 1])
%         imagesc(Try1)
%         cm=colormap;
%         cm(1,:)=[1,1,1];
%         cm(17,:)=cm(7,:);
%         cm(33,:)=cm(60,:);
%         cm(49,:)=[49,163,84]/255;
%         cm(end,:)=[222,45,38]/255;
%         colormap(cm)
%         title('T5 cluster/ROIs')
        
%         %% Plot for T4 (ON pixels)
%         
%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L2_ON)
%             
%             for i=1:length(ValidClusters_L2_ON(:,1))
%                 ClusterType=ValidClusters_L2_ON(i,3);
%                 Try1(ValidClusters_L2_ON(i,1),ValidClusters_L2_ON(i,2))=ClusterType;
%             end
%             masks_L2_ON=cell(1,max(ValidClusters_L2_ON(:,3)));
%             
%             for ii=1:max(ValidClusters_L2_ON(:,3))
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L2_ON{1,ii}=mask;
%             end
%             
%             figure('Color', [1 1 1])
%             subplot(2,2,2)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 2 - T4 ROIs')
%         else
%             figure('Color', [1 1 1])
%             subplot(2,2,2)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L2_ON=[];
%             
%             
%         end
%         %-----------------------------------------
%         
%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L4_ON)
%             for i=1:length(ValidClusters_L4_ON(:,1))
%                 ClusterType=ValidClusters_L4_ON(i,3);
%                 Try1(ValidClusters_L4_ON(i,1),ValidClusters_L4_ON(i,2))=ClusterType;
%             end
%             masks_L4_ON=cell(1,max(ValidClusters_L4_ON(:,3)));
%             
%             for ii=1:max(ValidClusters_L4_ON(:,3))
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L4_ON{1,ii}=mask;
%             end
%             
%             subplot(2,2,4)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 4 - T4 ROIs')
%         else
%             subplot(2,2,4)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L4_ON=[];
%             
%             
%         end
%         %-----------------------------------------
%         
%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L1_ON)
%             for i=1:length(ValidClusters_L1_ON(:,1))
%                 ClusterType=ValidClusters_L1_ON(i,3);
%                 Try1(ValidClusters_L1_ON(i,1),ValidClusters_L1_ON(i,2))=ClusterType;
%             end
%             
%             masks_L1_ON=cell(1,max(ValidClusters_L1_ON(:,3)));
%             for ii=1:max(ValidClusters_L1_ON(:,3))
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L1_ON{1,ii}=mask;
%             end
%             
%             subplot(2,2,1)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 1 - T4 ROIs')
%             
%         else
%             subplot(2,2,1)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L1_ON=[];
%             
%             
%         end
%         %-----------------------------------------
%         
%         
%         Try1=nan(size(DSI));
%         if ~isempty(ValidClusters_L3_ON)
%             for i=1:length(ValidClusters_L3_ON(:,1))
%                 ClusterType=ValidClusters_L3_ON(i,3);
%                 Try1(ValidClusters_L3_ON(i,1),ValidClusters_L3_ON(i,2))=ClusterType;
%             end
%             
%             masks_L3_ON=cell(1,max(ValidClusters_L3_ON(:,3)));
%             
%             for ii=1:max(ValidClusters_L3_ON(:,3))
%                 mask=zeros(size(DSI));
%                 mask(find(Try1==ii))=1;
%                 masks_L3_ON{1,ii}=mask;
%             end
%             
%             subplot(2,2,3)
%             imagesc(Try1)
%             cm=colormap('colorcube');
%             cm(1,:)=[1,1,1];
%             colormap(cm)
%             title('Layer 3 - T4 ROIs')
%             
%         else
%             subplot(2,2,3)
%             text(0.2, 0.5, 'No valid clusters found')
%             masks_L3_ON=[];
%             
%             
%         end
%         %-----------------------------------------
%         
%         Try1=ones(size(DSI));
%         
%         if ~isempty(ValidClusters_L2_ON)
%             for i=1:length(ValidClusters_L2_ON(:,1))
%                 ClusterType=ValidClusters_L2_ON(i,3);
%                 Try1(ValidClusters_L2_ON(i,1),ValidClusters_L2_ON(i,2))=2;
%             end
%         end
%         
%         if ~isempty(ValidClusters_L4_ON)
%             for i=1:length(ValidClusters_L4_ON(:,1))
%                 ClusterType=ValidClusters_L4_ON(i,3);
%                 Try1(ValidClusters_L4_ON(i,1),ValidClusters_L4_ON(i,2))=3;
%             end
%         end
%         
%         if ~isempty(ValidClusters_L1_ON)
%             for i=1:length(ValidClusters_L1_ON(:,1))
%                 ClusterType=ValidClusters_L1_ON(i,3);
%                 Try1(ValidClusters_L1_ON(i,1),ValidClusters_L1_ON(i,2))=4;
%             end
%         end
%         
%         if ~isempty(ValidClusters_L3_ON)
%             for i=1:length(ValidClusters_L3_ON(:,1))
%                 ClusterType=ValidClusters_L3_ON(i,3);
%                 Try1(ValidClusters_L3_ON(i,1),ValidClusters_L3_ON(i,2))=5;
%             end
%         end
%         
%         figure('Color', [1 1 1])
%         imagesc(Try1)
%         cm=colormap;
%         cm(1,:)=[1,1,1];
%         cm(17,:)=cm(7,:);
%         cm(33,:)=cm(60,:);
%         cm(49,:)=[49,163,84]/255;
%         cm(end,:)=[222,45,38]/255;
%         colormap(cm);
%         title('T4 cluster/ ROIs')
%         
%         
%         
%         % plot results
%         figure('Color', [1 1 1])
%         Pdir_im=Pdir;
%         Pdir_im(find(isnan(Pdir_im)==1))=0;
%         
%         h1=subplot(1,3,1)
%         imagesc(Pdir_im)
%         cm3=colormap;
%         cm3(1,:)=[1,1,1];
%         cm3(17,:)=[222,45,38]/255;
%         cm3(33,:)=[49,163,84]/255;
%         cm3(49,:)=cm3(60,:);
%         cm3(end,:)=cm3(7,:);
%         
%         title('Preferred Direction')
%         
%         h2=subplot(1,3,2)
%         imagesc(DSI)
%         title('DSI')
%         
%         
%         h3=subplot(1,3,3)
%         imagesc(CSI)
%         title('CSI')
%         
%         colormap(h1,cm3)
%         
%         %(range = (1+25):(length(t)-60); taken for dark edge presentation taken
%         %from movedges_figure_v1.m)
%         %Only take pixels with a DSI bigger than 0.3
%         
%         P_DSItoosmall=find(DSI<DSIthreshold);
%         Pdirthresh=Pdir_im;
%         Pdirthresh(P_DSItoosmall)=0;
%         figure('Color', [1 1 1])
%         imagesc(Pdirthresh)
%         title('Preferred Direction- after DSI threshold')
%         colormap(cm3);
%         
%         % Plot ONOFFmatrix
%         
%         ONorOFFpixel_im=ONorOFFpixel;
%         
%         ONorOFFpixel_im(find(isnan(ONorOFFpixel_im)==1))=3; % To have a different color for the background only for plotting!!!
%         figure('Color', [1 1 1]);
%         subplot(1,2,1)
%         imagesc(ONorOFFpixel_im);
%         title('Pixel assigned to ON (red) and OFF=(black) -- before thresholding ')
%         
%         P_DSItoosmall=find(DSI<DSIthreshold);
%         PONOFFthresh=ONorOFFpixel_im;
%         PONOFFthresh(P_DSItoosmall)=3;
%         
%         subplot(1,2,2)
%         imagesc(PONOFFthresh)
%         title(['-- after thresholding: DSI above ',num2str(DSIthreshold)])
%         
%         cm2=colormap;
%         cm2(1,:)=[0,0,0];
%         cm2(22,:)=[222,45,38]/255;
%         cm2(end,:)=[1,1,1];
%         colormap(cm2);
%         
%         
%         %Plot ONOFF matrix with additional CSI threshold
%         
%         figure('Color', [1 1 1]);
%         subplot(1,2,1)
%         imagesc(PONOFFthresh)
%         title(['Pixel assigned to ON (red) and OFF=(black)-- Before CSI thresholding: CSI above '])
%         
%         subplot(1,2,2)
%         
%         
%         P_CSItoosmall=find(CSI<CSIthreshold);
%         PONOFFthresh_CSI=PONOFFthresh;
%         PONOFFthresh_CSI(P_CSItoosmall)=3;
%         imagesc(PONOFFthresh_CSI)
%         title(['-- After CSI thresholding: CSI above ',num2str(CSIthreshold)])
%         colormap(cm2);
%         
%         
%         
        
        %% Put data together:
        
        nMasks=length(masks_L1_OFF)+length(masks_L2_OFF)+length(masks_L3_OFF)+length(masks_L4_OFF)+length(masks_L1_ON)+length(masks_L2_ON)+length(masks_L3_ON)+length(masks_L4_ON);
        Layer=ones(length(masks_L1_OFF),1);
        Layer=[Layer;2*ones(length(masks_L2_OFF),1)];
        Layer=[Layer;3*ones(length(masks_L3_OFF),1)];
        Layer=[Layer;4*ones(length(masks_L4_OFF),1)];
        
        T4_T5=5*ones(length(Layer),1);
        L=length(T4_T5);
        
        Layer=[Layer;1*ones(length(masks_L1_ON),1)];
        Layer=[Layer;2*ones(length(masks_L2_ON),1)];
        Layer=[Layer;3*ones(length(masks_L3_ON),1)];
        Layer=[Layer;4*ones(length(masks_L4_ON),1)];
        
        T4_T5=[T4_T5;4*ones(length(Layer)-L,1)];
        
        masked = zeros(str2double(in.xml.linesPerFrame),str2double(in.xml.pixelsPerLine));
        masks=[masks_L1_OFF,masks_L2_OFF,masks_L3_OFF,masks_L4_OFF,masks_L1_ON,masks_L2_ON,masks_L3_ON,masks_L4_ON]; % Contains the mask for each cluster
        %(which is a matrix with the dimensions of the recorded image, containing zeros for the background and ones for the position of the ROI
        smask = zeros(nMasks,1);
        for k = 1:nMasks   % Contains the number of Pixels for each ROI
            smask(k) = sum(sum(masks{k})); %Since masks contains zeros for backround and ones for pixels that are part of the ROI,
            % the sum of rows and columns will give the amount of pixels
            % describing the ROI
            masksi{k}= find(masks{k});  % will find the Position of ones in the masks for each ROI, so the x and y values of Pixels that are Part of the ROI
        end
        
        
             
        sNmask = sum(sum(Real_Background)); 
        Nmaski = find(Real_Background);
        avSignal1_CA = zeros(nMasks,nframes);
        dSignal1_CA = zeros(nMasks,nframes);
        
        for ind = 1:nframes  %%% collect data for raw responses (without averaged across stimulus epochs
            
            A = double(squeeze(ch1a_crop(:,:,ind)));
            B = double(squeeze(ch1a(:,:,ind)));

            for k = 1:nMasks
                
                masked = A(masksi{k});
                Nmasked = B((Nmaski));
                avSignal1_CA(k,ind) = sum(masked)./smask(k); % summed signal in a ROI, normalized by ROI size
                dSignal1_CA(k,ind) = avSignal1_CA(k,ind) - sum((Nmasked))./sNmask; % background subtraction (by signal in background normalized by background ROI size)
                
            end
        end
        
        %% Save important data in out structure :output of my function
        %New put all Information together into one structure
        % Save one struct. for Layers manually selected and one normal
        % cluster
        
        out = in;
        
        ClusterInfo.Layer=Layer;  %Contains the Layer assignment for each Cluster(ROI) with values between 1 and 4 according to 4 Layers in LP
        ClusterInfo.T4_T5=T4_T5;  %Contains the T4 or T5 assignment with values of 4 and 5 accordingly
        ClusterInfo.nMasks_CA=nMasks; %Number of masks/clusters/ROIs
        
        ClusterInfo.avSignal1_CA = avSignal1_CA;
        ClusterInfo.dSignal1_CA = dSignal1_CA;
        ClusterInfo.Epochs=ep_all;
        ClusterInfo.masks_CA=masks;
        

        
        ClusterInfo.DSI=DSIthreshold;
        ClusterInfo.CSI=CSIthreshold;
        ClusterInfo.timesSTD=timesSTD;
        ClusterInfo.Real_Background=Real_Background;

        
        % save masks of ROIs/individual clusters
        
        curDir = pwd;
        
        
        
        if ManuallySelect
            
            out.ClusterInfo_ManuallySelect=ClusterInfo;
             save('CA_information_ManuallySelect','masks', 'Layer', 'T4_T5','BW','BW2','DSI', 'CSI', 'timesSTD', 'crop','Real_Background');
        else
            out.ClusterInfo=ClusterInfo;
             save('CA_information','masks', 'Layer', 'T4_T5','BW','BW2','DSI', 'CSI', 'timesSTD','crop','Real_Background');
            
        end
        
        disp(('saved curMasks_CA'));
    else % If GoON=0
        
        disp(['The data are not recorded with the correct stimulus type. You need: DriftingEdge_LumDecLumInc_1period_20degPerSec_90deg_BlankRand_4Dirs! Your stimulus used is: ', in.stim_type])
        out=in;
        
    end %OF GoON case
    
%%    

elseif mode==2
    
    curDir = pwd;
    cd(curDir);
    dd1=dir('CA_information.mat');
    dd2=dir('CA_information_ManuallySelect.mat');
    dd=dir('CA_information*');
    if size(dd,1)==1
        load(dd.name);
        if size(dd1)
            x=1;
        elseif size(dd2)
            x=2;
        end
    elseif size(dd,1)>1
        
        x=inputdlg(['Two Information available: 1)"', dd(1).name, '" and 2)"', dd(2).name, '" --> Choose 1 or 2 :' ]);
        x=str2num(x{1,1});
        load(dd(x).name);
    end
    
    cd(curDir);
    % Generate ratio signals from all regions of interest - aligned data
    out = in;
    
    ch1a_crop=ch1a(crop(2):crop(2)+crop(4),crop(1):crop(1)+crop(3),:); 

    nMasks=length(masks);
    
    ClusterInfo.Layer=Layer;
    ClusterInfo.T4_T5=T4_T5;
    ClusterInfo.nMasks_CA=nMasks;
    ClusterInfo.avSignal1_CA = zeros(nMasks,nframes);
    ClusterInfo.dSignal1_CA = zeros(nMasks,nframes);
    ClusterInfo.DSI=DSI;
    ClusterInfo.CSI=CSI;
    ClusterInfo.timesSTD=timesSTD;
   
    smask = zeros(nMasks,1);
    
    for k = 1:nMasks
        smask(k) = sum(sum(masks{k}));
        masksi{k}= find(masks{k});
    end
    
    ClusterInfo.masks_CA=masks;
    
   
    
    sNmask = sum(sum(Real_Background)); 
    Nmaski = find(Real_Background);
    
    ClusterInfo.Real_Background=Real_Background;

    
    for ind = 1:nframes  %%% collect data for raw responses (without averaged across stimulus epochs
            
            A = double(squeeze(ch1a_crop(:,:,ind)));
            B = double(squeeze(ch1a(:,:,ind)));

            for k = 1:nMasks
                
                masked = A(masksi{k});
                Nmasked = B((Nmaski));
                ClusterInfo.avSignal1_CA(k,ind) = sum(masked)./smask(k); % summed signal in a ROI, normalized by ROI size
                ClusterInfo.dSignal1_CA(k,ind) = ClusterInfo.avSignal1_CA(k,ind) - sum((Nmasked))./sNmask; % background subtraction (by signal in background normalized by background ROI size)
                
            end
   end
    
    
    
    
    
    
    if x==1
        out.ClusterInfo=ClusterInfo;
    elseif x==2
        out.ClusterInfo_ManuallySelect=ClusterInfo;
    end
    
    
end  % if mode

    out.ch1a_crop=ch1a_crop;
    out.crop=crop; 
    
    cd(curDir)

end

