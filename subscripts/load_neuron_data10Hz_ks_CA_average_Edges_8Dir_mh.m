function out=load_neuron_data10Hz_ks_CA_average_Edges_8Dir_mh(in,locdir,DataType)

[fname,datetime,frames,zdepth,tottime,stimbouts,locnum,...
    activecells,inverted,stimcode,quality,driver,...
     moving,layer,wavelength,flyID, responsivefly,Layercheck,SingleFly,PTX_5_100,Washout5mM,Washout100mM]...
    =textread('MASTER_foldersummary_Miriam.txt',...
    '%s %s %f %f %f %f %f %s %f %s %f %s %f %s %f %f %f %f %f %f %f %f\r',...
    'headerlines',1,'delimiter','\t');

currd=pwd;
cd(locdir);


count=1;
for iFLY=1:length(in) % for each Fly
    
   FLYname=fname{in(iFLY)};

        
        x=load(FLYname);
        DATA=eval(['x.strct.',DataType]);
        
            % FIRST AVERAGE ALL CLUSTERS FOR T4/T5 EACH LAYER
            T5_A= find((DATA.T4_T5==5) .* (DATA.Layer==1));
            T5_B= find((DATA.T4_T5==5) .* (DATA.Layer==2));
            T5_C= find((DATA.T4_T5==5) .* (DATA.Layer==3));
            T5_D= find((DATA.T4_T5==5) .* (DATA.Layer==4));
            
            T4_A= find((DATA.T4_T5==4) .* (DATA.Layer==1));
            T4_B= find((DATA.T4_T5==4) .* (DATA.Layer==2));
            T4_C= find((DATA.T4_T5==4) .* (DATA.Layer==3));
            T4_D= find((DATA.T4_T5==4) .* (DATA.Layer==4));
            
            NROIS.T5_A=length(T5_A);
            NROIS.T5_B=length(T5_B);
            NROIS.T5_C=length(T5_C);
            NROIS.T5_D=length(T5_D);
            NROIS.T4_A=length(T4_A);
            NROIS.T4_B=length(T4_B);
            NROIS.T4_C=length(T4_C);
            NROIS.T4_D=length(T4_D);   
            
            Teta_deg=[90,45, 0, 315, 270, 225, 180, 135];
            Teta_rad=(Teta_deg*pi)/180;
            half=size(DATA.AV_ROIS_resp{1,1},2)/2;
            END=size(DATA.AV_ROIS_resp{1,1},2);
            %Get T5 Data
            L_dir_T5A=nan(1,length(T5_A));
            T5A_ALL_ROIS_resp=[];
            T5A_masks={};
            T5A_max={};
            for NROI=1:length(T5_A) 
             T5A_ALL_ROIS_resp=cat(3,T5A_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_A(NROI)});
             R_teta=max(T5A_ALL_ROIS_resp(:,1:half,NROI),[],2)';
             L_dir_T5A(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5A_masks{NROI}=DATA.masks_CA{T5_A(NROI)};
             [MAX,PD]=max(R_teta);
             T5A_max{NROI}.MAX=MAX;
             T5A_max{NROI}.PD=Teta_deg(PD);
            end  
            M_T5A=mean(T5A_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
            
            
            L_dir_T5B=nan(1,length(T5_B));
            T5B_ALL_ROIS_resp=[];
            T5B_masks={};
            T5B_max={};
            for NROI=1:length(T5_B)
             T5B_ALL_ROIS_resp=cat(3,T5B_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_B(NROI)});
             R_teta=max(T5B_ALL_ROIS_resp(:,1:half,NROI),[],2)';
             L_dir_T5B(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5B_masks{NROI}=DATA.masks_CA{T5_B(NROI)};
             [MAX,PD]=max(R_teta);
             T5B_max{NROI}.MAX=MAX;
             T5B_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T5B=mean(T5B_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            
            L_dir_T5C=nan(1,length(T5_C));
            T5C_ALL_ROIS_resp=[];
            T5C_masks={};
            T5C_max={};
            for NROI=1:length(T5_C)
             T5C_ALL_ROIS_resp=cat(3,T5C_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_C(NROI)});
             R_teta=max(T5C_ALL_ROIS_resp(:,1:half,NROI),[],2)';
             L_dir_T5C(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5C_masks{NROI}=DATA.masks_CA{T5_C(NROI)};
             [MAX,PD]=max(R_teta);
             T5C_max{NROI}.MAX=MAX;
             T5C_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T5C=mean(T5C_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            L_dir_T5D=nan(1,length(T5_D));
             T5D_ALL_ROIS_resp=[];
             T5D_masks={};
             T5D_max={};
            for NROI=1:length(T5_D)
             T5D_ALL_ROIS_resp=cat(3,T5D_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_D(NROI)});
             R_teta=max(T5D_ALL_ROIS_resp(:,1:half,NROI),[],2)';
             L_dir_T5D(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5D_masks{NROI}=DATA.masks_CA{T5_D(NROI)};
             [MAX,PD]=max(R_teta);
             T5D_max{NROI}.MAX=MAX;
             T5D_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T5D=mean(T5D_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
%            
%             figure; 
%             compass(L_dir_T5A,'g')
%             hold on 
%             compass(L_dir_T5B,'b')
%             compass(L_dir_T5C,'r')
%             compass(L_dir_T5D,'y')

            

            
            
            %Get T4 Data
            L_dir_T4A=nan(1,length(T4_A));
            T4A_ALL_ROIS_resp=[];
            T4A_masks={};
            T4A_max={};
            for NROI=1:length(T4_A) 
             T4A_ALL_ROIS_resp=cat(3,T4A_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_A(NROI)});
             R_teta=max(T4A_ALL_ROIS_resp(:,half+1:END,NROI),[],2)';
             L_dir_T4A(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4A_masks{NROI}=DATA.masks_CA{T4_A(NROI)};
             [MAX,PD]=max(R_teta);
             T4A_max{NROI}.MAX=MAX;
             T4A_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T4A=mean(T4A_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
            
            
            L_dir_T4B=nan(1,length(T4_B));
            T4B_ALL_ROIS_resp=[]; 
            T4B_masks={};
            T4B_max={};
            for NROI=1:length(T4_B)
             T4B_ALL_ROIS_resp=cat(3,T4B_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_B(NROI)});
             R_teta=max(T4B_ALL_ROIS_resp(:,half+1:END,NROI),[],2)';
             L_dir_T4B(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4B_masks{NROI}=DATA.masks_CA{T4_B(NROI)};
             [MAX,PD]=max(R_teta);
             T4B_max{NROI}.MAX=MAX;
             T4B_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T4B=mean(T4B_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
             
            L_dir_T4C=nan(1,length(T4_C));
            T4C_ALL_ROIS_resp=[]; 
            T4C_masks={};
            T4C_max={};
            for NROI=1:length(T4_C)
             T4C_ALL_ROIS_resp=cat(3,T4C_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_C(NROI)});
             R_teta=max(T4C_ALL_ROIS_resp(:,half+1:END,NROI),[],2)';
             L_dir_T4C(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4C_masks{NROI}=DATA.masks_CA{T4_C(NROI)};
             [MAX,PD]=max(R_teta);
             T4C_max{NROI}.MAX=MAX;
             T4C_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T4C=mean(T4C_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            
             L_dir_T4D=nan(1,length(T4_D));
             T4D_ALL_ROIS_resp=[]; 
             T4D_masks={};
             T4D_max={};
            for NROI=1:length(T4_D)
             T4D_ALL_ROIS_resp=cat(3,T4D_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_D(NROI)});
             R_teta=max(T4D_ALL_ROIS_resp(:,half+1:END,NROI),[],2)';
             L_dir_T4D(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4D_masks{NROI}=DATA.masks_CA{T4_D(NROI)};
             [MAX,PD]=max(R_teta);
             T4D_max{NROI}.MAX=MAX;
             T4D_max{NROI}.PD=Teta_deg(PD);
            end        
            M_T4D=mean(T4D_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            
            % Calc the vector average for each ROI using the formula from
            % Mazurek et al., 2014
%             figure; 
%             compass(L_dir_T4A,'g')
%             hold on 
%             compass(L_dir_T4B,'b')
%             compass(L_dir_T4C,'r')
%             compass(L_dir_T4D,'y')
            
            
%             figure; 
%              compass(mean(L_dir_T4C),'r')
%              hold on 
%              compass(mean(L_dir_T4A),'g')
%            
%             compass(mean(L_dir_T4B),'b')
%            
%             compass(mean(L_dir_T4D),'y')


            Z.T4A= L_dir_T4A;    Masks.T4A=T4A_masks;   MAXdeg.T4A= T4A_max; 
            Z.T4B= L_dir_T4B;    Masks.T4B=T4B_masks;   MAXdeg.T4B= T4B_max;   
            Z.T4C= L_dir_T4C;    Masks.T4C=T4C_masks;   MAXdeg.T4C= T4C_max;  
            Z.T4D= L_dir_T4D;    Masks.T4D=T4D_masks;   MAXdeg.T4D= T4D_max;  
            
            Z.T5A= L_dir_T5A;    Masks.T5A=T5A_masks;   MAXdeg.T5A= T5A_max;  
            Z.T5B= L_dir_T5B;    Masks.T5B=T5B_masks;   MAXdeg.T5B= T5B_max; 
            Z.T5C= L_dir_T5C;    Masks.T5C=T5C_masks;   MAXdeg.T5C= T5C_max; 
            Z.T5D= L_dir_T5D;    Masks.T5D=T5D_masks;   MAXdeg.T5D= T5D_max; 
            
           
           
            fps=x.strct.xml.framerate;
            t=[0:length(DATA.AV_ROIS_resp{1,1})-1]/fps;
            %ms: new timing vector for interpolation at 10Hz (0.1), starts at ends with 0.5/fps shift 
            
%              it=[0.5/fps:0.1:(length(DATA.AV_ROIS_resp{1,1})+0.5)/fps];
                        
%                it=[0:0.05:t(length(t))];            
            it=[0.5/fps:0.1:(0.5)/fps+7.9]; % here I change to 8 seconds
            
            iAV_ROI_Resp=[];
            for ddir=1:8
                if ~isempty(M_T5A)
                    iAV_ROI_Resp.iM_T5A(ddir,:)=interp1(t,M_T5A(ddir,:),it,'linear','extrap');
                else
                    iAV_ROI_Resp.iM_T5A=M_T5A; %If empty
                end
                
                if ~isempty(M_T5B)
                    iAV_ROI_Resp.iM_T5B(ddir,:)=interp1(t,M_T5B(ddir,:),it,'linear','extrap');
                else
                    iAV_ROI_Resp.iM_T5B=M_T5B;
                end 
                
                if ~isempty(M_T5C)
                    iAV_ROI_Resp.iM_T5C(ddir,:)=interp1(t,M_T5C(ddir,:),it,'linear','extrap');
                else
                    iAV_ROI_Resp.iM_T5C=M_T5C;
                end 
                
                if ~isempty(M_T5D) 
                    iAV_ROI_Resp.iM_T5D(ddir,:)=interp1(t,M_T5D(ddir,:),it,'linear','extrap');
                else 
                    iAV_ROI_Resp.iM_T5D=M_T5D;
                end
                
                
                if ~isempty(M_T4A)
                    iAV_ROI_Resp.iM_T4A(ddir,:)=interp1(t,M_T4A(ddir,:),it,'linear','extrap');
                else 
                    iAV_ROI_Resp.iM_T4A=M_T4A;
                end
                
                
                if ~isempty(M_T4B)
                    iAV_ROI_Resp.iM_T4B(ddir,:)=interp1(t,M_T4B(ddir,:),it,'linear','extrap');
                else 
                    iAV_ROI_Resp.iM_T4B=M_T4B;
                end 
                
                if ~isempty(M_T4C) 
                    iAV_ROI_Resp.iM_T4C(ddir,:)=interp1(t,M_T4C(ddir,:),it,'linear','extrap');
                else
                    iAV_ROI_Resp.iM_T4C=M_T4C;
                end 
                
                if ~isempty(M_T4D) 
                    iAV_ROI_Resp.iM_T4D(ddir,:)=interp1(t,M_T4D(ddir,:),it,'linear','extrap');
                else
                    iAV_ROI_Resp.iM_T4D=M_T4D;
                end
                
                
                
           end 




            out(count).iAV_ROI_Resp=iAV_ROI_Resp;
            out(count).Flyname=FLYname;   
            out(count).driver=driver{in(iFLY)};
            out(count).i_moving=moving(in(iFLY));
            out(count).i_responssivefly=responsivefly(in(iFLY));
            out(count).NROIs=NROIS;
            out(count).Z=Z;
            out(count).Masks=Masks;
            out(count).MAXdeg=MAXdeg;
            
            try
                AV = squeeze(sum(x.strct.ch1a,3))/size(x.strct.ch1a,3); % The average image
            catch 
                AV = squeeze(sum(x.strct.ch1,3))/size(x.strct.ch1,3); % The average image
            end 
            AV = im2double(AV);
            AV = AV./max(AV(:));
            out(count).AV=AV;   % Average intensity projection of the LP

            count=count+1;
        
%     end
end

cd(currd);