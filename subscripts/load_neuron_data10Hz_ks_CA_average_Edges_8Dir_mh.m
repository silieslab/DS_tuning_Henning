function out=load_neuron_data10Hz_ks_CA_average_Edges_8Dir_mh(in,locdir,DataType)

[fname,datetime,frames,zdepth,tottime,stimbouts,locnum,...
    activecells,inverted,stimcode,quality,driver,...
     moving,layer,wavelength,flyID, responsivefly,Layercheck,SingleFly,PTX_5_100,Washout5mM,Washout100mM]...
    =textread('MASTER_foldersummary_Miriam_sh.txt',...
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
            for NROI=1:length(T5_A) % for each ROI/cell:
             % Extract the response of the ROI to each direction   
             T5A_ALL_ROIS_resp=cat(3,T5A_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_A(NROI)});
             
             %Extract the maximal response of the ROI to each direction and
             %calculate the tuning vector (based on Mazurek et al., 2014)
             R_teta=max(T5A_ALL_ROIS_resp(:,1:half,NROI),[],2)';
             L_dir_T5A(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta)); % L_dir contains the tuning vector (later called Z) for each cell
             
             %Extracts the mask of the ROI
             T5A_masks{NROI}=DATA.masks_CA{T5_A(NROI)};
             
             %Extracts the preferred direction of the ROI, as the direction the neuron responded most to
             [MAX,PD]=max(R_teta);
             
             T5A_max{NROI}.MAX=MAX; % maximal response
             T5A_max{NROI}.PD=Teta_deg(PD); % direction the neuron responded most to
            end  
            
            %Repeat for each subtype:
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
           
            


            % SAVE the Data (
            Z.T4A= L_dir_T4A;    Masks.T4A=T4A_masks;   MAXdeg.T4A= T4A_max; 
            Z.T4B= L_dir_T4B;    Masks.T4B=T4B_masks;   MAXdeg.T4B= T4B_max;   
            Z.T4C= L_dir_T4C;    Masks.T4C=T4C_masks;   MAXdeg.T4C= T4C_max;  
            Z.T4D= L_dir_T4D;    Masks.T4D=T4D_masks;   MAXdeg.T4D= T4D_max;  
            
            Z.T5A= L_dir_T5A;    Masks.T5A=T5A_masks;   MAXdeg.T5A= T5A_max;  
            Z.T5B= L_dir_T5B;    Masks.T5B=T5B_masks;   MAXdeg.T5B= T5B_max; 
            Z.T5C= L_dir_T5C;    Masks.T5C=T5C_masks;   MAXdeg.T5C= T5C_max; 
            Z.T5D= L_dir_T5D;    Masks.T5D=T5D_masks;   MAXdeg.T5D= T5D_max; 
            


            %Save data in out.structure
            out(count).Flyname=FLYname;   
            out(count).NROIs=NROIS;
            out(count).Z=Z;
            out(count).Masks=Masks;
            out(count).MAXdeg=MAXdeg;
            out(count).PixelSize=x.strct.xml.micronsPerPixel.XAxis;
            
            
            AV = squeeze(sum(x.strct.ch1a_crop,3))/size(x.strct.ch1a_crop,3); % The average image
          
            AV = im2double(AV);
            AV = AV./max(AV(:));
            out(count).AV=AV;   % Average intensity projection of the LP

            count=count+1;
        
%     end
end

cd(currd);