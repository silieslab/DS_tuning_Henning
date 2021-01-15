function out=load_neuron_data10Hz_ks_CA_average_Stripes_mh(in,locdir,DataType)

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
        
            % FIRST FIND ALL CLUSTERS FOR T4/T5 EACH LAYER
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
            
            Teta_rad=([90,45, 0, 315, 270, 225, 180, 135]*pi)/180;
            
            %Get T5 Data
            L_dir_T5A=nan(1,length(T5_A));
            PDresp_T5A=nan(1,length(T5_A));
            NDresp_T5A=nan(1,length(T5_A));
            T5A_ALL_ROIS_resp=[];
            T5A_masks={};
            for NROI=1:length(T5_A) 
             T5A_ALL_ROIS_resp=cat(3,T5A_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_A(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T5_A(NROI)},[],2)';
             L_dir_T5A(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5A_masks{NROI}=DATA.masks_CA{T5_A(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T5A(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T5A(NROI)=R_teta(NDi);  
            end  
            M_T5A=mean(T5A_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
            
            L_dir_T5B=nan(1,length(T5_B));
            PDresp_T5B=nan(1,length(T5_B));
            NDresp_T5B=nan(1,length(T5_B));
            T5B_ALL_ROIS_resp=[];
            T5B_masks={};
            for NROI=1:length(T5_B)
             T5B_ALL_ROIS_resp=cat(3,T5B_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_B(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T5_B(NROI)},[],2)';
             L_dir_T5B(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5B_masks{NROI}=DATA.masks_CA{T5_B(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T5B(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T5B(NROI)=R_teta(NDi);
            end        
            M_T5B=mean(T5B_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            L_dir_T5C=nan(1,length(T5_C));
            PDresp_T5C=nan(1,length(T5_C));
            NDresp_T5C=nan(1,length(T5_C));
            T5C_ALL_ROIS_resp=[];
            T5C_masks={};
            for NROI=1:length(T5_C)
             T5C_ALL_ROIS_resp=cat(3,T5C_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_C(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T5_C(NROI)},[],2)';
             L_dir_T5C(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5C_masks{NROI}=DATA.masks_CA{T5_C(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T5C(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T5C(NROI)=R_teta(NDi);
            end        
            M_T5C=mean(T5C_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            
            L_dir_T5D=nan(1,length(T5_D));
            PDresp_T5D=nan(1,length(T5_D));
            NDresp_T5D=nan(1,length(T5_D));
             T5D_ALL_ROIS_resp=[];
             T5D_masks={};
            for NROI=1:length(T5_D)
             T5D_ALL_ROIS_resp=cat(3,T5D_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T5_D(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T5_D(NROI)},[],2)';
             L_dir_T5D(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T5D_masks{NROI}=DATA.masks_CA{T5_D(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T5D(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T5D(NROI)=R_teta(NDi);
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
            PDresp_T4A=nan(1,length(T4_A));
            NDresp_T4A=nan(1,length(T4_A));
            T4A_ALL_ROIS_resp=[];
            T4A_masks={};
            for NROI=1:length(T4_A) 
             T4A_ALL_ROIS_resp=cat(3,T4A_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_A(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T4_A(NROI)},[],2)';
             L_dir_T4A(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4A_masks{NROI}=DATA.masks_CA{T4_A(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T4A(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T4A(NROI)=R_teta(NDi);
            end        
            M_T4A=mean(T4A_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
            
            
            L_dir_T4B=nan(1,length(T4_B));
            PDresp_T4B=nan(1,length(T4_B));
            NDresp_T4B=nan(1,length(T4_B));
            T4B_ALL_ROIS_resp=[]; 
            T4B_masks={};
            for NROI=1:length(T4_B)
             T4B_ALL_ROIS_resp=cat(3,T4B_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_B(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T4_B(NROI)},[],2)';
             L_dir_T4B(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4B_masks{NROI}=DATA.masks_CA{T4_B(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T4B(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T4B(NROI)=R_teta(NDi);
            end        
            M_T4B=mean(T4B_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            
            L_dir_T4C=nan(1,length(T4_C));
            PDresp_T4C=nan(1,length(T4_C));
            NDresp_T4C=nan(1,length(T4_C));
            T4C_ALL_ROIS_resp=[]; 
            T4C_masks={};
            for NROI=1:length(T4_C)
             T4C_ALL_ROIS_resp=cat(3,T4C_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_C(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T4_C(NROI)},[],2)';
             L_dir_T4C(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4C_masks{NROI}=DATA.masks_CA{T4_C(NROI)};
             [PDr,PDi]=max(R_teta);
             PDresp_T4C(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T4C(NROI)=R_teta(NDi);
            end        
            M_T4C=mean(T4C_ALL_ROIS_resp,3); %Mean across Clusters for one FLy 
           
            
             L_dir_T4D=nan(1,length(T4_D));
             PDresp_T4D=nan(1,length(T4_D));
             NDresp_T4D=nan(1,length(T4_D));
             T4D_ALL_ROIS_resp=[]; 
             T4D_masks={};
            for NROI=1:length(T4_D)
             T4D_ALL_ROIS_resp=cat(3,T4D_ALL_ROIS_resp,DATA.AV_ROIS_resp{1,T4_D(NROI)});
             R_teta=max(DATA.AV_ROIS_resp{1,T4_D(NROI)},[],2)';
             L_dir_T4D(NROI)=sum(R_teta.*exp(1i*Teta_rad))/sum(abs(R_teta));
             T4D_masks{NROI}=DATA.masks_CA{T4_D(NROI)};
              [PDr,PDi]=max(R_teta);
             PDresp_T4D(NROI)=PDr;
             NDi=PDi+4; NDi(NDi>8)=NDi(NDi>8)-8;
             NDresp_T4D(NROI)=R_teta(NDi);
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


            Z.T4A= L_dir_T4A;    Masks.T4A=T4A_masks;   
            Z.T4B= L_dir_T4B;    Masks.T4B=T4B_masks;
            Z.T4C= L_dir_T4C;    Masks.T4C=T4C_masks;
            Z.T4D= L_dir_T4D;    Masks.T4D=T4D_masks;
            
            Z.T5A= L_dir_T5A;    Masks.T5A=T5A_masks;
            Z.T5B= L_dir_T5B;    Masks.T5B=T5B_masks;
            Z.T5C= L_dir_T5C;    Masks.T5C=T5C_masks;
            Z.T5D= L_dir_T5D;    Masks.T5D=T5D_masks;
            
           
            PD_resp.T4A=PDresp_T4A;     ND_resp.T4A=NDresp_T4A;
            PD_resp.T4B=PDresp_T4B;     ND_resp.T4B=NDresp_T4B;
            PD_resp.T4C=PDresp_T4C;     ND_resp.T4C=NDresp_T4C;
            PD_resp.T4D=PDresp_T4D;     ND_resp.T4D=NDresp_T4D;
            
            PD_resp.T5A=PDresp_T5A;     ND_resp.T5A=NDresp_T5A;
            PD_resp.T5B=PDresp_T5B;     ND_resp.T5B=NDresp_T5B;
            PD_resp.T5C=PDresp_T5C;     ND_resp.T5C=NDresp_T5C;
            PD_resp.T5D=PDresp_T5D;     ND_resp.T5D=NDresp_T5D;
            
            
            
            
            fps=x.strct.xml.framerate;
            t=[1:length(DATA.AV_ROIS_resp{1,1})]/fps;
            %ms: new timing vector for interpolation at 10Hz (0.1), starts at ends with 0.5/fps shift 
            
%              it=[0.5/fps:0.1:(length(DATA.AV_ROIS_resp{1,1})+0.5)/fps];
            it=[0.5/fps:0.1:(0.5)/fps+4.9]; % here I change to 8 seconds
            
            iAV_ROI_Resp=[];
            for ddir=1:8;
                try
                    iAV_ROI_Resp.iM_T5A(ddir,:)=interp1(t,M_T5A(ddir,:),it,'linear','extrap');
                catch 
                    iAV_ROI_Resp.iM_T5A=M_T5A; %If empty
                end
                
                try
                    iAV_ROI_Resp.iM_T5B(ddir,:)=interp1(t,M_T5B(ddir,:),it,'linear','extrap');
                catch
                    iAV_ROI_Resp.iM_T5B=M_T5B;
                end 
                
                try 
                    iAV_ROI_Resp.iM_T5C(ddir,:)=interp1(t,M_T5C(ddir,:),it,'linear','extrap');
                catch 
                    iAV_ROI_Resp.iM_T5C=M_T5C;
                end 
                
                try 
                    iAV_ROI_Resp.iM_T5D(ddir,:)=interp1(t,M_T5D(ddir,:),it,'linear','extrap');
                catch 
                    iAV_ROI_Resp.iM_T5D=M_T5D;
                end 
                
                try 
                    iAV_ROI_Resp.iM_T4A(ddir,:)=interp1(t,M_T4A(ddir,:),it,'linear','extrap');
                catch 
                    iAV_ROI_Resp.iM_T4A=M_T4A;
                end
                
                
                try
                    iAV_ROI_Resp.iM_T4B(ddir,:)=interp1(t,M_T4B(ddir,:),it,'linear','extrap');
                catch 
                    iAV_ROI_Resp.iM_T4B=M_T4B;
                end 
                
                try 
                    iAV_ROI_Resp.iM_T4C(ddir,:)=interp1(t,M_T4C(ddir,:),it,'linear','extrap');
                catch
                    iAV_ROI_Resp.iM_T4C=M_T4C;
                end 
                
                try 
                    iAV_ROI_Resp.iM_T4D(ddir,:)=interp1(t,M_T4D(ddir,:),it,'linear','extrap');
                catch
                    iAV_ROI_Resp.iM_T4D=M_T4D;
                end
                
                
                
            end 
            
           
            ROI_Resp.T5A=T5A_ALL_ROIS_resp;
            ROI_Resp.T5B=T5B_ALL_ROIS_resp;
            ROI_Resp.T5C=T5C_ALL_ROIS_resp;
            ROI_Resp.T5D=T5D_ALL_ROIS_resp;
            
            ROI_Resp.T4A=T4A_ALL_ROIS_resp;
            ROI_Resp.T4B=T4B_ALL_ROIS_resp;
            ROI_Resp.T4C=T4C_ALL_ROIS_resp;
            ROI_Resp.T4D=T4D_ALL_ROIS_resp;
           
            out(count).ROI_Resp=ROI_Resp; 
            out(count).iAV_ROI_Resp=iAV_ROI_Resp;
            out(count).Flyname=FLYname;   
            out(count).driver=driver{in(iFLY)};
            out(count).i_moving=moving(in(iFLY));
            out(count).i_responssivefly=responsivefly(in(iFLY));
            out(count).NROIs=NROIS;
            out(count).Z=Z;
            out(count).Masks=Masks;
            out(count).PD_resp=PD_resp;
            out(count).ND_resp=ND_resp;
            out(count).ZDepth=x.strct.xml.zdepth;
            out(count).Zoom=x.strct.xml.opticalZoom;
            
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