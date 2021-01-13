%%%%%%%%%%%%%%
% This script plots tuning maps shown in Figure4 and Extebded Data Fig. 4
% d-f 
%%%%%%%%%%%%%% 

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



Split_ONOFF=true;

load('MyColormap.mat')


 %%
  % Now Plot quiver plots 
 Quiver_T4A_all=[]; Quiver_T5A_all=[];
 Coord_T4A_all=[];  Coord_T5A_all=[];
 Quiver_T4B_all=[]; Quiver_T5B_all=[];
 Coord_T4B_all=[];  Coord_T5B_all=[];
 Quiver_T4C_all=[]; Quiver_T5C_all=[];
 Coord_T4C_all=[];  Coord_T5C_all=[];
 Quiver_T4D_all=[];  Quiver_T5D_all=[];
 Coord_T4D_all=[];  Coord_T5D_all=[];
 
 FlyID_T4A_all=[]; FlyID_T5A_all=[]; 
 FlyID_T4B_all=[]; FlyID_T5B_all=[]; 
 FlyID_T4C_all=[]; FlyID_T5C_all=[]; 
 FlyID_T4D_all=[]; FlyID_T5D_all=[]; 
 
 CellID_T4A_all=[]; CellID_T5A_all=[]; 
 CellID_T4B_all=[]; CellID_T5B_all=[]; 
 CellID_T4C_all=[]; CellID_T5C_all=[]; 
 CellID_T4D_all=[]; CellID_T5D_all=[]; 
 
 Field_T4A_all=[]; Field_T5A_all=[]; 
 Field_T4B_all=[]; Field_T5B_all=[]; 
 Field_T4C_all=[]; Field_T5C_all=[]; 
 Field_T4D_all=[]; Field_T5D_all=[]; 



 RANGE=[1:size(Control.T4T5_mb,2)];


FlyID=1; 
co=1;

name_before='no';
        
 for i=RANGE
     
     Quiver_T4A=[]; Quiver_T5A=[];
     Coord_T4A=[];  Coord_T5A=[];
     Quiver_T4B=[]; Quiver_T5B=[];
     Coord_T4B=[];  Coord_T5B=[];
     Quiver_T4C=[]; Quiver_T5C=[];
     Coord_T4C=[];  Coord_T5C=[];
     Quiver_T4D=[];  Quiver_T5D=[];
     Coord_T4D=[];  Coord_T5D=[];
 

    IFly=Control.T4T5_mb(i);
    iname=IFly.Flyname(1:11);
    if ~ strcmp(iname, name_before)
       co=co+3; 
       FlyID=FlyID+1; 
    end 
    
    Turn=IFly.turn;

    for ii=1:length(IFly.Masks.T4A)
        Quiver_T4A=[Quiver_T4A,IFly.Z.T4A(ii)*100];
        Y=IFly.RFCenter.T4A(1,ii);
        X=IFly.RFCenter.T4A(2,ii);
        Coord_T4A=[Coord_T4A,[-X;Y+Turn]];
    end
      
     for ii=1:length(IFly.Masks.T4B)
        Quiver_T4B=[Quiver_T4B,IFly.Z.T4B(ii)*100];
        Y=IFly.RFCenter.T4B(1,ii);
        X=IFly.RFCenter.T4B(2,ii);
        Coord_T4B=[Coord_T4B,[-X;Y+Turn]];
     end
    
      for ii=1:length(IFly.Masks.T4C)
        Quiver_T4C=[Quiver_T4C,IFly.Z.T4C(ii)*100];
        Y=IFly.RFCenter.T4C(1,ii);
        X=IFly.RFCenter.T4C(2,ii);
        Coord_T4C=[Coord_T4C,[-X;Y+Turn]];
      end
    
      for ii=1:length(IFly.Masks.T4D)
        Quiver_T4D=[Quiver_T4D,IFly.Z.T4D(ii)*100];
        Y=IFly.RFCenter.T4D(1,ii);
        X=IFly.RFCenter.T4D(2,ii);
        Coord_T4D=[Coord_T4D,[-X;Y+Turn]];
      end
    
      
      
    for ii=1:length(IFly.Masks.T5A)
        Quiver_T5A=[Quiver_T5A,IFly.Z.T5A(ii)*100];
        Y=IFly.RFCenter.T5A(1,ii);
        X=IFly.RFCenter.T5A(2,ii);
        Coord_T5A=[Coord_T5A,[-X;Y+Turn]];
    end
      
     for ii=1:length(IFly.Masks.T5B)
        Quiver_T5B=[Quiver_T5B,IFly.Z.T5B(ii)*100];
        Y=IFly.RFCenter.T5B(1,ii);
        X=IFly.RFCenter.T5B(2,ii);
        Coord_T5B=[Coord_T5B,[-X;Y+Turn]];
     end
    
      for ii=1:length(IFly.Masks.T5C)
        Quiver_T5C=[Quiver_T5C,IFly.Z.T5C(ii)*100];
        Y=IFly.RFCenter.T5C(1,ii);
        X=IFly.RFCenter.T5C(2,ii);
        Coord_T5C=[Coord_T5C,[-X;Y+Turn]];
        
      end
    
      for ii=1:length(IFly.Masks.T5D)
        Quiver_T5D=[Quiver_T5D,IFly.Z.T5D(ii)*100];
        Y=IFly.RFCenter.T5D(1,ii);
        X=IFly.RFCenter.T5D(2,ii);
        Coord_T5D=[Coord_T5D,[-X;Y+Turn]];
      end
      
 
    
     Quiver_T4A_all=[Quiver_T4A_all,Quiver_T4A]; Quiver_T5A_all=[Quiver_T5A_all,Quiver_T5A];
     Coord_T4A_all=[Coord_T4A_all,Coord_T4A];  Coord_T5A_all=[Coord_T5A_all,Coord_T5A];
     Quiver_T4B_all=[Quiver_T4B_all,Quiver_T4B]; Quiver_T5B_all=[Quiver_T5B_all,Quiver_T5B];
     Coord_T4B_all=[Coord_T4B_all,Coord_T4B];  Coord_T5B_all=[Coord_T5B_all,Coord_T5B];
     Quiver_T4C_all=[Quiver_T4C_all,Quiver_T4C]; Quiver_T5C_all=[Quiver_T5C_all,Quiver_T5C];
     Coord_T4C_all=[Coord_T4C_all,Coord_T4C];  Coord_T5C_all=[Coord_T5C_all,Coord_T5C];
     Quiver_T4D_all=[Quiver_T4D_all,Quiver_T4D];  Quiver_T5D_all=[Quiver_T5D_all,Quiver_T5D];
     Coord_T4D_all=[Coord_T4D_all,Coord_T4D];  Coord_T5D_all=[Coord_T5D_all,Coord_T5D];


     FlyID_T4A_all=[FlyID_T4A_all,ones(1, length(Quiver_T4A))*FlyID];
     FlyID_T4B_all=[FlyID_T4B_all,ones(1, length(Quiver_T4B))*FlyID];
     FlyID_T4C_all=[FlyID_T4C_all,ones(1, length(Quiver_T4C))*FlyID];
     FlyID_T4D_all=[FlyID_T4D_all,ones(1, length(Quiver_T4D))*FlyID]; 
     FlyID_T5A_all=[FlyID_T5A_all,ones(1, length(Quiver_T5A))*FlyID];
     FlyID_T5B_all=[FlyID_T5B_all,ones(1, length(Quiver_T5B))*FlyID];
     FlyID_T5C_all=[FlyID_T5C_all,ones(1, length(Quiver_T5C))*FlyID];
     FlyID_T5D_all=[FlyID_T5D_all,ones(1, length(Quiver_T5D))*FlyID];
     
     Field_T4A_all=[Field_T4A_all,ones(1, length(Quiver_T4A))*i];
     Field_T4B_all=[Field_T4B_all,ones(1, length(Quiver_T4B))*i];
     Field_T4C_all=[Field_T4C_all,ones(1, length(Quiver_T4C))*i];
     Field_T4D_all=[Field_T4D_all,ones(1, length(Quiver_T4D))*i]; 
     Field_T5A_all=[Field_T5A_all,ones(1, length(Quiver_T5A))*i];
     Field_T5B_all=[Field_T5B_all,ones(1, length(Quiver_T5B))*i];
     Field_T5C_all=[Field_T5C_all,ones(1, length(Quiver_T5C))*i];
     Field_T5D_all=[Field_T5D_all,ones(1, length(Quiver_T5D))*i];

     
     
     if ~isempty(length(IFly.CellID.T4A))
        CellID_T4A_all=[CellID_T4A_all, [1:length(IFly.CellID.T4A)]];
     end 
     if ~isempty(length(IFly.CellID.T4B))
        CellID_T4B_all=[CellID_T4B_all, [1:length(IFly.CellID.T4B)]]; 
     end 
     if ~isempty(length(IFly.CellID.T4C))
        CellID_T4C_all=[CellID_T4C_all, [1:length(IFly.CellID.T4C)]]; 
     end 
     if ~isempty(length(IFly.CellID.T4D))
        CellID_T4D_all=[CellID_T4D_all, [1:length(IFly.CellID.T4D)]]; 
     end 
     
     
      if ~isempty(length(IFly.CellID.T5A))
        CellID_T5A_all=[CellID_T5A_all, [1:length(IFly.CellID.T5A)]]; 
     end 
     if ~isempty(length(IFly.CellID.T5B))
        CellID_T5B_all=[CellID_T5B_all, [1:length(IFly.CellID.T5B)]]; 
     end 
     if ~isempty(length(IFly.CellID.T5C))
        CellID_T5C_all=[CellID_T5C_all, [1:length(IFly.CellID.T5C)]]; 
     end 
     if ~isempty(length(IFly.CellID.T5D))
        CellID_T5D_all=[CellID_T5D_all, [1:length(IFly.CellID.T5D)]]; 
     end 
     
   name_before=iname;
 end 
%%
 
 C=10; 
 Quiver_T4A_all_s=Quiver_T4A_all/C;
 Quiver_T4B_all_s=Quiver_T4B_all/C;
 Quiver_T4C_all_s=Quiver_T4C_all/C;
 Quiver_T4D_all_s=Quiver_T4D_all/C;
 Quiver_T5A_all_s=Quiver_T5A_all/C;
 Quiver_T5B_all_s=Quiver_T5B_all/C;
 Quiver_T5C_all_s=Quiver_T5C_all/C;
 Quiver_T5D_all_s=Quiver_T5D_all/C;
 Quiver_Try=10; 

%% Neurons in Layer A and B seem to have two subpopulation: Cluster Data 


% Use data from snob analysis

TA_T5=ClusterR.TA_T5;
TA_T4=ClusterR.TA_T4;

TB_T5=ClusterR.TB_T5;
TB_T4=ClusterR.TB_T4;

TC_T5=ClusterR.TC_T5;
TC_T4=ClusterR.TC_T4;

TD_T5=ClusterR.TD_T5;
TD_T4=ClusterR.TD_T4;



F9=figure('Position',[200 200 1850 500]);
subplot(2,3,1)
quiver(Coord_T5A_all(2,TA_T5==3)-34,Coord_T5A_all(1,TA_T5==3)+36,real(Quiver_T5A_all_s(TA_T5==3)),imag(Quiver_T5A_all_s(TA_T5==3)), 'AutoScale','off')
hold on 
quiver(Coord_T4A_all(2,TA_T4==3)-34,Coord_T4A_all(1,TA_T4==3)+36,real(Quiver_T4A_all_s(TA_T4==3)),imag(Quiver_T4A_all_s(TA_T4==3)), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['subtype A.I'])
hold on 
quiver(20,120,imag(Quiver_Try),real(Quiver_Try), 'AutoScale','off')

subplot(2,3,2)
quiver(Coord_T5A_all(2,TA_T5==1)-34,Coord_T5A_all(1,TA_T5==1)+36,real(Quiver_T5A_all_s(TA_T5==1)),imag(Quiver_T5A_all_s(TA_T5==1)), 'AutoScale','off')
hold on 
quiver(Coord_T4A_all(2,TA_T4==1)-34,Coord_T4A_all(1,TA_T4==1)+36,real(Quiver_T4A_all_s(TA_T4==1)),imag(Quiver_T4A_all_s(TA_T4==1)), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['subtype A.II'])




subplot(2,3,4)
quiver(Coord_T5B_all(2,TB_T5==1)-34,Coord_T5B_all(1,TB_T5==1)+36,real(Quiver_T5B_all_s(TB_T5==1)),imag(Quiver_T5B_all_s(TB_T5==1)), 'AutoScale','off')
hold on 
quiver(Coord_T4B_all(2,TB_T4==1)-34,Coord_T4B_all(1,TB_T4==1)+36,real(Quiver_T4B_all_s(TB_T4==1)),imag(Quiver_T4B_all_s(TB_T4==1)), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['subtype B.I'])


subplot(2,3,5)
quiver(Coord_T5B_all(2,TB_T5==2)-34,Coord_T5B_all(1,TB_T5==2)+36,real(Quiver_T5B_all_s(TB_T5==2)),imag(Quiver_T5B_all_s(TB_T5==2)), 'AutoScale','off')
hold on 
quiver(Coord_T4B_all(2,TB_T4==2)-34,Coord_T4B_all(1,TB_T4==2)+36,real(Quiver_T4B_all_s(TB_T4==2)),imag(Quiver_T4B_all_s(TB_T4==2)), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['subtype B.II'])


subplot(2,3,3)
quiver(Coord_T5C_all(2,~isnan(sum(Coord_T5C_all)))-34,Coord_T5C_all(1,~isnan(sum(Coord_T5C_all)))+36,real(Quiver_T5C_all_s(~isnan(sum(Coord_T5C_all)))),imag(Quiver_T5C_all_s(~isnan(sum(Coord_T5C_all)))), 'AutoScale','off')
hold on 
quiver(Coord_T4C_all(2,~isnan(sum(Coord_T4C_all)))-34,Coord_T4C_all(1,~isnan(sum(Coord_T4C_all)))+36,real(Quiver_T4C_all_s(~isnan(sum(Coord_T4C_all)))),imag(Quiver_T4C_all_s(~isnan(sum(Coord_T4C_all)))), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['subtype C'])

subplot(2,3,6)
quiver(Coord_T5D_all(2,~isnan(sum(Coord_T5D_all)))-34,Coord_T5D_all(1,~isnan(sum(Coord_T5D_all)))+36,real(Quiver_T5D_all_s(~isnan(sum(Coord_T5D_all)))),imag(Quiver_T5D_all_s(~isnan(sum(Coord_T5D_all)))), 'AutoScale','off')
hold on 
quiver(Coord_T4D_all(2,~isnan(sum(Coord_T4D_all)))-34,Coord_T4D_all(1,~isnan(sum(Coord_T4D_all)))+36,real(Quiver_T4D_all_s(~isnan(sum(Coord_T4D_all)))),imag(Quiver_T4D_all_s(~isnan(sum(Coord_T4D_all)))), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['subtype D'])



% subplot(2,3,3)
% quiver(Coord_T5C_all(2,~(TC_T5==2))-34,Coord_T5C_all(1,~(TC_T5==2))+36,real(Quiver_T5C_all_s(~(TC_T5==2))),imag(Quiver_T5C_all_s(~(TC_T5==2))), 'AutoScale','off')
% hold on 
% quiver(Coord_T4C_all(2,~(TC_T4==2))-34,Coord_T4C_all(1,~(TC_T4==2))+36,real(Quiver_T4C_all_s(~(TC_T4==2))),imag(Quiver_T4C_all_s(~(TC_T4==2))), 'AutoScale','off')
% axis('equal')
% set(gca,'XLim', [-34,44+2*45])
% set(gca,'YLim', [-17,36])
% title(['subtype C'])
% 
% subplot(2,3,6)
% quiver(Coord_T5D_all(2,TD_T5==1)-34,Coord_T5D_all(1,TD_T5==1)+36,real(Quiver_T5D_all_s(TD_T5==1)),imag(Quiver_T5D_all_s(TD_T5==1)), 'AutoScale','off')
% hold on 
% quiver(Coord_T4D_all(2,TD_T4==1)-34,Coord_T4D_all(1,TD_T4==1)+36,real(Quiver_T4D_all_s(TD_T4==1)),imag(Quiver_T4D_all_s(TD_T4==1)), 'AutoScale','off')
% axis('equal')
% set(gca,'XLim', [-34,44+2*45])
% set(gca,'YLim', [-17,36])
% title(['subtype D'])
% 

%% Average change of tuning: 
% horizontal: 
horizontal_steps=-20:10:75; 
AverTuning_AI=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5A_all(2,TA_T5==3)-34;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5A_all_s(TA_T5==3); 
    
    COORDi_T4=Coord_T4A_all(2,TA_T4==3)-34;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4A_all_s(TA_T4==3); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_AI(Ni)=Quiverav; 
    
end 

AverTuning_AII=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5A_all(2,TA_T5==1)-34;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5A_all_s(TA_T5==1); 
    
    COORDi_T4=Coord_T4A_all(2,TA_T4==1)-34;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4A_all_s(TA_T4==1); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_AII(Ni)=Quiverav; 
    
end 


figure('Position', [200 200, 800, 500])
subplot(2,3,1)
AverTuning_AI_NORM=AverTuning_AI./abs(AverTuning_AI);
Comp=compass(AverTuning_AI_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AI_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
LG{i}=num2str(horizontal_steps(i)+5); 

end
h=legend(LG); 
h.Position=[0.05    0.6739    0.0687    0.2370];
title('subtype A.I')

 
subplot(2,3,4)

AverTuning_AII_NORM=AverTuning_AII./abs(AverTuning_AII);
Comp=compass(AverTuning_AII_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AII_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
end

title('subtype A.II')




AverTuning_BI=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5B_all(2,TB_T5==1)-34;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5B_all_s(TB_T5==1); 
    
    COORDi_T4=Coord_T4B_all(2,TB_T4==1)-34;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4B_all_s(TB_T4==1); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_BI(Ni)=Quiverav; 
    
end 

AverTuning_BII=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5B_all(2,TB_T5==2)-34;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5B_all_s(TB_T5==2); 
    
    COORDi_T4=Coord_T4B_all(2,TB_T4==2)-34;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4B_all_s(TB_T4==2); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_BII(Ni)=Quiverav; 
    
end 


subplot(2,3,2)
AverTuning_BI_NORM=AverTuning_BI./abs(AverTuning_BI);
Comp=compass(AverTuning_BI_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AI_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
LG{i}=num2str(horizontal_steps(i)+5); 

end

title('subtype B.I')


 
subplot(2,3,5)

AverTuning_BII_NORM=AverTuning_BII./abs(AverTuning_BII);
Comp=compass(AverTuning_BII_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AII_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
end
 

title('subtype B.II')




TC_T5=ClusterR.TC_T5; 
TC_T4=ClusterR.TC_T4; 
TD_T5=ClusterR.TD_T5; 
TD_T4=ClusterR.TD_T4; 

AverTuning_C=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5C_all(2,~(TC_T5==2))-34;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5C_all_s(~(TC_T5==2)); 
    
    COORDi_T4=Coord_T4C_all(2,~(TC_T4==2))-34;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4C_all_s(~(TC_T4==2)); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_C(Ni)=Quiverav; 
    
end 

AverTuning_D=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5D_all(2,TD_T5==1)-34;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5D_all_s(TD_T5==1); 
    
    COORDi_T4=Coord_T4D_all(2,TD_T4==1)-34;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4D_all_s(TD_T4==1); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_D(Ni)=Quiverav; 
    
end 


subplot(2,3,6)

AverTuning_D_NORM=AverTuning_D./abs(AverTuning_D);
Comp=compass(AverTuning_D_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_D_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
LG{i}=num2str(horizontal_steps(i)+5); 

end

title('subtype D')


 
subplot(2,3,3)

AverTuning_C_NORM=AverTuning_C./abs(AverTuning_C);
Comp=compass(AverTuning_C_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AII_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
end
 
title('subtype C')



%%

% vertical: 
horizontal_steps=-10:10:30; 
AverTuning_AI=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5A_all(1,TA_T5==3)+36;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5A_all_s(TA_T5==3); 
    
    COORDi_T4=Coord_T4A_all(1,TA_T4==3)+36;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4A_all_s(TA_T4==3); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_AI(Ni)=Quiverav; 
    
end 

AverTuning_AII=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5A_all(1,TA_T5==1)+36;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5A_all_s(TA_T5==1); 
    
    COORDi_T4=Coord_T4A_all(1,TA_T4==1)+36;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4A_all_s(TA_T4==1); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_AII(Ni)=Quiverav; 
    
end 


figure('Position', [200 200, 800, 500])
subplot(2,3,1)
AverTuning_AI_NORM=AverTuning_AI./abs(AverTuning_AI);
Comp=compass(AverTuning_AI_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AI_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
LG{i}=num2str(horizontal_steps(i)+5); 

end
h=legend(LG); 
h.Position=[0.05    0.6739    0.0687    0.2370];

title('LayerA.I')

 
subplot(2,3,4)

AverTuning_AII_NORM=AverTuning_AII./abs(AverTuning_AII);
Comp=compass(AverTuning_AII_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AII_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
end
 
title('LayerA.II')




AverTuning_BI=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5B_all(1,TB_T5==1)+36;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5B_all_s(TB_T5==1); 
    
    COORDi_T4=Coord_T4B_all(1,TB_T4==1)+36;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4B_all_s(TB_T4==1); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_BI(Ni)=Quiverav; 
    
end 

AverTuning_BII=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5B_all(1,TB_T5==2)+36;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5B_all_s(TB_T5==2); 
    
    COORDi_T4=Coord_T4B_all(1,TB_T4==2)+36;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4B_all_s(TB_T4==2); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_BII(Ni)=Quiverav; 
    
end 


subplot(2,3,2)
AverTuning_BI_NORM=AverTuning_BI./abs(AverTuning_BI);
Comp=compass(AverTuning_BI_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AI_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
LG{i}=num2str(horizontal_steps(i)+5); 

end

title('LayerB.I')


 
subplot(2,3,5)

AverTuning_BII_NORM=AverTuning_BII./abs(AverTuning_BII);
Comp=compass(AverTuning_BII_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AII_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
end
 

title('LayerB.II')




TC_T5=ClusterR.TC_T5; 
TC_T4=ClusterR.TC_T4; 
TD_T5=ClusterR.TD_T5; 
TD_T4=ClusterR.TD_T4; 

AverTuning_C=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5C_all(1,~(TC_T5==2))+36;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5C_all_s(~(TC_T5==2)); 
    
    COORDi_T4=Coord_T4C_all(1,~(TC_T5==2))+36;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4C_all_s(~(TC_T5==2)); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_C(Ni)=Quiverav; 
    
end 

AverTuning_D=nan(1,length(horizontal_steps)); 

for Ni=1:length(horizontal_steps)
         
    N_hsteps = horizontal_steps(Ni); 
    COORDi_T5=Coord_T5D_all(1,TD_T5==1)+36;
    INDI_T5=find((COORDi_T5>= N_hsteps).* (COORDi_T5< N_hsteps+5)); 
    Quiveri_T5=Quiver_T5D_all_s(TD_T5==1); 
    
    COORDi_T4=Coord_T4D_all(1,TD_T4==1)+36;
    INDI_T4=find((COORDi_T4>= N_hsteps).* (COORDi_T4< N_hsteps+5)); 
    Quiveri_T4=Quiver_T4D_all_s(TD_T4==1); 
    
    Quiverav=nanmean([Quiveri_T5(INDI_T5),Quiveri_T4(INDI_T4)]); 
    
    AverTuning_D(Ni)=Quiverav; 
    
end 


subplot(2,3,6)

AverTuning_D_NORM=AverTuning_D./abs(AverTuning_D);
Comp=compass(AverTuning_D_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_D_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
LG{i}=num2str(horizontal_steps(i)+5); 

end

title('LayerD')


 
subplot(2,3,3)

AverTuning_C_NORM=AverTuning_C./abs(AverTuning_C);
Comp=compass(AverTuning_C_NORM);
cm=colormap('cool'); 
cmi=cm([1:floor(length(cm)/length(AverTuning_AII_NORM)):length(cm)],:)
for i=1:length(Comp)
  set(Comp(i),'color',cmi(i,:));
    set(Comp(i),'LineWidth',2);
end
 

title('LayerC')


%% Subgroups covering screen for each fly 
%Extended Data Fig. 4d 

ID=5; % ID=5 is the FlyID shown in the manuscript, any other exmaple can be plotted as well by changing this to a value between 2 and 15
F9=figure;
subplot(2,1,1)
Ind=find((TA_T5==3).*(FlyID_T5A_all==ID));
Ind2=find((TA_T4==3).*(FlyID_T4A_all==ID)); 
P1=plot(Coord_T5A_all(2,Ind)-34,Coord_T5A_all(1,Ind)+36, 'kx'); 
hold on 
plot(Coord_T4A_all(2,Ind2)-34,Coord_T4A_all(1,Ind2)+36, 'kx')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerA'])

Ind=find((TA_T5==1).*(FlyID_T5A_all==ID));
Ind2=find((TA_T4==1).*(FlyID_T4A_all==ID)); 
P2=plot(Coord_T5A_all(2,Ind)-34,Coord_T5A_all(1,Ind)+36, 'rx'); 
hold on 
plot(Coord_T4A_all(2,Ind2)-34,Coord_T4A_all(1,Ind2)+36, 'rx')
legend([P1(1), P2(1)],'A.I', 'A.II')

subplot(2,1,2)
Ind=find((TB_T5==1).*(FlyID_T5B_all==ID));
Ind2=find((TB_T4==1).*(FlyID_T4B_all==ID)); 
P3=plot(Coord_T5B_all(2,Ind)-34,Coord_T5B_all(1,Ind)+36, 'kx'); 
hold on 
plot(Coord_T4B_all(2,Ind2)-34,Coord_T4B_all(1,Ind2)+36, 'kx')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerB'])

Ind=find((TB_T5==2).*(FlyID_T5B_all==ID));
Ind2=find((TB_T4==2).*(FlyID_T4B_all==ID)); 
P4=plot(Coord_T5B_all(2,Ind)-34,Coord_T5B_all(1,Ind)+36, 'rx');
hold on 
plot(Coord_T4B_all(2,Ind2)-34,Coord_T4B_all(1,Ind2)+36, 'rx')
legend([P3(1), P4(1)], 'B.I', 'B.II')

Savename=['RF_Center_Analysis_FlyID:', num2str(ID)];
% set(F9,'PaperSize', [50 20])
                            

    



