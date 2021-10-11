
function Z= averageDirectionVectors(T4T5_mb)
%First concatinate all numbers of each cell for each Fly 
% I can then take an average of all Cells for one certain type of cell 
FlyIDs= [ones(1,9), ones(1,9)*2,ones(1,8)*3,ones(1,12)*4,ones(1,10)*5,ones(1,8)*6,ones(1,8)*7,ones(1,6)*8,ones(1,9)*9,ones(1,5)*10,ones(1,4)*11,ones(1,5)*12,ones(1,10)*13,ones(1,11)*14];

Z_T5A_ALL=[];Z_T5B_ALL=[];Z_T5C_ALL=[];Z_T5D_ALL=[];
Z_T4A_ALL=[];Z_T4B_ALL=[];Z_T4C_ALL=[];Z_T4D_ALL=[];
FlIDs_T4A=[]; FlIDs_T4B=[];FlIDs_T4C=[];FlIDs_T4D=[];
FlIDs_T5A=[]; FlIDs_T5B=[];FlIDs_T5C=[];FlIDs_T5D=[];

for NFlies= 1:length(T4T5_mb)
    
     % T5-------------------------------------------------
     Z_T5A= eval(['T4T5_mb(', num2str(NFlies),').Z.T5A']);
     Z_T5A_ALL=[Z_T5A_ALL,Z_T5A];
     Z_T5A_mean(NFlies)=mean(Z_T5A);
     Z_T5A_std(NFlies)=std(Z_T5A);
     Ang=angle(Z_T5A); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T5A_circstd(NFlies)=circ_std(Ang');
     Z_T5A_circste(NFlies)=circ_std(Ang')/sqrt(length(Ang));
   
     FlIDs_T5A=[FlIDs_T5A,ones(1,length(Z_T5A))*FlyIDs(NFlies)];
     
     
     Z_T5B= eval(['T4T5_mb(', num2str(NFlies),').Z.T5B']);
     Z_T5B_ALL=[Z_T5B_ALL,Z_T5B];
     Z_T5B_mean(NFlies)=mean(Z_T5B);
     Z_T5B_std(NFlies)=std(Z_T5B);
      Ang=angle(Z_T5B); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T5B_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T5B_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T5B=[FlIDs_T5B,ones(1,length(Z_T5B))*FlyIDs(NFlies)];

     
     Z_T5C= eval(['T4T5_mb(', num2str(NFlies),').Z.T5C']);
     Z_T5C_ALL=[Z_T5C_ALL,Z_T5C];
     Z_T5C_mean(NFlies)=mean(Z_T5C);
     Z_T5C_std(NFlies)=std(Z_T5C);
      Ang=angle(Z_T5C); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T5C_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T5C_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T5C=[FlIDs_T5C,ones(1,length(Z_T5C))*FlyIDs(NFlies)];

     
     
     
     Z_T5D= eval(['T4T5_mb(', num2str(NFlies),').Z.T5D']);
     Z_T5D_ALL=[Z_T5D_ALL,Z_T5D];
     Z_T5D_mean(NFlies)=mean(Z_T5D);
     Z_T5D_std(NFlies)=std(Z_T5D);
      Ang=angle(Z_T5D); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T5D_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T5D_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T5D=[FlIDs_T5D,ones(1,length(Z_T5D))*FlyIDs(NFlies)];

     
     % T4-------------------------------------------------
     Z_T4A= eval(['T4T5_mb(', num2str(NFlies),').Z.T4A']);
     Z_T4A_ALL=[Z_T4A_ALL,Z_T4A];
     Z_T4A_mean(NFlies)=mean(Z_T4A);
     Z_T4A_std(NFlies)=std(Z_T4A);
      Ang=angle(Z_T4A); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4A_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T4A_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T4A=[FlIDs_T4A,ones(1,length(Z_T4A))*FlyIDs(NFlies)];

     
     Z_T4B= eval(['T4T5_mb(', num2str(NFlies),').Z.T4B']);
     Z_T4B_ALL=[Z_T4B_ALL,Z_T4B];
     Z_T4B_mean(NFlies)=mean(Z_T4B);
     Z_T4B_std(NFlies)=std(Z_T4B);
      Ang=angle(Z_T4B); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4B_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T4B_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T4B=[FlIDs_T4B,ones(1,length(Z_T4B))*FlyIDs(NFlies)];

     
     Z_T4C= eval(['T4T5_mb(', num2str(NFlies),').Z.T4C']);
     Z_T4C_ALL=[Z_T4C_ALL,Z_T4C];
     Z_T4C_mean(NFlies)=mean(Z_T4C);
     Z_T4C_std(NFlies)=std(Z_T4C);
      Ang=angle(Z_T4C); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4C_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T4C_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T4C=[FlIDs_T4C,ones(1,length(Z_T4C))*FlyIDs(NFlies)];

     
     
     Z_T4D= eval(['T4T5_mb(', num2str(NFlies),').Z.T4D']);
     Z_T4D_ALL=[Z_T4D_ALL,Z_T4D];
     Z_T4D_mean(NFlies)=mean(Z_T4D);
     Z_T4D_std(NFlies)=std(Z_T4D);
      Ang=angle(Z_T4D); 
%      ng=find(Ang)<=0; 
%      Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4D_circstd(NFlies)=circ_std(angle(Ang)');
     Z_T4D_circste(NFlies)=circ_std(angle(Ang)')/sqrt(length(Ang));
     FlIDs_T4D=[FlIDs_T4D,ones(1,length(Z_T4D))*FlyIDs(NFlies)];

end 
    
    Z.T4A.ALL=Z_T4A_ALL;
    Z.T4A.M=Z_T4A_mean;
    Z.T4A.STD=Z_T4A_std;
    Z.T4A.circSTD=Z_T4A_circstd; 
    Z.T4A.circSTE=Z_T4A_circste; 
    Z.T4A.FlyIDs=FlIDs_T4A;
    
    
    Z.T4B.ALL=Z_T4B_ALL;
    Z.T4B.M=Z_T4B_mean;
    Z.T4B.STD=Z_T4B_std;
    Z.T4B.circSTD=Z_T4B_circstd; 
    Z.T4B.circSTE=Z_T4B_circste; 
    Z.T4B.FlyIDs=FlIDs_T4B;
    
    Z.T4C.ALL=Z_T4C_ALL;
    Z.T4C.M=Z_T4C_mean;
    Z.T4C.STD=Z_T4C_std;
    Z.T4C.circSTD=Z_T4C_circstd; 
    Z.T4C.circSTE=Z_T4C_circste; 
    Z.T4C.FlyIDs=FlIDs_T4C;
    
    Z.T4D.ALL=Z_T4D_ALL;
    Z.T4D.M=Z_T4D_mean;
    Z.T4D.STD=Z_T4D_std;
    Z.T4D.circSTD=Z_T4D_circstd; 
    Z.T4D.circSTE=Z_T4D_circste; 
    Z.T4D.FlyIDs=FlIDs_T4D;
    
    
    Z.T5A.ALL=Z_T5A_ALL;
    Z.T5A.M=Z_T5A_mean;
    Z.T5A.STD=Z_T5A_std;
    Z.T5A.circSTD=Z_T5A_circstd; 
    Z.T5A.circSTE=Z_T5A_circste; 
    Z.T5A.FlyIDs=FlIDs_T5A;
    
    Z.T5B.ALL=Z_T5B_ALL;
    Z.T5B.M=Z_T5B_mean;
    Z.T5B.STD=Z_T5B_std;
    Z.T5B.circSTD=Z_T5B_circstd; 
    Z.T5B.circSTE=Z_T5B_circste; 
    Z.T5B.FlyIDs=FlIDs_T5B;
    
    Z.T5C.ALL=Z_T5C_ALL;
    Z.T5C.M=Z_T5C_mean;
    Z.T5C.STD=Z_T5C_std;
    Z.T5C.circSTD=Z_T5C_circstd; 
    Z.T5C.circSTE=Z_T5C_circste; 
    Z.T5C.FlyIDs=FlIDs_T5C;
    
    Z.T5D.ALL=Z_T5D_ALL;
    Z.T5D.M=Z_T5D_mean;
    Z.T5D.STD=Z_T5D_std;
    Z.T5D.circSTD=Z_T5D_circstd; 
    Z.T5D.circSTE=Z_T5D_circste; 
    Z.T5D.FlyIDs=FlIDs_T5D;

    
end 