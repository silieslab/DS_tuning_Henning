
function Z= circ_distribution(T4T5_mb_Bright,T4T5_mb_Dark)
%First concatinate all numbers of each cell for each Fly 
% I can then take an average of all Cells for one certain type of cell 
Z_T5A_ALL=[];Z_T5B_ALL=[];Z_T5C_ALL=[];Z_T5D_ALL=[];
Z_T4A_ALL=[];Z_T4B_ALL=[];Z_T4C_ALL=[];Z_T4D_ALL=[];

for NFlies= 1:length(T4T5_mb_Bright)
    
     % T5-------------------------------------------------
     Z_T5A= eval(['T4T5_mb_Dark(', num2str(NFlies),').Z.T5A']);
     Z_T4A= eval(['T4T5_mb_Bright(', num2str(NFlies),').Z.T4A']);

     Ang=[angle(Z_T5A),angle(Z_T4A)]; 
     ng=find(Ang)<=0; 
     Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4T5A_circstd(NFlies)=circ_std(Ang');
     Z_T4T5A_circste(NFlies)=circ_std(Ang')/sqrt(length(Ang));
     
     
     Z_T5B= eval(['T4T5_mb_Dark(', num2str(NFlies),').Z.T5B']);
     Z_T4B= eval(['T4T5_mb_Bright(', num2str(NFlies),').Z.T4B']);

     Ang=[angle(Z_T5B),angle(Z_T4B)]; 
     ng=find(Ang)<=0; 
     Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4T5B_circstd(NFlies)=circ_std(Ang');
     Z_T4T5B_circste(NFlies)=circ_std(Ang')/sqrt(length(Ang));
     
     
     Z_T5C= eval(['T4T5_mb_Dark(', num2str(NFlies),').Z.T5C']);
     Z_T4C= eval(['T4T5_mb_Bright(', num2str(NFlies),').Z.T4C']);

     Ang=[angle(Z_T5C),angle(Z_T4C)]; 
     ng=find(Ang)<=0; 
     Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4T5C_circstd(NFlies)=circ_std(Ang');
     Z_T4T5C_circste(NFlies)=circ_std(Ang')/sqrt(length(Ang));
     
     
     
     
     Z_T5D= eval(['T4T5_mb_Dark(', num2str(NFlies),').Z.T5D']);
     Z_T4D= eval(['T4T5_mb_Bright(', num2str(NFlies),').Z.T4D']);

     Ang=[angle(Z_T5D),angle(Z_T4D)]; 
     ng=find(Ang)<=0; 
     Ang(ng)=Ang(ng)+2*pi; %Convert from -pi->pi scale to 0->2pi scale 
     Z_T4T5D_circstd(NFlies)=circ_std(Ang');
     Z_T4T5D_circste(NFlies)=circ_std(Ang')/sqrt(length(Ang));
     
     
     
end 
    
    Z.T4T5A_circstd=Z_T4T5A_circstd;
    Z.T4T5B_circstd=Z_T4T5B_circstd;
    Z.T4T5C_circstd=Z_T4T5C_circstd;
    Z.T4T5D_circstd=Z_T4T5D_circstd;

    Z.T4T5A_circste=Z_T4T5A_circste;
    Z.T4T5B_circste=Z_T4T5B_circste;
    Z.T4T5C_circste=Z_T4T5C_circste;
    Z.T4T5D_circste=Z_T4T5D_circste;
    
end 