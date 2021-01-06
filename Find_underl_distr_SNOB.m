% Find underlying distributions: Snob Analysis

%% Linearize data 
Data_A=angle([Z.T5A.ALL,Z.T4A.ALL])'; % here I can just use the scale from -pi to pi, because most neurons in that Layer will fall into that scale
Data_B=convert_angle(angle([Z.T5B.ALL,Z.T4B.ALL]),'rad')';   % here I need to convert to a scale of 0 to 2pi 
Data_C=angle([Z.T5C.ALL,Z.T4C.ALL])'; %same as in A
Data_D=angle([Z.T5D.ALL,Z.T4D.ALL])'; %same as in A

% run the snob analysis:
mm_A = snob(Data_A, {'norm',1});
mm_B = snob(Data_B, {'norm',1});
mm_C = snob(Data_C, {'norm',1});
mm_D = snob(Data_D, {'norm',1});

% get Info from Snob:
mm_Summary(mm_A);
mm_Summary(mm_B);
mm_Summary(mm_C);
mm_Summary(mm_D);

% Plot the snob results:
figure
mm_PlotModel1d(mm_A,Data_A,1);

figure
mm_PlotModel1d(mm_B,Data_B,1);

figure
mm_PlotModel1d(mm_C,Data_C,1);

figure
mm_PlotModel1d(mm_D,Data_D,1);


%find the class with highest probability for each neuron 
[~, M_A]=max(mm_A.r'); 
[~, M_B]=max(mm_B.r'); 
[~, M_C]=max(mm_C.r'); 
[~, M_D]=max(mm_D.r'); 


TA_T5=M_A(1:length(Z.T5A.ALL));
TA_T4=M_A(length(Z.T5A.ALL)+1:end);

TB_T5=M_B(1:length(Z.T5B.ALL));
TB_T4=M_B(length(Z.T5B.ALL)+1:end);

TC_T5=M_C(1:length(Z.T5C.ALL));
TC_T4=M_C(length(Z.T5C.ALL)+1:end);

TD_T5=M_D(1:length(Z.T5D.ALL));
TD_T4=M_D(length(Z.T5D.ALL)+1:end);


%% 
% each time you run the analysis the assignment of class1, 2 or 3 will
% differ, therefore I saved the outcome of the analysis in Snob_Cluster_Info.mat


ClusterR.TA_T4=TA_T4;
ClusterR.TA_T5=TA_T5;
ClusterR.TB_T4=TB_T4;
ClusterR.TB_T5=TB_T5;
ClusterR.TC_T4=TC_T4;
ClusterR.TC_T5=TC_T5;
ClusterR.TD_T4=TD_T4;
ClusterR.TD_T5=TD_T5;
ClusterR.mm_A=mm_A; 
ClusterR.mm_B=mm_B;
ClusterR.mm_C=mm_C; 
ClusterR.mm_D=mm_D; 

%%


