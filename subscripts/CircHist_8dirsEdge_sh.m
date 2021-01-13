function [F2, F3]=CircHist_8dirsEdge_sh(cond)
%%
%% Plot the circular histogram (Figure 1b):
nbins= 72;

Z=averageDirectionVectors(cond.T4T5_mb);
Z_All_A=[Z.T4A.ALL,Z.T5A.ALL];
Z_All_B=[Z.T4B.ALL,Z.T5B.ALL];
Z_All_C=[Z.T4C.ALL,Z.T5C.ALL];
Z_All_D=[Z.T4D.ALL,Z.T5D.ALL];
Z_All=[Z_All_A,Z_All_B,Z_All_C,Z_All_D];

Colors=[23,106,0;0,0,154;240,0,0;214,214,0];

% Plot Figure 1b 
F2=figure('Position', [200 200 1200 900]);
subAx1 = subplot(2, 2, 1, polaraxes);
obj1=CircHist(convert_angle(angle(Z_All_A)),nbins,'parent', subAx1);
title('T4\T5 tuning distribution- layerA', 'FontSize', 14)
subAx2 = subplot(2, 2, 2, polaraxes);
obj2=CircHist(convert_angle(angle(Z_All_B)),nbins,'parent', subAx2);
title('layerB', 'FontSize', 14)
subAx3 = subplot(2, 2, 3, polaraxes);
obj3=CircHist(convert_angle(angle(Z_All_C)),nbins,'parent', subAx3);
title('layerC', 'FontSize', 14)
subAx4 = subplot(2, 2, 4, polaraxes);
obj4=CircHist(convert_angle(angle(Z_All_D)),nbins,'parent', subAx4);
title('layerD', 'FontSize', 14)

for ax=1:4
currobj=eval(['obj', num2str(ax)]);
% Adjust appearance:
currobj.colorBar = Colors(ax,:)/255;
 % change color of bars
currobj.avgAngH.LineStyle = '--'; % make average-angle line dashed
currobj.avgAngH.LineWidth = 1; % make average-angle line thinner
currobj.colorAvgAng = [.5 .5 .5]; % change average-angle line color
% remove offset between bars and plot-center
rl = rlim; % get current limits
currobj.setRLim([-20, 100]); % set lower limit to 0
% draw circle at r == 0.5 (where r == 1 would be the outer plot edge)
rl = rlim;
% obj1.drawCirc((rl(2) - rl(1)) /2, '--b', 'LineWidth', 2)
currobj.scaleBarSide = 'right'; % draw rho-axis on the right side of the plot
currobj.polarAxs.ThetaZeroLocation = 'right'; % rotate the plot to have 0� on the right side
% draw resultant vector r as arrow
delete(currobj.rH)
currobj.drawArrow(currobj.avgAng, currobj.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'k')

end 

%% Split Layer A and B into subtypes based on Snob analysis (Finite gaussian mixture model)
Colorsn=[ 23,106,0; 151,217,0; 0,0,154; 117,0,154; 240,0,0; 214,214,0];

%load snob cluster Info 
load('Data/Data_Edges/Snob_Cluster_Info.mat')

TA_T5=ClusterR.TA_T5; 
TA_T4=ClusterR.TA_T4; 

TB_T5=ClusterR.TB_T5; 
TB_T4=ClusterR.TB_T4;

TC_T5=ClusterR.TC_T5; 
TC_T4=ClusterR.TC_T4; 

TD_T5=ClusterR.TD_T5; 
TD_T4=ClusterR.TD_T4; 
%

%% Plot the circular histogram for each subtype (Figure 1c):
% all subplots were later merged in illustrator 

F3=figure('Position', [200 200 900 900]);
subAx1 = subplot(3, 2, 1, polaraxes);
obj1=CircHist(convert_angle(angle([Z.T5A.ALL(TA_T5==3),Z.T4A.ALL(TA_T4==3)])),nbins,'parent', subAx1);
title('subtype A.I', 'FontSize', 14)

subAx2 = subplot(3, 2, 2, polaraxes);
obj2=CircHist(convert_angle(angle([Z.T5A.ALL(TA_T5==1),Z.T4A.ALL(TA_T4==1)])),nbins,'parent', subAx2);
title('subtype A.II', 'FontSize', 14)

subAx3 = subplot(3, 2, 3, polaraxes);
obj3=CircHist(convert_angle(angle([Z.T5B.ALL(TB_T5==1),Z.T4B.ALL(TB_T4==1)])),nbins,'parent', subAx3);
title('subtype B.I', 'FontSize', 14)

subAx4 = subplot(3, 2, 4, polaraxes);
obj4=CircHist(convert_angle(angle([Z.T5B.ALL(TB_T5==2),Z.T4B.ALL(TB_T4==2)])),nbins,'parent', subAx4);
title('subtype B.II', 'FontSize', 14)

subAx5 = subplot(3, 2, 5, polaraxes);
obj5=CircHist(convert_angle(angle([Z.T5C.ALL(~(TC_T5==2)),Z.T4C.ALL(~(TC_T4==2))])),nbins,'parent', subAx5);
title('subtype C', 'FontSize', 14)

subAx6 = subplot(3, 2, 6, polaraxes);
obj6=CircHist(convert_angle(angle([Z.T5D.ALL(TD_T5==1),Z.T4D.ALL(TD_T4==1)])),nbins,'parent', subAx6);
title('subtype D', 'FontSize', 14)


Colorsi=[Colorsn(1,:);Colorsn(2,:);Colorsn(3,:);Colorsn(4,:);Colorsn(5,:); Colorsn(6,:) ];

for ax=1:6
currobj=eval(['obj', num2str(ax)]);
% Adjust appearance:
currobj.colorBar = Colorsi(ax,:)/255;
 % change color of bars
currobj.avgAngH.LineStyle = '--'; % make average-angle line dashed
currobj.avgAngH.LineWidth = 1; % make average-angle line thinner
currobj.colorAvgAng = [.5 .5 .5]; % change average-angle line color
% remove offset between bars and plot-center
rl = rlim; % get current limits
currobj.setRLim([-20, 100]); % set lower limit to 0
% draw circle at r == 0.5 (where r == 1 would be the outer plot edge)
% rl = rlim;
% obj1.drawCirc((rl(2) - rl(1)) /2, '--b', 'LineWidth', 2)
currobj.scaleBarSide = 'right'; % draw rho-axis on the right side of the plot
currobj.polarAxs.ThetaZeroLocation = 'right'; % rotate the plot to have 0� on the right side
% draw resultant vector r as arrow
delete(currobj.rH)
currobj.drawArrow(currobj.avgAng, currobj.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'k')

end 

F3.Renderer='Painters';
%% Compass plots for each subtype (Extended Data figure d)
% separat for T4 and T5

F6=figure('Position', [200,200,1000,500])
subplot(2,4,1)
Zd1=compass(Z.T5A.ALL(TA_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(2,:)/255)
end  
hold on 
Zd2=compass(Z.T5A.ALL(TA_T5==3));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(1,:)/255)
end 
Zd3=compass(Z.T5A.ALL(TA_T5==2));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
L=legend([Zd1(1),Zd2(1),Zd3(1)], {'subtype A.II', 'subtype A.I', 'noise'});
L.Position(1)=0.095;
L.Position(2)=0.85;
title('T5A')
% 
subplot(2,4,5)
Zd1=compass(Z.T4A.ALL(TA_T4==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(2,:)/255)
end  
hold on 
Zd2=compass(Z.T4A.ALL(TA_T4==3));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(1,:)/255)
end 
Zd3=compass(Z.T4A.ALL(TA_T4==2));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
% % 
 title('T4A')
% % 
% % 
subplot(2,4,2)
Zd1=compass(Z.T5B.ALL(TB_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(3,:)/255)
end  
hold on 
Zd2=compass(Z.T5B.ALL(TB_T5==2));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(4,:)/255)
end 
Zd3=compass(Z.T5B.ALL(TB_T5==3));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
L=legend([Zd1(1),Zd2(1),Zd3(1)], {'subtype B.I', 'subtype B.II', 'noise'});
L.Position(1)=0.3;
L.Position(2)=0.85;
title('T5B')
% 

subplot(2,4,6)
Zd=compass(Z.T4B.ALL(TB_T4==1));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(3,:)/255)
end 
hold on 
Zd=compass(Z.T4B.ALL(TB_T4==2));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(4,:)/255)
end 
Zd=compass(Z.T4B.ALL(TB_T4==3));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[0.5,0.5,0.5])
end 
 title('T4B')
% 

subplot(2,4,3)
Zd1=compass(Z.T5C.ALL(TC_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',[133 0 0]/255)
end 
hold on 
Zd2=compass(Z.T5C.ALL(TC_T5==3));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(5,:)/255)
end 
Zd3=compass(Z.T5C.ALL(TC_T5==2))
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
L=legend([Zd1(1),Zd2(1),Zd3(1)], {'subtype C', 'subtype C', 'noise'});
L.Position(1)=0.51;
L.Position(2)=0.85;
title('T5C')
% 
subplot(2,4,7)
Zd=compass(Z.T4C.ALL(TC_T4==1));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[133 0 0]/255)
end
hold on 
Zd=compass(Z.T4C.ALL(TC_T4==3));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(5,:)/255)
end 
Zd=compass(Z.T4C.ALL(TC_T4==2));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[0.5,0.5,0.5])
end 
 title('T4C')

 
subplot(2,4,4)
Zd1=compass(Z.T5D.ALL(TD_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(6,:)/255)
end 
hold on 
Zd2=compass(Z.T5D.ALL(TD_T5==2));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',[0.5,0.5,0.5])
end 
L=legend([Zd1(1),Zd2(1)], {'subtype D', 'noise'});
L.Position(1)=0.72;
L.Position(2)=0.87;
title('T5D')


subplot(2,4,8)
Zd=compass(Z.T4D.ALL(TD_T4==1));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(6,:)/255)
end 
hold on 
Zd=compass(Z.T4D.ALL(TD_T4==2));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[0.5,0.5,0.5])
end 

title('T4D')
% 
F6.Renderer='Painters';


%% T4/T5 plotted together (Figure 1d)

F7=figure('Position', [200,200,1000,300]);
subplot(1,4,1)
Zd1=compass(Z.T5A.ALL(TA_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(2,:)/255)
end  
hold on 
Zd1=compass(Z.T4A.ALL(TA_T4==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(2,:)/255)
end  
compass(mean([Z.T4A.ALL(TA_T4==1),Z.T5A.ALL(TA_T5==1)]), 'k')

Zd2=compass(Z.T5A.ALL(TA_T5==3));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(1,:)/255)
end 
Zd2=compass(Z.T4A.ALL(TA_T4==3));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(1,:)/255)
end 
compass(mean([Z.T4A.ALL(TA_T4==3),Z.T5A.ALL(TA_T5==3)]), 'k')

Zd3=compass(Z.T5A.ALL(TA_T5==2));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
Zd3=compass(Z.T4A.ALL(TA_T4==2));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
compass(mean([Z.T4A.ALL(TA_T4==2),Z.T5A.ALL(TA_T5==2)]), 'k')

L=legend([Zd1(1),Zd2(1),Zd3(1)], {'subtype A.II', 'subtype A.I', 'noise'});
L.Position(1)=0.095;
L.Position(2)=0.85;
text(-2,-1.5,'T4/T5 tuning distribution of all subtypes', 'Fontsize', 12)

% 

subplot(1,4,2)
Zd1=compass(Z.T5B.ALL(TB_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(3,:)/255)
end  
hold on 
Zd=compass(Z.T4B.ALL(TB_T4==1));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(3,:)/255)
end 
compass(mean([Z.T4B.ALL(TB_T4==1),Z.T5B.ALL(TB_T5==1)]), 'k')

Zd2=compass(Z.T5B.ALL(TB_T5==2));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(4,:)/255)
end 
Zd=compass(Z.T4B.ALL(TB_T4==2));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(4,:)/255)
end 
compass(mean([Z.T4B.ALL(TB_T4==2),Z.T5B.ALL(TB_T5==2)]), 'k')

Zd3=compass(Z.T5B.ALL(TB_T5==3));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
Zd=compass(Z.T4B.ALL(TB_T4==3));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[0.5,0.5,0.5])
end 
compass(mean([Z.T4B.ALL(TB_T4==3),Z.T5B.ALL(TB_T5==3)]), 'k')

L=legend([Zd1(1),Zd2(1),Zd3(1)], {'subtype B.I', 'subtype B.II', 'noise'});
L.Position(1)=0.3;
L.Position(2)=0.85;
% 

subplot(1,4,3)
Zd1=compass(Z.T5C.ALL(TC_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',[133 0 0]/255)
end 
hold on 
Zd=compass(Z.T4C.ALL(TC_T4==1));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[133 0 0]/255)
end
compass(mean([Z.T4C.ALL(TC_T4==1),Z.T5C.ALL(TC_T5==1)]), 'k')

Zd2=compass(Z.T5C.ALL(TC_T5==3));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',Colorsn(5,:)/255)
end 
Zd=compass(Z.T4C.ALL(TC_T4==3));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(5,:)/255)
end 
compass(mean([Z.T4C.ALL(TC_T4==3),Z.T5C.ALL(TC_T5==3)]), 'k')

Zd3=compass(Z.T5C.ALL(TC_T5==2));
for Zdi=1:length(Zd3)
set(Zd3(Zdi),'Color',[0.5,0.5,0.5])
end 
Zd=compass(Z.T4C.ALL(TC_T4==2));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[0.5,0.5,0.5])
end 
compass(mean([Z.T4C.ALL(TC_T4==2),Z.T5C.ALL(TC_T5==2)]), 'k')

L=legend([Zd1(1),Zd2(1),Zd3(1)], {'subtype C', 'subtype C', 'noise'});
L.Position(1)=0.51;
L.Position(2)=0.85;


subplot(1,4,4)
Zd1=compass(Z.T5D.ALL(TD_T5==1));
for Zdi=1:length(Zd1)
set(Zd1(Zdi),'Color',Colorsn(6,:)/255)
end 
hold on 
Zd=compass(Z.T4D.ALL(TD_T4==1));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',Colorsn(6,:)/255)
end 
compass(mean([Z.T4D.ALL(TD_T4==1),Z.T5D.ALL(TD_T5==1)]), 'k')

Zd2=compass(Z.T5D.ALL(TD_T5==2));
for Zdi=1:length(Zd2)
set(Zd2(Zdi),'Color',[0.5,0.5,0.5])
end 
Zd=compass(Z.T4D.ALL(TD_T4==2));
for Zdi=1:length(Zd)
set(Zd(Zdi),'Color',[0.5,0.5,0.5])
end 
compass(mean([Z.T4D.ALL(TD_T4==2),Z.T5D.ALL(TD_T5==2)]), 'k')

L=legend([Zd1(1),Zd2(1)], {'subtype D', 'noise'});
L.Position(1)=0.72;
L.Position(2)=0.87;


F7.Renderer='Painters';


end 