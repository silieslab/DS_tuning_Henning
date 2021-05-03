
%%Save Data for Luis

% Coord_T5A_I=Coord_T5A_all(:,TA_T5==1); 
% Coord_T5A_II=Coord_T5A_all(:,TA_T5==3); 
% 
% Coord_T4A_I=Coord_T4A_all(:,TA_T4==1); 
% Coord_T4A_II=Coord_T4A_all(:,TA_T4==3); 
% 
% Coord_T5B_I=Coord_T5B_all(:,TB_T5==1); 
% Coord_T5B_II=Coord_T5B_all(:,TB_T5==2); 
% 
% Coord_T4B_I=Coord_T4B_all(:,TB_T4==1); 
% Coord_T4B_II=Coord_T4B_all(:,TB_T4==2); 
% 
% Coord_T5C=Coord_T5C_all; 
% Coord_T4C=Coord_T4C_all; 
% 
% Coord_T5D=Coord_T5D_all; 
% Coord_T4D=Coord_T4D_all; 
% 
% 
% 
% Quiver_T5A_I=Quiver_T5A_all_s(TA_T5==1); 
% Quiver_T5A_II=Quiver_T5A_all_s(TA_T5==3); 
% 
% Quiver_T4A_I=Quiver_T4A_all_s(TA_T4==1); 
% Quiver_T4A_II=Quiver_T4A_all_s(TA_T4==3); 
% 
% Quiver_T5B_I=Quiver_T5B_all_s(TB_T5==1); 
% Quiver_T5B_II=Quiver_T5B_all_s(TB_T5==2); 
% 
% Quiver_T4B_I=Quiver_T4B_all_s(TB_T4==1); 
% Quiver_T4B_II=Quiver_T4B_all_s(TB_T4==2); 
% 
% Quiver_T5C=Quiver_T5C_all_s; 
% Quiver_T4C=Quiver_T4C_all_s; 
% 
% Quiver_T5D=Quiver_T5D_all_s; 
% Quiver_T4D=Quiver_T4D_all_s; 



F9=figure('Position',[200 200 1466 814]);
subplot(3,2,1)
quiver(Coord_T5A_I(2,:)-34,Coord_T5A_I(1,:)+36,real(Quiver_T5A_I),imag(Quiver_T5A_I), 'AutoScale','off')
hold on 
quiver(Coord_T4A_I(2,:)-34,Coord_T4A_I(1,:)+36,real(Quiver_T4A_I),imag(Quiver_T4A_I), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerA.I'])


subplot(3,2,2)
quiver(Coord_T5A_II(2,:)-34,Coord_T5A_II(1,:)+36,real(Quiver_T5A_II),imag(Quiver_T5A_II), 'AutoScale','off')
hold on 
quiver(Coord_T4A_II(2,:)-34,Coord_T4A_II(1,:)+36,real(Quiver_T4A_II),imag(Quiver_T4A_II), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerA.II'])


subplot(3,2,3)
quiver(Coord_T5B_I(2,:)-34,Coord_T5B_I(1,:)+36,real(Quiver_T5B_I),imag(Quiver_T5B_I), 'AutoScale','off')
hold on 
quiver(Coord_T4B_I(2,:)-34,Coord_T4B_I(1,:)+36,real(Quiver_T4B_I),imag(Quiver_T4B_I), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerB.I'])



subplot(3,2,4)
quiver(Coord_T5B_II(2,:)-34,Coord_T5B_II(1,:)+36,real(Quiver_T5B_II),imag(Quiver_T5B_II), 'AutoScale','off')
hold on 
quiver(Coord_T4B_II(2,:)-34,Coord_T4B_II(1,:)+36,real(Quiver_T4B_II),imag(Quiver_T4B_II), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerB.II'])


subplot(3,2,5)
quiver(Coord_T5C(2,:)-34,Coord_T5C(1,:)+36,real(Quiver_T5C),imag(Quiver_T5C), 'AutoScale','off')
hold on 
quiver(Coord_T4C(2,:)-34,Coord_T4C(1,:)+36,real(Quiver_T4C),imag(Quiver_T4C), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerC'])


subplot(3,2,6)
quiver(Coord_T5D(2,:)-34,Coord_T5D(1,:)+36,real(Quiver_T5D),imag(Quiver_T5D), 'AutoScale','off')
hold on 
quiver(Coord_T4D(2,:)-34,Coord_T4D(1,:)+36,real(Quiver_T4D),imag(Quiver_T4D), 'AutoScale','off')
axis('equal')
set(gca,'XLim', [-34,44+2*45])
set(gca,'YLim', [-17,36])
title(['LayerD'])



% save('Flow_field_Data_new', 'Coord_T5A_I', 'Coord_T5A_II', 'Coord_T4A_I', 'Coord_T4A_II',...
%                             'Coord_T5B_I', 'Coord_T5B_II', 'Coord_T4B_I', 'Coord_T4B_II',...
%                             'Coord_T5C', 'Coord_T4C', 'Coord_T5D', 'Coord_T4D', ...
%                             'Quiver_T5A_I', 'Quiver_T5A_II', 'Quiver_T4A_I', 'Quiver_T4A_II', ...
%                             'Quiver_T5B_I', 'Quiver_T5B_II', 'Quiver_T4B_I', 'Quiver_T4B_II',...
%                             'Quiver_T5C','Quiver_T4C','Quiver_T5D','Quiver_T4D')
% 


