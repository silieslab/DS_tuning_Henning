
function [F1,F2]= CompassPlot_Zdepth_8dirsEdge(cond, con_name, Average)

      
    Z_depth=[30,35,45,50,60,65,75,80];
    
    Z_A_color=[128,255,78;78,235,38;78,205,0;23,164,0;0,153,0;0,124,0;0,104,0;0,76,0]/255;
    Z_B_color=[190,210,255;141,158,255;101,116,255;71,86,255;0,0,232;0,0,154;0,0,119;0,0,78]/255;
    Z_C_color=[255,181,149;255,141,109;255,98,88;240,0,0;220,0,0;177,0,0;134,0,0;94,0,0]/255;
    Z_D_color=[254,254,159;254,254,43;254,244,0;234,234,0;214,214,0;194,194,0;174,174,0;144,144,0]/255;
    
          
%     Z_depth=[30,45,60,75];

%     Z_A_color=[78,235,38;78,205,0;23,164,0;0,153,0;0,124,0]/255;
%     Z_B_color=[71,86,255;0,0,232;0,0,154;0,0,119;0,0,78]/255;
%     Z_C_color=[255,141,109;255,88,78;240,0,0;177,0,0;134,0,0]/255;
%     Z_D_color=[254,254,0;234,234,0;214,214,0;194,194,0;174,174,0]/255;
    
    
    NCells_T4A=0;      NCells_T5A=0;    
    NCells_T4B=0;      NCells_T5B=0;
    NCells_T4C=0;      NCells_T5C=0;
    NCells_T4D=0;      NCells_T5D=0;


    %% Initiate compass plots 
    
    theta=[90, 45, 0, 315, 270, 225, 180, 130, 90];
    onedeg=2*pi/360; %in rad
    theta=theta*onedeg;
        
    R_teta=[1,0,0,0,0,0,0,0];
    L=sum(R_teta.*exp(1i*theta(1:end-1)))/sum(R_teta); %Einheitsvector length of one
        
   
    F1=figure(1);                     
    
    subplot(2,2,1)                      
    P=compass(L);                       
    set(P, 'Visible', 'off')
    hold on;
    
    subplot(2,2,2)
    P=compass(L);
    set(P, 'Visible', 'off')
    hold on
    
    subplot(2,2,3)
    P=compass(L);
    set(P, 'Visible', 'off')
    hold on
    
    subplot(2,2,4)
    P=compass(L);
    set(P, 'Visible', 'off')
    hold on
    
    
%     F2=figure(2);
%     
%     subplot(2,2,1)                      
%     P=compass(L);                       
%     set(P, 'Visible', 'off')
%     hold on;
%     
%     subplot(2,2,2)
%     P=compass(L);
%     set(P, 'Visible', 'off')
%     hold on
%     
%     subplot(2,2,3)
%     P=compass(L);
%     set(P, 'Visible', 'off')
%     hold on
%     
%     subplot(2,2,4)
%     P=compass(L);
%     set(P, 'Visible', 'off')
%     hold on
    
    
 %%   
    
    for N_depth=1:length(Z_depth)
        
        Z_cond= averageDirectionVectors_Zdepth_8dirEdge(cond.T4T5_mb,Z_depth(N_depth));
        
        
        
        if ~ Average
            Z=Z_cond;
            
            %Bright Stripe - T4
            figure(1)
            subplot(2,2,1)
            Comp=compass(Z.T4A.ALL);
            for i=1:length(Comp)
                set(Comp(i),'color',Z_A_color(N_depth,:));
            end
            Comp=compass(mean(Z.T4A.ALL),'k');
            NCells_T4A=NCells_T4A+length(Z.T4A.ALL); %Count cells
            
            subplot(2,2,2);
            Comp=compass(Z.T4B.ALL, 'b');
            for i=1:length(Comp)
                set(Comp(i),'color',Z_B_color(N_depth,:));
            end
            Comp=compass(mean(Z.T4B.ALL),'k');
            NCells_T4B=NCells_T4B+length(Z.T4B.ALL);

           
            subplot(2,2,3);
            Comp=compass(Z.T4C.ALL);
            for i=1:length(Comp)
                set(Comp(i),'color',Z_C_color(N_depth,:));
            end
            Comp=compass(mean(Z.T4C.ALL),'k');
            NCells_T4C=NCells_T4C+length(Z.T4C.ALL);  
            
            
            subplot(2,2,4);
            Comp=compass(Z.T4D.ALL);
            for i=1:length(Comp)
                set(Comp(i),'color', Z_D_color(N_depth,:));
            end
            Comp=compass(mean(Z.T4D.ALL),'k');
            NCells_T4D=NCells_T4D+length(Z.T4D.ALL);
          
            
            %Dark Stripe - T5
            Z=Z_cond;
%              figure(2)
            subplot(2,2,1)
            Comp=compass(Z.T5A.ALL);
            for i=1:length(Comp)
                set(Comp(i),'color', Z_A_color(N_depth,:));
            end
            Comp=compass(mean(Z.T5A.ALL),'k');
            NCells_T5A=NCells_T5A+length(Z.T5A.ALL);
            
            subplot(2,2,2);
            Comp=compass(Z.T5B.ALL, 'b');
            for i=1:length(Comp)
                set(Comp(i),'color', Z_B_color(N_depth,:));
            end
            Comp=compass(mean(Z.T5B.ALL),'k');
            NCells_T5B=NCells_T5B+length(Z.T5B.ALL);
            
            subplot(2,2,3);
            Comp=compass(Z.T5C.ALL);
            for i=1:length(Comp)
                set(Comp(i),'color', Z_C_color(N_depth,:));
            end
            Comp=compass(mean(Z.T5C.ALL),'k');
            NCells_T5C=NCells_T5C+length(Z.T5C.ALL);
            
            subplot(2,2,4);
            Comp=compass(Z.T5D.ALL);
            for i=1:length(Comp)
                set(Comp(i),'color', Z_D_color(N_depth,:));
            end
            Comp=compass(mean(Z.T5D.ALL),'k');
            NCells_T5D=NCells_T5D+length(Z.T5D.ALL);
            
        elseif Average
            
            Z=Z_cond;
            
            %Bright Stripe - T4
            figure(1)
            subplot(2,2,2)        
            Comp=compass(Z.T4A.M);
            for i=1:length(Comp)
                set(Comp(i),'color',Z_A_color(N_depth,:));
            end 
            title([con_name,'  T4:']);
            
            subplot(2,2,2);
            Comp=compass(Z.T4B.M, 'b');
            for i=1:length(Comp)
                set(Comp(i),'color',Z_B_color(N_depth,:));
            end
           
            subplot(2,2,3);
            Comp=compass(Z.T4C.M);
            for i=1:length(Comp)
                set(Comp(i),'color',Z_C_color(N_depth,:));
            end
            
            subplot(2,2,4);
            Comp=compass(Z.T4D.M);
            for i=1:length(Comp)
                set(Comp(i),'color', Z_D_color(N_depth,:));
            end
           
            
            
            %Dark Stripe - T5
            Z=Z_cond_OFF;
            figure(2)
            subplot(2,2,1)
            Comp=compass(Z.T5A.M);
            for i=1:length(Comp)
                set(Comp(i),'color',Z_A_color(N_depth,:));
            end
            title([con_name,' N Flies= ',num2str(length(Z.T5A.M)), '   T5:']);
            
            subplot(2,2,2);
            Comp=compass(Z.T5B.M, 'b');
            for i=1:length(Comp)
                set(Comp(i),'color',Z_B_color(N_depth,:));
            end
            
            
            subplot(2,2,3);
            Comp=compass(Z.T5C.M);
            for i=1:length(Comp)
                set(Comp(i),'color',Z_C_color(N_depth,:));
            end


            subplot(2,2,4);
            Comp=compass(Z.T5D.M);
            for i=1:length(Comp)
                set(Comp(i),'color', Z_D_color(N_depth,:));
            end

        end
        
    end
    
    
 if ~Average   
     
 figure(1)
 subplot(2,2,1)    
 title([con_name,'  T4+T5: LayerA - N Cells= ',num2str(NCells_T4A + NCells_T5A)]);

 subplot(2,2,2)   
 title(['LayerB - N Cells= ',num2str(NCells_T4B + NCells_T5B)]);

 subplot(2,2,3)
 title(['LayerC - N Cells= ',num2str(NCells_T4C + NCells_T5C)]);

 subplot(2,2,4)    
 title(['LayerD - N Cells= ',num2str(NCells_T4D + NCells_T5D)]);
 
 
 %% If you plot T4 and T5 seperate:
 
%  figure(1)
%  subplot(2,2,1)    
%  title([con_name,'  T4: LayerA - N Cells= ',num2str(NCells_T4A)]);
% 
%  subplot(2,2,2)   
%  title(['LayerB - N Cells= ',num2str(NCells_T4B)]);
% 
%  subplot(2,2,3)
%  title(['LayerC - N Cells= ',num2str(NCells_T4C)]);
% 
%  subplot(2,2,4)    
%  title(['LayerD - N Cells= ',num2str(NCells_T4D)]);
%  
%  figure(2)
%  subplot(2,2,1)    
%  title([con_name,'  T5: LayerA - N Cells= ',num2str(NCells_T5A)]);
% 
%  subplot(2,2,2)   
%  title(['LayerB - N Cells= ',num2str(NCells_T5B)]);
% 
%  subplot(2,2,3)
%  title(['LayerC - N Cells= ',num2str(NCells_T5C)]);
% 
%  subplot(2,2,4)    
%  title(['LayerD - N Cells= ',num2str(NCells_T5D)]);
%  
 
 end 

            
end