
function [F1,F2]= CompassPlot(cond_ON, cond_OFF,con_name, Average)

if nargin==4
    
    
    
    Z_cond_OFF= averageDirectionVectors(cond_OFF.T4T5_mb);
    Z_cond_ON= averageDirectionVectors(cond_ON.T4T5_mb);
    
    theta=[90, 45, 0, 315, 270, 225, 180, 130, 90];
    onedeg=2*pi/360; %in rad
    theta=theta*onedeg;
    
    R_teta=[1,0,0,0,0,0,0,0];
    L=sum(R_teta.*exp(1i*theta(1:end-1)))/sum(R_teta);
    
    if ~ Average
        Z=Z_cond_ON;
        
        %Bright Stripe - T4
        F1=figure ;
        subplot(2,2,1)
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        
        Comp=compass(Z.T4A.ALL);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        Comp=compass(mean(Z.T4A.ALL),'k');
        title([con_name,'  T4: LayerA - N Cells= ',num2str(length(Z.T4A.ALL))]);
        
        
        subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4B.ALL, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        Comp=compass(mean(Z.T4B.ALL),'k');
        title(['LayerB - N Cells= ',num2str(length(Z.T4B.ALL))]);
        
        
        subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4C.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        Comp=compass(mean(Z.T4C.ALL),'k');
        title(['LayerC - N Cells= ',num2str(length(Z.T4C.ALL))]);
        
        
        subplot(2,2,4);
        P=compass(L);set(P, 'Visible', 'off');hold on
        Comp=compass(Z.T4D.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        Comp=compass(mean(Z.T4D.ALL),'k');
        title(['LayerD - N Cells= ',num2str(length(Z.T4D.ALL))]);
        
        
        
        %Dark Stripe - T5
        Z=Z_cond_OFF;
        F2=figure;
        subplot(2,2,1)
        P=compass(L); set(P, 'Visible', 'off');hold on ;
        Comp=compass(Z.T5A.ALL);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        Comp=compass(mean(Z.T5A.ALL),'k');
        title([con_name,'  T5: LayerA - N Cells= ',num2str(length(Z.T5A.ALL))]);
        
        subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5B.ALL, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        Comp=compass(mean(Z.T5B.ALL),'k');
        title(['LayerB - N Cells= ',num2str(length(Z.T5B.ALL))]);
        
        
        subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5C.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        Comp=compass(mean(Z.T5C.ALL),'k');
        title(['LayerC - N Cells= ',num2str(length(Z.T5C.ALL))]);
        
        
        subplot(2,2,4);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5D.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        Comp=compass(mean(Z.T5D.ALL),'k');
        title(['LayerD - N Cells= ',num2str(length(Z.T5D.ALL))]);
        
    elseif Average
        
        Z=Z_cond_ON;
        
        %Bright Stripe - T4
        F1=figure ;
        subplot(1,2,2)
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        
        Comp=compass(Z.T4A.M);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        %     Comp=compass(mean(Z.T4A.M),'k');
        %     title([con_name,'  T4: N Flies= ',num2str(length(Z.T4A.M))]);
        title([con_name,'  T4:']);
        
        %     subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4B.M, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        %     Comp=compass(mean(Z.T4B.M),'k');
        %     title(['LayerB - N Flies= ',num2str(length(Z.T4B.M))]);
        
        
        %     subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4C.M);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        %     Comp=compass(mean(Z.T4C.M),'k');
        %     title(['LayerC - N Flies= ',num2str(length(Z.T4C.M))]);
        
        
        %     subplot(2,2,4);
        P=compass(L);set(P, 'Visible', 'off');hold on
        Comp=compass(Z.T4D.M);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        %     Comp=compass(mean(Z.T4D.M),'k');
        %     title(['LayerD - N Flies= ',num2str(length(Z.T4D.M))]);
        
        
        
        %Dark Stripe - T5
        Z=Z_cond_OFF;
        %     F2=figure;
        subplot(1,2,1)
        P=compass(L); set(P, 'Visible', 'off');hold on ;
        Comp=compass(Z.T5A.M);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        %     Comp=compass(mean(Z.T5A.M),'k');
        title([con_name,' N Flies= ',num2str(length(Z.T5A.M)), '   T5:']);
        
        %     subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5B.M, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        %     Comp=compass(mean(Z.T5B.M),'k');
        %     title(['LayerB - N Flies= ',num2str(length(Z.T5B.M))]);
        
        
        %     subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5C.M);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        %     Comp=compass(mean(Z.T5C.M),'k');
        %     title(['LayerC - N Flies= ',num2str(length(Z.T5C.M))]);
        
        
        %     subplot(2,2,4);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5D.M);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        %     Comp=compass(mean(Z.T5D.M),'k');
        %     title(['LayerD - N Flies= ',num2str(length(Z.T5D.M))]);
        
        F2=figure;
    end
elseif nargin<4
    
    Condition=cond_ON;
    Average=con_name;
    con_name=cond_OFF;
    
    Z_cond_ON= averageDirectionVectors(Condition.T4T5_mb);
    Z_cond_OFF= averageDirectionVectors(Condition.T4T5_mb);
    
    theta=[90, 67.5, 45, 22.5, 0, 337.5,315,292.5, 270,247.5, 225, 202.5, 180,157.5, 135, 112.5,90];
    onedeg=2*pi/360; %in rad
    theta=theta*onedeg;
    
    R_teta=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    L=sum(R_teta.*exp(1i*theta(1:end-1)))/sum(R_teta);
    
    if ~ Average
        Z=Z_cond_ON;
        
        %Bright Stripe - T4
        F1=figure ;
        subplot(2,2,1)
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        
        Comp=compass(Z.T4A.ALL);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        Comp=compass(mean(Z.T4A.ALL),'k');
        title([con_name,'  T4: LayerA - N Cells= ',num2str(length(Z.T4A.ALL))]);
        
        
        subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4B.ALL, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        Comp=compass(mean(Z.T4B.ALL),'k');
        title(['LayerB - N Cells= ',num2str(length(Z.T4B.ALL))]);
        
        
        subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4C.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        Comp=compass(mean(Z.T4C.ALL),'k');
        title(['LayerC - N Cells= ',num2str(length(Z.T4C.ALL))]);
        
        
        subplot(2,2,4);
        P=compass(L);set(P, 'Visible', 'off');hold on
        Comp=compass(Z.T4D.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        Comp=compass(mean(Z.T4D.ALL),'k');
        title(['LayerD - N Cells= ',num2str(length(Z.T4D.ALL))]);
        
        
        
        %Dark Stripe - T5
        Z=Z_cond_OFF;
        F2=figure;
        subplot(2,2,1)
        P=compass(L); set(P, 'Visible', 'off');hold on ;
        Comp=compass(Z.T5A.ALL);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        Comp=compass(mean(Z.T5A.ALL),'k');
        title([con_name,'  T5: LayerA - N Cells= ',num2str(length(Z.T5A.ALL))]);
        
        subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5B.ALL, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        Comp=compass(mean(Z.T5B.ALL),'k');
        title(['LayerB - N Cells= ',num2str(length(Z.T5B.ALL))]);
        
        
        subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5C.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        Comp=compass(mean(Z.T5C.ALL),'k');
        title(['LayerC - N Cells= ',num2str(length(Z.T5C.ALL))]);
        
        
        subplot(2,2,4);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5D.ALL);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        Comp=compass(mean(Z.T5D.ALL),'k');
        title(['LayerD - N Cells= ',num2str(length(Z.T5D.ALL))]);
        
        
        
    elseif Average
        
        Z=Z_cond_ON;
        
        %Bright Stripe - T4
        F1=figure ;
        subplot(1,2,2)
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        
        Comp=compass(Z.T4A.M);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        %     Comp=compass(mean(Z.T4A.M),'k');
        %     title([con_name,'  T4: N Flies= ',num2str(length(Z.T4A.M))]);
        title([con_name,'  T4:']);
        
        %     subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4B.M, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        %     Comp=compass(mean(Z.T4B.M),'k');
        %     title(['LayerB - N Flies= ',num2str(length(Z.T4B.M))]);
        
        
        %     subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T4C.M);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        %     Comp=compass(mean(Z.T4C.M),'k');
        %     title(['LayerC - N Flies= ',num2str(length(Z.T4C.M))]);
        
        
        %     subplot(2,2,4);
        P=compass(L);set(P, 'Visible', 'off');hold on
        Comp=compass(Z.T4D.M);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        %     Comp=compass(mean(Z.T4D.M),'k');
        %     title(['LayerD - N Flies= ',num2str(length(Z.T4D.M))]);
        
        
        
        %Dark Stripe - T5
        Z=Z_cond_OFF;
        %     F2=figure;
        subplot(1,2,1)
        P=compass(L); set(P, 'Visible', 'off');hold on ;
        Comp=compass(Z.T5A.M);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0]);
        end
        %     Comp=compass(mean(Z.T5A.M),'k');
        title([con_name,' N Flies= ',num2str(length(Z.T5A.M)), '   T5:']);
        
        %     subplot(2,2,2);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5B.M, 'b');
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]);
        end
        %     Comp=compass(mean(Z.T5B.M),'k');
        %     title(['LayerB - N Flies= ',num2str(length(Z.T5B.M))]);
        
        
        %     subplot(2,2,3);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5C.M);
        for i=1:length(Comp);
            set(Comp(i),'color',[1 0 0]);
        end
        %     Comp=compass(mean(Z.T5C.M),'k');
        %     title(['LayerC - N Flies= ',num2str(length(Z.T5C.M))]);
        
        
        %     subplot(2,2,4);
        P=compass(L);
        set(P, 'Visible', 'off')
        hold on
        Comp=compass(Z.T5D.M);
        for i=1:length(Comp);
            set(Comp(i),'color', [.7 .7 0]);
        end
        %     Comp=compass(mean(Z.T5D.M),'k');
        %     title(['LayerD - N Flies= ',num2str(length(Z.T5D.M))]);
        
        F2=figure;
        
    end
    
    
    
    
end