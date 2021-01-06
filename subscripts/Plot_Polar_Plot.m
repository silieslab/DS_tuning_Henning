function F=Plot_Polar_Plot(cond_ON,cond_OFF, av_Vector)



COLOR_OF_PLOT_GREEN = [0 .5 0];% GREEN
% COLOR_CLOUD =[.5 .7 .5];% GREEN

COLOR_OF_PLOT_RED = [1 0 0]; %RED
% COLOR_CLOUD = [1 .5 .5]; %RED cloud

COLOR_OF_PLOT_YELLOW = [.7 .7 0];% YELLOW
COLOR_CLOUD = [.9 .9 .5];%YELLOW cloud

COLOR_OF_PLOT_BLUE = [0 0 1]; %BLUE;
% COLOR_CLOUD = [.5 .5 1]; %blue cloud
Color=[COLOR_OF_PLOT_GREEN;COLOR_OF_PLOT_BLUE;COLOR_OF_PLOT_RED;COLOR_OF_PLOT_YELLOW];


if nargin>2
    
    theta=[90, 45, 0, 315, 270, 225, 180, 130, 90];
    onedeg=2*pi/360; %in rad
    theta=theta*onedeg;
    F=figure('Color', [1 1 1],'Position', [400 400 800 400]);
    subplot(1,2,1)
    
    P = polar(theta, 4 * ones(size(theta)));
    set(P, 'Visible', 'off')
    hold on
    for i=1:4 %four Layers
        LAYERS=['A','B','C','D'];
        
        tuningPdirPFly=squeeze(max(eval(['cond_OFF.ALL_Flies.T5',LAYERS(i)]),[],2)); %tuning (strength of response) per direction and per fly
        mtuningPdir=nanmean(tuningPdirPFly,2);
        stuningPdir=nanstd(tuningPdirPFly')/sqrt(length(tuningPdirPFly));
        rho=[mtuningPdir;mtuningPdir(1)]';
        PPp(i)=polar(theta,rho);%'Color',Color(i,:));
        PPp(i).Color=Color(i,:);
        hold on
        PolarPlot_std_new
    end
    [PPp.LineWidth] = deal(1.5);
    LegendName={'LayerA', 'LayerB', 'LayerC', 'LayerD'};
    legend(PPp,LegendName,'Location', 'southwest');
    
    
    set(gca,'YTickLabel','')
    set(gca,'FontSize',13)
    
    if av_Vector
        
        Z_Cond_OFF= averageDirectionVectors(cond_OFF.T4T5_mb);
        
        Comp=compass(mean(Z_Cond_OFF.T5A.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0],'LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_OFF.T5B.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_OFF.T5C.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[1 0 0]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_OFF.T5D.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color', [.7 .7 0]','LineWidth',2);
        end
        
    end
    
    
    
    
    subplot(1,2,2)
    
    
    P = polar(theta, 4 * ones(size(theta)));
    set(P, 'Visible', 'off')
    hold on
    for i=1:4 %four Layers
        LAYERS=['A','B','C','D'];
        
        tuningPdirPFly=squeeze(max(eval(['cond_ON.ALL_Flies.T4',LAYERS(i)]),[],2)); %tuning (strength of response) per direction and per fly
        mtuningPdir=nanmean(tuningPdirPFly,2);
        stuningPdir=nanstd(tuningPdirPFly')/sqrt(length(tuningPdirPFly));
        rho=[mtuningPdir;mtuningPdir(1)]';
        PPp(i)=polar(theta,rho);
        PPp(i).Color=Color(i,:);
        set(PPp(i),'LineWidth', 1.5)
        
        % rlim([-5 15])
        %set(PPp,'FontSize',24)
        hold on
        
        PolarPlot_std_new
        
    end
    title('T4 responses to Bright Stripes');
    set(gca,'YTickLabel','')
    set(gca,'FontSize',13)
    set(F,'PaperSize', [8 8])
    
    
    
    if av_Vector
        
        Z_Cond_ON= averageDirectionVectors(cond_OFF.T4T5_mb);
        
        Comp=compass(mean(Z_Cond_ON.T4A.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0],'LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_ON.T4B.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_ON.T4C.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[1 0 0]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_ON.T4D.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color', [.7 .7 0]','LineWidth',2);
        end
        
    end
    
    
else
    Condition=cond_ON;
    av_Vector=cond_OFF;
    Ndirs=size(squeeze(max(eval(['Condition.ALL_Flies.T5A']),[],2)),1); %tuning (strength of response) per direction and per fly

    if Ndirs==8
        theta=[90, 45, 0, 315, 270, 225, 180, 130, 90];
    elseif Ndirs==16
        theta=[90, 67.5, 45, 22.5, 0, 337.5,315,292.5, 270,247.5, 225, 202.5, 180,157.5, 135, 112.5,90];
    end
        
    onedeg=2*pi/360; %in rad
    theta=theta*onedeg;
    F=figure('Color', [1 1 1],'Position', [400 400 800 400]);
    subplot(1,2,1)
    
    P = polar(theta, 4 * ones(size(theta)));
    set(P, 'Visible', 'off')
    hold on
    for i=1:4 %four Layers
        LAYERS=['A','B','C','D'];
        tuningPdirPFly=squeeze(max(eval(['Condition.ALL_Flies.T5',LAYERS(i)]),[],2)); %tuning (strength of response) per direction and per fly
        mtuningPdir=nanmean(tuningPdirPFly,2);
        stuningPdir=nanstd(tuningPdirPFly')/sqrt(length(tuningPdirPFly));
        rho=[mtuningPdir;mtuningPdir(1)]';
        PPp(i)=polar(theta,rho);%'Color',Color(i,:));
        PPp(i).Color=Color(i,:);
        hold on
        PolarPlot_std_new
    end
    [PPp.LineWidth] = deal(1.5);
    LegendName={'LayerA', 'LayerB', 'LayerC', 'LayerD'};
    legend(PPp,LegendName,'Location', 'southwest');
    
    
    set(gca,'YTickLabel','')
    set(gca,'FontSize',13)
    
    if av_Vector
        
        Z_Cond_OFF= averageDirectionVectors(Condition.T4T5_mb);
        
        Comp=compass(mean(Z_Cond_OFF.T5A.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0],'LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_OFF.T5B.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_OFF.T5C.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[1 0 0]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_OFF.T5D.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color', [.7 .7 0]','LineWidth',2);
        end
        
    end
    
    
    
    
    subplot(1,2,2)
    
    
    P = polar(theta, 4 * ones(size(theta)));
    set(P, 'Visible', 'off')
    hold on
    for i=1:4 %four Layers
        LAYERS=['A','B','C','D'];
        
        tuningPdirPFly=squeeze(max(eval(['Condition.ALL_Flies.T4',LAYERS(i)]),[],2)); %tuning (strength of response) per direction and per fly
        mtuningPdir=nanmean(tuningPdirPFly,2);
        stuningPdir=nanstd(tuningPdirPFly')/sqrt(length(tuningPdirPFly));
        rho=[mtuningPdir;mtuningPdir(1)]';
        PPp(i)=polar(theta,rho);
        PPp(i).Color=Color(i,:);
        set(PPp(i),'LineWidth', 1.5)
        
        % rlim([-5 15])
        %set(PPp,'FontSize',24)
        hold on
        
        PolarPlot_std_new
        
    end
    title('T4 responses to Bright Stripes');
    set(gca,'YTickLabel','')
    set(gca,'FontSize',13)
    set(F,'PaperSize', [8 8])
    
    
    
    if av_Vector
        
        Z_Cond_ON= averageDirectionVectors(Condition.T4T5_mb);
        
        Comp=compass(mean(Z_Cond_ON.T4A.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 .5 0],'LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_ON.T4B.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[0 0 1]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_ON.T4C.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color',[1 0 0]','LineWidth',2);
        end
        Comp=compass(mean(Z_Cond_ON.T4D.ALL)*5);
        for i=1:length(Comp)
            set(Comp(i),'color', [.7 .7 0]','LineWidth',2);
        end
        
    end
end

end