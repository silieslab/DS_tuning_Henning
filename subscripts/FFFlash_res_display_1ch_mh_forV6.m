function [h2,h3,h4,out] = FFFlash_res_display_1ch_mh_forV6(in,SelectLayers)
%% Find stimulus timing - modify after clean-up script is written

out = in;
nframes = in.xml.frames;
fps = in.xml.framerate;

thresh = 0.5;
mask = out.fstimval>thresh;
dmask = mask(2:end)-mask(1:end-1);
dmask = [0; dmask];

if isfield(in, 'ch1a_crop')
    
    ch1a=in.ch1a_crop;
end
%% Plot ratio signals from ROIs together with stimulus timing
try
    nMasks = length(in.masks_CA);
    T='Layers calculated';
catch
    if SelectLayers
        
        in2=in.ClusterInfo_ManuallySelect;
        if isfield(in2, 'masks_CA_m')
            Masks=in2.masks_CA_m;
            nMasks = length(Masks);
            avSignal1_CA=in2.avSignal1_CA_m;
            dSignal1_CA=in2.dSignal1_CA_m;
            T='Layers selected manually - merged ROIs';
        else
            Masks=in.ClusterInfo_ManuallySelect.masks_CA;
            nMasks = length(Masks);
            T='Layers selected manually';
            avSignal1_CA=in2.avSignal1_CA;
            dSignal1_CA=in2.dSignal1_CA;
            
        end
        
    else
        
        in2=in.ClusterInfo;
        
        if isfield(in2, 'masks_CA_m')
            Masks=in2.masks_CA_m;
            nMasks = length(Masks);
            avSignal1_CA=in2.avSignal1_CA_m;
            dSignal1_CA=in2.dSignal1_CA_m;
            T='Layers calculated - merged ROIs';
        else
            Masks=in2.masks_CA;
            nMasks = length(Masks);
            T='Layers calculated';
            avSignal1_CA=in2.avSignal1_CA;
            dSignal1_CA=in2.dSignal1_CA;
            
        end
    end
end
%cm = [colormap('lines');colormap('lines');colormap('lines')];
L1ON=[59,235,53]/255;
L1OFF=[37,134,33]/255;
L2ON=[61,94,255]/255;
L2OFF=[39,59,159]/255;
L3ON=[250,68,68]/255;
L3OFF=[156,37,37]/255;
L4ON=[249,236,62]/255;
L4OFF=[191,181,50]/255;
cm=[];

if isfield(in2, 'Layer_m')
    Layer=in2.Layer_m;
    T4_T5=in2.T4_T5_m;
else
    Layer=in2.Layer;
    T4_T5=in2.T4_T5;
end


for i=1:length(Layer)
    if     Layer(i) == 1 && T4_T5(i)==4
        cm=[cm;L1ON];
    elseif Layer(i) == 1 && T4_T5(i)==5
        cm=[cm;L1OFF];
    elseif Layer(i)== 2 && T4_T5(i)==4
        cm=[cm;L2ON];
    elseif Layer(i)== 2 && T4_T5(i)==5
        cm=[cm;L2OFF];
    elseif Layer(i)== 3 && T4_T5(i)==4
        cm=[cm;L3ON];
    elseif Layer(i)== 3 && T4_T5(i)==5
        cm=[cm;L3OFF];
    elseif Layer(i)== 4 && T4_T5(i)==4
        cm=[cm;L4ON];
    elseif Layer(i)== 4 && T4_T5(i)==5
        cm=[cm;L4OFF];
    end
end



h1=figure;
for i = 1:nMasks
    plot((1:nframes)/fps,(avSignal1_CA(i,:)/mean(avSignal1_CA(i,:))+2*(i-1)),'color',cm(i,:));hold on;
end
xlabel('time (sec)');
title('Signal in region of interest - before background substraction, aligned data');
hold off


h2=figure;
for i = 1:nMasks
    plot((1:nframes)/fps,(dSignal1_CA(i,:)/mean(dSignal1_CA(i,:)))+4*(i-1),'color',cm(i,:),'linewidth',2);hold on;
end

plot((1:nframes)/fps,out.fstimval*3,'LineWidth',2);
axis([0 nframes/fps 0 i*4+1]);
xlabel('time (sec)');


inds = find(dmask~=0);
for k = 1:length(inds)
    if(dmask(inds(k))>0)
        line([inds(k)/fps inds(k)/fps],[0 nMasks+1],'color',cm(i,:),'LineWidth',1,'LineStyle','-');
    else
        line([inds(k)/fps inds(k)/fps],[0 nMasks+1],'color',cm(i,:),'LineWidth',1,'LineStyle','--');
    end
end
xlabel('time (sec)');
title('Ratio, YFP, CFP signals - aligned data');
hold off



% Create a colored map of ROIs
if exist('in.xml.linesperframe','var') && exist('in.xml.pixperline','var')
    CMask = zeros(in.xml.linesperframe, in.xml.pixperline, 3); % Luis 13.11.2015
else
    CMask = zeros(size(Masks{1,1},1), size(Masks{1,1},2), 3);% Miri 09.11.2018
    %adapted to the size of the cropped image!!
end

for i = 1:nMasks
    curColor = cm(i,:);
    %curMask = cat(3,curColor(1).*flipud(in.masks{i}),curColor(2).*flipud(in.masks{i}),curColor(3).*flipud(in.masks{i}));
    curMask = cat(3,curColor(1).*Masks{i},curColor(2).*Masks{i},curColor(3).*Masks{i});
    
    CMask = CMask + curMask;
end

if(~isfield(in,'AV1'))
    AV = squeeze(sum(ch1a,3))/nframes; % The average image
    AV = im2double(AV);
    AV = AV./max(AV(:));
else
    AV = in.AV1;
end

h3=figure;
imshow(AV,[],'InitialMagnification',600);
hold(h3.Children,'on');
hcllus=imshow(CMask,'InitialMagnification',600, 'Parent', h3.Children);
set(hcllus,'AlphaData',0.5);
h3.Children.Title.String=[sprintf('average aligned image at z-depth %0.5g microns',in.xml.zdepth),T];
hold(h3.Children,'off')
hold off 
% again with clusters

if exist('in.xml.linesperframe') && exist('in.xml.pixperline')
    CMask = zeros(in.xml.linesperframe, in.xml.pixperline, 3); % Luis 13.11.2015
else
    CMask = zeros(size(Masks{1,1},1), size(Masks{1,1},2), 3);% Miri 09.11.2018
    %adapted to the size of the cropped image!!
end


Fi=figure;
cm2=colormap('ColorCube');
close(Fi)
cm2=[cm2; cm2; cm2];

for i = 1:nMasks
    curColor = cm2(i,:);
    %curMask = cat(3,curColor(1).*flipud(in.masks{i}),curColor(2).*flipud(in.masks{i}),curColor(3).*flipud(in.masks{i}));
    curMask = cat(3,curColor(1).*Masks{i},curColor(2).*Masks{i},curColor(3).*Masks{i});
    
    CMask = CMask + curMask;
end


h4=figure;
h41=imshow(AV,[],'InitialMagnification',600);
hold(h4.Children,'on');h = imshow(CMask,'InitialMagnification',600, 'Parent', h4.Children);
set(h,'AlphaData',0.5);
h4.Children.Title.String=[sprintf('average aligned image at z-depth %0.5g microns',in.xml.zdepth),T];
hold(h4.Children,'off')

end
% %% Z-stack analysis
%
% if 0
% curDir = pwd;
% cd(in.fileloc);
% [fname,pathname]=uigetfile('*','Choose the corresponding z-stack...');
% cd(curDir);
% if fname ~= 0
%     z=read_z_stack(pathname);
% else
%     return;
% end
%
% %here
% zd=z.xml.zdepth;
% [dum,choose]=min(abs(zd-in.xml.zdepth));
%
% figure; hold on;
% imagesc(squeeze(z.ims(:,:,choose))); colormap('gray');
% plot(20+[0 50/z.xml.xres],(size(z.ims,1)-25)*[1 1],'color',ones(1,3),'linewidth',3);
% title(['Z = ' num2str(zd(choose)) ' microns, ' num2str(zd(choose)-zd(1)) ' microns from the top']);
% set(gca,'data',[1 1 1],'xtick',[],'ytick',[],'xlim',[0 size(z.ims,2)],'ylim',[0 size(z.ims,1)],'ydir','normal');
%
% [xin,yin]=ginput;
% y1=round(min(yin)); y2=round(max(yin));
% x1=round(min(xin)); x2=round(max(xin));
% maxproj=squeeze(max(z.ims(y1:y2,x1:x2,:),[],1))';
% figure; hold on;
% imagesc(maxproj); colormap('gray');
% plot([1 size(maxproj,2)],choose*[1 1],'color',ones(1,3),'linewidth',2);
% plot(20+[0 50/z.xml.xres],(size(maxproj,1)-5)*[1 1],'color',ones(1,3),'linewidth',3);
% set(gca,'xtick',[],'ytick',[],'xlim',[1 size(maxproj,2)],'ylim',[1 size(maxproj,1)],'ydir','reverse');
% set(gca,'data',[abs(zd(2)-zd(1)),z.xml.xres, 1]);
%
% end
