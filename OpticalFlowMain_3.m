%% preprocess images - imagej
% run normalize local contrast
% 40 40 3 stretched center
% run guassian blur 2.0 pixel


%% some constants for making figures

lightgreencolor = [102/255 102/255 102/255];
darkgreencolor = [153/255 153/255 153/255];
greycolor = [209/255 209/255 209/255];

%% The actual optical is done in this section, the rest of the code is just data visualization

clear all
close all

file = uigetfile('*.tif');
filename=['realoutput/' file];
opticFlow = opticalFlowFarneback('NeighborhoodSize', 5,'FilterSize',15);
%opticFlow = opticalFlowLK('NoiseThreshold',0.05);

info = imfinfo(file);
num_images = numel(info);
videocontainer = zeros(info(1).Height,info(1).Width,num_images);

for k = 1:num_images
    videocontainer(:,:,k) = imread(file, k);
end
videocontainer = videocontainer - min(videocontainer(:))+1;
videocontainer=videocontainer/max(videocontainer(:));
for k = 1:num_images

    frameGray = videocontainer(:,:,k);
    flowAll{k}= estimateFlow(opticFlow,frameGray); 
end

outputFileName=char(strcat([file '.mat']));
%save(outputFileName,'flowAll')
%save(outputFileName,'flowAll','videocontainer')
save(outputFileName,'-v7.3')    


%% let's makes some graphs for the circle data


load('circle1.tif.mat')
flowAll1 = flowAll;
load('circle2.tif.mat')
flowAll2 = flowAll;
load('circle3.tif.mat')
flowAll3 = flowAll;
%%
averagecircleu = (flowAll1{1}.Vx+ flowAll2{1}.Vx+flowAll3{1}.Vx);    
averagecirclev = (flowAll1{1}.Vy+ flowAll2{1}.Vy+flowAll3{1}.Vy);
for n=2:21
averagecircleu = averagecircleu+ (flowAll1{n}.Vx+ flowAll2{n}.Vx+flowAll3{n}.Vx);    
averagecirclev = averagecirclev+ (flowAll1{n}.Vy+ flowAll2{n}.Vy+flowAll3{n}.Vy);
end
averagecircleu = averagecircleu./3;
averagecirclev = averagecirclev./3;

figure
rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
hold on
rectangle('Position',[600,600,400,400],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor,'Curvature',[1 1])  % choose circle color here
scaling = 160;
blocksize = 160;
%scaling = 32;
%blocksize = 32;

x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);

uoffset = mean(mean(averagecircleu(600:1000,600:1000)));
voffset = mean(mean(averagecirclev(600:1000,600:1000)));
 u3= imresize(imgaussfilt(averagecircleu-uoffset,blocksize),1/scaling ,'bicubic'); % remove drift
 v3= imresize(imgaussfilt(averagecirclev-voffset,blocksize),1/scaling ,'bicubic');  %remove drift
 
q =quiver(x,y,u3.*8,v3.*8,0); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q =quiver(x,y,u3,v3,0); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q.Color = 'white';
q.Color = 'black';

% final adjustments made to the plot
set(gca,'Ydir','reverse')
axis equal
axis off
%title(['displacement vectors - all circles'])
%set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'color','none')

tempy = 1:1600;
tempy = repmat(tempy',[1 1600]);
tempx = 1:1600;
tempx = repmat(tempx,[1600 1]);
centerx = 800;
centery = 800;
distance = sqrt((tempy-centery).^2 + (tempx-centerx).^2)*211.2/1600; %distance in um
sumflow = sqrt( (averagecircleu-uoffset).^2+(averagecirclev-voffset).^2)*211.2/1600;
counter1 = 1;
binsize = 5;
lindistance = distance(:);
linsumflow = sumflow(:);
bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:145
    bindistance(counter1) = distancerange+binsize/2;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
figure
rectangle('Position',[0*211.2/1600 0 200*211.2/1600 3],'FaceColor','w','EdgeColor','k','LineStyle', '--')
hold on
fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
xlim([0  1130*211.2/1600])
%title('Distance from center vs estimated flow')
xlabel('Distance from center of circle')
ylabel('estimated optical flow in \mum (mean and SD)')


%% over time 1

h = figure;
filename = 'circleperframe-autoscaled.gif';

for n = 2:81
    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 .*80;
    averagecirclev = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 .*80;
    
    rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
    hold on
    rectangle('Position',[600,600,400,400],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor,'Curvature',[1 1])  % choose circle color here
    scaling = 80;
    blocksize = 80;
    %scaling = 32;
    %blocksize = 32;
    % mean(averagecircleu(:))
    % mean(averagecirclev(:))
    x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
    y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);

    uoffset = mean(mean(averagecircleu(600:1000,600:1000)));
    voffset = mean(mean(averagecirclev(600:1000,600:1000)));
     u3= imresize(imgaussfilt(averagecircleu-uoffset,blocksize),1/scaling ,'bicubic'); % remove drift
     v3= imresize(imgaussfilt(averagecirclev-voffset,blocksize),1/scaling ,'bicubic');  %remove drift

    q =quiver(x,y,u3,v3,1); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
    %q.Color = 'white';
    q.Color = 'black';

    % final adjustments made to the plot
    set(gca,'Ydir','reverse')
    axis equal
    axis off
    %title(['displacement vectors - all circles'])
    %set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'color','none')
    axis tight manual
  
    
    
    
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    
    
    
    hold off
    
    
    
    
end


%% over time 2
h = figure;
filename = 'circle2perframe.gif';

for n = 2:81
  
    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 ;
    averagecirclev = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 ;
    
    tempy = 1:1600;
    tempy = repmat(tempy',[1 1600]);
    tempx = 1:1600;
    tempx = repmat(tempx,[1600 1]);
    centerx = 800;
    centery = 800;
    distance = sqrt((tempy-centery).^2 + (tempx-centerx).^2)*211.2/1600; %distance in um
    sumflow = sqrt( (averagecircleu-mean(averagecircleu(:))).^2+(averagecirclev-mean(averagecirclev(:))).^2)*211.2/1600;
    counter1 = 1;
    binsize = 5;
    lindistance = distance(:);
    linsumflow = sumflow(:);
    bindistance = [];
    meanflowmag = [];
    stdflowmag = [];
    I= [];
    for distancerange =0:binsize:145
        bindistance(counter1) = distancerange+binsize/2;
        I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
        meanflowmag(counter1)= mean(linsumflow(I));
        stdflowmag(counter1)=std(linsumflow(I));
        counter1 = counter1+1;
    end
    rectangle('Position',[0*211.2/1600 0 1130*211.2/1600 0.3],'FaceColor','w','EdgeColor','w')
    rectangle('Position',[0*211.2/1600 0 200*211.2/1600 0.3],'FaceColor','w','EdgeColor','k','LineStyle', '--')
    hold on
    fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
    plot(bindistance,meanflowmag,'-k','LineWidth',2)
    xlim([0  1130*211.2/1600])
    %title('Distance from center vs estimated flow')
    xlabel('Distance from center of circle')
    ylabel('estimated optical flow in \mum (mean and SD)')
  
    
    
    
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    
    
    
    hold off
end

%%  100 rectangles 

%manually grab all the u and vs
%   allrect100u(:,:,3) = u;  
%   allrect100v(:,:,3) = v;
%  
%  save('allrect100uv','allrect100u','allrect100v')
%  save('allrect50uv','allrect100u','allrect100v')

%load('allrect100uv.mat')

averagecircleu = mean(allrect100u,3);
averagecirclev = mean(allrect100v,3);

scaling = 160;
blocksize = 160;
%scaling = 32;          % for sup
%blocksize = 32;        % for sup

figure
rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
hold on
rectangle('Position',[1,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
rectangle('Position',[1200,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);
u3= imresize(imgaussfilt(averagecircleu-mean(averagecircleu(:)),blocksize),1/scaling,'bicubic'); % remove drift
v3= imresize(imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize),1/scaling,'bicubic');  %remove drift

q =quiver(x,y,u3.*8,v3.*8,0);                   %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q =quiver(x,y,u3,v3,0);      %for sup          %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q.Color = 'white';
q.Color = 'black';

% final adjustments made to the plot
set(gca,'Ydir','reverse')
axis equal
axis off
%title(['displacement vectors - all rectangles 100 /mum apart'])
%set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'color','none')

distance = [1:1600] *211.2/1600;
counter1 = 1;
binsize = 5;
% lindistance = repmat(distance(:)', [1600 1]);
% linsumflow = (averagecircleu+10.5875)*211.2/1600;
%just the middle
lindistance = repmat(distance(:)', [400 1]);
linsumflow = (averagecircleu(600:1000,:)+10.5875)*211.2/1600;
bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:210
    bindistance(counter1) = distancerange+binsize/8;
    %bindistance(counter1) = distancerange;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    %quartileflowmag(counter1,:)= quantile(linsumflow(I),[0.05,0.25,0.5,0.75,0.95]);
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
figure
rectangle('Position',[0 -5 220 10],'FaceColor','w','EdgeColor','w')
hold on
rectangle('Position',[0*211.2/1600 -5 400*211.2/1600 10],'FaceColor','w','EdgeColor','k','LineStyle', '--')
rectangle('Position',[1200*211.2/1600 -5 400*211.2/1600 10],'FaceColor','w','EdgeColor','k','LineStyle', '--')
%plot(bindistance,quartileflowmag(:,3),'-r','LineWidth',2)
%fill([bindistance, fliplr(bindistance)], [quartileflowmag(:,2)', fliplr(quartileflowmag(:,4)')],'r','EdgeColor','none','FaceAlpha',0.5)
%fill([bindistance, fliplr(bindistance)], [quartileflowmag(:,1)', fliplr(quartileflowmag(:,5)')],'r','EdgeColor','none','FaceAlpha',0.2)
fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
plot([1,211],[0,0],'-k')

%title('horizontal flow during crosslinking of rectangles 100\mum apart')
xlabel('distance along the horizontal direction (\mum)')
ylabel('estimated optical flow in the horizontal direction (\mum)')
xlim([0 211.2])


%% over time 


 load('rect100-1.tif.mat')
 flowAll1 = flowAll;
 load('rect100-2.tif.mat')
 flowAll2 = flowAll;
 load('rect100-3.tif.mat')
 flowAll3 = flowAll;

%%
h = figure;
filename = 'r100perframe-manualscaled.gif';

for n = 2:21
    rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
    hold on
    rectangle('Position',[1,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
    rectangle('Position',[1200,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
    scaling = 80;
    blocksize = 80;
    % mean(averagecircleu(:))
    % mean(averagecirclev(:))
    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 .*20;
    averagecirclev = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 .*20;
    x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
    y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);
u3= imresize(imgaussfilt(averagecircleu-mean(averagecircleu(:)),blocksize),1/scaling,'bicubic'); % remove drift
v3= imresize(imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize),1/scaling,'bicubic');  %remove drift
    q =quiver(x,y,u3,v3,0); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
    q.Color = 'black';
    % final adjustments made to the plot
    set(gca,'Ydir','reverse')
    %title(['displacement vectors - all circles'])
    %set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'color','none')
    axis tight manual equal off
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    hold off
end

%%

h = figure;
filename = 'r100perframe2.gif';

for n = 2:21

    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 ;
    %averagesinglerectv = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 ;

distance = [1:1600] *211.2/1600;
counter1 = 1;
binsize = 5;

lindistance = repmat(distance(:)', [400 1]);
linsumflow = (averagecircleu(600:1000,:)-mean(averagecircleu(:)))*211.2/1600;
bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:210
    bindistance(counter1) = distancerange+binsize/8;
    %bindistance(counter1) = distancerange;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    %quartileflowmag(counter1,:)= quantile(linsumflow(I),[0.05,0.25,0.5,0.75,0.95]);
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
rectangle('Position',[0 -0.6 220 1.2],'FaceColor','w','EdgeColor','w')
hold on
rectangle('Position',[0*211.2/1600 -0.6 400*211.2/1600 1.2],'FaceColor','w','EdgeColor','k','LineStyle', '--')
rectangle('Position',[1200*211.2/1600 -0.6 400*211.2/1600 1.2],'FaceColor','w','EdgeColor','k','LineStyle', '--')
fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
plot([1,211],[0,0],'-k')
xlabel('distance along the horizontal direction (\mum)')
ylabel('estimated optical flow in the horizontal direction (\mum)')
xlim([0 211.2])
    axis tight manual 
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    hold off
end





%% 50um rect
%load('allrect50uv.mat')
averagecircleu = mean(allrect100u,3);
averagecirclev = mean(allrect100v,3);

scaling = 160;
blocksize = 160;
scaling = 32;          % for sup
blocksize = 32;        % for sup
figure
rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
hold on
rectangle('Position',[200,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
rectangle('Position',[1000,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);
 u3= imresize(imgaussfilt(averagecircleu-mean(averagecircleu(:)),blocksize),1/scaling,'bicubic'); % remove drift
 v3= imresize(imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize),1/scaling,'bicubic');  %remove drift
% imagesc(sqrt(imgaussfilt(averagecircleu-mean(averagecircleu(:)),blocksize).^2+imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize).^2)*211.2/1600)  %remove drift
% colorbar
% caxis([0 5])
% colormap jet

%q =quiver(x,y,u3.*8,v3.*8,0);                   %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
q =quiver(x,y,u3,v3,0);      %for sup          %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q.Color = 'white';
q.Color = 'black';

% final adjustments made to the plot
set(gca,'Ydir','reverse')
axis equal
axis off
%title(['displacement vectors - all rectangles 100 /mum apart'])
%title(['displacement vectors - all rectangles 50 /mum apart'])
%set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'color','none')




distance = [1:1600] *211.2/1600;
counter1 = 1;
binsize = 5;
%lindistance = repmat(distance(:)', [1600 1]);
%linsumflow = (averagecircleu+10.5875)*211.2/1600;
%just the middle
lindistance = repmat(distance(:)', [400 1]);
linsumflow = (averagecircleu(601:1000,:)+10.5875)*211.2/1600;


bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:210
    bindistance(counter1) = distancerange+binsize/8;
    %bindistance(counter1) = distancerange;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    %quartileflowmag(counter1,:)= quantile(linsumflow(I),[0.05,0.25,0.5,0.75,0.95]);
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
figure
rectangle('Position',[0 -5 220 10],'FaceColor','w','EdgeColor','w')
hold on
rectangle('Position',[200*211.2/1600 -5 400*211.2/1600 10],'FaceColor','w','EdgeColor','k','LineStyle', '--')
rectangle('Position',[1000*211.2/1600 -5 400*211.2/1600 10],'FaceColor','w','EdgeColor','k','LineStyle', '--')
% plot(bindistance,quartileflowmag(:,3),'-r','LineWidth',2)
% fill([bindistance, fliplr(bindistance)], [quartileflowmag(:,2)', fliplr(quartileflowmag(:,4)')],'r','EdgeColor','none','FaceAlpha',0.5)
% fill([bindistance, fliplr(bindistance)], [quartileflowmag(:,1)', fliplr(quartileflowmag(:,5)')],'r','EdgeColor','none','FaceAlpha',0.2)

fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
plot([1,211],[0,0],'-k')

%title('horizontal flow during crosslinking of rectangles 50\mum apart')
xlabel('distance along the horizontal direction (\mum)')
ylabel('estimated optical flow in the horizontal direction (\mum)')
xlim([0 211.2])


%% over time 


 load('rect50-1.tif.mat')
 flowAll1 = flowAll;
 load('rect50-2.tif.mat')
 flowAll2 = flowAll;
 load('rect50-3.tif.mat')
 flowAll3 = flowAll;

%%
h = figure;
filename = 'r50perframe-autoscaled.gif';

for n = 2:21
rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
hold on
rectangle('Position',[200,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
rectangle('Position',[1000,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
    scaling = 80;
    blocksize = 80;
    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 .*20;
    averagecirclev = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 .*20;
    
    x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
    y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);

u3= imresize(imgaussfilt(averagecircleu-mean(averagecircleu(:)),blocksize),1/scaling,'bicubic'); % remove drift
v3= imresize(imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize),1/scaling,'bicubic');  %remove drift

    q =quiver(x,y,u3,v3,1); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
    q.Color = 'black';

    % final adjustments made to the plot
    set(gca,'Ydir','reverse')
    %title(['displacement vectors - all circles'])
    %set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'color','none')
    axis tight manual equal off

     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    hold off

end
%%

h = figure;
filename = 'r50perframe2.gif';

for n = 2:21

    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 ;
    %averagesinglerectv = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 ;

distance = [1:1600] *211.2/1600;
counter1 = 1;
binsize = 5;

lindistance = repmat(distance(:)', [400 1]);
linsumflow = (averagecircleu(600:1000,:)-mean(averagecircleu(:)))*211.2/1600;
bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:210
    bindistance(counter1) = distancerange+binsize/8;
    %bindistance(counter1) = distancerange;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    %quartileflowmag(counter1,:)= quantile(linsumflow(I),[0.05,0.25,0.5,0.75,0.95]);
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
rectangle('Position',[0 -0.6 220 1.2],'FaceColor','w','EdgeColor','w')
hold on
rectangle('Position',[200*211.2/1600 -0.6 400*211.2/1600 1.2],'FaceColor','w','EdgeColor','k','LineStyle', '--')
rectangle('Position',[1000*211.2/1600 -0.6 400*211.2/1600 1.2],'FaceColor','w','EdgeColor','k','LineStyle', '--')
fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
plot([1,211],[0,0],'-k')
xlabel('distance along the horizontal direction (\mum)')
ylabel('estimated optical flow in the horizontal direction (\mum)')
xlim([0 211.2])
    axis tight manual 
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    hold off
end




%% single rectangles

%manually grab all the u and vs
% singlerectu(:,:,3) = u;  
% singlerectv(:,:,3) = v;
% 
% save('singlerectuv','singlerectu','singlerectv')

%load('singlerectuv.mat')

averagesinglerectu = mean(singlerectu,3);
averagesinglerectv = mean(singlerectv,3);

scaling = 160;
blocksize = 160;
%scaling = 32;          % for sup
%blocksize = 32;        % for sup
figure
rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
hold on
rectangle('Position',[1,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here
x = [1:scaling:size(averagesinglerectu,1)] + round(scaling/2); %set up the grid for where the arrows go
y = [1:scaling:size(averagesinglerectv,2)] + round(scaling/2);
 u3= imresize(imgaussfilt(averagesinglerectu,blocksize),1/scaling,'bicubic'); % remove drift
 v3= imresize(imgaussfilt(averagesinglerectv-mean(averagesinglerectv(:)),blocksize),1/scaling,'bicubic');  %remove drift
% imagesc(sqrt(imgaussfilt(averagecircleu-mean(averagecircleu(:)),blocksize).^2+imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize).^2)*211.2/1600)  %remove drift
% colorbar
% caxis([0 5])
% colormap jet
q =quiver(x,y,u3.*8,v3.*8,0);                   %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q =quiver(x,y,u3,v3,0);      %for sup          %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
%q.Color = 'white';
q.Color = 'black';
% final adjustments made to the plot
set(gca,'Ydir','reverse')
axis equal
axis off
%set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'color','none')


distance = [1:1600] *211.2/1600;
counter1 = 1;
binsize = 5;
%lindistance = repmat(distance(:)', [1600 1]);
%linsumflow = averagesinglerectu*211.2/1600;
%just the middle
lindistance = repmat(distance(:)', [400 1]);
linsumflow = averagesinglerectu(600:1000,:)*211.2/1600;
bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:210
    bindistance(counter1) = distancerange+binsize/8;
    %bindistance(counter1) = distancerange;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    %quartileflowmag(counter1,:)= quantile(linsumflow(I),[0.05,0.25,0.5,0.75,0.95]);
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
figure
rectangle('Position',[0 -5 220 10],'FaceColor','w','EdgeColor','w')
hold on
rectangle('Position',[1*211.2/1600 -5 400*211.2/1600 10],'FaceColor','w','EdgeColor','k','LineStyle', '--')
% plot(bindistance,quartileflowmag(:,3),'-r','LineWidth',2)
% fill([bindistance, fliplr(bindistance)], [quartileflowmag(:,2)', fliplr(quartileflowmag(:,4)')],'r','EdgeColor','none','FaceAlpha',0.5)
% fill([bindistance, fliplr(bindistance)], [quartileflowmag(:,1)', fliplr(quartileflowmag(:,5)')],'r','EdgeColor','none','FaceAlpha',0.2)
% 
fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
plot([1,211],[0,0],'-k')

%title('horizontal flow during crosslinking of single rectangles')
xlabel('distance along the horizontal direction (\mum)')
ylabel('estimated optical flow in the horizontal direction (\mum)')
xlim([0 211.2])

%% over time 


 load('rect1.tif.mat')
 flowAll1 = flowAll;
 load('rect2.tif.mat')
 flowAll2 = flowAll;
 load('rect3.tif.mat')
 flowAll3 = flowAll;

%%
h = figure;
filename = 'rperframe-manualscaled.gif';

for n = 2:21
rectangle('Position',[1,1,1600,1600],'FaceColor',darkgreencolor,'EdgeColor',darkgreencolor,'LineWidth',3) % choose background color here
hold on
rectangle('Position',[1,1,400,1600],'FaceColor',lightgreencolor,'EdgeColor',lightgreencolor)  % choose circle color here

    scaling = 80;
    blocksize = 80;
    % mean(averagecircleu(:))
    % mean(averagecirclev(:))
    
    averagecircleu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 .*20;
    averagecirclev = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 .*20;
    
    x = [1:scaling:size(averagecircleu,1)] + round(scaling/2); %set up the grid for where the arrows go
    y = [1:scaling:size(averagecirclev,2)] + round(scaling/2);

    u3= imresize(imgaussfilt(averagecircleu,blocksize),1/scaling,'bicubic'); % remove drift
    v3= imresize(imgaussfilt(averagecirclev-mean(averagecirclev(:)),blocksize),1/scaling,'bicubic');  %remove drift
    
    q =quiver(x,y,u3,v3,0); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
    q.Color = 'black';

    % final adjustments made to the plot
    set(gca,'Ydir','reverse')
    %title(['displacement vectors - all circles'])
    %set(gcf, 'Position', get(0, 'Screensize'));
    set(gca,'color','none')
    axis tight manual equal off
  
    
    
    
    
    
    
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    
    
    
    hold off
    
    
    
    
end



%%

h = figure;
filename = 'rperframe2.gif';

for n = 2:21

    averagesinglerectu = (flowAll1{n}.Vx + flowAll2{n}.Vx + flowAll3{n}.Vx) ./3 ;
    averagesinglerectv = (flowAll1{n}.Vy + flowAll2{n}.Vy + flowAll3{n}.Vy) ./3 ;

distance = [1:1600] *211.2/1600;
counter1 = 1;
binsize = 5;
%lindistance = repmat(distance(:)', [1600 1]);
%linsumflow = averagesinglerectu*211.2/1600;
%just the middle
lindistance = repmat(distance(:)', [400 1]);
linsumflow = averagesinglerectu(600:1000,:)*211.2/1600;
bindistance = [];
meanflowmag = [];
stdflowmag = [];
I= [];
for distancerange =0:binsize:210
    bindistance(counter1) = distancerange+binsize/8;
    %bindistance(counter1) = distancerange;
    I = lindistance<(distancerange+binsize/2) & lindistance>distancerange;
    %quartileflowmag(counter1,:)= quantile(linsumflow(I),[0.05,0.25,0.5,0.75,0.95]);
    meanflowmag(counter1)= mean(linsumflow(I));
    stdflowmag(counter1)=std(linsumflow(I));
    counter1 = counter1+1;
end
rectangle('Position',[0 -0.5 220 1],'FaceColor','w','EdgeColor','w')
hold on
rectangle('Position',[1*211.2/1600 -0.5 400*211.2/1600 1],'FaceColor','w','EdgeColor','k','LineStyle', '--')
fill([bindistance, fliplr(bindistance)], [meanflowmag-stdflowmag, fliplr(meanflowmag+stdflowmag)],greycolor,'EdgeColor','none')
plot(bindistance,meanflowmag,'-k','LineWidth',2)
plot([1,211],[0,0],'-k')
%title('horizontal flow during crosslinking of single rectangles')
xlabel('distance along the horizontal direction (\mum)')
ylabel('estimated optical flow in the horizontal direction (\mum)')
xlim([0 211.2])


    axis tight manual 
  
    
    
    
    
    
    
     drawnow 
     % Capture the plot as an image 
     frame = getframe(h); 
     im = frame2im(frame); 
     [imind,cm] = rgb2ind(im,256); 
     % Write to the GIF File 
     if n == 2 
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
     else 
         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
     end 
    
    
    
    hold off
    
    
    
end

