clear
close all

% updated 2021/5/20 
% Written by Qingda Hu
% Based on the formulas from
% https://www.files.ethz.ch/structuralgeology/sms/NumModRocDef/Strain_Tensors.pdf
% and https://ocw.mit.edu/courses/mechanical-engineering/2-080j-structural-mechanics-fall-2013/course-notes/MIT2_080JF13_Lecture2.pdf

 
%% preprocessing of data
% we want to load the the Ux and Uy (displacement, in units of pixels)
% we also want to keep track of the dx (in units of pixels)

%if you want to just load one: choose from all the .mat files that have
%flowAll objects in them
load('circle3.tif.mat')     %load .mat with flowAll object

Ux = flowAll{2}.Vx;         % load Ux between first and second frame 
Uy = flowAll{2}.Vy;         % load Uy between first and second frame 
for n=3:21                  % loop to 21 if you want 500 scans, 81 if you want 2000 scans (where it exists)
Ux = Ux+ flowAll{n}.Vx;     % collect all the Ux 
Uy = Uy+ flowAll{n}.Vy;
end


% % we may want a more complicated processing process involving some
% % complicated filtering to account for the poresize of the fiber network 
% Ux = zeros(size(flowAll{2}.Vx));         
% Uy = zeros(size(flowAll{2}.Vx));         
% for n=2:21                  % loop to 21 if you want 500 scans, 81 if you want 2000 scans (where it exists)
% Ux = Ux+ ordfilt2( flowAll{n}.Vx, 6, true(3));     % collect all the Ux 
% Uy = Uy+ ordfilt2( flowAll{n}.Vy, 6, true(3));
% end



%if you want to remove drift
Ux= Ux - mean(Ux(:)); % global drift removal in x
Uy= Uy - mean(Uy(:)); % global drift removal in y
% you may want to consider local drift remove vs global

%you may want to smooth or resize the images
scale = 0.05;
Ux = imresize( imgaussfilt(  Ux, 5/scale), scale,'bicubic');
Uy = imresize( imgaussfilt(  Uy, 5/scale), scale,'bicubic');

% if the image has been scaled down, we need to account for that. If not, h needs to be 1 or else code will error 
try
h = 1/scale;
catch 
    h=1;
end


%if you want to visualize Ux and Uy at this stage before moving forward
figure
x = ((1:size(Ux,1)) + round(1/2))*h; %set up the grid for where the arrows go
y = ((1:size(Ux,2)) + round(1/2))*h;
q =quiver(x,y,Ux,Uy,0); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
q.Color = 'black';
set(gca,'Ydir','reverse')
axis equal
axis off
set(gca,'color','none')
figure
imagesc( sqrt(Ux.^2 + Uy.^2) /7.5758)
colorbar


%% Calculating displacement gradient tensor
clear estrain wrotation Vs Ds  % clearning these variable helps to avoid errors due to changes in scaling between runs

[Uxx,Uxy] = gradient(Ux,h,h); % and we are done!
[Uyx,Uyy] = gradient(Uy,h,h);

% We may want to rearrange H so that we can call the nummbers easier
for i=1:size(Ux,1)
    for j=1:size(Ux,2)
        Hdisplacementgradient(i,j,:,:) =  [Uxx(i,j) Uxy(i,j); Uyx(i,j) Uyy(i,j)];
    end
end

% if we want engineering strain for small deformations
for i=1:size(Ux,1)
    for j=1:size(Ux,2)
        estrain(i,j,:,:) = (squeeze(Hdisplacementgradient(i,j,:,:)) + squeeze(Hdisplacementgradient(i,j,:,:))')./2 ;      % strain
        wrotation(i,j,:,:)  =(squeeze(Hdisplacementgradient(i,j,:,:)) - squeeze(Hdisplacementgradient(i,j,:,:))')./2  ;   % rotation
        [V,D] = eig(squeeze(estrain(i,j,:,:)));                                      % V are the eigenvectors and D are the eigenvalues
        Vs(i,j,:,:) = V;                                                    % I think its [x1,x2;y1,y2] (i think)
        Ds(i,j,:,:) = D;
    end
end        

%if you want Cauchy-Green
for i=1:size(Ux,1)
    for j=1:size(Ux,2)
        F = squeeze(Hdisplacementgradient(i,j,:,:)) + eye(2);
        estrain(i,j,:,:) = F'*F ;      % Cauchy-Green strain tensor
    end
end      


%% We actually need to go about finding the change in anisotropy due to the displacement (with affine displacement assumption)
% can we directly find the change in anisotropy without going to strain and
% is that faster/better or not?
% can we just insert text vectors and ask how they are transformed based on
% matrix F?
for i=1:size(Ux,1)
    for j=1:size(Ux,2)
        F = Hdisplacementgradient(i,j,:,:) + eye(2);
        alignmentdistribution = [0.1:0.1:360];
        testvectors = [sind(alignmentdistribution) ; cosd(alignmentdistribution)];
    end
end    



%Using engineering strain (code is probably correct but need more testing)
alignmentdistributions = [];
for i=1:size(Vs,1)
    for j=1:size(Vs,2)
        alignmentdistribution = [0.1:0.1:360];
        alignmentdistribution = atand(   (sind(alignmentdistribution)*(1+Ds(i,j,1,1)) )  ./   (cosd(alignmentdistribution)*(1+Ds(i,j,2,2)) )); 
        alignmentdistribution = alignmentdistribution + atand(Vs(i,j,1,1)/Vs(i,j,2,1));
        alignmentdistribution(alignmentdistribution<-90) = alignmentdistribution(alignmentdistribution<-90) + 180;
        alignmentdistribution(alignmentdistribution>90) = alignmentdistribution(alignmentdistribution>90) - 180;
        alignmentdistributions(i,j,:) = alignmentdistribution;
    end
end














%% old code from here down

figure
%quiver(squeeze(estrain(:,:,1,1)),squeeze((estrain(:,:,2,2)))) ;
subplot(2,2,1)
imagesc(estrain(:,:,1,1))
title('e11')
axis equal tight
colorbar
subplot(2,2,2)
imagesc(estrain(:,:,1,2))
title('e12')
axis equal tight
colorbar
subplot(2,2,3)
imagesc(estrain(:,:,2,1))
title('e21')
axis equal tight
colorbar
subplot(2,2,4)
imagesc(estrain(:,:,2,2))
title('e22')
axis equal tight
colorbar


figure
subplot(2,2,1)
imagesc(wrotation(:,:,1,1))
title('w11')
axis equal tight
colorbar
subplot(2,2,2)
imagesc(wrotation(:,:,1,2))
title('w12')
axis equal tight
colorbar
subplot(2,2,3)
imagesc(wrotation(:,:,2,1))
title('w21')
axis equal tight
colorbar
subplot(2,2,4)
imagesc(wrotation(:,:,2,2))
title('w22')
axis equal tight
colorbar

figure 
quiver(Vs(:,:,1,1).*Ds(:,:,1,1),Vs(:,:,2,1).*Ds(:,:,1,1),1)
hold on
quiver(Vs(:,:,1,2).*Ds(:,:,2,2),Vs(:,:,2,2).*Ds(:,:,2,2),1)
hold off
axis equal tight




figure
subplot(2,2,1)
imagesc(Ds(:,:,1,1))
title('\lambda1')
axis equal tight
colorbar
subplot(2,2,2)
imagesc(Ds(:,:,2,2))
title('\lambda2')
axis equal tight
colorbar



% 
% figure 
% quiver(Vs(3:end-2,3:end-2,1,2).*Ds(3:end-2,3:end-2,1,1)+ Vs(3:end-2,3:end-2,1,2).*Ds(3:end-2,3:end-2,2,2)   , Vs(3:end-2,3:end-2,2,2).*Ds(3:end-2,3:end-2,1,1)+ Vs(3:end-2,3:end-2,2,2).*Ds(3:end-2,3:end-2,2,2)  ,1)
%%


%% first attempt at visualization

figure
histogram(alignmentdistributions(:))

temp1 = alignmentdistributions(25:32,2:10,:);
figure
subplot(2,2,1)
histogram(temp1(:) )
temp1 = alignmentdistributions(12:17, 12:17,:);
subplot(2,2,2)
histogram(temp1(:) )
temp1 = alignmentdistributions(2:8, 2:8,:);
subplot(2,2,3)
histogram(temp1(:) )
temp1 = alignmentdistributions(20:25, 20:25,:);
subplot(2,2,4)
histogram(temp1(:) )


temp1 = alignmentdistributions([3:9,25:31],5:25,:);
figure
subplot(2,2,1)
histogram(temp1(:) )
temp1 = alignmentdistributions(11:21, 5:25,:);
subplot(2,2,2)
histogram( temp1(:))