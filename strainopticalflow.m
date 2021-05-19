clear
close all

load('allcircleuv.mat')
% load('allrect100uv.mat')

 
%%
clear estrain wrotation Vs Ds
scale = 1;
% Ux= mean(allcircleu ,3) -  mean(mean(mean(allcircleu(600:1000,600:1000,:))));
% Uy= mean(allcirclev ,3) -  mean(mean(mean(allcirclev(600:1000,600:1000,:))));

Ux= allcircleu(:,:,1) -  mean(mean(mean(allcircleu(600:1000,600:1000,1))));
Uy= allcirclev(:,:,1) -  mean(mean(mean(allcirclev(600:1000,600:1000,1))));

% Ux= mean(allrect100u ,3) - mean(allrect100u(:));
% Uy= mean(allrect100v ,3) - mean(allrect100v(:));

%Ux = imresize( imgaussfilt(  Ux, 1/scale), scale,'bicubic');
%Uy = imresize( imgaussfilt(  Uy, 1/scale), scale,'bicubic');



h = 1/scale;

figure
x = [1:h:size(allcircleu,1)] + round(h/2); %set up the grid for where the arrows go
y = [1:h:size(allcircleu,2)] + round(h/2);
% x = [1:h:size(allrect100u,1)] + round(h/2); %set up the grid for where the arrows go
% y = [1:h:size(allrect100u,2)] + round(h/2);
q =quiver(x,y,Ux,Uy,1); %draw arrows, arrows are exact length and not scaled by any factors if scale=0[var5]. otherwise they are scaled.
q.Color = 'black';
set(gca,'Ydir','reverse')
axis equal
axis off
set(gca,'color','none')









[Uxx,Uxy] = gradient(Ux,h,h);
[Uyx,Uyy] = gradient(Uy,h,h);


for i=1:size(Ux,1)
    for j=1:size(Ux,2)
        Fdisplacementgradient =  [Uxx(i,j) Uxy(i,j); Uyx(i,j) Uyy(i,j)] ;
        [V,D] = eig((Fdisplacementgradient + Fdisplacementgradient')./2);
        estrain(i,j,:,:) = (Fdisplacementgradient + Fdisplacementgradient')./2 ;
        wrotation(i,j,:,:)  =(Fdisplacementgradient - Fdisplacementgradient')./2  ; 
        Vs(i,j,:,:) = V;
        Ds(i,j,:,:) = D;
    end
end


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
alignmentdistributions = [];
for i=1:size(Vs,1)
    for j=1:size(Vs,2)
        alignmentdistribution = [0.1:0.1:360];
        alignmentdistribution = atand(   (sind(alignmentdistribution)*(1+Ds(i,j,1,1)) )  ./   (cosd(alignmentdistribution)*(1+Ds(i,j,2,2)) ));
        %alignmentdistribution = atand(   (sind(alignmentdistribution)*(1-0.1 ))  ./   (cosd(alignmentdistribution)*(1+0 ) ) );
        alignmentdistribution = alignmentdistribution - atand(Vs(i,j,2,1)/Vs(i,j,1,1));
        alignmentdistribution(alignmentdistribution<-90) = alignmentdistribution(alignmentdistribution<-90) + 180;
        alignmentdistribution(alignmentdistribution>90) = alignmentdistribution(alignmentdistribution>90) - 180;
        alignmentdistributions(i,j,:) = alignmentdistribution;
    end
end

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


temp1 = alignmentdistributions([3:9,25:31,],5:25,:);
figure
subplot(2,2,1)
histogram(temp1(:) )
temp1 = alignmentdistributions([11:21], 5:25,:);
subplot(2,2,2)
histogram( temp1(:))