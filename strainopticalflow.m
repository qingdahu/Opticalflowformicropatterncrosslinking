clear
close all

load('allcircleuv.mat')
%load('allrect100uv.mat')

 
%%
clear estrain wrotation Vs Ds
scale = 0.02;
Ux= mean(allcircleu ,3) - mean(allcircleu(:));
Uy= mean(allcirclev ,3) - mean(allcirclev(:));
% Ux= mean(allrect100u ,3) - mean(allrect100u(:));
% Uy= mean(allrect100v ,3) - mean(allrect100v(:));

Ux = imresize( imgaussfilt(  Ux, 1/scale), scale,'bicubic');
Uy = imresize( imgaussfilt(  Uy, 1/scale), scale,'bicubic');

h = 1/scale;




[Uxy,Uxx] = gradient(Ux,h,h);
[Uyy,Uyx] = gradient(Uy,h,h);


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
