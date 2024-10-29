%% This code is modified from CLFM_fish_Reconstruction_CPU.m file
%Reference: Zhang, Z., Bai, L., Cong, L. et al. Imaging volumetric dynamics at high speed in mouse and zebrafish brain with confocal light field microscopy. Nat Biotechnol 39, 74–83 (2021).

%% Volumetric reconstruction for freely-swimming zebrafish
% Scanning wide-field tomography for high-speed, mesoscale-volumetric imaging of biodynamics in vivo
% SYSTEM REQUIREMENTS
% Memory: 128 GB RAM
% Graphics: 24 GB RAM
% MATLAB: R2021a

%% This code involving several functional modules: 
%1. 3D PSFs simulation for random access
%2. Pixel-wise alignment
%3. 3D reconstruction by RL deconvolution

close all
clear

%% set workpath to this code
local_address=mfilename('fullpath');
[pathstr,namestr]=fileparts(local_address);
cd(pathstr);
addpath(pathstr);

%% Load wide-field PSFs for calculating 3D PSFs in random access VoIs

load(['psf_wide_field.mat']); % load wide field PSF
group1_info=[1:4]';
group2_info=[1:4]';
load('delta.mat')
structured_mode=1;
  
psfWAVE_STACK=gpuArray(psfWAVE_STACK);

%% Load raw images and Segment image now.

load(['data.mat']);
load('psf_part.mat', 'fitresult_V_vs_move1')
load('psf_part.mat', 'fitresult_V_vs_move2')
clear move
move1(1,:)=fitresult_V_vs_move1(V(2,:),V(1,:));
move1(2,:)=fitresult_V_vs_move2(V(2,:),V(1,:));
move1=move1';
move1(:,1)=move1(:,1)-0.08;
move1(:,2)=move1(:,2)+0.087;

load(['seg_parameter_raw_data.mat']);
center=[1317,1085];
images=loadtiff('demo_data\demo_data_1second.tif');


range=1:36;
frame_num=size(range,2);
move1=move1(range,:);

[views1,views2]=image_seg_free_zebrafish ( V(:,range),images(:,:,range),center,fit_group1,fit_group2,frame_num);% image segmentations with PWA

stack_seg_wide_group1=views1(:,:,:,:);

clear views1 views2 

times=size(stack_seg_wide_group1,4);
disp(['total fov frame is ' num2str(times)]);
image_seg=stack_seg_wide_group1(:,:,:,1);

ObjRECON_t_group1=zeros(size(image_seg,1),size(image_seg,2),81,times,'single'); 


for time=1:times
    
    
savepath=[ pathstr '\demo_reconstruction\']; 

 if exist(savepath)==7,
   ; 
else
   mkdir(savepath);   
 end
 
 
disp(['This is time : ' num2str(time)]);
tic;

if time<1.1     %calculate 3D PSFs
for ii=1:size(psfWAVE_STACK,3)
      
       depthi=psfWAVE_STACK(:,:,ii);
       depthi=padarray(depthi,[3169,3169]);
       frequency=fftshift(fft2(depthi));
       
   for views=1:4
      
       delta1=squeeze(delta_mm_g1(views,:,:));
       center=mean(delta1,2)*6.5;
       m_center=move1(time,:)*1000;
        center=center+m_center';
       center=center/(150/33)+[3302,3302]';
       center=round(center);
       lens=zeros(6603,6603);
       
       if views>1.5
       for xcen=center(1)-620:center(1)+620
           for ycen=center(2)-620:center(2)+620
               if (xcen-center(1)).^2+(ycen-center(2)).^2<380000
                   lens(xcen,ycen)=1;
               end
           end
       end
       else
       for xcen=center(1)-440:center(1)+440
           for ycen=center(2)-440:center(2)+440
               if (xcen-center(1)).^2+(ycen-center(2)).^2<193600
                   lens(xcen,ycen)=1;
               end
           end
       end
       end
       psf_domain=ifft2(lens.*frequency);
       psf_l=imresize(abs(psf_domain(3151:3453,3154:3450)).^2,[101,99],'cubic');
       psf_lens(:,:,ii,views)=psf_l';
   end
   
end

psf_lens_write=psf_lens./(1e24);

end

if time>1  
    if norm(V(:,range(1)+time-1)-V(:,range(1)-1+time-1))>0.1
        clear psf_lens  psf_domain  psf_l delta2 center
     for ii=1:size(psfWAVE_STACK,3)
      
       depthi=psfWAVE_STACK(:,:,ii);
       depthi=padarray(depthi,[3169,3169]);
       frequency=fftshift(fft2(depthi));
       
   for views=1:4
       delta1=squeeze(delta_mm_g1(views,:,:));
       center=mean(delta1,2)*6.5+move1(time,:)*1000;
       center=center/(150/33)+[3302,3302];
       center=round(center);
       lens=zeros(6603,6603);
       
       if views>1.5
       for xcen=center(1)-620:center(1)+620
           for ycen=center(2)-620:center(2)+620
               if (xcen-center(1)).^2+(ycen-center(2)).^2<380000
                   lens(xcen,ycen)=1;
               end
           end
       end
       else
       for xcen=center(1)-440:center(1)+440
           for ycen=center(2)-440:center(2)+440
               if (xcen-center(1)).^2+(ycen-center(2)).^2<193600
                   lens(xcen,ycen)=1;
               end
           end
       end
       end
       psf_domain=ifft2(lens.*frequency);
       psf_l=imresize(abs(psf_domain(3151:3453,3154:3450)).^2,[101,99],'cubic');
       psf_lens(:,:,ii,views)=psf_l';
    end
       
     end
psf_lens_write=psf_lens./(1e24);
  
   
    else
        psf_lens=psf_lens;
    end
end
psf_lens_group1=gather(psf_lens);

toc;

image_seg1=stack_seg_wide_group1(:,:,:,time);
%image_seg2=stack_seg_wide_group2(:,:,:,time);
Img_group1=image_seg1;
%Img_group2=image_seg2;
Img_group1=single(Img_group1);
%Img_group2=single(Img_group2);

ximg=size(Img_group1,1);
yimg=size(Img_group1,2);

factor=1;

Img_group1= imresize(Img_group1,[factor*ximg factor*yimg],'cubic');
%Img_group2= imresize(Img_group2,[factor*ximg factor*yimg],'cubic');



 Info1=[1:4]';
Info2=[1:4]'; %%如果不想要某个透镜，可以删除它

%savename='recon' num2str(numberi) '.tif';
% reconstruction parameters
ReconROIx=size(image_seg,1)/2-0.5; 
ReconROIy=size(image_seg,2)/2-0.5;% half-width of object
FitROIx=size(image_seg,1)/2-0.5; 
FitROIy=size(image_seg,2)/2-0.5; 

% half-width of calibration area (>=ReconROI)
ItN=15; % number of iterations
SNR=200; % estimated signal noise ratio
show_reconresult=1; 
PSFROIx=49;
PSFROIy=50;

psf_lensn_group1=zeros(factor*99,factor*101,size(psf_lens_group1,3),size(group1_info,1),'single');

for ii=1:size(Info1,1)
    
    for zz=1:size(psf_lens_group1,3)
         psf_lensn_group1(:,:,zz,ii)=imresize(psf_lens_group1(:,:,zz,Info1(ii)),[factor*99 factor*101],'cubic');
    end
end


psf_lensn_group1=gpuArray(psf_lensn_group1);

psf_lensn_group1=65535.*psf_lensn_group1./(max(psf_lens_group1(:)));

%%
% prepare for reconstruction

for viewx=1:4
    PSF_power_z1=sum(sum(squeeze(psf_lensn_group1(:,:,:,viewx)),1),2);
    PSF_power_z1=squeeze(PSF_power_z1);
    for zd=1:81
        psf_lensn_group1(:,:,zd,viewx)=psf_lensn_group1(:,:,zd,viewx)./PSF_power_z1(zd)*PSF_power_z1(10);
    end
end

Nz=size(psf_lensn_group1,3);
PSF_power1=sum(psf_lensn_group1(:));
PSF_power_zn1=sum(sum(sum(psf_lensn_group1,1),2),4);
PSF_power_zn1=single(squeeze(repmat(PSF_power_zn1,2*ReconROIx+1,2*ReconROIy+1,1,1)));
PSF_power_zn1=gather(PSF_power_zn1);

% prepare for calibration
x=[-FitROIx:FitROIx];
y=[-FitROIy:FitROIy];
[x y]=meshgrid(y,x);
x=single(x);
y=single(y);
LensN1=size(Info1,1);

%% Group1 Reconstruction
disp('reconstruction for group 2');
for ii=1:LensN1

     x_f_shift(:,:,ii)=x;
     y_f_shift(:,:,ii)=y;
     x_ff_shift(:,:,ii)=x;
     y_ff_shift(:,:,ii)=y;

end


ImgExp1=gpuArray(single(Img_group1));

%%
psf_lensn_group1=single(psf_lensn_group1);
ObjRecon1=ones(2*ReconROIx+1,2*ReconROIy+1,Nz,'single','gpuArray'); % estimated object    
ImgEstROI1=zeros(2*FitROIx+1,2*FitROIy+1,LensN1,'single','gpuArray'); % estimated sub-image
RatioROI1=zeros(2*ReconROIx+1,2*ReconROIy+1,LensN1,'single'); % sub-ratio
RatioAvg1=ObjRecon1*0; 
ImgEst1=Img_group1*0;
Ratio1=ImgEst1+1;
RatioAvg1=gather(RatioAvg1);
RatioROI_lens1=zeros(2*ReconROIx+1,2*ReconROIy+1,LensN1,'single');
%%
for iti=1:ItN
   
     psf_lensn_group1=gpuArray(psf_lensn_group1);
    ImgEstROI1=gpuArray(ImgEstROI1);
     ObjRecon1=gpuArray(ObjRecon1);
   for jj=1:LensN1
        ImgEstROI1(:,:,jj)=sum(real(ifft2(fft2(ifftshift(ifftshift(...
           padarray(psf_lensn_group1(:,:,:,jj),[FitROIx-PSFROIx FitROIy-PSFROIy 0],0,'both')...
          ,1),2)).*fft2(padarray(ObjRecon1,[FitROIx-ReconROIx FitROIy-ReconROIy],0,'both')))),3);
    end 
   
    ImgEst1=ImgEst1*0;
   
    for ii=1:LensN1
        Img_lens=interp2(x,y,ImgEstROI1(:,:,ii),x_ff_shift(:,:,ii),y_ff_shift(:,:,ii),'linear',0);
          Img_lens=gather(Img_lens);
           ImgEst1(:,:,ii)=  ImgEst1(:,:,ii)+Img_lens(FitROIx-ReconROIx+1:FitROIx+ReconROIx+1,FitROIy-ReconROIy+1:FitROIy+ReconROIy+1);
    end
    ImgEst1=ImgEst1/(PSF_power1/Nz);
    Ratio1=ImgExp1./(ImgEst1+mean(ImgEst1(:))/SNR);
    Ratio_exp=Ratio1; 
    RatioAvg1=gpuArray(RatioAvg1);
     RatioROI1=gpuArray(RatioROI1);
   RatioROI_lens1=gpuArray(RatioROI_lens1);
    for ii=1:LensN1
        RatioROI_lens1(:,:,ii)=interp2(x,y, Ratio_exp(:,:,ii)...
             ,x_f_shift(:,:,ii),y_f_shift(:,:,ii),'linear',0);  
        RatioROI1(:,:,ii)=RatioROI_lens1(FitROIx-ReconROIx+1:FitROIx+ReconROIx+1,FitROIy-ReconROIy+1:FitROIy+ReconROIy+1,ii);
    end
    RatioAvg1=RatioAvg1*0;
   Ratio_exp=gather(Ratio_exp);
    RatioROI_lens1=gather(RatioROI_lens1);
    Ratio1=gather(Ratio1);
    ImgEst1=gather(ImgEst1);
    ImgEstROI1=gather(ImgEstROI1);
    ImgExp1=gather(ImgExp1);
    ObjRecon1=gather(ObjRecon1);

    for ii=1:LensN1
       RatioAvg1=RatioAvg1+max(real(ifft2(fft2(repmat(RatioROI1(:,:,ii),1,1,Nz))...
            .*conj(fft2(ifftshift(ifftshift(...
           padarray(psf_lensn_group1(:,:,:,ii),[ReconROIx-PSFROIx ReconROIy-PSFROIy 0],0,'both')...
            ,1),2))))),0)./PSF_power_zn1;
    end    

    ObjRecon1=ObjRecon1.*RatioAvg1;
 
    % show MIP of reconstruction result
    if show_reconresult==1
        figure(1000);
        subplot(1,3,1);
        imagesc(squeeze(max(ObjRecon1(:,:,:),[],3)));
        title (['iteration ' num2str(iti) ' xy MIP']);
        xlabel ('x');
        ylabel ('y');
        axis equal;
        subplot(1,3,2);
        imagesc(squeeze(max(ObjRecon1(:,:,:),[],2)));
        title(['iteration ' num2str(iti) ' yz MIP']);
        xlabel('z');
        ylabel('y');
        axis equal;
        subplot(1,3,3);
        imagesc(squeeze(max(ObjRecon1(:,:,:),[],1)));
        title(['iteration ' num2str(iti) ' xz MIP']);
        xlabel('z');
        ylabel('x');
        axis equal;
        drawnow
    else
    end
    RatioROI1=gather(RatioROI1);
    RatioAvg1=gather(RatioAvg1);
    psf_lensn_group1=gather(psf_lensn_group1);
    
end

%%
tic;
frame=max(ObjRecon1,[],3);
frame(frame<mean(frame(:))+2)=0;
L=bwlabel(frame);

 stat=regionprops('table',L,'Area','Centroid',...
   'MajorAxisLength',"Orientation");
Ar=cat(1,stat.Area);
ind=find(Ar==max(Ar));
frame(find(L~=ind))=0;
mask=frame;
mask(mask>0)=1;
nums=find(mask>0);
nums=size(nums,1);
for znn=1:size(ObjRecon1,3)
    ObjRecon1x(:,:,znn)=ObjRecon1(:,:,znn).*mask;
    zobj=squeeze(ObjRecon1x(:,:,znn));
    image_power_z1(znn,1)=sum(zobj(:))./nums;
end

 image_power_z1mean=mean(image_power_z1,1);
  for znn=1:size(ObjRecon1,3) 
     ObjRecon1(:,:,znn)=ObjRecon1(:,:,znn)./image_power_z1(znn,1).*image_power_z1mean;
  end
imagepower1(:,time)=gather(image_power_z1);
toc;
ObjRECON_t_group1(:,:,:,time)=gather(ObjRecon1);

end

ObjRECON_t_group1=uint16(round(65535.*ObjRECON_t_group1./max(ObjRECON_t_group1(:))));
 save ([savepath, 'ObjRecon_parts1.mat'], 'ObjRECON_t_group1','imagepower1', '-v7.3');
 
