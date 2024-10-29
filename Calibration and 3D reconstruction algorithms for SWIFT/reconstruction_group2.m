%% This code is modified from CLFM_fish_Reconstruction_CPU.m file
%Reference: Zhang, Z., Bai, L., Cong, L. et al. Imaging volumetric dynamics at high speed in mouse and zebrafish brain with confocal light field microscopy. Nat Biotechnol 39, 74–83 (2021).

%% multi-subVoVs volumetric reconstruction
% Scanning wide-field tomography for high-speed, mesoscale-volumetric imaging of biodynamics in vivo
% SYSTEM REQUIREMENTS
% Memory: 128 GB RAM
% Graphics: 24 GB RAM
% MATLAB: R2021a

%% This code involving several functional modules: 
%1. Computational optical sectionning by SIM
%2. Pixel-wise alignment for multi-subVoVs
%3. Dual-focus reconstruction
%4. 3D reconstruction by RL deconvolution

%%  Parameters
% data structure:    
%                    PSF (psfx, psfy, zstacks, views, fovs)
%       for example       (79,   81,   15,      4,     9   )
%                    Image_Segment (x,    y,    views, structured modes, fovs, times)
%       for example                (1331, 1339, 4,     1/3,              9,    100  )
%       structrue modes : 
%       1:Wide illumination: 1-wide      3:Structured illumination: 1-sim1 2-sim2 3-sim3

clear;      
close all

%% set workpath to this code
load('supporting_data.mat');
local_address=mfilename('fullpath');
[pathstr,namestr]=fileparts(local_address);
cd(pathstr);
addpath(pathstr);

%% Get image name list

[namelist,pathname] = get_imagelist;
image_nums=size(namelist,1);
VoVS=9;  %The numbers


 subVoV1=4:9:image_nums;  %% subVoV1 list
 subVoV2=3:9:image_nums;  %% subVoV2
 subVoV3=2:9:image_nums;  %% subVoV3
 subVoV4=5:9:image_nums;  %% subVoV4
 subVoV5=9:9:image_nums;  %% subVoV5
 subVoV6=1:9:image_nums;  %% subVoV6
 subVoV7=6:9:image_nums;  %% subVoV7
 subVoV8=7:9:image_nums;  %% subVoV8 
 subVoV9=8:9:image_nums;  %% subVoV9 list
 VoV=[subVoV1;subVoV2;subVoV3;subVoV4;subVoV5;subVoV6;subVoV7;subVoV8;subVoV9]; % Full VoV lists
 
%% load psf files for reconstruction

structured_mode=1;
group2_info=[1,5,6,7]; % The index of views in group1 
psf_group2=single(psf_group2);


%% Load raw images and Segment image now.
if VoVS>8
    VoVnum=[1,2,3,4,5,6,7,8,9];
end


for fovs=1:VoVS  % or size(image_part1,5)
    
   [stack_thisfov] = readraw1(5120,4000,'uint8',namelist(VoV(fovs,:)),pathname);% give the list of this fov 
   [image_part]=imageparts_BFOV_demo(stack_thisfov,psf_part01,structured_mode,views_range); % Imageseg: Raw images to multi-views images
 
   % Parameters now:
   % image_part: (x, y, views, structured_mode, t)
   % structured_mode : defined in the begining of this code
   % centers_part_psf_fov: (views, x and y, fovs)
   % views_range : the edge of each views
   clear stack_thisfov;
   
if structured_mode==1 
   stack_seg_wide_group2=squeeze(image_part(:,:,group2_info,1,:)); %(x, y, views, t)
end

if structured_mode==3 
   [SIM1,uniform1]=SIM_demo(image_part);
   stack_seg_wide_group2=squeeze(SIM1(:,:,group2_info,:)); %(x, y, views, t)
end
%Wide field mode and sturctured illumination mode are switched here.

psf_lens_group1=squeeze(psf_group1(:,:,:,group1_info,VoVnum(fovs)));
psf_lens_group2=squeeze(psf_group2(:,:,:,[1,2,3,4],VoVnum(fovs)));
%Info=[1,2,3,4]'; % If you want to delete this view, delete its number from Info.

times=size(stack_seg_wide_group2,4);
disp(['total fov frame is ' num2str(times)]);
image_seg=stack_seg_wide_group2(:,:,:,1);

ObjRECON_t_group2=zeros(size(image_seg,1),size(image_seg,2),size(psf_lens_group2,3),times,'single'); 
% (x,y,z,t) for each sub-VoV

for time=1:times
disp(['This is time : ' num2str(time)]);
image_seg2=stack_seg_wide_group2(:,:,:,time);
Img_group2=image_seg2;
Img_group2=single(Img_group2);

ximg=size(Img_group2,1);
yimg=size(Img_group2,2);

factor=1;

Img_group2= imresize(Img_group2,[factor*ximg factor*yimg],'cubic');

% This part is for Pixel-wise alignment:
for views=5:7
    fitresult1339=fitresults{views-1,1,VoVnum(fovs)};
    fitresult1331=fitresults{views-1,2,VoVnum(fovs)};
    % The fitresults are early calibrated by pixel-wise alignment.
    [x y]=meshgrid(1:1331,1:1339);
    [b]=fitresult1339(x,y);
    y=y+b;
    [a]=fitresult1331(x,y);
    x=x+a;
    [xq,yq]=meshgrid(1:1331,1:1339);
    Img_group2(:,:,views)=interp2(xq,yq,Img_group2(:,:,views-3),x(:,:),y(:,:),'cubic',0);
    % Img_group2 is an aberration corrected 3D image (x,y,views)in group2.
end


% To conpensate the non-uniform illumination by laser.
for views=2:4
    Img_group2(:,:,views)=Img_group2(:,:,views).*intensity_inverse_group2(:,:,views);
end

savepath=[pathstr '\demo_Reconstruction\fov' num2str(fovs)];

if exist(savepath)==7,
   ; 
else
   mkdir(savepath);   
end
 
 Info1=[1:4]';

% reconstruction parameters
ReconROIx=size(image_seg,1)/2-0.5; 
ReconROIy=size(image_seg,2)/2-0.5;
FitROIx=size(image_seg,1)/2-0.5; 
FitROIy=size(image_seg,2)/2-0.5; 

% half-width of calibration area (>=ReconROI)
ItN=10; % number of iterations
SNR=200; % estimated signal noise ratio
show_reconresult=1; % show max intensity projection of reconstruction result every iteration(1) or not(0)
PSFROIx=49;
PSFROIy=50;

psf_lensn_group2=zeros(factor*99,factor*101,size(psf_lens_group2,3),size(group1_info,1),'single');

for ii=1:size(Info1,1)
psf_lensn_group2(:,:,:,ii)=psf_lens_group2(:,:,:,Info1(ii));
end

psf_lensn_group2=gpuArray(psf_lensn_group2);

%%
% prepare for reconstruction
Nz=size(psf_lensn_group2,3);
PSF_power2=sum(psf_lensn_group2(:));
PSF_power_zn2=sum(sum(sum(psf_lensn_group2,1),2),4);
PSF_power_zn2=single(squeeze(repmat(PSF_power_zn2,2*ReconROIx+1,2*ReconROIy+1,1,1)));
PSF_power_zn2=gather(PSF_power_zn2);
x=[-FitROIx:FitROIx];
y=[-FitROIy:FitROIy];
[x y]=meshgrid(y,x);
x=single(x);
y=single(y);
LensN2=size(Info1,1);

%% Group1 Reconstruction
disp('reconstruction for group 1');
for ii=1:LensN2
    x_f_shift(:,:,ii)=x;
    y_f_shift(:,:,ii)=y;
    x_ff_shift(:,:,ii)=x;
    y_ff_shift(:,:,ii)=y;
 end

ImgExp2=gpuArray(single(Img_group2));

psf_lensn_group2=single(psf_lensn_group2);
ObjRecon2=ones(2*ReconROIx+1,2*ReconROIy+1,Nz,'single','gpuArray'); 
ImgEstROI2=zeros(2*FitROIx+1,2*FitROIy+1,LensN2,'single','gpuArray'); 
RatioROI2=zeros(2*ReconROIx+1,2*ReconROIy+1,LensN2,'single'); 
RatioAvg2=ObjRecon2*0; 
ImgEst2=Img_group2*0;
Ratio2=ImgEst2+1;
RatioAvg2=gather(RatioAvg2);
RatioROI_lens2=zeros(2*ReconROIx+1,2*ReconROIy+1,LensN2,'single');

%%
for iti=1:ItN
 
     psf_lensn_group2=gpuArray(psf_lensn_group2);
    ImgEstROI2=gpuArray(ImgEstROI2);
     ObjRecon2=gpuArray(ObjRecon2);
   for jj=1:LensN2
        ImgEstROI2(:,:,jj)=sum(real(ifft2(fft2(ifftshift(ifftshift(...
           padarray(psf_lensn_group2(:,:,:,jj),[FitROIx-PSFROIx FitROIy-PSFROIy 0],0,'both')...
          ,1),2)).*fft2(padarray(ObjRecon2,[FitROIx-ReconROIx FitROIy-ReconROIy],0,'both')))),3);
    end  
   
    ImgEst2=ImgEst2*0;
   
    for ii=1:LensN2
        Img_lens=interp2(x,y,ImgEstROI2(:,:,ii),x_ff_shift(:,:,ii),y_ff_shift(:,:,ii),'linear',0);
          Img_lens=gather(Img_lens);
           ImgEst2(:,:,ii)=  ImgEst2(:,:,ii)+Img_lens(FitROIx-ReconROIx+1:FitROIx+ReconROIx+1,FitROIy-ReconROIy+1:FitROIy+ReconROIy+1);
    end
    ImgEst2=ImgEst2/(PSF_power2/Nz);
    Ratio2=ImgExp2./(ImgEst2+mean(ImgEst2(:))/SNR);
    Ratio_exp=Ratio2;  
    RatioAvg2=gpuArray(RatioAvg2);
     RatioROI2=gpuArray(RatioROI2);
   RatioROI_lens2=gpuArray(RatioROI_lens2);
    for ii=1:LensN2
        RatioROI_lens2(:,:,ii)=interp2(x,y, Ratio_exp(:,:,ii)...
             ,x_f_shift(:,:,ii),y_f_shift(:,:,ii),'linear',0);
        RatioROI2(:,:,ii)=RatioROI_lens2(FitROIx-ReconROIx+1:FitROIx+ReconROIx+1,FitROIy-ReconROIy+1:FitROIy+ReconROIy+1,ii);
    end
    RatioAvg2=RatioAvg2*0;

   Ratio_exp=gather(Ratio_exp);
    RatioROI_lens2=gather(RatioROI_lens2);
    Ratio2=gather(Ratio2);
    ImgEst2=gather(ImgEst2);
    ImgEstROI2=gather(ImgEstROI2);
    ImgExp2=gather(ImgExp2);
    ObjRecon2=gather(ObjRecon2);

    for ii=1:LensN2
       RatioAvg2=RatioAvg2+max(real(ifft2(fft2(repmat(RatioROI2(:,:,ii),1,1,Nz))...
            .*conj(fft2(ifftshift(ifftshift(...
           padarray(psf_lensn_group2(:,:,:,ii),[ReconROIx-PSFROIx ReconROIy-PSFROIy 0],0,'both')...
            ,1),2))))),0)./PSF_power_zn2;
    end    
  
    ObjRecon2=ObjRecon2.*RatioAvg2;
 
    % show MIP of reconstruction result
    if show_reconresult==1
        figure(1000);
        subplot(1,3,1);
        imagesc(squeeze(max(ObjRecon2(:,:,:),[],3)));
        title (['iteration ' num2str(iti) ' xy MIP']);
        xlabel ('x');
        ylabel ('y');
        axis equal;
        subplot(1,3,2);
        imagesc(squeeze(max(ObjRecon2(:,:,:),[],2)));
        title(['iteration ' num2str(iti) ' yz MIP']);
        xlabel('z');
        ylabel('y');
        axis equal;
        subplot(1,3,3);
        imagesc(squeeze(max(ObjRecon2(:,:,:),[],1)));
        title(['iteration ' num2str(iti) ' xz MIP']);
        xlabel('z');
        ylabel('x');
        axis equal;
        drawnow
    else
    end
    RatioROI2=gather(RatioROI2);
    RatioAvg2=gather(RatioAvg2);
    psf_lensn_group2=gather(psf_lensn_group2);
    
end

ObjRECON_t_group2(:,:,:,time)=gather(ObjRecon2);

end

ObjRECON_t_group2=uint16(round(65535.*ObjRECON_t_group2./max(ObjRECON_t_group2(:))));
 save ([savepath, '\ObjRecon_parts2.mat'], 'ObjRECON_t_group2', '-v7.3');
 
end