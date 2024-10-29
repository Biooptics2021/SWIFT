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
group1_info=[1,2,3,4]; % The index of views in group1 
psf_group1=single(psf_group1);


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
   stack_seg_wide_group1=squeeze(image_part(:,:,group1_info,1,:)); %(x, y, views, t)
end

if structured_mode==3 
   [SIM1,uniform1]=SIM_demo(image_part);
   stack_seg_wide_group1=squeeze(SIM1(:,:,group1_info,:)); %(x, y, views, t)
end
%Wide field mode and sturctured illumination mode are switched here.

psf_lens_group1=squeeze(psf_group1(:,:,:,group1_info,VoVnum(fovs)));
psf_lens_group2=squeeze(psf_group2(:,:,:,group1_info,VoVnum(fovs)));
%Info=[1,2,3,4]'; % If you want to delete this view, delete its number from Info.

times=size(stack_seg_wide_group1,4);
disp(['total fov frame is ' num2str(times)]);
image_seg=stack_seg_wide_group1(:,:,:,1);

ObjRECON_t_group1=zeros(size(image_seg,1),size(image_seg,2),size(psf_lens_group1,3),times,'single'); 
% (x,y,z,t) for each sub-VoV

for time=1:times
disp(['This is time : ' num2str(time)]);
image_seg1=stack_seg_wide_group1(:,:,:,time);
Img_group1=image_seg1;
Img_group1=single(Img_group1);

ximg=size(Img_group1,1);
yimg=size(Img_group1,2);

factor=1;

Img_group1= imresize(Img_group1,[factor*ximg factor*yimg],'cubic');

% This part is for Pixel-wise alignment:
for views=2:4
    fitresult1339=fitresults{views-1,1,VoVnum(fovs)};
    fitresult1331=fitresults{views-1,2,VoVnum(fovs)};
    % The fitresults are early calibrated by pixel-wise alignment.
    [x y]=meshgrid(1:1331,1:1339);
    [b]=fitresult1339(x,y);
    y=y+b;
    [a]=fitresult1331(x,y);
    x=x+a;
    [xq,yq]=meshgrid(1:1331,1:1339);
    Img_group1(:,:,views)=interp2(xq,yq,Img_group1(:,:,views),x(:,:),y(:,:),'cubic',0);
    % Img_group1 is an aberration corrected 3D image (x,y,views).
end


% To conpensate the non-uniform illumination by laser.
for views=2:4
    Img_group1(:,:,views)=Img_group1(:,:,views).*intensity_inverse_group1(:,:,views);
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

psf_lensn_group1=zeros(factor*99,factor*101,size(psf_lens_group1,3),size(group1_info,1),'single');

for ii=1:size(Info1,1)
psf_lensn_group1(:,:,:,ii)=psf_lens_group1(:,:,:,Info1(ii));
end

psf_lensn_group1=gpuArray(psf_lensn_group1);

%%
% prepare for reconstruction
Nz=size(psf_lensn_group1,3);
PSF_power1=sum(psf_lensn_group1(:));
PSF_power_zn1=sum(sum(sum(psf_lensn_group1,1),2),4);
PSF_power_zn1=single(squeeze(repmat(PSF_power_zn1,2*ReconROIx+1,2*ReconROIy+1,1,1)));
PSF_power_zn1=gather(PSF_power_zn1);
x=[-FitROIx:FitROIx];
y=[-FitROIy:FitROIy];
[x y]=meshgrid(y,x);
x=single(x);
y=single(y);
LensN1=size(Info1,1);

%% Group1 Reconstruction
disp('reconstruction for group 1');
for ii=1:LensN1
    x_f_shift(:,:,ii)=x;
    y_f_shift(:,:,ii)=y;
    x_ff_shift(:,:,ii)=x;
    y_ff_shift(:,:,ii)=y;
 end

ImgExp1=gpuArray(single(Img_group1));

psf_lensn_group1=single(psf_lensn_group1);
ObjRecon1=ones(2*ReconROIx+1,2*ReconROIy+1,Nz,'single','gpuArray'); 
ImgEstROI1=zeros(2*FitROIx+1,2*FitROIy+1,LensN1,'single','gpuArray'); 
RatioROI1=zeros(2*ReconROIx+1,2*ReconROIy+1,LensN1,'single'); 
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

ObjRECON_t_group1(:,:,:,time)=gather(ObjRecon1);

end

ObjRECON_t_group1=uint16(round(65535.*ObjRECON_t_group1./max(ObjRECON_t_group1(:))));
 save ([savepath, '\ObjRecon_parts1.mat'], 'ObjRECON_t_group1', '-v7.3');
 
end