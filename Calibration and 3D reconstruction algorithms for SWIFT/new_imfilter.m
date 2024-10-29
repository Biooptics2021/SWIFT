function [image_fuse]=new_imfilter(SIM1,uniform1,sigmaW)

alpha=0.2;
waveletGaussiansRatio = 2;
sigmaLo = 2.86 * sigmaW * (1 + waveletGaussiansRatio/2);

for tt=1:size(uniform1,4)
    
    for views=1:size(uniform1,3)
        low_uniform(:,:,views,tt)= imgaussfilt(uniform1(:,:,views,tt),sigmaLo);
        high_uniform(:,:,views,tt)=uniform1(:,:,views,tt)-low_uniform(:,:,views,tt);
        low_sim(:,:,views,tt)=imgaussfilt(SIM1(:,:,views,tt),sigmaLo);
        image_fuse(:,:,views,tt)=alpha.*low_sim(:,:,views,tt)+high_uniform(:,:,views,tt);
         image_fuse(:,:,views,tt)= image_fuse(:,:,views,tt);
    end
    
    image_fuse(image_fuse<0)=0;
    
end