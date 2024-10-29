function [SIM1,uniform1]=SIM_demo(image_part)

display(['SIM_demo ' ]);

tic;

for views=1:size(image_part,3)
    for i=1:size(image_part,5)  
        SIM1(:,:,views,i)=((image_part(:,:,views,1,i)-image_part(:,:,views,2,i)).^2+(image_part(:,:,views,1,i)-image_part(:,:,views,3,i)).^2+(image_part(:,:,views,2,i)-image_part(:,:,views,3,i)).^2).^0.5;
        uniform1(:,:,views,i)=(image_part(:,:,views,1,i)+image_part(:,:,views,2,i)+image_part(:,:,views,3,i))./3;
    end
end

disp('Using RSIM process');
SIM1=new_imfilter(SIM1,uniform1,20);

toc;