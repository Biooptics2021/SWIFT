function [image_part]=imageparts_BFOV_demo(stack_image,centers_part,structured_mode,views_range)

%单次仅对特定视场做分割。

image_part=zeros(1339,1331,7,structured_mode,floor(size(stack_image,3)/structured_mode),'single');  
%                 x,  y,  views, structured_mode, time
sub1=0;
sub2=0;

for views=1:7
    
 x0=views_range(views,1);
    y0=views_range(views,3);
    
sub1=sub1+(centers_part(views,1)-x0)/7;
sub2=sub2+(centers_part(views,2)-y0)/7;

end
sub1=round(sub1)-10;
sub2=round(sub2)-9;

for t=1:floor(size(stack_image,3)/structured_mode)
    
    tic;
   
    
      
        for structure_image=1:structured_mode  %此时先写个最简单的wide模式
            
image=squeeze(stack_image(:,:,structured_mode*t-structured_mode+structure_image));


image1=zeros(5120,4000,'single');
image1(:,:)=image;
clear image;
image=image1;
clear image1;

for views=1:7
    
   
image_part(:,:,views,structure_image,t)=image(centers_part(views,1)-sub1:centers_part(views,1)-sub1+1338,centers_part(views,2)-sub2:centers_part(views,2)-sub2+1330); 

end

        end
   
    
toc;
end
