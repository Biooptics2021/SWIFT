function [views1,views2]=image_seg_free_zebrafish ( V,stacks,center,fit_group1,fit_group2,frame_num)

 load(['laser_profile.mat']);

center=center+200;

for ii=1:frame_num
    
    v(1)=V(2,ii);
    v(2)=V(1,ii);
    
    image=stacks(:,:,ii);
    
    image=padarray(image,[200 200],0, 'both');
    
    views1(:,:,1,ii)=image(center(2)-265:center(2)+265,center(1)-270:center(1)+270);
    
    for group1=2:4
        
        fit_result1=fit_group1{group1-1,1};
        fit_result2=fit_group1{group1-1,2};
        
        delta2560=fit_result1(v(1),v(2));
        delta2160=fit_result2(v(1),v(2));
        views1(:,:,group1,ii)=image(round(center(2)+delta2160-265):round(center(2)+delta2160+265),round(center(1)+delta2560-270):round(center(1)+delta2560+270));
        clear  delta2560 delta2160 fit_result1 fit_result2
    end
    
    views2(:,:,1,ii)=image(center(2)-265:center(2)+265,center(1)-270:center(1)+270);
    
    for group2=2:4
        
        fit_result1=fit_group2{group2-1,1};
        fit_result2=fit_group2{group2-1,2};
        
        delta2560=fit_result1(v(1),v(2));
        delta2160=fit_result2(v(1),v(2));
        views2(:,:,group2,ii)=image(round(center(2)+delta2160-265):round(center(2)+delta2160+265),round(center(1)+delta2560-270):round(center(1)+delta2560+270));
        clear  delta2560 delta2160 fit_result1 fit_result2
    end
    
end
1;
 for ii=1:size(views2,4)
     for viewsn=1:size(views2,3)
         views2(:,:,viewsn,ii)=uint16(single(squeeze(views2(:,:,viewsn,ii))).*laser_views2(:,:,viewsn));
         views1(:,:,viewsn,ii)=uint16(single(squeeze(views1(:,:,viewsn,ii))).*laser_views1(:,:,viewsn));
     end
 end
 
 for ii=1:size(views2,4)
 temp=medfilt2(squeeze(views2(:,:,1,ii)),[3 3]);
 views2(430:end,:,3,ii)=temp(430:end,:)*1.78;
 views2(1:78,:,2,ii)=temp(1:78,:)*1.78;
 end
 
 for ii=1:size(views1,4)
 temp=medfilt2(squeeze(views1(:,:,1,ii)),[3 3]);
 views1(:,486:end,2,ii)=temp(:,486:end)*1.78;
 views1(464:end,:,3,ii)=temp(464:end,:)*1.78;
 views1(1:120,:,3,ii)=temp(1:120,:)*1.78;
 end
    
            
       
        
    