function [namelist,pathname] = get_imagelist

[filename, pathname] = uigetfile('*.raw', '读入图像');
filename =[pathname,'\*.raw'];
namelist = dir(filename);%Get the namelist of all raw images
img_num = length(namelist);%Get the total number of all raw images
