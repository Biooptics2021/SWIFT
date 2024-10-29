function [stack1] = readraw1(M,N,type,namelist,pathname)
% type is 8bit/16bit 'uint8'/'uint16'

%M=5120;
%N=5120;
i1=0;
img_num=size(namelist,1)
stack1=zeros(M,N,img_num,'single');
for number=1:img_num
    tic;
    i1=i1+1;

filename=namelist(number).name;
f1 = fopen([pathname, filename], 'r');
data = fread(f1, type);
fclose(f1);
len = length(data);
%M=5120;
%N=5120;
k = len/(M*N);
image = reshape(data, M, N, k);
stack1(:,:,i1)=image;
toc;
end