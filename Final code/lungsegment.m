function [sout]=lungsegment(Input)

%%%%% To resize an input image file
Input = imresize(Input,[128,128]);
[r c p]=size(Input);

if p==3
   Input =rgb2gray(Input);
end

Input = medfilt2(Input);
Input   = double(Input);
[r c]   = size(Input);
Length  = r*c; 
wd1=r;
wd2=c;
Dataset = reshape(Input,[Length,1]);

No = inputdlg('Input the No of Cluster(3/4):');
No = cell2mat(No);
No = str2num(No);
   
if No == 3
[class]=FSeg(Dataset,No);
AA1=reshape(class(1,:),wd1,wd2); %reshape class 1 into an image
AA2=reshape(class(2,:),wd1,wd2); %reshape class 1 into an image
AA3=reshape(class(3,:),wd1,wd2); %reshape class 1 into an image

cd Clusim
imwrite(AA1,'seg1.bmp');
imwrite(AA2,'seg2.bmp');
imwrite(AA3,'seg3.bmp');
cd ..

figure('Name','Segmented Results','MenuBar','none');
subplot(2,2,1);imshow(AA1,[]);
title('Region1');
subplot(2,2,2);imshow(AA2,[]);
title('Region2');
subplot(2,2,3);imshow(AA3,[]);
title('Region3');

elseif No==4         

[class]=FSeg(Dataset,No); %try with 4 classes
AA1=reshape(class(1,:),wd1,wd2); %reshape class 1 into an image
AA2=reshape(class(2,:),wd1,wd2); %reshape class 1 into an image
AA3=reshape(class(3,:),wd1,wd2); %reshape class 1 into an image
AA4=reshape(class(4,:),wd1,wd2); %reshape class 1 into an image

cd Clusim
imwrite(AA1,'seg1.bmp');
imwrite(AA2,'seg2.bmp');
imwrite(AA3,'seg3.bmp');
imwrite(AA4,'seg4.bmp');
cd ..


figure('Name','Segmented Results','MenuBar','none');
subplot(2,2,1);imshow(AA1,[]);
subplot(2,2,2);imshow(AA2,[]);
subplot(2,2,3);imshow(AA3,[]);
subplot(2,2,4);imshow(AA4,[]);

else
    
    warndlg('Enter the Valid Choice');

end

%%%%%%%%% Input the Segmeted Image into Post processing
cd Clusim
   [file path] = uigetfile('*.bmp','Pick a Segmented Image File');
   I = imread(file);
cd ..
Bw = im2bw(I);

figure('Name','PostProcessing','MenuBar','none');
subplot(2,2,1);
imshow(Bw), title('Input Image');

dim = size(Bw);
row = round(dim(1)/2);
col = min(find(Bw(row,:)));

%%%%Tracing image
boundary = bwtraceboundary(Bw,[row, col],'W');

subplot(2,2,2);
imshow(Bw), title('Traced boundary');

hold on;
%%%%Display traced boundarytrack
plot(boundary(:,2),boundary(:,1),'g','LineWidth',2);
hold off

%%% Filling inner region
Out_image=imfill(Bw,'holes');
subplot(2,2,3);
imshow(Out_image);
title('Region Filling');

%%% Lung Region Extraction
rr=abs(imsubtract(Bw,Out_image));

subplot(2,2,4);
imshow(rr);
title('Lung parenchyma region');

sout=imcomplement(rr);

% figure('Name','Segmented Image','MenuBar','none');
% imshow(sout,[]);

