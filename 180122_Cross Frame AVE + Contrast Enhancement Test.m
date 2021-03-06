clear all
%
Folder_Path='C:\TuanShu\180122_Dynamic Focusing Image Guiding Test\Bilinear_Extract step 2_Gain0_it10_2\';
Center=[500 600];
Number_of_Frame_AVE=100;
Half_Width=100;
ROI=[Center(1)-Half_Width Center(1)+Half_Width Center(2)-Half_Width Center(2)+Half_Width];
Normalization='count';  %'count' (default) | 'probability' | 'countdensity' | 'pdf' | 'cumcount' | 'cdf'
%

Data_Save_Folder=[Folder_Path '\Processed Data\'];
if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end
%
File_list=dir(Folder_Path);
Original_Length_file_list=length(File_list);
for p=1:Original_Length_file_list
    [pathstr,name,ext] = fileparts(File_list(Original_Length_file_list+1-p).name);
    if File_list(Original_Length_file_list+1-p).isdir ~= 0
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(File_list(Original_Length_file_list+1-p).name,'.') == 1
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(File_list(Original_Length_file_list+1-p).name,'..') == 1
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(ext,'.PNG')
        File_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(File_list(Original_Length_file_list+1-p).name,'Processed Data') == 1
        File_list(Original_Length_file_list+1-p)=[];
    end
end
Length_file_list=length(File_list);
%
Image_Temp=0;
for p=1:min(Length_file_list,Number_of_Frame_AVE)
    Image_Temp=Image_Temp+double(imread([Folder_Path '\' File_list(p).name],'tif'));
    disp(p);
end
Image_Temp=(Image_Temp/min(Length_file_list,Number_of_Frame_AVE));
%%
NumberofTilt=16;
subplot(1,1,1)
imagesc(tonemap(Image_Temp,'AdjustLightness',[0.1 1],'AdjustSaturation',1,'NumberOfTiles',[NumberofTilt NumberofTilt]));
axis equal
axis off
%%
imwrite(uint8(Image_Temp),[Data_Save_Folder '\' File_list(p).name sprintf('_AVE%g_ori contrast.tif',Number_of_Frame_AVE)],'tiff');

%% Contrast Enhancement
Bit=8;  
Brightness_Level=128;
Contrast_Degree=65;
A=Brightness_Level;
B=(2^Bit)/tan(Contrast_Degree/180*pi);
Image_Enhanced=((Image_Temp)/B*(2^Bit)+A-2^(Bit-1));
subplot(1,2,1)
imagesc(tonemap(Image_Enhanced,'AdjustLightness',[0.01 1],'AdjustSaturation',2,'NumberOfTiles',[16 16]));
axis equal
axis off
subplot(1,2,2)
imagesc(tonemap(Image_Temp,'AdjustLightness',[0.01 1],'AdjustSaturation',2,'NumberOfTiles',[16 16]));
axis equal
axis off
%histogram(Image_Enhanced);
%% HSV CLAHE
Image_HSV=rgb2hsv(Image_Enhanced);
subplot(1,1,1)
imagesc(Image_HSV(:,:,1));
%%
Value=Image_HSV(:,:,3);
Value(Image_HSV(:,:,3)<0)=0;
Image_HSV(:,:,3)=Value;
imagesc(Image_HSV(:,:,3));
%%
Image_HE=adapthisteq(Image_HSV(:,:,3));
imagesc(Image_HE./max(Image_HE(:)));
caxis([0 0.03]);
colormap(gray);
%imwrite(Image_Temp,[Data_Save_Folder '\' 'AVE.tif'],'tiff');
