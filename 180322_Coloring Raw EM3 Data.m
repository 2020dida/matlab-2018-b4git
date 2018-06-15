clear all
fclose('all')
%% File Location
Folder_Path='C:\TuanShu\180209_EM3 Data for Processing\3D Data after PostAVE\';
File_Name='20171018181034_1x_3.5it_F_Back to 8.3mW_3D.raw';
%%
File_Path=[Folder_Path File_Name];
%% Parameter: Data
If_Need_Reverse=1;  
Known_Slide_Thickness=0.56; %micron
PostAVE_2=3;
FFT_Length=8;   %Unit: frame after PostAVE_2
FFTBin_For_Color=2;
Phase=0;        %如果Phase=0, 則自動取目前Frame前後的影像算Color
Color_Sigma=2;
%% Parameters: File
Row=1024;
Column=1024;
Data_Format='double';
if strcmp(Data_Format,'double')
    Byte_per_Pixel=8;
end
Frame_Skipped=100;
Maximum_Number_of_Frame_to_Read=500;
%% Parameters: File (Derived)
FileInfo=dir(File_Path);
Number_of_Frame_Total=FileInfo.bytes/Row/Column/Byte_per_Pixel;
Number_of_Frame_to_Read=min(Number_of_Frame_Total-Frame_Skipped,Maximum_Number_of_Frame_to_Read);
%% File Loading
Image_Stack=zeros(Row,Column,Number_of_Frame_to_Read);
fin=fopen(File_Path);
fseek(fin, Row*Column*Byte_per_Pixel*Frame_Skipped, 'bof'); 
for p=1:Number_of_Frame_to_Read
    Image_Stack(:,:,p)=fread(fin,[Row,Column],Data_Format);
    disp(p);
end
fclose(fin);

if If_Need_Reverse ==1
    Image_Stack=Image_Stack(:,:,end:-1:1);
end
%% 2nd PostAVE (有點像以下計算中會用到之最小厚度, 跟Color和成像都有關
Reduced_Stack=TSBinning(Image_Stack,3,PostAVE_2);
Reduced_Length=size(Reduced_Stack,3);
%% Look the data
Frame_Number=17;
ABTH=0.995;
CMin=0.05;
imagesc(TSAutoBrightness(Reduced_Stack(:,:,Frame_Number),ABTH,CMin));
colormap(gray)
axis equal
%% To Generate the Color Stack
Color_Stack=zeros(Row,Column,Reduced_Length);
for p=1:Reduced_Length
    Index_Start=max(1,p-floor(FFT_Length/2));
    Index_End=min(Reduced_Length,p+floor(FFT_Length/2)-1);
    Temp_Stack=Reduced_Stack(:,:,Index_Start:Index_End);
    Temp_Mean=mean(Temp_Stack,3);
    FFT_Stack=fft(Temp_Stack-repmat(Temp_Mean,[1 1 size(Temp_Stack,3)]),[],3);
    Temp_Color=imgaussfilt(abs(FFT_Stack(:,:,FFTBin_For_Color))./Temp_Mean,Color_Sigma);
    Color_Stack(:,:,p)=Temp_Color/max(Temp_Color(:));
    disp(p);
end
%% Look the Color
Frame_Number=157;
ABTH=0.995;
CMin=0.05;
imagesc(TSAutoBrightness(Color_Stack(:,:,Frame_Number),ABTH,CMin));
colormap(gray)
axis equal
%% Test
Frame_Number=30;
Frame_Thickness=8;
V_min=20;
V_max=90;

Cmin_adj=-0.06;
Cmax_adj=1.7;
Color_Map_Used=(Color_Stack(:,:,Frame_Number+Phase)-Cmin_adj)/(Cmax_adj-Cmin_adj);
RGB_Image=TSHSVcoloring(mean(Reduced_Stack(:,:,(Frame_Number-floor(Frame_Thickness/2)):(Frame_Number+floor(Frame_Thickness/2)-1)),3),Color_Map_Used,[V_min V_max]);
imagesc(RGB_Image);
%xlim([300 700]);
%ylim([300 700]);

axis equal
imwrite(RGB_Image,[Folder_Path sprintf('%s_FN%g_FT%g.tif',File_Name,Frame_Number,Frame_Thickness)],'tiff');
%% Normalized_Stack
Max_Array=max(max(Reduced_Stack,[],1),[],2);
Norm_Stack=Reduced_Stack./repmat(Max_Array,[Row Column 1]);
%% 
Frame_Number=80;
Frame_Thickness=8;
V_min=15;
V_max=80;

Cmin_adj=0.01;
Cmax_adj=0.55;
False_Color_Map=(Norm_Stack(:,:,Frame_Number+Phase)-Cmin_adj)/(Cmax_adj-Cmin_adj);

% subplot(1,2,1)
% histogram(Color_Stack(:,:,Frame_Number+Phase));
% xlim([0 1]);
% subplot(1,2,2)
% histogram(False_Color_Map);
% xlim([0 1]);
% %%
subplot(1,1,1)
RGB_Image_False=TSHSVcoloring(mean(Reduced_Stack(:,:,(Frame_Number-floor(Frame_Thickness/2)):(Frame_Number+floor(Frame_Thickness/2)-1)),3),False_Color_Map,[V_min V_max]);
imagesc(RGB_Image_False);
axis equal
imwrite(RGB_Image_False,[Folder_Path sprintf('%s_FN%g_FT%g_False.tif',File_Name,Frame_Number,Frame_Thickness)],'tiff');

%% False with Axial Blur
Blur_Window=5;
Norm_Blur_Stack=zeros(Row,Column,Reduced_Length);
for p=1:Reduced_Length
    Index_Start=max(1,p-floor(Blur_Window/2));
    Index_End=min(Reduced_Length,p+floor(Blur_Window/2)-1);
    Temp_Stack=Reduced_Stack(:,:,Index_Start:Index_End);
    Temp_Mean=imgaussfilt(mean(Temp_Stack,3),Color_Sigma);
    Norm_Blur_Stack(:,:,p)=Temp_Mean./max(Temp_Mean(:));
    disp(p);
end

%% False Color with blur map

Frame_Number=72.72727;
Frame_Thickness=8;
V_min=15;
V_max=80;
%
Cmin_adj=0.05;
Cmax_adj=0.86;

False_Color_Map=(Norm_Blur_Stack(:,:,Frame_Number+Phase)-Cmin_adj)/(Cmax_adj-Cmin_adj);

% subplot(1,2,1)
% histogram(Color_Stack(:,:,Frame_Number+Phase));
% xlim([0 1]);
% subplot(1,2,2)
% histogram(False_Color_Map);
% xlim([0 1]);
%
subplot(1,1,1)
RGB_Image_False_Blur=TSHSVcoloring(mean(Reduced_Stack(:,:,(Frame_Number-floor(Frame_Thickness/2)):(Frame_Number+floor(Frame_Thickness/2)-1)),3),False_Color_Map,[V_min V_max]);
imagesc(RGB_Image_False_Blur);
axis equal
imwrite(RGB_Image_False_Blur,[Folder_Path sprintf('%s_FN%g_FT%g_False_Blur.tif',File_Name,Frame_Number,Frame_Thickness)],'tiff');
