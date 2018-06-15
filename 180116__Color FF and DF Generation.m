clear all
%%
If_Mask=1;
If_Normalzied_FF=1;
Sigma=25;
Diameter=1200;
Folder_Path='C:\TuanShu\180116_New FF_R1.7_B2.1\DF0_Gain0_b128_c45_it6_R1.7_B2.1_3\';
Row=1224;
Colomn=1024;

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

Frame=length(File_list);

   %% File Loading
Image_R=zeros(Colomn,Row);
Image_G=zeros(Colomn,Row);
Image_B=zeros(Colomn,Row);
for p=1:length(File_list)
    Image_Temp=imread([Folder_Path '\' File_list(p).name],'tif');
    Image_R=Image_R+double(Image_Temp(:,:,1));
    Image_G=Image_G+double(Image_Temp(:,:,2));
    Image_B=Image_B+double(Image_Temp(:,:,3));
    disp(p);
end

if If_Normalzied_FF ==0
    Image_R=Image_R/length(File_list);
    Image_G=Image_G/length(File_list);
    Image_B=Image_B/length(File_list);
elseif If_Normalzied_FF ==1
    R_Max=max(Image_R(:))/length(File_list);
    G_Max=max(Image_G(:))/length(File_list);
    B_Max=max(Image_B(:))/length(File_list);
    Image_R=Image_R/max(Image_R(:));
    Image_G=Image_G/max(Image_G(:));
    Image_B=Image_B/max(Image_B(:));
end
    

if If_Mask == 1 %計算Mask外之數值並以之
    Image_M=(Image_R+Image_G+Image_B)/3;
    Image_M=Image_M./max(Image_M(:));
    X_Image_Array=mean(Image_M,2)';
    Y_Image_Array=mean(Image_M,1);
    X_Grid_1D=1:length(X_Image_Array);
    Y_Grid_1D=1:length(Y_Image_Array);
    X_Grid=repmat((1:length(X_Image_Array))',[1 length(Y_Image_Array)]);
    Y_Grid=repmat(1:length(Y_Image_Array),[length(X_Image_Array) 1]);
    X_Center=sum(X_Image_Array.*X_Grid_1D)/(sum(X_Image_Array));
    Y_Center=sum(Y_Image_Array.*Y_Grid_1D)/(sum(Y_Image_Array));
    Mask=zeros(size(Image_R,1),size(Image_R,2));
    Mask(((X_Grid-X_Center).^2+(Y_Grid-Y_Center).^2)>(Diameter/2)^2)=1;
    Mask=imgaussfilt(Mask,Sigma);
  Image_R=Image_R./(Image_R.^Mask)*R_Max;
  Image_G=Image_G./(Image_G.^Mask)*G_Max;
  Image_B=Image_B./(Image_B.^Mask)*B_Max;
end


Image_Flat_Field=zeros(Colomn,Row,3);
Image_Flat_Field(:,:,1)=Image_R;
Image_Flat_Field(:,:,2)=Image_G;
Image_Flat_Field(:,:,3)=Image_B;

Image_Flat_Field=uint8(Image_Flat_Field);
imagesc(Image_Flat_Field);
mkdir([Folder_Path '\Result\']);
imwrite(Image_Flat_Field,[Folder_Path '\Result\FF.tiff'],'tiff');
imwrite(Image_Flat_Field/2,[Folder_Path '\Result\DF2.tiff'],'tiff');
imwrite(Image_Flat_Field/4,[Folder_Path '\Result\DF4.tiff'],'tiff');
imwrite(Image_Flat_Field/8,[Folder_Path '\Result\DF8.tiff'],'tiff');
%imwrite(Mask,[Folder_Path '\Result\Mask.tiff'],'tiff');

%%

A=[1 2; 1 2];

EX=[2 2; 1 1];

B=A.^EX;