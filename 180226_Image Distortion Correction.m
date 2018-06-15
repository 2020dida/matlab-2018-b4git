clear all
%%
Folder_Path='C:\TuanShu\180226_IG of standard target\Single Image Test\';
Row=1224;
Colomn=1024;
IntrinsicMatrix = [700 0 0; 0 700 0; 612 512 1];
radialDistortion = [0 0]; 
cameraParams = cameraParameters('IntrinsicMatrix',IntrinsicMatrix,'RadialDistortion',radialDistortion); 

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

   % File Loading

for p=1:Frame
    Image_Ori = imread([Folder_Path '\' File_list(p).name],'tif');
    Image_Cor = undistortImage(Image_Ori,cameraParams);
    subplot(2,1,1)
    imagesc(Image_Ori);
    axis equal
    axis off
    subplot(1,1,1)
    imagesc(Image_Cor);
    axis equal
    axis off
    disp(p);
end