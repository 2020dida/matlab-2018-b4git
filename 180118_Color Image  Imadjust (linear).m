clear all
%%
Folder_Path='C:\TuanShu\180117_Histogram Analysis\Change R at Contrast 70\';
Center=[500 600];
Half_Width=100;
ROI=[Center(1)-Half_Width Center(1)+Half_Width Center(2)-Half_Width Center(2)+Half_Width];
Normalization='count';  %'count' (default) | 'probability' | 'countdensity' | 'pdf' | 'cumcount' | 'cdf'
%%

Data_Save_Folder=[Folder_Path '\Processed Data\'];
if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end
%%
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
%%
for p=1:Length_file_list
    Image_Temp=imadjustimread([Folder_Path '\' File_list(p).name],'tif');
    Image_Temp(:,:,3)=imadjust(Image_Temp(:,:,3));
    imwrite(hsv2rgb(Image_Temp),[Data_Save_Folder '\' File_list(p).name '_HEQ.tif'],'tiff');
    disp(p);
end
histogram(Image_Temp)