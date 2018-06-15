clear all
%%
root_folder_path='C:\TuanShu\';
last_folder_name='180330_B100 data\';
Number=[1];

folder_path=[root_folder_path last_folder_name];

Data_Save_Folder=[folder_path '\Processed Data\'];

if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end



file_list=dir(folder_path);
Original_Length_file_list=length(file_list);
for p=1:Original_Length_file_list
    [pathstr,name,ext] = fileparts(file_list(Original_Length_file_list+1-p).name);
    if file_list(Original_Length_file_list+1-p).isdir ~= 0
        file_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(file_list(Original_Length_file_list+1-p).name,'.') == 1
        file_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(file_list(Original_Length_file_list+1-p).name,'..') == 1
        file_list(Original_Length_file_list+1-p)=[];
    elseif strcmp(file_list(Original_Length_file_list+1-p).name,'Processed Data') == 1
        file_list(Original_Length_file_list+1-p)=[];
    end
end

if isempty(Number)==0
    file_list=file_list(Number);
end


   %% File Loading
for p=1:length(file_list)
    %% Read Header First
    clear Header
    fin=fopen([folder_path file_list(p).name]);

    fseek(fin,	0	, 'bof');	Header.	FileHeaderVersion	=	fread(fin,	1	,	'long'	);
    fseek(fin,	4	, 'bof');	Header.	Fileheaderlength	=	fread(fin,	1	,	'uint64'	);
    fseek(fin,	12	, 'bof');	Header.	MachineSN	=	fread(fin,	1	,	'uint64'	);
    fseek(fin,	20	, 'bof');	Header.	ImageTypeStringLength	=	fread(fin,	1	,	'long'	);
    fseek(fin,	57	, 'bof');	Header.	ImageProviderStringLength	=	fread(fin,	1	,	'long'	);
    fseek(fin,	82	, 'bof');	Header.	ImageDataFormat	=	fread(fin,	1	,	'uint32'	);
    fseek(fin,	86	, 'bof');	Header.	DataTimeStringLength	=	fread(fin,	1	,	'long'	);
    fseek(fin,	104	, 'bof');	Header.	DataBodyLength	=	fread(fin,	1	,	'uint64'	);
    fseek(fin,	112	, 'bof');	Header.	DataType	=	fread(fin,	1	,	'long'	);
    fseek(fin,	116	, 'bof');	Header.	FOVRowSize	=	fread(fin,	1	,	'long'	);
    fseek(fin,	120	, 'bof');	Header.	FOVColumnSize	=	fread(fin,	1	,	'long'	);
    fseek(fin,	124	, 'bof');	Header.	FOVVolumeSize	=	fread(fin,	1	,	'long'	);
    fseek(fin,	128	, 'bof');	Header.	PixelRowSize	=	fread(fin,	1	,	'long'	);
    fseek(fin,	132	, 'bof');	Header.	PixelColumnSize	=	fread(fin,	1	,	'long'	);
    fseek(fin,	136	, 'bof');	Header.	FirstFOVRowIndex	=	fread(fin,	1	,	'long'	);
    fseek(fin,	140	, 'bof');	Header.	FirstFOVColumnIndex	=	fread(fin,	1	,	'long'	);
    fseek(fin,	144	, 'bof');	Header.	Scanspeed	=	fread(fin,	1	,	'long'	);
    fseek(fin,	148	, 'bof');	Header.	Startdepth	=	fread(fin,	1	,	'single'	);
    fseek(fin,	152	, 'bof');	Header.	Enddepth	=	fread(fin,	1	,	'single'	);
    fseek(fin,	156	, 'bof');	Header.	FOVImageThickness	=	fread(fin,	1	,	'short'	);
    fseek(fin,	160	, 'bof');	Header.	ExposureTime	=	fread(fin,	1	,	'uint32'	);
    fseek(fin,	164	, 'bof');	Header.	Bit	=	fread(fin,	1	,	'uint32'	);
    fseek(fin,	168	, 'bof');	Header.	UserCommentStringLenghth	=	fread(fin,	1	,	'long'	);
    fseek(fin,	1196	, 'bof');	Header.	PixelHorizontalRatio	=	fread(fin,	1	,	'short'	);
    fseek(fin,	1200	, 'bof');	Header.	PixelVerticalRatio	=	fread(fin,	1	,	'short'	);

    %fclose(fin);
    %% Read Image
    %fin=fopen([folder_path file_list(p).name]);
    fseek(fin, 1204, 'bof');
    Image_Cell{p}=fread(fin,[Header.FOVColumnSize*Header.PixelColumnSize,Header.FOVRowSize*Header.PixelRowSize],'single');
    fclose('all');
    disp(p);
end
%%
imagesc(Image_Cell{p})
colormap(gray)
caxis([0 100])
%%
  
   if If_Save_Raw3D == 1
       Processed_Data_Path=[Data_Save_Folder '\' sprintf('%s_BeforeDemodulation.raw',folder_list(QQQ).name)];
        mkdir(Data_Save_Folder);
        fid = fopen(Processed_Data_Path, 'w+');
        fwrite(fid, Image_Stack, 'single');
        %fwrite(fid, Phase_Stack, 'double');
        
        fclose(fid); 
        fclose('all');
   end

   
   
Header_Size=10240;
folder_path='I:\iio Data\160721_Stitched Data\';

If_Read_Dark_Frame=0;

If_Read_FF=0;

Data=1;
%Circ_Offset=0;
Test_Number=1;

File_list={'20160721152619_Pork_after C4 1st Z-com_S133_E1137_1X_light increased_degreased 3.bin'}; 
    %'20160721152619_Pork_after C4 1st Z-com_S133_E1137_1X_light increased_degreased 3.bin'; 
    %'20160721153111_Pork_after C4 1st Z-com_S133_E1137_1X_light increased_degreased 3_4 blocks.bin'};

file_path=[folder_path File_list{Test_Number}];


fin=fopen(file_path);
Header=fread(fin,Header_Size,'ulong');

Colomn_per_FOV=648;
Row_per_FOV=488;

Row=Header(2);
Colomn=Header(1);

Number_of_Row_FOV=Row/Row_per_FOV;
Number_of_Colomn_FOV=Colomn/Colomn_per_FOV;


fin=fopen(file_path);
fseek(fin, Header_Size*4, 'bof');
Image=fread(fin,[Row,Colomn],'single');
%Image=circshift(Image,Circ_Offset);
%
imagesc(Image);

colormap(gray);
%%
Min=4.5;
Max=10;    %因目前沒有平均, 所以訊號越強, 範圍就要越大



Image_Norm=(Image-Min)/((Max-Min));
Image_Norm(Image_Norm>1)=1;

Image_Norm(Image_Norm<0)=0;


imagesc(Image_Norm);
axis equal
fclose('all');
xlim([0 Row]);
ylim([0 Colomn]);
colormap(gray);

axis off
fclose('all');
imwrite(Image_Norm,[folder_path,File_list{Test_Number},'_Image_Norm','.png']);

%%
q=8;
r=5;

Current_Image=Image_Norm(((q-1)*Row_per_FOV+1):q*Row_per_FOV,((r-1)*Colomn_per_FOV+1):r*Colomn_per_FOV);

imagesc(Current_Image);

axis equal
fclose('all');
xlim([0 Row_per_FOV]);
ylim([0 Colomn_per_FOV]);
colormap(gray);

axis off
%% try to find the Dark frame & flat-field from multiple wide scan images
if If_Read_Dark_Frame == 0
    TH=3;
    TH_map=ones(Row_per_FOV,Colomn_per_FOV)*TH;
    DF_map=ones(Row_per_FOV,Colomn_per_FOV)*TH*100;


    for p=1:length(File_list)

        file_path=[folder_path File_list{p}];


        fin=fopen(file_path);
        Header=fread(fin,Header_Size,'ulong');

        Colomn_per_FOV=648;
        Row_per_FOV=488;

        Row=Header(2);
        Colomn=Header(1);

        Number_of_Row_FOV=Row/Row_per_FOV;
        Number_of_Colomn_FOV=Colomn/Colomn_per_FOV;


        fin=fopen(file_path);
        fseek(fin, Header_Size*4, 'bof');
        Image=fread(fin,[Row,Colomn],'single');

        for q=3:Number_of_Row_FOV
            for r=1:Number_of_Colomn_FOV
                Current_Image=Image(((q-1)*Row_per_FOV+1):q*Row_per_FOV,((r-1)*Colomn_per_FOV+1):r*Colomn_per_FOV);
                Current_Image(Current_Image<TH)=TH*100;
                DF_map=max(TH_map,min(DF_map,Current_Image));
            end
        end

        disp(p);

    end

    imagesc(DF_map);

    Times=4;
    Search_Window_Half_Size=10;
    Search_Window_Weight=fspecial('gaussian', Search_Window_Half_Size*2+1,Search_Window_Half_Size/2);    %gaussian mask
    Search_Window_Weight=Search_Window_Weight/max(Search_Window_Weight(:)); %因為是要取max而非平均, 所以filter的最大值應該是1
    Dark_frame_Temp=DF_map;
    Dark_frame_Out=DF_map;
    for w=1:Times
        for p=1:size(Dark_frame_Temp,1)
            for q=1:size(Dark_frame_Temp,2)
                x_min=max((p-Search_Window_Half_Size+1),1);
                x_max=min((p+Search_Window_Half_Size),size(Dark_frame_Temp,1));
                y_min=max((q-Search_Window_Half_Size+1),1);
                y_max=min((q+Search_Window_Half_Size),size(Dark_frame_Temp,2));
                Window_Temp=Search_Window_Weight((x_min-p+Search_Window_Half_Size+1):(x_max-(p+Search_Window_Half_Size)+Search_Window_Half_Size*2+1),(y_min-q+Search_Window_Half_Size+1):(y_max-(q+Search_Window_Half_Size)+Search_Window_Half_Size*2+1)).*Dark_frame_Temp(x_min:x_max,y_min:y_max);
                Dark_frame_Out(p,q)=max(Window_Temp(:));
            end
        end
        Dark_frame_Temp=Dark_frame_Out;
        disp(w);
    end

    subplot(1,2,1)

    imagesc(DF_map);
        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(DF_map,2)]);
        ylim([0 size(DF_map,1)]);

    subplot(1,2,2)

    imagesc(Dark_frame_Out);
        colormap('gray');
        %caxis([0 1]);
        axis equal
        xlim([0 size(Dark_frame_Out,2)]);
        ylim([0 size(Dark_frame_Out,1)]);
    Dark_frame=Dark_frame_Out;


else
    Dark_frame=dlmread([folder_path,'Dark_frame.txt']);
end

%% Search for 2 value
% 
% DF_map_New=DF_map;
% 
% DF_map_New(:,1:Colomn_per_FOV/2)=1.2*ones(Row_per_FOV,Colomn_per_FOV/2);
% DF_map_New(:,(Colomn_per_FOV/2+1):end)=1.05*ones(Row_per_FOV,Colomn_per_FOV/2);
% % 
% % 

if If_Read_FF == 0
    %%
    Total_Image=zeros(Row_per_FOV,Colomn_per_FOV);
    Total_Number_of_Frame=0;
    for p=1:length(File_list)

        file_path=[folder_path File_list{p}];


        fin=fopen(file_path);
        Header=fread(fin,Header_Size,'ulong');

        Colomn_per_FOV=648;
        Row_per_FOV=488;

        Row=Header(2);
        Colomn=Header(1);

        Number_of_Row_FOV=Row/Row_per_FOV;
        Number_of_Colomn_FOV=Colomn/Colomn_per_FOV;


        fin=fopen(file_path);
        fseek(fin, Header_Size*4, 'bof');
        Image=fread(fin,[Row,Colomn],'single');

        for q=1:Number_of_Row_FOV
            for r=1:Number_of_Colomn_FOV
                Current_Image=Image(((q-1)*Row_per_FOV+1):q*Row_per_FOV,((r-1)*Colomn_per_FOV+1):r*Colomn_per_FOV)-Dark_frame;
                Total_Image=Total_Image+Current_Image;
            end
        end
        Total_Number_of_Frame=Total_Number_of_Frame+Number_of_Row_FOV*Number_of_Colomn_FOV;
        disp(p);

    end
    Averaged_Image=Total_Image/Total_Number_of_Frame;
    FF_map=1./(Averaged_Image/mean(Averaged_Image(:)));
    imagesc(Averaged_Image);
    axis equal
    fclose('all');
    xlim([0 Row_per_FOV]);
    ylim([0 Colomn_per_FOV]);
    colormap(gray);

    axis off

    imagesc(FF_map);

    axis equal
    fclose('all');
    xlim([0 Row_per_FOV]);
    ylim([0 Colomn_per_FOV]);
    colormap(gray);

    axis off
else
    FF_map=dlmread([folder_path,'FF.txt']);
end
%%
Test_Number=1;
Factor=1.2;
file_path=[folder_path File_list{Test_Number}];


fin=fopen(file_path);
Header=fread(fin,Header_Size,'ulong');

Colomn_per_FOV=648;
Row_per_FOV=488;

Row=Header(2);
Colomn=Header(1);

Number_of_Row_FOV=Row/Row_per_FOV;
Number_of_Colomn_FOV=Colomn/Colomn_per_FOV;


DF_map_large=repmat(Dark_frame,[Number_of_Row_FOV Number_of_Colomn_FOV]);
FF_map_large=repmat(FF_map,[Number_of_Row_FOV Number_of_Colomn_FOV]);

fin=fopen(file_path);
fseek(fin, Header_Size*4, 'bof');
Image=fread(fin,[Row,Colomn],'single');
%Image=circshift(Image,Circ_Offset);
%
imagesc(Image);

Image_Correlcted=(Image-DF_map_large*Factor).*FF_map_large;


imagesc(Image_Correlcted);
axis equal
fclose('all');
xlim([0 Row]);
ylim([0 Colomn]);
colormap(gray);

axis off
%%

Min=0;
Max=50;%4; 


Image_Correlcted_Norm=(Image_Correlcted-Min)/(Max-Min);

Image_Correlcted_Norm(Image_Correlcted_Norm>1)=1;

Image_Correlcted_Norm(Image_Correlcted_Norm<0)=0;


imagesc(Image_Correlcted_Norm);
axis equal
fclose('all');
xlim([0 Row]);
ylim([0 Colomn]);
colormap(gray);

axis off

imwrite(Image_Correlcted_Norm,[folder_path,File_list{Test_Number},'_stiched_image','.png']);
%dlmwrite([folder_path,'Dark_frame.txt'],Dark_frame,'delimiter','\t','newline','pc','precision', '%.6f');
%dlmwrite([folder_path,'FF.txt'],FF_map,'delimiter','\t','newline','pc','precision', '%.6f');
