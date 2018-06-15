clear all
cd('C:\TuanShu\MATLAB\TSLib');
%%
DF_Constant=1.3;
root_folder_path='C:\TuanShu\';
last_folder_name='180330_B100 data\';
Number=[1];

folder_path=[root_folder_path last_folder_name];

Data_Save_Folder=[root_folder_path '\Processed Data\'];

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
    fseek(fin, Header.	Fileheaderlength, 'bof');
    %Image_Cell{p}=fread(fin,[Header.FOVColumnSize*Header.PixelColumnSize,Header.FOVRowSize*Header.PixelRowSize],'single');
    Stitched_Image=fread(fin,[Header.FOVColumnSize*Header.PixelColumnSize,Header.FOVRowSize*Header.PixelRowSize],'single');
    
    fclose('all');
    disp(p);
end
%%
imagesc(Stitched_Image)
colormap(gray)
caxis([0 100])
%% FF&DF finder
DF_Temp_Image_Stack=99999*ones(Header.PixelColumnSize,Header.PixelRowSize,2);
FFaddDF_Temp_Image=zeros(Header.PixelColumnSize,Header.PixelRowSize);

for p=1:Header.FOVColumnSize
    for q=1:Header.FOVRowSize
        Column_StartIndex=(p-1)*Header.PixelColumnSize+1;
        Column_EndIndex=p*Header.PixelColumnSize;

        Row_StartIndex=(q-1)*Header.PixelRowSize+1;
        Row_EndIndex=q*Header.PixelRowSize;
        %DF_Temp_Image_Stack(:,:,2)=Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex);
        %DF_Temp_Image_Stack(:,:,1)=min(DF_Temp_Image_Stack,[],3);
        FFaddDF_Temp_Image=FFaddDF_Temp_Image+Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex);
    end
    disp(p);
end
FFaddDF_Temp_Image=FFaddDF_Temp_Image/p/q;


DF_Image=((FFaddDF_Temp_Image).^0.5)*DF_Constant;
FF_Image=FFaddDF_Temp_Image-DF_Image;

%% Correction
Corrected_Stitched_Image=Stitched_Image;
for p=1:Header.FOVColumnSize
    for q=1:Header.FOVRowSize
        Column_StartIndex=(p-1)*Header.PixelColumnSize+1;
        Column_EndIndex=p*Header.PixelColumnSize;

        Row_StartIndex=(q-1)*Header.PixelRowSize+1;
        Row_EndIndex=q*Header.PixelRowSize;
        %DF_Temp_Image_Stack(:,:,2)=Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex);
        %DF_Temp_Image_Stack(:,:,1)=min(DF_Temp_Image_Stack,[],3);
        Corrected_Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex)=(Corrected_Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex)-DF_Image)./FF_Image;
;
    end
    disp(p);
end
%%
ABTH=0.995;
Cmin=0.2;
Reduced_Image=TSAutoBrightness(Corrected_Stitched_Image,ABTH,Cmin);
imagesc(Reduced_Image)
colormap(gray)
imwrite(Reduced_Image,[Data_Save_Folder last_folder_name 'Thumbnail.png']);

