clear all
cd('C:\TuanShu\MATLAB\TSLib');
%% Data-related Parameters
root_folder_path='C:\TuanShu\180402_MSK Contrast Agent Test\';
Set_Number=[2];
Measurement_Number=[2];
%% Image-related Parameters
Glass_TH=2;   %unit: mean of the whole, image
If_FFDF=1;
If_Glass_Identification=1;
If_Background_Identification=1;
DF_Constant=1;
%folder_path=[root_folder_path last_folder_name];

%% Make Folder
Data_Save_Folder=[root_folder_path '\Processed Data\'];

if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end
%% Set Listing
folder_list=dir(root_folder_path);
Original_Length_folder_list=length(folder_list);
for p=1:Original_Length_folder_list
    [pathstr,name,ext] = fileparts(folder_list(Original_Length_folder_list+1-p).name);
    if folder_list(Original_Length_folder_list+1-p).isdir == 0
        folder_list(Original_Length_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Length_folder_list+1-p).name,'.') == 1
        folder_list(Original_Length_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Length_folder_list+1-p).name,'..') == 1
        folder_list(Original_Length_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Length_folder_list+1-p).name,'Processed Data') == 1
        folder_list(Original_Length_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Length_folder_list+1-p).name,'Others') == 1
        folder_list(Original_Length_folder_list+1-p)=[];
    end
end

if isempty(Set_Number)==0
    folder_list=folder_list(Set_Number);
end
%%

for q=1:length(folder_list)
    %% Measurement Listing
    file_list=dir([root_folder_path folder_list(q).name]);
    Original_Length_file_list=length(file_list);
    for p=1:Original_Length_file_list
        [pathstr,name,ext] = fileparts(file_list(Original_Length_file_list+1-p).name);
        if file_list(Original_Length_file_list+1-p).isdir ~= 0
            file_list(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list(Original_Length_file_list+1-p).name,'.') == 1
            file_list(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list(Original_Length_file_list+1-p).name,'..') == 1
            file_list(Original_Length_file_list+1-p)=[];
        end
    end

    if isempty(Measurement_Number)==0
        file_list=file_list(Measurement_Number);
    end


           %% File Loading
    for p=1:length(file_list)
            %% Read Header First
            clear Header
            fin=fopen([root_folder_path folder_list(q).name '\' file_list(p).name]);

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
            fseek(fin,	168	, 'bof');	Header.	UserCommentStringLength	=	fread(fin,	1	,	'long'	);
            fseek(fin,	168+4+Header.	UserCommentStringLength	, 'bof');	Header.	PixelHorizontalRatio	=	fread(fin,	1	,	'single'	);
            fseek(fin,	168+4+Header.	UserCommentStringLength+4	, 'bof');	Header.	PixelVerticalRatio	=	fread(fin,	1	,	'single'	);

            %fclose(fin);
            %% Read Image
            %fin=fopen([folder_path file_list(p).name]);
            fseek(fin, Header.	Fileheaderlength, 'bof');
            %Image_Cell{p}=fread(fin,[Header.FOVColumnSize*Header.PixelColumnSize,Header.FOVRowSize*Header.PixelRowSize],'single');
            Stitched_Image=fread(fin,[Header.FOVColumnSize*Header.PixelColumnSize,Header.FOVRowSize*Header.PixelRowSize],'single');

            fclose('all');
            disp(p);
            %% Glass Identification (try to exclude Glass FOV)
        if If_Glass_Identification == 1
            Small_Image=imresize(Stitched_Image,[Header.FOVColumnSize Header.FOVRowSize]);
            imagesc(Small_Image)
            Mask_Glass=Small_Image>Glass_TH*median(Small_Image(:));
            Mask_Glass_Full=imresize(Mask_Glass,[Header.FOVColumnSize*Header.PixelColumnSize Header.FOVRowSize*Header.PixelRowSize]);
            Not_Glass_Count=Header.FOVColumnSize*Header.FOVRowSize-sum(Mask_Glass(:));
        end
       %%
        
            %%
        if If_FFDF ==1
            %% FF&DF finder
            %DF_Temp_Image_Stack=99999*ones(Header.PixelColumnSize,Header.PixelRowSize,2);
            FFaddDF_Temp_Image=zeros(Header.PixelColumnSize,Header.PixelRowSize);

            for Ci=1:Header.FOVColumnSize
                for Ri=1:Header.FOVRowSize
                    Column_StartIndex=(Ci-1)*Header.PixelColumnSize+1;
                    Column_EndIndex=Ci*Header.PixelColumnSize;

                    Row_StartIndex=(Ri-1)*Header.PixelRowSize+1;
                    Row_EndIndex=Ri*Header.PixelRowSize;
                    %DF_Temp_Image_Stack(:,:,2)=Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex);
                    %DF_Temp_Image_Stack(:,:,1)=min(DF_Temp_Image_Stack,[],3);
                    Temp_Image=Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex);
                    Temp_Image=Temp_Image./max(Temp_Image(:));
                    if If_Glass_Identification == 1
                        FFaddDF_Temp_Image=FFaddDF_Temp_Image+Mask_Glass(Ci,Ri)*Temp_Image;
                    else
                        FFaddDF_Temp_Image=FFaddDF_Temp_Image+Temp_Image;
                    end
                end
                disp(p);
            end
            if If_Glass_Identification == 1
                FFaddDF_Temp_Image=FFaddDF_Temp_Image/Not_Glass_Count;
            else
                FFaddDF_Temp_Image=FFaddDF_Temp_Image/Header.FOVColumnSize/Header.FOVRowSize;
            end

            DF_Image=((FFaddDF_Temp_Image).^0.5)*DF_Constant;
            FF_Image=FFaddDF_Temp_Image-DF_Image;
            FF_Image=FF_Image/max(FF_Image(:));
            %% Correction
            Corrected_Stitched_Image=Stitched_Image;
            for Ci=1:Header.FOVColumnSize
                for Ri=1:Header.FOVRowSize
                    Column_StartIndex=(Ci-1)*Header.PixelColumnSize+1;
                    Column_EndIndex=Ci*Header.PixelColumnSize;

                    Row_StartIndex=(Ri-1)*Header.PixelRowSize+1;
                    Row_EndIndex=Ri*Header.PixelRowSize;
                    %DF_Temp_Image_Stack(:,:,2)=Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex);
                    %DF_Temp_Image_Stack(:,:,1)=min(DF_Temp_Image_Stack,[],3);
                    Corrected_Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex)=(Corrected_Stitched_Image(Column_StartIndex:Column_EndIndex,Row_StartIndex:Row_EndIndex)-DF_Image)./FF_Image;
                end
                disp(Ci);
            end
        end
     %% 
        %%
        ABTH=0.97;
        Cmin=0.4;
        if If_FFDF == 1
            if If_Glass_Identification == 1
                Reduced_Image=TSAutoBrightness(Corrected_Stitched_Image.*(~Mask_Glass_Full),ABTH,Cmin);
            else
            Reduced_Image=TSAutoBrightness(Corrected_Stitched_Image,ABTH,Cmin);
                %Reduced_Image=adapthisteq(Corrected_Stitched_Image/max(Corrected_Stitched_Image(:)),'Distribution','rayleigh','ClipLimit',0.1,'NumTiles',[Header.FOVColumnSize,Header.FOVRowSize]);
            end
        else
            if If_Glass_Identification == 1
                Reduced_Image=TSAutoBrightness(Stitched_Image.*(~Mask_Glass_Full),ABTH,Cmin);
            else
                Reduced_Image=TSAutoBrightness(Stitched_Image,ABTH,Cmin);
            end
            
            %Reduced_Image=adapthisteq(Stitched_Image/max(Stitched_Image(:)),'Distribution','rayleigh','ClipLimit',0.2,'NumTiles',[Header.FOVColumnSize,Header.FOVRowSize]);
            %Reduced_Image=adapthisteq(Stitched_Image,'NumTiles',[8,8]);
            %Reduced_Image=adapthisteq(Stitched_Image/max(Stitched_Image(:)),'ClipLimit',0.1);
            %Reduced_Image=adapthisteq(Stitched_Image/max(Stitched_Image(:)));
            
        end

        %%
        subplot(1,1,1)
        imagesc(Reduced_Image)
        
%         subplot(1,2,2)
%         histogram(Corrected_Stitched_Image)
% 
%         xlim([0 10])
        %%
        colormap(gray)
        %caxis([0.05 0.3])
        imwrite(Reduced_Image,[Data_Save_Folder file_list(p).name '_Thumbnail.png']);
            
    end
    disp(q);
end