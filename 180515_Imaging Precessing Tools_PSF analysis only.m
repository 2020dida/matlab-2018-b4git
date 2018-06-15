clear all
addpath('C:\TuanShu\MATLAB\TSLib');
%% PSF related Parameters
If_Plot_PSF=1;
PSF_XROI=[512 512];
PSF_YROI=[5 5];
Z_WindowSize=4000;  %frame
Group_Index_Array=[1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3];
%%
Row=1024;
Colomn=11;%50 
Colomn_Read=11;
Number_of_Row_Skipped=0;
Byte_Skip=2*Row*Number_of_Row_Skipped;

If_New_Folder_Structure=1;
If_Bandpass=1;
%%
Frame_Rate=4000;%2000
Stitched_Y_ROI=[];
%% Some Scanning Info (Mainly for FFT)
Scanning_Speed=140/8.4286*8; %micron/sec
Wavelength_SP=0.5;  %micron
Wavelength_LP=1.3;  %micron
%%
root_folder_path='C:\TuanShu\';
last_folder_name='180515_PSF comp\Selected Data\';
Number=[];
Wavelength_of_Interest_Array=0.4:0.01:1.2; %micron
RI=1.406;    %Refractive Index
Predicted_Central_Wavelength=0.78;

%% For Last File Number Estimation
Bit_per_Pixel=16;
Frame_Size=Row*Colomn*Bit_per_Pixel/8;  %Byte


parent_folder_path=[root_folder_path last_folder_name];
folder_list=dir(parent_folder_path);
Original_Lnegth_folder_list=length(folder_list);

for p=1:Original_Lnegth_folder_list
    if folder_list(Original_Lnegth_folder_list+1-p).isdir == 0
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'.') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'..') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'Processed Data') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    elseif strcmp(folder_list(Original_Lnegth_folder_list+1-p).name,'backup') == 1
        folder_list(Original_Lnegth_folder_list+1-p)=[];
    end
end

Data_Save_Folder=[parent_folder_path '\Processed Data\'];

if exist(Data_Save_Folder)==0
    mkdir(Data_Save_Folder);
end

if isempty(Number)==0
    folder_list=folder_list(Number);
end

PSF_Matrix=zeros(Z_WindowSize,length(folder_list)*(PSF_XROI(2)-PSF_XROI(1))*(PSF_YROI(2)-PSF_YROI(1)));
PSF_Env_Matrix=zeros(Z_WindowSize,length(folder_list)*(PSF_XROI(2)-PSF_XROI(1))*(PSF_YROI(2)-PSF_YROI(1)));
Spectrum_Matrix=zeros(1316,length(folder_list));
for QQQ=1:length(folder_list)
    %% File Listing
    if If_New_Folder_Structure ==1
        folder_path=[parent_folder_path '\' folder_list(QQQ).name '\Raw Image\'];
    else
        folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];
    end

    file_list_ori=dir(folder_path);
    Original_Length_file_list=length(file_list_ori);

    for p=1:Original_Length_file_list
        [pathstr,name,ext] = fileparts(file_list_ori(Original_Length_file_list+1-p).name);
        if file_list_ori(Original_Length_file_list+1-p).isdir ~= 0
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'.') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'..') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(ext,'.PNG')
            file_list_ori(Original_Length_file_list+1-p)=[];
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'Processed Data') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        end
    end
    
    %% To Generate the real File List and Frame_Number_in_File_List Based on Data Size Estimation
    file_list=[];
    Frame_Number_in_File_List=[];
    for p=1:length(file_list_ori)
        Frame_Number_Calculated_Based_on_Size=file_list_ori(p).bytes/Frame_Size;
        file_list=[file_list;repmat(file_list_ori(p),[Frame_Number_Calculated_Based_on_Size 1])]; 
        Frame_Number_in_File_List=[Frame_Number_in_File_List;[1:Frame_Number_Calculated_Based_on_Size]']; 
    end

    Frame=length(file_list);
    
       %% File Loading
    Image_Stack=zeros(Row,Colomn_Read,Frame);
    for p=1:length(file_list)
        fin=fopen([folder_path file_list(p).name]);
        fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(p)-1), 'bof');
        Image_Stack(:,:,p)=fread(fin,[Row,Colomn_Read],'uint16');

        fclose('all');
        disp(p);
    end
   
    %% FFT to Frequency Domain
    if (If_Bandpass==1)
        Image_Stack_MeanMap=mean(Image_Stack,3);
        Axial_SR_Raw=Scanning_Speed/Frame_Rate; %micron
        Optical_Full_BW=3E8/(Axial_SR_Raw*1E-6);
        Frequency_Resolution=Optical_Full_BW/size(Image_Stack,3);
        Frequency_SP=3E8/(Wavelength_LP*(1E-6)/RI/2);
        Frequency_LP=3E8/(Wavelength_SP*(1E-6)/RI/2);
        Index_SP=round(Frequency_SP/Frequency_Resolution);
        Index_LP=round(Frequency_LP/Frequency_Resolution);
        Image_Stack_FFT=fft(Image_Stack-repmat(Image_Stack_MeanMap,[1 1 size(Image_Stack,3)]),[],3);
        %Image_Stack_FFT=fft(Image_Stack,[],3);
        %%
         Image_Stack_FFT(:,:,1:Index_SP)=0;
         Image_Stack_FFT(:,:,Index_LP:end)=0;
         %Image_Stack_FFT(:,:,1:(end-Index_LP))=0;
         %Image_Stack_FFT(:,:,(end-Index_SP):end)=0;
         
         subplot(1,1,1)
         plot(squeeze(mean(mean(abs(Image_Stack_FFT),1),2)))
         AVE_Spectrum=squeeze(mean(mean(abs(Image_Stack_FFT),1),2));
         
         xlim([Index_SP Index_LP]);
         Spectrum_Matrix(:,QQQ)=AVE_Spectrum(Index_SP:Index_LP);
         Spectrum_Matrix(:,QQQ)=Spectrum_Matrix(:,QQQ)/max(Spectrum_Matrix(:,QQQ));
    end
    %% If_Bandpass
    if If_Bandpass ==1
            Image_Stack=2*real(ifft(Image_Stack_FFT,[],3));
            Image_Stack_Env=2*abs(ifft(Image_Stack_FFT,[],3));
    end
%% PSF relater
    if If_Plot_PSF == 1
        Image_Stack_ROI=Image_Stack(PSF_XROI(1):PSF_XROI(2),PSF_YROI(1):PSF_YROI(2),:);
        Image_Stack_Env_ROI=Image_Stack_Env(PSF_XROI(1):PSF_XROI(2),PSF_YROI(1):PSF_YROI(2),:);
        [maxvalue Max_Index_Map]=max(Image_Stack_Env_ROI,[],3);
        Image_Stack_ROI_Zwindowed=zeros(size(Image_Stack_ROI,1),size(Image_Stack_ROI,2),Z_WindowSize);
        Image_Stack_Env_ROI_Zwindowed=zeros(size(Image_Stack_ROI,1),size(Image_Stack_ROI,2),Z_WindowSize);
        for p=1:size(Image_Stack_ROI,1)
            for q=1:size(Image_Stack_ROI,2)
                Total_Index=(QQQ-1)*size(Image_Stack_ROI,1)*size(Image_Stack_ROI,2)+(p-1)*size(Image_Stack_ROI,2)+q;
                PSF_Matrix(:,Total_Index)=squeeze(Image_Stack_ROI(p,q,(Max_Index_Map(p,q)-floor(Z_WindowSize/2)):(Max_Index_Map(p,q)+floor(Z_WindowSize/2)-1)));
                PSF_Env_Matrix(:,Total_Index)=squeeze(Image_Stack_Env_ROI(p,q,(Max_Index_Map(p,q)-floor(Z_WindowSize/2)):(Max_Index_Map(p,q)+floor(Z_WindowSize/2)-1)));
            end
        end
    end
   
end
%%
Calibrated_SR=Predicted_Central_Wavelength/RI/2/8.4286;
%SR=Scanning_Speed/Frame_Rate;
SR=Calibrated_SR;
X=[SR:SR:SR*Z_WindowSize]-SR*Z_WindowSize/2;
plot(X,PSF_Matrix);
plot(X,PSF_Env_Matrix);

%% Normalization and FWHM calc
for p=1:size(PSF_Matrix,2)
    PSF_Matrix_Norm(:,p)=(PSF_Matrix(:,p)-min(PSF_Env_Matrix(:,p)))/(max(PSF_Env_Matrix(:,p))-min(PSF_Env_Matrix(:,p)));
    PSF_Env_Matrix_Norm(:,p)=(PSF_Env_Matrix(:,p)-min(PSF_Env_Matrix(:,p)))/(max(PSF_Env_Matrix(:,p))-min(PSF_Env_Matrix(:,p)));
end

plot(X,PSF_Env_Matrix_Norm);

for p=1:size(PSF_Matrix,2)
    Left_Index=find(PSF_Env_Matrix_Norm(:,p)>0.5,1,'first');
    Right_Index=find(PSF_Env_Matrix_Norm(:,p)>0.5,1,'last');
    FWHM_Array(p)=(Right_Index-Left_Index)*Scanning_Speed/Frame_Rate;
end

%% AVE Spectrum and Envelope by Group
Number_of_Group=max(Group_Index_Array);
Group_Count=[10 10 7];
Env_by_Group=zeros(size(PSF_Matrix,1),Number_of_Group);
Spectrum_by_Group=zeros(1316,Number_of_Group);
for p=1:size(PSF_Matrix,2)
    Env_by_Group(:,Group_Index_Array(p))=Env_by_Group(:,Group_Index_Array(p))+PSF_Env_Matrix_Norm(:,p);
    Spectrum_by_Group(:,Group_Index_Array(p))=Spectrum_by_Group(:,Group_Index_Array(p))+Spectrum_Matrix(:,p);
end

for p=1:Number_of_Group
    Env_by_Group(:,p)=(Env_by_Group(:,p)-min(Env_by_Group(:,p)))/(max(Env_by_Group(:,p))-min(Env_by_Group(:,p)));
    Spectrum_by_Group(:,p)=(Spectrum_by_Group(:,p)-min(Spectrum_by_Group(:,p)))/(max(Spectrum_by_Group(:,p))-min(Spectrum_by_Group(:,p)));
    
end
%%
Index=14;
plot(Spectrum_Matrix(:,Index))
%%
plot(X,Env_by_Group)
legend('SMirau','MLOBJ-same illu','MLOBJ-adj illu')