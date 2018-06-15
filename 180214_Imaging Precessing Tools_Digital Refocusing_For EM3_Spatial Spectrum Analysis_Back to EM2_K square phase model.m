clear all
cd('C:\TuanShu\MATLAB\TSLib');
%%
If_Save_3DXSpectrum=0;
%%
    Lateral_Scanning_Speed=0.96;   %mm/sec
    Frame_Rate=4000;%2000

    Row=1024;
    Colomn=11;%50 
    Stitched_Y_ROI=[];
%% Line-scan Parameter
FF_Length=7;
If_Demodulation=1;          %0 for no demodulation, 1 for frame shifting after demodulation, 2 for frame shifting before demodulation
If_Demosaic=1;
Order=1;                    %Up or Down, 1 for Up scanning, 2 for Down Scanning, ignored for decimation
Overlapping_Mode=2;         %1 for decimation, 2 for averaging
Lateral_SR=0.48;            %micron
Y_Pixel_Phase=1;            %ignored if average, unit: pixel, i.e. 1 for uppest line, 6 for central line (11 pixel total)


%% Some Scanning Info (Mainly for FFT)
Scanning_Speed=140; %micron/sec
Wavelength_SP=0.3;  %micron
Wavelength_LP=1.2;  %micron
%%
root_folder_path='C:\TuanShu\';
last_folder_name='20180208 plastic GP3\';
Number=[3];
%%
Task='Refocusing';    %CCDZ or PZT
If_New_Folder_Structure=1;
%% Basic Operations
If_Reverse=0;
Normalization_ROI=[400 600];
%% Partial Spectrum
if strcmp(Task,'CCDZ')
    If_Speckle=0;
    If_normalization=1;
    If_DF=1;
    If_FF=0;
    If_Bandpass=1;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=1;   %Output the abs (env)
elseif strcmp(Task,'PZT')
    If_Speckle=0;
    If_normalization=0;
    If_DF=1;
    If_FF=0;
    If_Bandpass=0;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=0;   %Output the abs (env)
elseif strcmp(Task,'Speckle')
    If_Speckle=1;
    If_normalization=1;
    If_DF=1;
    If_FF=0;
    If_Bandpass=1;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=2;   %2 for both
elseif strcmp(Task,'LineScan')
    If_Speckle=0;
    If_normalization=1;
    If_DF=0;
    If_FF=0;
    If_Bandpass=0;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=1;   %2 for both
elseif strcmp(Task,'Refocusing')
    If_Speckle=0;
    If_normalization=1;
    If_DF=0;
    If_FF=0;
    If_Bandpass=1;
    BW=0.32;        %Guassian sig
    Offset=0;    %percent
    If_Hilbert=0;   %2 for both
end
%% Image Sharpness Estimation
if strcmp(Task,'CCDZ')
    If_Sharpness_Est=1;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 225];%[100 350];    %micron
elseif strcmp(Task,'PZT')
    If_Sharpness_Est=0;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 350];    %micron
    
elseif strcmp(Task,'Speckle')
    If_Sharpness_Est=0;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 350];    %micron
elseif strcmp(Task,'LineScan')
    If_Sharpness_Est=0;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 350];    %micron
elseif strcmp(Task,'Refocusing')
    If_Sharpness_Est=0;
    If_Lateral_Sharpness_Only=1;
    AVE_postHilbert=4*2*2;
    Axial_Resolution=0.56;  %
    Axial_ROI=[100 350];    %micron
end


%% General STFT
Ave_PreSpectrogram=1;
Time_Window=200/Ave_PreSpectrogram;
nBin=1024;
Spectrogram_X_ROI=[1 1024];%[200 800];
Spectrogram_Y_ROI=[1 Colomn];
%% STFT 1D (PZT Speed Analysis)
if strcmp(Task,'CCDZ')
    If_Refocusing=0;
    If_LineScan=0;
    If_STFT_1D=0;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
    If_Save_Raw=0;

elseif strcmp(Task,'PZT')
    If_Refocusing=0;
    If_LineScan=0;
    If_STFT_1D=1;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.6:0.01:1.2;%0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
    If_Save_Raw=0;

elseif strcmp(Task,'Speckle')
    If_Refocusing=0;
    If_LineScan=0;
    If_STFT_1D=0;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.5:0.01:1.2;%0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
    If_Save_Raw=1;

elseif strcmp(Task,'LineScan')
    If_Refocusing=0;
    If_LineScan=1;
    If_STFT_1D=0;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.5:0.01:1.2;%0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
    If_Save_Raw=0;
elseif strcmp(Task,'Refocusing')
    If_Refocusing=1;
    Refocusing_Coef=1;
    Refocusing_XROI=[1 1024];
    Refocusing_YROI=[1 11];
    Refocusing_ZeroPaddingRatio=1;
    If_LineScan=0;
    If_STFT_1D=0;   % only use the central 1D map for STFT, reducing data szie
    Wavelength_of_Interest_Array=0.5:0.01:1.2;%0.6:0.01:1.2; %micron
    RI=1.406;    %Refractive Index
    Predicted_Central_Wavelength=0.78;
    If_Save_Raw=0;
end



%% N-point Calculation
If_Npoint=0; %Along Z, ~= along Time
N=4;
PreAVE=2;
PostAVE=2;%16/AVE_preDecorr/N;
Auto_Brightness_TH=0.99;
If_Enhace_Deep_Signal=1;

BI=-0.1;
BF=1;

%% Save Raw File

%%
N=4;

% Data format related
Byte_Skip=0;
% Processing related
Column_Binning_Factor=1;
Row_Binning_Factor=1;
Lateral_Ave_Factor=1;   %New Param 170413
Axial_Ave_Factor=2*N; %8*8=64, *N for N-point

%Y_Ave_Factor=8;
PreAVE=2;%[4*4*4];
Product_of_Axial_Decimation_Factor_and_Ave_Factor=PreAVE;    %should be 64 for 0.28/4, use 8 here
PostAVE=2;%16/AVE_preDecorr/N;

Phase_Shift=0;

Number_of_Frame_per_File=2;


% For Last File Number Estimation
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
    file_list=downsample(file_list,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/PreAVE),Phase_Shift);
    Frame_Number_in_File_List=downsample(Frame_Number_in_File_List,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/PreAVE),Phase_Shift);
    %
    File_Reading_Range=[1 11000];
    file_list=file_list(File_Reading_Range(1):File_Reading_Range(2));
    Frame_Number_in_File_List=Frame_Number_in_File_List(File_Reading_Range(1):File_Reading_Range(2));
    Frame=length(file_list);
    
       %% File Loading
    Image_Stack=zeros(Row,Colomn,Frame);
    for p=1:length(file_list)
        fin=fopen([folder_path file_list(p).name]);
        fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(p)-1), 'bof');
        Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint16');
        fclose('all');
        disp(p);
    end
        
    Ave_Array=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
   %% Basic Operations 
   if If_Reverse == 1
       Image_Stack=Image_Stack(:,:,size(Image_Stack,3):-1:1);
   end
   Ave_Array=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
   
   
   if If_normalization ==1  %如果DF是只減min, 那DF就可以在normalization前
        Image_Stack=Image_Stack./repmat(Ave_Array,[size(Image_Stack,1),size(Image_Stack,2),1])*mean(Ave_Array);
    elseif If_normalization ==2
        Image_Stack=Image_Stack./Image_Stack;
    end    
    Med_Image=median(Image_Stack,3);

      
    if If_FF ==1    %即使DF只減min, 先減DF的話FF還是很可能fail, 因為主要的FF都被減掉, 只剩carrier和移動中pattern誤差(後者可能較大)
        Image_Stack=Image_Stack./repmat(Med_Image,[1,1,size(Image_Stack,3)])*mean(Med_Image(:));
    end
    
   Min_Image=min(Image_Stack,[],3);
   
   if If_DF ==1
        Image_Stack=Image_Stack-repmat(Min_Image,[1,1,size(Image_Stack,3)]);
    end
    Ave_Array_After_DF=mean(mean(Image_Stack(Normalization_ROI(1):Normalization_ROI(2),:,:),1),2);
  

  
   
   


    %% FFT to Frequency Domain
    if (If_Bandpass==1)

        Axial_SR_Raw=Scanning_Speed/Frame_Rate; %micron
        Optical_Full_BW=3E8/(Axial_SR_Raw*1E-6);
        Frequency_Resolution=Optical_Full_BW/size(Image_Stack,3);
        Frequency_SP=3E8/(Wavelength_LP*(1E-6)/RI/2);
        Frequency_LP=3E8/(Wavelength_SP*(1E-6)/RI/2);
        
        if If_Refocusing == 0

            Index_SP=round(Frequency_SP/Frequency_Resolution);
            Index_LP=round(Frequency_LP/Frequency_Resolution);
            Image_Stack_FFT=fft(Image_Stack,[],3);
        elseif If_Refocusing == 1

            Index_SP=round(Frequency_SP/Frequency_Resolution*Refocusing_ZeroPaddingRatio);
            Index_LP=round(Frequency_LP/Frequency_Resolution*Refocusing_ZeroPaddingRatio);
            Image_Stack_Refocusing_ROI=Image_Stack(Refocusing_XROI(1):Refocusing_XROI(2),Refocusing_YROI(1):Refocusing_YROI(2),:);
            
            BND_Image=mean(Image_Stack_Refocusing_ROI,3);
            %BND_Array=mean(mean(Image_Stack_Refocusing_ROI,1),2);
            %BND_Array=BND_Array./mean(BND_Array);
            Image_Stack_Refocusing_ROI=Image_Stack_Refocusing_ROI-repmat(BND_Image,[1 1 size(Image_Stack_Refocusing_ROI,3)]);%.*repmat(BND_Array,[size(Image_Stack_Refocusing_ROI,1) size(Image_Stack_Refocusing_ROI,2) 1]);
            Image_Stack_Refocusing_ROI(:,:,(size(Image_Stack,3)+1):(size(Image_Stack,3)*Refocusing_ZeroPaddingRatio))=0;
            Image_Stack_FFT=fft(Image_Stack_Refocusing_ROI,[],3);
        end
        %%
        %imagesc(Image_Stack_Refocusing_ROI(:,:,191));
        plot(squeeze(mean(mean(Image_Stack_Refocusing_ROI,1),2)))
        %%
         subplot(1,1,1)
         plot(squeeze(mean(mean(abs(Image_Stack_FFT),1),2)))
         xlim([Index_SP Index_LP]);
         Image_Stack_FFT(:,:,1:Index_SP)=0;
         Image_Stack_FFT(:,:,Index_LP:end)=0;
    end
    %% If_Bandpass
    if If_Bandpass ==1
        %Mean_Spectrum=squeeze(mean(mean(abs(Image_Stack_FFT),1),2));    %%only used for center finding
        %Mean_Spectrum(Mean_Spectrum<max(Mean_Spectrum)*0.2)=0;          %only >20% peak are considered
        %Peak_Index=sum([1:length(Mean_Spectrum)]'.*Mean_Spectrum)/sum(Mean_Spectrum);
        %Filter_Function=zeros(1,1,length(Mean_Spectrum));
        %Filter_Function(1,1,:)=gaussmf([1:length(Mean_Spectrum)]',[BW*length(Mean_Spectrum) Peak_Index+round(size(Image_Stack_FFT,3)*Offset)]);
        %plot(squeeze(mean(mean(abs(Image_Stack_FFT),1),2)).*squeeze(Filter_Function));
        %Image_Stack_Ori=Image_Stack;
        %Image_Stack_Ori=Image_Stack_Ori(:,:,size(Image_Stack_Ori,3):-1:1);
        if If_Hilbert == 0
            %Image_Stack=2*real(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            %Phase_Stack=angle(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Image_Stack=2*real(ifft(Image_Stack_FFT,[],3));
            Phase_Stack=angle(ifft(Image_Stack_FFT,[],3));

        elseif If_Hilbert == 1
            %Image_Stack=2*abs(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            %Phase_Stack=angle(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Image_Stack=2*abs(ifft(Image_Stack_FFT,[],3));
            Phase_Stack=angle(ifft(Image_Stack_FFT,[],3));
         elseif If_Hilbert == 2
            %Image_Stack=2*abs(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            %Phase_Stack=angle(ifft(Image_Stack_FFT.*repmat(Filter_Function,[size(Image_Stack_FFT,1) size(Image_Stack_FFT,2) 1]),[],3));
            Image_Stack=2*real(ifft(Image_Stack_FFT,[],3));
            Phase_Stack=angle(ifft(Image_Stack_FFT,[],3));     
            Image_Stack_Env=2*abs(ifft(Image_Stack_FFT,[],3));

            
        end

        %%
                %imagesc(squeeze(Image_Stack(:,6,:))')

    end
%% Refocusing
if If_Refocusing == 1
    %%
    Image_Stack_Complex=ifft(Image_Stack_FFT,[],3);

    Image_Stack_FFTX=fft(Image_Stack_Complex,[],1);
      %% Unit Conversion (lps/mm) and Tick Label
    KXTick_Separation=250;    %cycle/mm
    ZTick_Separation=50;    %micron
    X_Spatial_Frequency_Sampling_Resolution=1/(Lateral_SR*size(Image_Stack_FFTX,1)*1E-3);   %%cycle/mm
    Axial_Spatial_Sampling_Resolution=Predicted_Central_Wavelength/2/RI/N/PreAVE*Ave_PreSpectrogram;    %%micron
    Z_Total_Depth=Axial_Spatial_Sampling_Resolution.*size(Image_Stack_FFTX,3);
    X_Spatial_Frequency_Max=X_Spatial_Frequency_Sampling_Resolution*size(Image_Stack_FFTX,1)/2;    %cycle/mm,
    KXTick_Array=(-1*round(X_Spatial_Frequency_Max/KXTick_Separation):round(X_Spatial_Frequency_Max/KXTick_Separation))*KXTick_Separation;
    KXTick_Index_Array=round((KXTick_Array+X_Spatial_Frequency_Max)/X_Spatial_Frequency_Sampling_Resolution);
    for p=1:length(KXTick_Index_Array)
        KXTickLabel_Array{p}=KXTick_Array(p);
    end
    
    ZTick_Array=0:ZTick_Separation:round(Z_Total_Depth/ZTick_Separation)*ZTick_Separation;
    ZTick_Index_Array=round(ZTick_Array/Axial_Spatial_Sampling_Resolution);
    for p=1:length(ZTick_Index_Array)
        ZTickLabel_Array{p}=ZTick_Array(p);
    end 
  %%
    Near_Center_Supression_Ratio=3;
    Image_Stack_FFTX(1:2,:,:)=Image_Stack_FFTX(1:2,:,:)/Near_Center_Supression_Ratio;
    Image_Stack_FFTX(size(Image_Stack_FFTX,1),:,:)=Image_Stack_FFTX(size(Image_Stack_FFTX,1),:,:)/Near_Center_Supression_Ratio;

    subplot(1,1,1)
    imagesc(squeeze(mean(fftshift(abs(Image_Stack_FFTX),1),2))');

    set(gca,'XTick',KXTick_Index_Array)
    set(gca,'XTickLabel',KXTickLabel_Array)

    set(gca,'YTick',ZTick_Index_Array)
    set(gca,'YTickLabel',ZTickLabel_Array)

    xlabel('Spatial Frequency (cycle/mm)')
    ylabel('Depth (\mum)')

    caxis([0 500])
%%
   Phase_Stack_FFTX=angle(Image_Stack_FFTX);
   
   KX_Grid=repmat(reshape([1:size(Image_Stack_FFTX,1)],[size(Image_Stack_FFTX,1) 1 1]),[1 size(Image_Stack_FFTX,2) size(Image_Stack_FFTX,3)]);
   KX_Grid((round(size(Image_Stack_FFTX,1)/2)+1):end,:,:)=KX_Grid(round(size(Image_Stack_FFTX,1)/2):-1:1,:,:);

   Z_Grid=repmat(reshape([1:size(Image_Stack_FFTX,3)],[1 1 size(Image_Stack_FFTX,3)]),[size(Image_Stack_FFTX,1) size(Image_Stack_FFTX,2) 1]);
% Phase Model: exp(i*C*(Z-Z0)*(KX-KX0)^2), 3 const model
    %C=0:1:10; %time 1E-9

    %Phase_Corr_Factor=exp(i*C*(1E-9).*(Z_Grid-Z0).*((KX_Grid-KX0).^2));
    % Try to Find the Opt C for Specific depth
%%
    C=25; %time 1E-9
    Z0_Pixel=0;
    Z0=Z0_Pixel*PostAVE*PreAVE*N;
    KX0=0;
    Z_Test_Pixel=200;
    Z_Test=Z_Test_Pixel*PostAVE*PreAVE*N+1;
    %Normalized_Sharpness_Array=zeros([length(C) 1]);
    Entropy_Array=zeros([length(C) 1]);
    for p=1:length(C)
        Image_Stack_FFTX_Image_Corrected=abs(ifft((Image_Stack_FFTX(:,:,(Z_Test:Z_Test+PostAVE*PreAVE*N)).*exp(i*C(p)*(1E-9).*(Z_Grid(:,:,(Z_Test:Z_Test+PostAVE*PreAVE*N))-Z0).*((KX_Grid(:,:,(Z_Test:Z_Test+PostAVE*PreAVE*N))-KX0).^2))),[],1));
        %Normalized_Sharpness_Array(p)=sum(abs(diff(sum(sum(Image_Stack_FFTX_Image_Corrected,3),2),1,1)))/mean(mean(mean(Image_Stack_FFTX_Image_Corrected)));
        Entropy_Array(p)=entropy(squeeze(mean(Image_Stack_FFTX_Image_Corrected,3)));
        disp(p);
    end
%     imagesc(squeeze(mean(Image_Stack_FFTX_Image_Corrected,2)));
%     caxis([0 5]);
%     plot(C,Entropy_Array/max(Entropy_Array));
    %
    [value index_Best]=max(Entropy_Array);
    C_Best=C(index_Best);
    %
    Image_Stack_FFTX_Corrected=Image_Stack_FFTX.*exp(i*C_Best*(1E-9).*(Z_Grid-Z0).*((KX_Grid-KX0).^2));
    %%
    subplot(2,1,1)
    imagesc(squeeze(mean(fftshift(abs( Image_Stack_FFTX),1),2))');
    title('Spatial Frequency Spectrum (Amplitude) Before Correction')

    set(gca,'XTick',KXTick_Index_Array)
    set(gca,'XTickLabel',KXTickLabel_Array)

    set(gca,'YTick',ZTick_Index_Array)
    set(gca,'YTickLabel',ZTickLabel_Array)

    xlabel('Spatial Frequency (cycle/mm)')
    ylabel('Depth (\mum)')
    
    caxis([0 500])

    subplot(2,1,2)
    imagesc(squeeze(mean(fftshift(abs(Image_Stack_FFTX_Corrected),1),2))');
    title('Spatial Frequency Spectrum (Amplitude) After Correction')

    set(gca,'XTick',KXTick_Index_Array)
    set(gca,'XTickLabel',KXTickLabel_Array)

    set(gca,'YTick',ZTick_Index_Array)
    set(gca,'YTickLabel',ZTickLabel_Array)

    xlabel('Spatial Frequency (cycle/mm)')
    ylabel('Depth (\mum)')
    caxis([0 500])
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    saveas(gcf,[Data_Save_Folder sprintf('%s_Spatial Frequency Spectrum.png',folder_list(QQQ).name)]);
    %%
    Image_Stack_Corrected=abs(ifft(Image_Stack_FFTX_Corrected,[],1));

    Reduced_Image_Corrected=squeeze(TSBinning(mean(Image_Stack_Corrected,2),3,PostAVE*PreAVE*N))'; %2 for PostAVE*PreAVE
%
    Normalized_Sharpness_atZTest=sum(abs(diff(Reduced_Image_Corrected(Z_Test_Pixel,:),1,2)))/mean(Reduced_Image_Corrected(Z_Test_Pixel,:))
%
    % Auto-brightness
    Auto_Brightness_TH=0.995;   %%0.995
    nbin=250;
    Histogram=hist(Reduced_Image_Corrected(:),nbin);
    Total=sum(Histogram);
    Int_Hist = cumsum(Histogram)/Total;

    C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Reduced_Image_Corrected(:));
    C_min=C_max*0.1;            %%0.1
    Normalized_Image=(Reduced_Image_Corrected-C_min)/(C_max-C_min);
    Normalized_Image(Normalized_Image<0)=0;
    Normalized_Image(Normalized_Image>1)=1;
    %
    subplot(1,1,1)
    imagesc(Normalized_Image);
    colormap(gray);
    % ylim([300 400])
    imwrite(Normalized_Image,[Data_Save_Folder sprintf('%s_DAO_C%g_Z%g.tif',folder_list(QQQ).name,C_Best,Z0_Pixel)],'tiff');

        %%
%     R=0.03125;
%     FrameN=2200;
%     %TEMPTEMP=unwrap(angle(fftshift(Image_Stack_FFTX(:,6,FrameN),1)));
%     TEMPTEMP=abs((fftshift(Image_Stack_FFTX(:,6,FrameN),1)));
%     plot(TEMPTEMP);
% 
% %     imagesc(real(fftshift(fftshift(Image_Stack_FFTXY(:,:,FrameN),1),2)));
%      caxis([0 R*max(max(TEMPTEMP))])
     %%
     if If_Save_3DXSpectrum == 1
        Processed_Data_Path=[Data_Save_Folder '\' sprintf('%s_FFTX_ABS.raw',folder_list(QQQ).name)];
        mkdir(Data_Save_Folder);
        fid = fopen(Processed_Data_Path, 'w+');
        fwrite(fid, fftshift(abs(Image_Stack_FFTX),1), 'single');
        %fwrite(fid, Phase_Stack, 'double');
        
        fclose(fid); 
        fclose('all');
     end
end
%% Line Scan
if If_LineScan == 1
    if If_Demodulation == 0
        Image_Stack_Used=Image_Stack;
        Lateral_Scanning_Speed_Pixel=Lateral_Scanning_Speed*1000/Lateral_SR;
        Lateral_Pixel_Spacing_per_Frame=round(Lateral_Scanning_Speed_Pixel/Frame_Rate);
    elseif If_Demodulation == 1
        N=4;
        Temp=0;
        Temp_2=0;
        Length_Used_Pre=floor(size(Image_Stack,3)/N)*N;
        for p=1:N
            Temp=Temp+Image_Stack(:,:,(N-(p-1)):N:(Length_Used_Pre-(p-1)));
            Temp_2=Temp_2+Image_Stack(:,:,(N-(p-1)):N:(Length_Used_Pre-(p-1))).^2;
            disp(p);
        end
        Image_Stack_Used=(N*Temp_2-Temp.^2).^0.5*(2^0.5)/N; 
        Lateral_Scanning_Speed_Pixel=Lateral_Scanning_Speed*1000/Lateral_SR;
        Lateral_Pixel_Spacing_per_Frame=round(Lateral_Scanning_Speed_Pixel/(Frame_Rate/N));  
        
        % 在Z方向demoulation完後, 再拼接
        
        Stitched_Image=zeros(size(Image_Stack_Used,1),size(Image_Stack_Used,3)*Lateral_Pixel_Spacing_per_Frame);

       %
        %Stitched_Image=reshape(Image_Stack,[size(Image_Stack,1) size(Image_Stack,2)*size(Image_Stack,3)]);
        %Stitched_Image=Stitched_Image/max(Stitched_Image(:));
        %if Overlapping_Mode ==1

        for p=1:size(Image_Stack_Used,3)
            if Overlapping_Mode ==1

                Stitched_Image(:,((p-1)*Lateral_Pixel_Spacing_per_Frame+1):(p*Lateral_Pixel_Spacing_per_Frame))=Image_Stack_Used(:,Y_Pixel_Phase:(Y_Pixel_Phase+Lateral_Pixel_Spacing_per_Frame-1),p);
            elseif Overlapping_Mode ==2
                if ((p-1)*Lateral_Pixel_Spacing_per_Frame+size(Image_Stack_Used,2))<=size(Stitched_Image,2)
                    if Order == 1
                        Stitched_Image(:,((p-1)*Lateral_Pixel_Spacing_per_Frame+1):((p-1)*Lateral_Pixel_Spacing_per_Frame+size(Image_Stack_Used,2)))=Stitched_Image(:,((p-1)*Lateral_Pixel_Spacing_per_Frame+1):((p-1)*Lateral_Pixel_Spacing_per_Frame+size(Image_Stack_Used,2)))+Image_Stack_Used(:,end:-1:1,p);
                    elseif Order == 2
                        Stitched_Image(:,((p-1)*Lateral_Pixel_Spacing_per_Frame+1):((p-1)*Lateral_Pixel_Spacing_per_Frame+size(Image_Stack_Used,2)))=Stitched_Image(:,((p-1)*Lateral_Pixel_Spacing_per_Frame+1):((p-1)*Lateral_Pixel_Spacing_per_Frame+size(Image_Stack_Used,2)))+Image_Stack_Used(:,:,p);

                    end
                else
                    continue;
                end
            end
         end
        
        
        
    elseif If_Demodulation == 2
        %% 概念是, 直接沿著Camera Y方向做demodulation, 但必須先在Frame之間加上offset, 令Camera Y方向恰為Z方向
        N=4;
        Temp=0;
        Temp_2=0;
        Image_Stack_Temp=Image_Stack;
        %% 只做Y方向之FF, 且找local之flat field
        %DF_beforeDemodulation=median(Image_Stack_Temp,3);
        FF_3D=smooth3(Image_Stack_Temp,'box',[1 1 FF_Length]);
        %FF_beforeDemodulation=median(Image_Stack_Temp,3);

        Real_XFF=mean(mean(FF_3D,2),3); %因為不想在X方向做FF, 所以再乘回來
        %Image_Stack_Temp=(Image_Stack_Temp)./repmat(FF_beforeDemodulation,[1 1 size(Image_Stack_Temp,3)]);
        Image_Stack_Temp=(Image_Stack_Temp)./FF_3D;
        Image_Stack_Temp=Image_Stack_Temp.*repmat(Real_XFF,[1 size(Image_Stack_Temp,2) size(Image_Stack_Temp,3)]);
                %%

        Length_Used_Pre=floor(size(Image_Stack_Temp,2)/N)*N; %Note, 這裡是在Y方向做N-point
        Lateral_Scanning_Speed_Pixel=Lateral_Scanning_Speed*1000/Lateral_SR;
        Frame_Spacing_per_Lateral_Pixel_PreNPoint=round(Frame_Rate/Lateral_Scanning_Speed_Pixel); 

        for p=1:size(Image_Stack_Temp,2)
            if Order == 1
                Image_Stack_Temp(:,size(Image_Stack_Temp,2)-p+1,:)=circshift(Image_Stack_Temp(:,size(Image_Stack_Temp,2)-p+1,:),p*Frame_Spacing_per_Lateral_Pixel_PreNPoint,3);
            elseif Order == 2
                Image_Stack_Temp(:,p,:)=circshift(Image_Stack_Temp(:,p,:),p*Frame_Spacing_per_Lateral_Pixel_PreNPoint,3);
            end
            disp(p);

        end
        % 在shift完後, N-point前, 再沿要N-point的方向做一次normalization看看
        %D_PreN_Norm_Array=mean(mean(Image_Stack_Temp,1),3);
        %Image_Stack_Temp=Image_Stack_Temp./repmat(D_PreN_Norm_Array,[size(Image_Stack_Temp,1) 1 size(Image_Stack_Temp,3)]);
        for p=1:N
            Temp=Temp+Image_Stack_Temp(:,(N-(p-1)):N:(Length_Used_Pre-(p-1)),:);
            Temp_2=Temp_2+Image_Stack_Temp(:,(N-(p-1)):N:(Length_Used_Pre-(p-1)),:).^2;
            disp(p);
        end
        Image_Stack_Used=(N*Temp_2-Temp.^2).^0.5*(2^0.5)/N; 
        
        Stitched_Image=squeeze(mean(Image_Stack_Used,2));
        imagesc(Stitched_Image);
        % 因為直接在Y方向demodulation了, 所以不用另外拚, 反而因為實際上在時間的offset為subpixel
        % (例如如果掃描速度為0.24 mm/sec, frame rate為2000fps,
        % 則在時間方向取樣密度為0.12micron), 所以要在橫向做平均, 如果剛好取到臨界值就沒這問題
        % 甚至, 在橫向AVE前, 還要先做normalization (之類的)
        
        % 在這裡也做DF
        Stitched_LineImage=median(Stitched_Image,2);
        Stitched_Array=median(Stitched_Image,1);
        %Stitched_Image=Stitched_Image-repmat(Stitched_Array,[size(Stitched_Image,1) 1]);
        Stitched_Image=TSBinning(Stitched_Image,2,Frame_Spacing_per_Lateral_Pixel_PreNPoint);
        
%%
    imagesc(Stitched_Image);
    %%
%%
    end

%     elseif Overlapping_Mode ==2
%         if Order == 1
%             for p=1:size(Image_Stack,2)
%                 Stitched_Image=Stitched_Image+circshift(squeeze(Image_Stack(:,p,:)),-1*Lateral_Pixel_Spacing_per_Frame*(p-1),2);
%             end
%         elseif Order == 2
%             for p=1:size(Image_Stack,2)
%                 Stitched_Image=Stitched_Image+circshift(squeeze(Image_Stack(:,p,:)),1*Lateral_Pixel_Spacing_per_Frame*(p-1),2);
%             end
%         elseif Order == 3
%             for p=1:size(Image_Stack,2)
%                 Stitched_Image=Stitched_Image+circshift(squeeze(Image_Stack(:,p,:)),-1*10*(p-1),2);
%             end
%         end
%     end

    Stitched_Image=Stitched_Image/max(Stitched_Image(:));
    
    % If_Demosaic
    if If_Demosaic == 1
        Stitched_AVE_Array=mean(Stitched_Image,1);
        Stitched_MED_Array=median(Stitched_Image,1);
        Stitched_MED_LineImage=median(Stitched_Image,2);
        Stitched_MED_LineImage_Norm=Stitched_MED_LineImage;

        Stitched_Image=Stitched_Image./repmat(Stitched_MED_Array,[size(Stitched_Image,1) 1])./repmat(Stitched_MED_LineImage_Norm,[1 size(Stitched_Image,2)]);

    end
    
    
    %

    Auto_Brightness_TH=0.995;  %0.995 0.998 for for 0.24 with frame shifting b4 dem
    nbin=250;
    Histogram=hist(Stitched_Image(:),nbin);
    Total=sum(Histogram);
    Int_Hist = cumsum(Histogram)/Total;

    C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Stitched_Image(:));
    C_min=C_max*0.1;%0.1 0.3 for 0.24 with frame shifting b4 dem
    Stitched_Image=(Stitched_Image-C_min)/(C_max-C_min);
    Stitched_Image(Stitched_Image<0)=0;
    Stitched_Image(Stitched_Image>1)=1;
    
    subplot(1,1,1)
    imagesc(Stitched_Image);
    axis equal
    colormap(gray);
    
    %
    if Overlapping_Mode ==1
        imwrite(Stitched_Image,[Data_Save_Folder '\' folder_list(QQQ).name sprintf('_Line_Scan_Decimated_Y%g_DMode%g.png',Y_Pixel_Phase,If_Demodulation)],'png');
    elseif Overlapping_Mode ==2
        imwrite(Stitched_Image,[Data_Save_Folder '\' folder_list(QQQ).name sprintf('_Line_Scan_Averaged_O%g_DMode%g_FFlength%g.png',Order,If_Demodulation,FF_Length)],'png');
    end

end
%% Speckle Analysis
    if If_Speckle == 1
        Reduced_Image=(TSBinning(TSBinning(Image_Stack_Env,3,16),2,11));
        Reduced_Image_Stack=zeros(size(Reduced_Image,1),size(Reduced_Image,2)*11,size(Reduced_Image,3)*16);
        for p=1:16
            Reduced_Image_Stack(:,:,p:16:((size(Reduced_Image,3)-1)*16+p))=repmat(Reduced_Image,[1 11 1]);
            disp(p);
        end
        for p=1:11
            %Fused_Image{p} = imfuse(squeeze(Image_Stack_Env(:,p,:))'*10,squeeze(Image_Stack(:,p,:))'*3,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
            Fused_Image{p}(:,:,2) = squeeze(Reduced_Image_Stack(:,p,:))'*3;
            Fused_Image{p}(:,:,3) = squeeze(Image_Stack_Env(:,p,:))'*3;
            Fused_Image{p}(:,:,1) = squeeze(Image_Stack(:,p,:))'*2; 
                disp(p);
            imwrite(Fused_Image{p},[Data_Save_Folder sprintf('%s_Fused_%g.png',folder_list(QQQ).name,p)],'png');
        end
        imshow(Fused_Image{1});
        %%
        
        Image_to_be_Filtered=squeeze(TSBinning(Image_Stack_Env(:,5,:),3,16))';
        Image_Filtered=AdaptiveMedian(Image_to_be_Filtered,3,3,0.3,2,3,3);
        %
        subplot(1,2,1)
        imagesc(Image_to_be_Filtered)
        caxis([ 0 100])
        xlim([400 600]);
        ylim([0 400]);

        axis equal
        axis off

        subplot(1,2,2)
        imagesc(Image_Filtered)
        caxis([ 0 100])
        xlim([400 600]);
        ylim([0 400]);

        axis equal
        axis off

        %%
            Fused_Image{1}(:,:,2) = squeeze(Reduced_Image_Stack(:,p,:))'*3;
            Fused_Image{1}(:,:,3) = squeeze(Image_Stack_Env(:,p,:))'*3;
            Fused_Image{1}(:,:,1) = squeeze(Image_Stack(:,p,:))'*2 ;
        imagesc(Fused_Image{1});
    end
        %% Image Sharpness Estimation_Est == 1

    if If_Sharpness_Est
        
        %% Post-FFT Axial Ave
        Reduced_Image=squeeze(mean(Image_Stack,2))';
       
        Temp=0;
        Length_Used_Post=floor(size(Reduced_Image,1)/AVE_postHilbert)*AVE_postHilbert;
        for p=1:AVE_postHilbert
           Temp=Temp+Reduced_Image((AVE_postHilbert-(p-1)):AVE_postHilbert:(Length_Used_Post-(p-1)),:);
        end
        Reduced_Image=Temp/AVE_postHilbert;   
        %% Generating Axial Profile for Reduced Image
        Z_Profile=Axial_Resolution*([0:(size(Reduced_Image,1)-1)]);
       Index_Range=[find(Z_Profile>Axial_ROI(1),1,'first') find(Z_Profile>Axial_ROI(2),1,'first')];
        
        
        %% Basic Index        
        [Gx, Gz]=gradient(Reduced_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz=0;
        end
        S=sqrt(Gx.^2+Gz.^2)/mean(mean(Ave_Array))*4096;  %!!!!! 用最原始intensity norm, 有點像NN
        sharpness=mean(mean(S(Index_Range(1):Index_Range(2),:)));
        %% Intensity Profile
        intensity=mean(mean(Reduced_Image));
        I_Profile=mean(Reduced_Image,2);
        %% No normalization at all        
        S_NoNorm=S*mean(mean(Ave_Array))/4096;
        sharpness_NoNorm=mean(mean(S_NoNorm(Index_Range(1):Index_Range(2),:)));
        %% Depth Enhacement

        Linear_Enhace_Coef=[BI:(BF-BI)/(size(Reduced_Image,1)-1):BF]';
        Enhanced_Image=Reduced_Image.*repmat(Linear_Enhace_Coef,[1 size(Reduced_Image,2)]);
        Enhanced_Axial_Profile=mean(Enhanced_Image,2);
        %
        %
        %% Auto-brightness
        Auto_Brightness_TH=0.995;   %%0.995
        nbin=250;
        Histogram=hist(Enhanced_Image(:),nbin);
        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Enhanced_Image(:));
        C_min=C_max*0.1;            %%0.1
        Normalized_Image=(Enhanced_Image-C_min)/(C_max-C_min);
        Normalized_Image(Normalized_Image<0)=0;
        Normalized_Image(Normalized_Image>1)=1;
        %
        subplot(1,1,1)
        imagesc(Normalized_Image);
        colormap(gray);
        %% Processed Index        
        [Gx_Processed, Gz_Processed]=gradient(Normalized_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz_Processed=0;
        end
        S_Processed=sqrt(Gx_Processed.^2+Gz_Processed.^2);
        sharpness_Processed=mean(mean(S_Processed(Index_Range(1):Index_Range(2),:)));
   
        %% Normalized Gradient
        [Gx_Norm_temp, Gz_Norm_temp]=gradient(Reduced_Image);
        if If_Lateral_Sharpness_Only ==1
            Gz_Norm_temp=0;
        end
        Gx_Norm=Gx_Norm_temp./Reduced_Image;
        Gz_Norm=Gz_Norm_temp./Reduced_Image;
        Gx_Norm(Reduced_Image==0)=0;
        Gz_Norm(Reduced_Image==0)=0;
        
        S_Norm=sqrt(Gx_Norm.^2+Gz_Norm.^2);    %在此normalized S下, noise之貢獻被高估
        %imagesc(abs(S_Norm));
        
        sharpness_Norm=mean(mean((S_Norm(Index_Range(1):Index_Range(2),:))));
        S_Profile=mean(S,2);
        S_Processed_Profile=mean(S_Processed,2);
        S_Norm_Profile=mean(S_Norm,2);
        
        subplot(4,1,1)
        plot(Z_Profile,S_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 10])
        xlabel('Depth (micron)')
        ylabel('Sharpness')
        title(sprintf('Sharpness Profile, S=%.04f',sharpness))

        subplot(4,1,2)
        plot(Z_Profile,I_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 30])
        xlabel('Depth (micron)')
        ylabel('Intensity')
        title(sprintf('Intensity Profile, I=%.04f',intensity))

        subplot(4,1,3)
        plot(Z_Profile,S_Processed_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 0.15])
        xlabel('Depth (micron)')
        ylabel('Processed Sharpness')
        title(sprintf('Processed Sharpness Profile, S_p_r_o_c=%.04f',sharpness_Processed))
        
        subplot(4,1,4)
        plot(Z_Profile,S_Norm_Profile,'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([0 0.3])
        xlabel('Depth (micron)')
        ylabel('Normalized Sharpness')
        title(sprintf('Normalized Sharpness Profile, S_n_o_r_m=%.04f',sharpness_Norm))
        saveas(gcf,[Data_Save_Folder sprintf('%s_Sharpness.png',folder_list(QQQ).name)]);

        fprintf('\nSharpness Index: %g',sharpness);
        fprintf('\nProcessed Sharpness Index: %g',sharpness_Processed);
        fprintf('\nPixel-wise Normalized Sharpness Index: %g\n\n',sharpness_Norm);

        imwrite(Normalized_Image,[Data_Save_Folder sprintf('%s_Sharpness=%04g_Processed Sharpness=%04g_Pixel-wise Normalized Sharpness=%04g.png',folder_list(QQQ).name,sharpness,sharpness_Processed,sharpness_Norm)],'png');
        dlmwrite([Data_Save_Folder 'Sharpness.txt'],[QQQ sharpness sharpness_Processed sharpness_Norm sharpness_NoNorm intensity] ,'delimiter','\t','newline','pc','precision', '%.6f','-append');

    end
  
       %% Spectrogram for STFT
    if (If_STFT_1D == 1)
        Image_Virtual_B_Scan=reshape(Image_Stack(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),Spectrogram_Y_ROI(1):Spectrogram_Y_ROI(2),:),(Spectrogram_X_ROI(2)-Spectrogram_X_ROI(1)+1)*(Spectrogram_Y_ROI(2)-Spectrogram_Y_ROI(1)+1),[]);
        Temp=0;
        Length_Used_PreSpectrogram=floor(size(Image_Virtual_B_Scan,2)/Ave_PreSpectrogram)*Ave_PreSpectrogram;
        for p=1:Ave_PreSpectrogram
           Temp=Temp+Image_Virtual_B_Scan(:,(Ave_PreSpectrogram-(p-1)):Ave_PreSpectrogram:(Length_Used_PreSpectrogram-(p-1)));
        end
        Image_Virtual_B_Scan_Aved=Temp/Ave_PreSpectrogram;
        
        Image_Spectrogram=zeros(size(Image_Virtual_B_Scan_Aved,1),nBin/2+1,floor(size(Image_Virtual_B_Scan_Aved,2)/Time_Window));
        for p=1:size(Image_Virtual_B_Scan_Aved,1)
            Image_Spectrogram(p,:,:)=spectrogram(Image_Virtual_B_Scan_Aved(p,:),Time_Window,0,nBin);  %1st 0 for zero overlapping, nBin for DFT bins, it will outout only half spectrum (512+1)
            disp(p);
        end
    end
       
       %% STFT 1D (PZT Speed Analysis)
    if If_STFT_1D
        Ave_Spectrogram=squeeze(mean(abs(Image_Spectrogram(Spectrogram_X_ROI(1):Spectrogram_X_ROI(2),:,:)),1));
        imagesc(abs(Ave_Spectrogram));  %max frequency = 0.5 * sampling frequency
        caxis([0 500]);
        Axial_Spatial_Sampling_Resolution=Predicted_Central_Wavelength/2/RI/N/PreAVE*Ave_PreSpectrogram;                    %micron
        Axial_Spatial_Sampling_Resolution_OPD=Predicted_Central_Wavelength/2/N/PreAVE*Ave_PreSpectrogram;                    %micron
        Axial_Spatial_Sampling_Frequency=1/Axial_Spatial_Sampling_Resolution_OPD;   %Sample/micron-OPD
        Half_Axial_Spatial_Sampling_Frequency=0.5*Axial_Spatial_Sampling_Frequency; %*0.5 FOR DFT
        Predicted_Central_Frequency=1/(Predicted_Central_Wavelength/2); %based on OPD
        Predicted_Central_Bin=nBin/2/Half_Axial_Spatial_Sampling_Frequency*Predicted_Central_Frequency;
        Max_Spectral_Frequency=3E8*Half_Axial_Spatial_Sampling_Frequency/2*1E6;   %/2 for roundtrip; *1E6 for micron to mm, TiSa ~=3.8E14
        Predicted_Spectral_Frequency=3E8/(Predicted_Central_Wavelength)*1E6;
        
        Spectral_Frequency_Array=0:(Max_Spectral_Frequency/(nBin/2)):Max_Spectral_Frequency;
        Spectral_Frequency_of_Interst_Array=3E8./Wavelength_of_Interest_Array*1E6;
        Ave_Spectrogram_Based_on_Wavelength=interp2(repmat(1:size(Ave_Spectrogram,2),[size(Ave_Spectrogram,1) 1]),repmat(Spectral_Frequency_Array',[1 size(Ave_Spectrogram,2)]),Ave_Spectrogram,repmat(1:size(Ave_Spectrogram,2),[length(Spectral_Frequency_of_Interst_Array) 1]),repmat(Spectral_Frequency_of_Interst_Array',[1 size(Ave_Spectrogram,2)]));
        [Max_Value Max_Index]=max(Ave_Spectrogram_Based_on_Wavelength,[],1);
        Signal_Frequency_Array=(1:size(Ave_Spectrogram,1))*(Frame_Rate/Ave_PreSpectrogram/2)/size(Ave_Spectrogram,1);
        %%  17/11/12 More Data Analysis
        Spectrograqm_Contrast_Index=max(Ave_Spectrogram_Based_on_Wavelength,[],1)./min(Ave_Spectrogram_Based_on_Wavelength,[],1);
        plot(Spectrograqm_Contrast_Index);
        %%
        Spectrogram_Norm_LB=100;
        Ave_Spectrogram_Based_on_Wavelength_Norm=Ave_Spectrogram_Based_on_Wavelength./repmat(Max_Value,[size(Ave_Spectrogram_Based_on_Wavelength,1) 1]);
        [Max_Value_2 Max_Index_2]=max(Ave_Spectrogram(Spectrogram_Norm_LB:end,:),[],1);
        Max_Index_2=Max_Index_2+Spectrogram_Norm_LB-1;
        Ave_Spectrogram_Norm=Ave_Spectrogram./repmat(Max_Value_2,[size(Ave_Spectrogram,1) 1]);
        Ave_Spectrogram_Norm(Ave_Spectrogram_Norm>1)=1;
        Bin_FPC=[(1:size(Ave_Spectrogram_Based_on_Wavelength_Norm,2))' (Wavelength_of_Interest_Array(Max_Index)/Predicted_Central_Wavelength*N*PreAVE)'];
        
        dlmwrite([Data_Save_Folder sprintf('%s_Bin_FPC.txt',folder_list(QQQ).name)],Bin_FPC,'delimiter','\t','newline','pc','precision', '%.6f');
        %% FPC
        subplot(2,1,1)
        imagesc(Ave_Spectrogram_Based_on_Wavelength_Norm);
        colormap(parula); %parula
        hold on
        plot(1:length(Max_Index),Max_Index,'ro');
        hold off
        set(gca,'YTick',(1:8:length(Wavelength_of_Interest_Array)),'YTickLabel',round(downsample(Wavelength_of_Interest_Array/Predicted_Central_Wavelength*N*PreAVE,8)*10)/10,'FontSize', 14);
        xlabel('Temporal (Axial) Bins');
        ylabel('Frame per Carrier');
        pbaspect([2 1 1]);
        saveas(gcf,[Data_Save_Folder sprintf('%s_FPC.png',folder_list(QQQ).name)]);
        %% Signal Frequency
        subplot(2,1,2)
        imagesc(Ave_Spectrogram_Norm);
        colormap(parula); %parula
        hold on
        plot(1:length(Max_Index_2),Max_Index_2,'ro');
        hold off
        set(gca,'YTick',(1:40:length(Signal_Frequency_Array)),'YTickLabel',round(downsample(Signal_Frequency_Array,40)),'FontSize', 14);
        xlabel('Temporal (Axial) Bins');
        ylabel('Signal Frequency (Hz)');
        pbaspect([2 1 1]);
        saveas(gcf,[Data_Save_Folder sprintf('%s_Time Frequency.png',folder_list(QQQ).name)]);
    end
    %% N-point
    if If_Npoint == 1
        Temp=0;
        Length_Used_Pre=floor(size(Image_Stack,3)/AVE_preNp)*AVE_preNp;
        for p=1:AVE_preNp
           Temp=Temp+Image_Stack(:,:,(AVE_preNp-(p-1)):AVE_preNp:(Length_Used_Pre-(p-1)));
        end
        Image_Stack=Temp/AVE_preNp;
        
        Temp=0;
        Temp_2=0;
        Length_Used_Pre=floor(size(Image_Stack,3)/N)*N;
        for p=1:N
            Temp=Temp+Image_Stack(:,:,(N-(p-1)):N:(Length_Used_Pre-(p-1)));
            Temp_2=Temp_2+Image_Stack(:,:,(N-(p-1)):N:(Length_Used_Pre-(p-1))).^2;
            disp(p);
        end
        Image_Stack=(N*Temp_2-Temp.^2).^0.5*(2^0.5)/N; 
        
        Temp=0;
        Length_Used_Post=floor(size(Image_Stack,3)/AVE_postNp)*AVE_postNp;
        for p=1:AVE_postNp
           Temp=Temp+Image_Stack(:,:,(AVE_postNp-(p-1)):AVE_postNp:(Length_Used_Post-(p-1)));
        end
        Image_Stack=Temp/AVE_postNp;   
        Image_Stack=Image_Stack(:,:,size(Image_Stack,3):-1:1);

        
%%
        if If_Enhace_Deep_Signal == 1

        BI=-0.1;
        BF=1;
        Linear_Enhace_Coef=reshape([BI:(BF-BI)/(size(Image_Stack,3)-1):BF],1,1,[]);

            Image_Stack=Image_Stack.*repmat(Linear_Enhace_Coef,[size(Image_Stack,1) size(Image_Stack,2) 1]);

        end   
        
        %%
        Npoint_Image=squeeze(mean(Image_Stack,2))';

     
 
        nbin=250;
        Histogram=hist(Npoint_Image(:),nbin);
        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        C_max=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Npoint_Image(:));
        C_min=C_max*0.05;

        Npoint_Image=(Npoint_Image-C_min)/(C_max-C_min);
        Npoint_Image(Npoint_Image>1)=1;
        Npoint_Image(Npoint_Image<0)=0;
        imagesc(Npoint_Image);
        colormap(gray);
        
        %%
        if If_Bandpass == 0
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint.png',folder_list(QQQ).name)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint.tif',folder_list(QQQ).name)],'tiff');

        else
            imwrite(Npoint_Image,[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.png',folder_list(QQQ).name,BW,Offset)],'png');
            imwrite(uint16(Npoint_Image*2^16),[Data_Save_Folder sprintf('%s_Npoint_Filtered_BW%g_Offset%g.tif',folder_list(QQQ).name,BW,Offset)],'tiff');
        end

    end
    
    %% Post N-point Speckle Axial
    
    %% Save Raw File
    if If_Save_Raw == 1
        Processed_Data_Path=[Data_Save_Folder '\' sprintf('%s.raw',folder_list(QQQ).name)];
        mkdir(Data_Save_Folder);
        fid = fopen(Processed_Data_Path, 'w+');
        fwrite(fid, Image_Stack, 'single');
        %fwrite(fid, Phase_Stack, 'double');
        
        fclose(fid); 
        fclose('all');
    end
end



