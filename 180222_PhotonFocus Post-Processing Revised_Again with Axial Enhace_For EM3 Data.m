clearvars 
%%
% Old Style: If_Enhace_Deep_Signal=2,
% Ratio_of_DF_Subtraction_before_Enhance=0
% New Style: If_Enhace_Deep_Signal=1,
% Ratio_of_DF_Subtraction_before_Enhance=0.8

Number=[1:4];%9
Trail_Number=35;
System='780';
Camera='MV1';    %MV1 or B0620;
Set='EM3';
If_Reshape_Array=0;
if strcmp(Set,'EM3')
    root_folder_path='C:\TuanShu\20180214-EM3 data\';
    last_folder_name='40x Forearm_Raw';
elseif strcmp(Set,'Viper')
    root_folder_path='C:\TuanShu\181122_with and without normalization test\';
    last_folder_name='180117';
elseif strcmp(Set,'EM2')
    root_folder_path='C:\TuanShu\180104_Pork Test\0.4A_3Hz\0.5mmsec\';
    last_folder_name='5Hz_Different Amp\';
    
end


If_New_Folder_Structure=0;
If_Reverse=1;
If_Auto_Brightness=1;   %0: Do nothing, 1: Histogram Auto-brightness, 2: CLAHE
Auto_Brightness_TH=0.99;    % Param for Histogram Auto-brightness
If_CLAHE=0;
ClipLimit=0.005;             % Param for CLAHE, 因為CLAHE不做normalization, 所以最好還是在Autobrightness後做
if strcmp(Set,'EM3')
    Imaging_Depth_Range=[0 342];
    Imaging_Depth_Range=[0 400];
    X_ROI=[1 1024];
    NN_ROI=[1 1024];%[350 650];
    NN_Y_Range=4;       %這行怪怪的, 好像應該是[3 7] 類似這樣 
elseif strcmp(Set,'Viper')
    Imaging_Depth_Range=[0 342];
    Imaging_Depth_Range=[0 400];
    X_ROI=[1 1024];
    NN_ROI=[1 1024];%[350 650];
    NN_Y_Range=4;       %這行怪怪的, 好像應該是[3 7] 類似這樣 
elseif strcmp(Set,'EM2')
    Imaging_Depth_Range=[0 342];
    Imaging_Depth_Range=[0 400];
    X_ROI=[1 1024];
    NN_ROI=[1 1024];%[350 650];
    NN_Y_Range=4;       %這行怪怪的, 好像應該是[3 7] 類似這樣 
end

If_Enhace_Deep_Signal=0;    %2 for 170725 style
Ratio_of_DF_Subtraction_before_Enhance=0.8;%0.8;    1for low DF case 171012
Exponantial_Coef=0.05;   % dB/micron    0.075
BI=0;%-0.1;
BF=1;

If_FolderNaming=1;
Naxial=4;        %the N along axial
N=4;      %the N for N-point calculation
If_Calculate_NN=0;
If_Search_Glass=0;
If_Save_Raw=0;

If_Save_3D=1;

Assumed_Glass_Interface_Index=1450;  %711;
Axial_Sampling_Resolution=0.56;
If_Norm=1;

If_Axial_Profile=0;
If_Subtract_DF_for_Axial_Profile=1;

% Data format related
Row=1024;
Colomn=1024;%1024;
Y_Offset=0;
% 
if strcmp(Camera,'MV1')
    Number_of_Frame_per_File=1000;  %2
elseif strcmp(Camera,'B0620')
    Number_of_Frame_per_File=1;
end
Byte_Skip=Y_Offset*Row*2;

% Processing related
Column_Binning_Factor=1;
Row_Binning_Factor=1;

Lateral_Ave_Factor=1;   %New Param 170413
Axial_ave_Factor=2;%2;
Product_of_Axial_Decimation_Factor_and_Ave_Factor=2;%round(2);
PreAVE=2;

Phase_Shift=0;

NN_Window_Size=60-30;  %100
Signal_Window_Size=20;
Separate_Size=25;   %25

% For Last File Number Estimation
Bit_per_Pixel=16;
Frame_Size=Row*Colomn*Bit_per_Pixel/8;  %Byte

% Scan Parent Folder
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
Matrix_Record=[];
Axial_Profile_Matrix=[];
Axial_Profile_10Perc_Matrix=[];
Noise_Array=[];
for QQQ=1:length(folder_list)
    if If_New_Folder_Structure == 1
        folder_path=[parent_folder_path '\' folder_list(QQQ).name '\Raw Image\' ];
    else
        folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];
    end
    %cd(folder_path);
    Max_Number_of_Frame=3000000;

    %%
    file_list_ori=dir(folder_path);
    Original_Length_file_list=length(file_list_ori);
    if Original_Length_file_list <3
        continue;
    end
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
        elseif strcmp(ext,'.db')
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

    %%
    file_list=downsample(file_list,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/PreAVE),Phase_Shift);
    Frame_Number_in_File_List=downsample(Frame_Number_in_File_List,floor(Product_of_Axial_Decimation_Factor_and_Ave_Factor/PreAVE),Phase_Shift);
    Frame=length(file_list);
    %%
    Ave_Temp=zeros(Row,Colomn,PreAVE);

    After_Npoint_Frame_Length=floor(Frame/N/PreAVE);
    After_Npoint_Image_Stack=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));

    if If_Calculate_NN ==1
        Intensity_Stack=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,min(Max_Number_of_Frame,After_Npoint_Frame_Length));
    end
    

    X=[1:Frame];
    Frame_Ave=zeros([min(Max_Number_of_Frame,After_Npoint_Frame_Length)*Naxial 1]);
    for p=1:min(Max_Number_of_Frame,After_Npoint_Frame_Length)
        Npoint_Temp=zeros(Row/Row_Binning_Factor,Colomn/Column_Binning_Factor,Naxial);
        for q=1:Naxial
            for r=1:PreAVE
                %%      
                Current_Count=(p-1)*Naxial*PreAVE+(q-1)*PreAVE+r;
                file_path=[folder_path file_list(Current_Count).name];

                fin=fopen(file_path);

                fseek(fin, Byte_Skip+Row*Colomn*2*(Frame_Number_in_File_List(Current_Count)-1), 'bof');
                %QQ=fread(fin,[Row,Colomn],'uint16');   
                Ave_Temp(:,:,r)=fread(fin,[Row,Colomn],'uint16'); %*Frame   不知為何, 看起來就像是要除16
                %%
                %fclose(fin);
                fclose('all');
            end
            Mean=mean(Ave_Temp,3);

            % Row Binning
            Temp=0;
            Reduced_Width=size(Mean,1)/Row_Binning_Factor;
            for tt=1:Row_Binning_Factor
                Temp=Temp+Mean((Row_Binning_Factor-(tt-1)):Row_Binning_Factor:(Row_Binning_Factor*Reduced_Width)-(tt-1),:);
            end
            Temp=Temp/Row_Binning_Factor;

            Reduced_Length=size(Mean,2)/Column_Binning_Factor;
            for tt=1:Column_Binning_Factor
                Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)+Temp(:,(Column_Binning_Factor-(tt-1)):Column_Binning_Factor:(Column_Binning_Factor*Reduced_Length)-(tt-1));
            end
            Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Column_Binning_Factor;


            Frame_Ave((p-1)*Naxial+q)=mean(mean(Npoint_Temp(NN_ROI(1):NN_ROI(2),:,q)));    %Always calculated

            if If_Norm ==1
                Npoint_Temp(:,:,q)=Npoint_Temp(:,:,q)/Frame_Ave((p-1)*Naxial+q);
            end
        end
        After_Npoint_Image_Stack(:,:,p)=((size(Npoint_Temp,3)*sum(Npoint_Temp.^2,3)-sum(Npoint_Temp,3).^2).^0.5)*(2^0.5)/size(Npoint_Temp,3);


        if If_Calculate_NN ==1
            Intensity_Stack(:,:,p)=mean(Npoint_Temp,3);
        end

        disp(p);
    end

    if If_Norm ==1
        After_Npoint_Image_Stack=After_Npoint_Image_Stack*mean(Frame_Ave);
        if If_Calculate_NN ==1
            Intensity_Stack=Intensity_Stack*mean(Frame_Ave);
        end
    end
    
        %% Glass interface detection

    if If_Search_Glass ==1
        [temp_value Glass_Interface_Index]=max(squeeze(mean(mean(After_Npoint_Image_Stack,1),2)));
        if Glass_Interface_Index+(NN_Window_Size+Separate_Size)>size(After_Npoint_Image_Stack,3)
            Glass_Interface_Index=Assumed_Glass_Interface_Index;
        end
    else
        Glass_Interface_Index=Assumed_Glass_Interface_Index;
    end
    
%% NN
    if If_Calculate_NN ==1

        if If_Reverse ==1
            NN_Range=[Glass_Interface_Index+Separate_Size+1 Glass_Interface_Index+Separate_Size+NN_Window_Size];
            Signal_Range=[Glass_Interface_Index-Signal_Window_Size/2+1 Glass_Interface_Index+Signal_Window_Size/2];
        else
            NN_Range=[Glass_Interface_Index-Separate_Size-NN_Window_Size Glass_Interface_Index-Separate_Size-1];
            Signal_Range=[Glass_Interface_Index-Signal_Window_Size/2+1 Glass_Interface_Index+Signal_Window_Size/2];
        end
        NN_Map=std(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2))./Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),0,3);
        IE_Map=max(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2))./Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2)),[],3);
        Signal_Map=max(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2)),[],3);
        Noise_Map=std(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),0,3);
        DF_Map=mean(After_Npoint_Image_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),3);
        Intensity_Map_Noise=mean(Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)),3);
        Intensity_Map_Sig=mean(Intensity_Stack(NN_ROI(1):NN_ROI(2),NN_Y_Range,Signal_Range(1):Signal_Range(2)),3);
        NN=mean(NN_Map(:));
        IE=max(IE_Map(:));
        DF=mean(DF_Map(:));
        Signal=max(Signal_Map(:));
        Noise=mean(Noise_Map(:));
        Intensity=mean(Intensity_Map_Noise(:));
        IE_Sub_Map=((Signal_Map.^2-DF.^2).^0.5)./Intensity_Map_Sig;
        IE_Sub=max(IE_Sub_Map(:));
    else
        DF=0;
        Noise=0;
    end
    
    
    disp(QQQ);
    
%     %% If Subtract DF
%     if If_Subtract_DF == 1
%         After_Npoint_Image_Stack_DF_Sub=real((After_Npoint_Image_Stack.^2-DF.^2).^0.5);
%         After_Npoint_Image_Stack_DF_Sub(isnan(After_Npoint_Image_Stack_DF_Sub))=0;
%     end
    %%
    Maximum_Axial_Frame=round(Imaging_Depth_Range(2)/Axial_Sampling_Resolution);
    Minimum_Axial_Frame=round(Imaging_Depth_Range(1)/Axial_Sampling_Resolution);
    Temp=0;
    Axial_Length_Original=size(After_Npoint_Image_Stack,3);

    Axial_Length_Used=floor(Axial_Length_Original/Axial_ave_Factor)*Axial_ave_Factor;
    Reduced_Length=Axial_Length_Used/Axial_ave_Factor;
    for p=1:Axial_ave_Factor
        Temp=Temp+After_Npoint_Image_Stack(:,:,(Axial_ave_Factor-(p-1)):Axial_ave_Factor:(Axial_ave_Factor*Reduced_Length)-(p-1));

    end
    Reduced_Stack=Temp/Axial_ave_Factor;
    Glass_Interface_Index=round(Glass_Interface_Index/Axial_ave_Factor);
    %Reduced_Image=squeeze(Reduced_Stack(:,size(Reduced_Stack,2)/2,:))';
    %Reduced_Image=squeeze(mean(Reduced_Stack(:,1:Y_Ave_Factor/Column_Binning_Factor,:),2))';
    Reduced_Image=squeeze(mean(Reduced_Stack(:,:,:),2))';

    

    if If_Reverse ==1
        Reduced_Image=Reduced_Image(size(Reduced_Image,1):-1:1,:);
        Glass_Interface_Index=size(Reduced_Image,1)-Glass_Interface_Index;
    end   

        Reduced_Image=Reduced_Image(max(1,Minimum_Axial_Frame):min(size(Reduced_Image,1),Maximum_Axial_Frame),X_ROI(1):X_ROI(2));

   %% 
    if If_Axial_Profile == 1
        Axial_Profile_Folder_Name=[Data_Save_Folder '\Axial Profiles\'];
        if exist(Axial_Profile_Folder_Name)==0
            mkdir(Axial_Profile_Folder_Name);
        end
        
        Z_Profile=Axial_Sampling_Resolution*([0:(size(Reduced_Image,1)-1)]);
        Reduced_Image_ROI_Sort=sort(Reduced_Image(:,NN_ROI(1):NN_ROI(2)),2,'descend');
        if If_Subtract_DF_for_Axial_Profile == 0
            Axial_Profile=real(mean(Reduced_Image(:,NN_ROI(1):NN_ROI(2)),2));
            Axial_Profile_Top10Perc=mean(Reduced_Image_ROI_Sort(:,1:floor(size(Reduced_Image_ROI_Sort,2)*0.1)),2);

        elseif If_Subtract_DF_for_Axial_Profile == 1
            Axial_Profile=real((mean(Reduced_Image(:,NN_ROI(1):NN_ROI(2)),2).^2-DF^2).^0.5);
            Axial_Profile(isnan(Axial_Profile))=0;
            Axial_Profile_Top10Perc=real((mean(Reduced_Image_ROI_Sort(:,1:floor(size(Reduced_Image_ROI_Sort,2)*0.1)),2).^2-DF^2).^0.5);
            Axial_Profile_Top10Perc(isnan(Axial_Profile_Top10Perc))=0;
        end

        plot(Z_Profile,20*log10(Axial_Profile),Z_Profile,20*log10(Axial_Profile_Top10Perc),'LineWidth',2)
        xlim([0 max(Z_Profile)])
        ylim([-5 60])
        xlabel('Depth (micron)')
        ylabel('OCT Signal (dB)')
        if If_Calculate_NN == 1
            hold on
            plot([0 max(Z_Profile)],[20*log10(Noise) 20*log10(Noise)],'r','LineWidth',2);
            hold off
        end
            legend('Average','Top 10% Average','Noise Floor')

        saveas(gcf,[Axial_Profile_Folder_Name '\' sprintf('%s_Axial_Profile_%gAVE.png',folder_list(QQQ).name,PreAVE)]);
    end

    %% Axial Enhacement
    if If_Enhace_Deep_Signal ==1
        Z_Profile=Axial_Sampling_Resolution*([0:(size(Reduced_Image,1)-1)]);
        Assume_Axial_Profile=10.^((-1*Exponantial_Coef*(Z_Profile-Axial_Sampling_Resolution*Glass_Interface_Index))/20);
        Assume_Axial_Profile(Z_Profile<=(Axial_Sampling_Resolution*Glass_Interface_Index))=1;
        Exponential_Enhace_Coef=1./Assume_Axial_Profile;
        Reduced_Image=((Reduced_Image.^2-(DF*Ratio_of_DF_Subtraction_before_Enhance)^2).^0.5).*repmat(Exponential_Enhace_Coef',[1 size(Reduced_Image,2)]);
        
    elseif If_Enhace_Deep_Signal == 2

        Linear_Enhace_Coef=[BI:(BF-BI)/(size(Reduced_Image,1)-1):BF]';

        
        Reduced_Image=real((Reduced_Image.^2-(DF*Ratio_of_DF_Subtraction_before_Enhance)^2).^0.5).*repmat(Linear_Enhace_Coef,[1 size(Reduced_Image,2)]);
        Reduced_Image(isnan(Reduced_Image))=0;
    end    
    
    
%%        
    %clear After_Npoint_Image_Stack
     
    
    if If_Auto_Brightness == 1
        nbin=250;
        Histogram=hist(Reduced_Image(:),nbin);

        Total=sum(Histogram);
        Int_Hist = cumsum(Histogram)/Total;

        plot(Int_Hist);

        C_max_Auto=find(Int_Hist>Auto_Brightness_TH,1,'first')/nbin*max(Reduced_Image(:));
        C_max=C_max_Auto;
        C_min=C_max*0.05;
    else

        C_max=30;
        C_min=1.5;
    end

    Reduced_Image_normalized=(Reduced_Image-C_min)/(C_max-C_min);
    Reduced_Image_normalized(Reduced_Image_normalized<0)=0;
    Reduced_Image_normalized(Reduced_Image_normalized>1)=1;
   % CLAHE
   if If_CLAHE == 1
        %Reduced_Image_normalized = adapthisteq(single(Reduced_Image_normalized),'ClipLimit',ClipLimit,'NumTiles',[36,72],'Distribution','rayleigh','Range','original','alpha',0.7);
       Reduced_Image_normalized = adapthisteq(single(Reduced_Image_normalized),'ClipLimit',ClipLimit,'NumTiles',[12,38],'Distribution','rayleigh');

   end
    
    % Lateral Ave
    Temp=0;
    Lateral_Length_Original=size(Reduced_Image_normalized,2);
    Reduced_Length=Lateral_Length_Original/Lateral_Ave_Factor;
    for p=1:Lateral_Ave_Factor
       Temp=Temp+Reduced_Image_normalized(:,(Lateral_Ave_Factor-(p-1)):Lateral_Ave_Factor:(Lateral_Ave_Factor*Reduced_Length)-(p-1));
    end
    Lateral_Binned_Image=Temp/Lateral_Ave_Factor;
 
    % Image Output

    subplot(1,1,1)
    imagesc(Lateral_Binned_Image);
    axis off
    axis equal
    caxis([0 1]);
    colormap(gray);
    if If_Norm == 1
        imwrite(Lateral_Binned_Image,[Data_Save_Folder folder_list(QQQ).name '_with norm.png'],'png');
    else
        imwrite(Lateral_Binned_Image,[Data_Save_Folder folder_list(QQQ).name '_without norm.png'],'png');
    end
    %%
    subplot(3,1,1)
    plot(Frame_Ave,'r');
    xlabel('Frame Index (After PreAVE)');
    ylabel('Signal (DN)');
    xlim([1 length(Frame_Ave)]);
    subplot(3,1,2)
    plot(Frame_Ave,'r');
    xlabel('Frame Index (After PreAVE)');
    ylabel('Signal (DN)');
    xlim([1 500]);
    subplot(3,1,3)
    plot(Frame_Ave,'.r');
    xlabel('Frame Index (After PreAVE)');
    ylabel('Signal (DN)');
    xlim([1 100]);
    saveas(gcf,[Data_Save_Folder '\' folder_list(QQQ).name '_Power Stability.png']);
    %%
    if (If_Calculate_NN ==1) && length(NN_Y_Range)==1
        if Glass_Interface_Index+(NN_Window_Size+Separate_Size)>size(After_Npoint_Image_Stack,3)
            Test_Image=squeeze(After_Npoint_Image_Stack(X_ROI(1):X_ROI(2),NN_Y_Range,:))';

        else
            Test_Image=squeeze(After_Npoint_Image_Stack(X_ROI(1):X_ROI(2),NN_Y_Range,NN_Range(1):NN_Range(2)))';
        end
        Test_Image=Test_Image/20;
        dlmwrite([Data_Save_Folder 'Inten_NN Record.txt'],[QQQ Intensity NN],'delimiter','\t','newline','pc','-append');
        dlmwrite([Data_Save_Folder 'IE_Sub Record.txt'],[QQQ IE_Sub],'delimiter','\t','newline','pc','-append');
        imwrite(Test_Image,[Data_Save_Folder sprintf('NN_%03d_%s.png',QQQ,System)],'png');
        Matrix_Record=[Matrix_Record;[QQQ,Intensity,NN]];
    end
    
    if If_Save_Raw == 1
        fid = fopen([Data_Save_Folder sprintf('%03d_%s.raw',QQQ,System)], 'w+');
        fwrite(fid, Lateral_Binned_Image', 'single');
        fclose(fid);
    end

    
    if If_Save_3D == 1
        fid = fopen([Data_Save_Folder sprintf('%s_3D.raw',folder_list(QQQ).name)], 'w+');
        fwrite(fid, Reduced_Stack, 'double');
        fclose(fid);
    end
    
    if If_Axial_Profile == 1
        Axial_Profile_Matrix=[Axial_Profile_Matrix; Axial_Profile'];
        Axial_Profile_10Perc_Matrix=[Axial_Profile_10Perc_Matrix; Axial_Profile_Top10Perc'];
        Noise_Array=[Noise_Array Noise];
    end
    
    
end
 

if If_Axial_Profile == 1
    [maxvalue maxindex]=max(Axial_Profile_Matrix,[],2);;
    for p=1:size(Axial_Profile_Matrix,1)
        Axial_Profile_Matrix_Aligned(p,:)=circshift(Axial_Profile_Matrix(p,:),[1, -(round(maxindex(p)-70))]);
        Axial_Profile_10Perc_Matrix_Aligned(p,:)=circshift(Axial_Profile_10Perc_Matrix(p,:),[1, -(round(maxindex(p)-70))]);

    end
    Axial_Profile_Mean=mean(Axial_Profile_Matrix_Aligned,1);
    Axial_Profile_10Perc_Mean=mean(Axial_Profile_10Perc_Matrix_Aligned,1);
    Noise_Mean=mean(Noise_Array);
    Axial_Profile_Mean_STD=std(Axial_Profile_Matrix_Aligned,1,1);
    Axial_Profile_10Perc_Mean_STD=std(Axial_Profile_10Perc_Matrix_Aligned,1,1);
    
    plot(Z_Profile,20*log10(Axial_Profile_Mean),Z_Profile,20*log10(Axial_Profile_10Perc_Mean),'b','LineWidth',2)
    
%     hold on
%     errorbar(downsample(Z_Profile,20),20*log10(downsample(Axial_Profile_Mean,20)),20*log20(downsample(Axial_Profile_Mean_STD,20)),'.m','LineWidth',2)
%     errorbar(downsample(Z_Profile,20),downsample(Axial_Profile_10Perc_Mean,20),downsample(Axial_Profile_10Perc_Mean_STD,20),'.m','LineWidth',2)
%     hold off
    
    xlim([0 max(Z_Profile)])
    ylim([0 50])
    xlabel('Depth (micron)')
    ylabel('Mean OCT Signal (dB)')
    
    if If_Calculate_NN == 1
        hold on
        plot([0 max(Z_Profile)],[20*log10(Noise_Mean) 20*log10(Noise_Mean)],'r','LineWidth',2);
        hold off
    end
        legend(sprintf('Average of %g Measurements',size(Axial_Profile_Matrix_Aligned,1)),sprintf('Average of %g Measurements (Top 10 percent)',size(Axial_Profile_Matrix_Aligned,1)),'Noise Floor')

    fig = gcf;
    fig.PaperPositionMode = 'auto';
    saveas(gcf,[Data_Save_Folder '\' 'Axial_Profile_Mean.png']);
end

%%
fclose('all');
