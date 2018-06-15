clearvars 
%%
cd('C:\TuanShu\MATLAB\TSLib');
%%
Max_Number_of_Frame=200;
X_ROI=[1201 1300];
Y_ROI=[1001 1100];
Number=[];%9
Set='Basler';
Green_Handling='Decimate2'; % Decimate2 or Decimate1 or Average or Add
%Color='G';
%%
if (mod(X_ROI(2)-X_ROI(1)+1,2) == 1) || (mod(Y_ROI(2)-Y_ROI(1)+1,2) == 1)
    error('ROI sizes not Even');
end;
%% Generating Bayer Mask
X_Offset=1;
Y_Offset=1;
% Mask=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1);
% X_Grid=repmat([1:(X_ROI(2)-X_ROI(1)+1)]',[1 Y_ROI(2)-Y_ROI(1)+1])+X_Offset;
% Y_Grid=repmat([1:(Y_ROI(2)-Y_ROI(1)+1)],[X_ROI(2)-X_ROI(1)+1 1])+X_Offset;
% R_Mask=(mod(X_Grid,2)==0) & (mod(Y_Grid,2)==0);
% B_Mask=(mod(X_Grid+1,2)==0) & (mod(Y_Grid+1,2)==0);
% G1_Mask=(mod(X_Grid+1,2)==0) & (mod(Y_Grid,2)==0);
% G2_Mask=(mod(X_Grid,2)==0) & (mod(Y_Grid+1,2)==0);
% 
% 
% if strcmp(Green_Handling,'Decimate2') == 1
% 	G_Mask=G1_Mask;
% elseif strcmp(Green_Handling,'Decimate1') == 1
% 	G_Mask=G2_Mask;
% elseif strcmp(Green_Handling,'Average') == 1
% 	G_Mask=G2_Mask + G1_Mask;
% 
% elseif strcmp(Green_Handling,'Add') == 1
% 	G_Mask=G2_Mask + G1_Mask;
% 
% end
% imagesc(G_Mask);
% 
% R_Mask_1D=R_Mask(:);
% G_Mask_1D=G_Mask(:);
% B_Mask_1D=B_Mask(:);
% 
% Mask_1D=double(Mask_1D);
% Mask_1D(Mask_1D==0)=NaN;
%% PTC related param
Noise_Removal_Method=1; %1: do nothing, 2: normalization, 3: filtering, 4: 4-point

% Feature for Filtering mode
Filtering_Ratio=0.02;       %�N��O, �YArray���׬�A, �h�N�~��(�C�W)��A*Filtering_Ratio��pixels�]��0
% Feature for 4-point mode
Averaging_Factor=1;
N=4;          %for N-point
Unbias_Estimator=0.9213;

if strcmp(Set,'Basler')
    root_folder_path='C:\TuanShu\180226_no demosaic (comp elast and oil)\';
    last_folder_name='Oil_R3';
    Row=2448;
    Colomn=2048;
    IfColor=1;
end
Bit_per_Pixel=24;
Frame_Size=Row*Colomn*Bit_per_Pixel/8;
   %Header=fread(fin,Header_Size,'ulong');

%fseek(fin, 256*4, 'bof');


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
Small_X=(X_ROI(2)-X_ROI(1)+1)/2;
Small_Y=(Y_ROI(2)-Y_ROI(1)+1)/2;
Mean_Array_Map_R=zeros(Small_X,Small_Y,length(folder_list));
Mean_Array_Map_G=zeros(Small_X,Small_Y,length(folder_list));
Mean_Array_Map_B=zeros(Small_X,Small_Y,length(folder_list));
STD_Array_Map_R=zeros(Small_X,Small_Y,length(folder_list));
STD_Array_Map_G=zeros(Small_X,Small_Y,length(folder_list));
STD_Array_Map_B=zeros(Small_X,Small_Y,length(folder_list));
for QQQ=1:length(folder_list)

    folder_path=[parent_folder_path '\' folder_list(QQQ).name '\'];

    %cd(folder_path);

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
        elseif strcmp(file_list_ori(Original_Length_file_list+1-p).name,'Processed Data') == 1
            file_list_ori(Original_Length_file_list+1-p)=[];
        end
    end
    %% To Generate the real File List and Frame_Number_in_File_List Based on Data Size Estimation
    file_list=[];
    Frame_Number_in_File_List=[];
    for p=1:length(file_list_ori)
        Frame_Number_Calculated_Based_on_Size=(file_list_ori(p).bytes)/Frame_Size;
        if Frame_Number_Calculated_Based_on_Size>1.1
            file_list=[file_list;repmat(file_list_ori(p),[Frame_Number_Calculated_Based_on_Size 1])]; 
        else
            file_list=[file_list file_list_ori(p)];
        end
        Frame_Number_in_File_List=[Frame_Number_in_File_List;[1:Frame_Number_Calculated_Based_on_Size]']; 
    end

    %%
    Frame=min(length(file_list),Max_Number_of_Frame);
    %%
    X=[1:Frame];
    %Image_Stack=zeros(X_ROI(2)-X_ROI(1)+1,Y_ROI(2)-Y_ROI(1)+1,Frame);
    Image_Stack_R=zeros(Small_X,Small_Y,Frame);
    Image_Stack_G=zeros(Small_X,Small_Y,Frame);
    Image_Stack_B=zeros(Small_X,Small_Y,Frame);
    
    for p=1:Frame
        file_path=[folder_path file_list(p).name];
        %QQ=fread(fin,[Row,Colomn],'uint16');
        if  IfColor == 0
            fin=fopen(file_path);
            fseek(fin, HeaderSize+Frame_Size*(Frame_Number_in_File_List(p)-1), 'bof');
            Image_Stack(:,:,p)=fread(fin,[Row,Colomn],'uint8');
            
        elseif IfColor == 1
            Image_Temp=double(imread(file_path,'tiff'));
            
            %R_Mask=(mod(X_Grid,2)==0) & (mod(Y_Grid,2)==0);
            %B_Mask=(mod(X_Grid+1,2)==0) & (mod(Y_Grid+1,2)==0);
            %G1_Mask=(mod(X_Grid+1,2)==0) & (mod(Y_Grid,2)==0);
            %G2_Mask=(mod(X_Grid,2)==0) & (mod(Y_Grid+1,2)==0);
            Image_Stack_R(:,:,p)=Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset,(1:2:(2*Small_X))+X_ROI(1)+X_Offset)';
            if strcmp(Green_Handling,'Decimate2') == 1
                Image_Stack_G(:,:,p)=Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset+1,(1:2:(2*Small_X))+X_ROI(1)+X_Offset)';
            elseif strcmp(Green_Handling,'Decimate1') == 1
                Image_Stack_G(:,:,p)=Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset,(1:2:(2*Small_X))+X_ROI(1)+X_Offset+1)';
            elseif strcmp(Green_Handling,'Average') == 1
                Image_Stack_G(:,:,p)=(Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset+1,(1:2:(2*Small_X))+X_ROI(1)+X_Offset)'+Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset,(1:2:(2*Small_X))+X_ROI(1)+X_Offset+1)')/2;
            elseif strcmp(Green_Handling,'Add') == 1
                Image_Stack_G(:,:,p)=(Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset+1,(1:2:(2*Small_X))+X_ROI(1)+X_Offset)'+Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset,(1:2:(2*Small_X))+X_ROI(1)+X_Offset+1)');

            end
            Image_Stack_B(:,:,p)=Image_Temp((1:2:(2*Small_Y))+Y_ROI(1)+Y_Offset+1,(1:2:(2*Small_X))+X_ROI(1)+X_Offset+1)';
            
        else
            continue;
        end
        fclose('all');
            %%
        disp(p);
    end
%     %%
%     subplot(1,1,1)
%     imagesc(Image_Stack_R(:,:,1));
    %% PTC related Calcuation (for R, G, B, repectivily)
    for p=1:3   %1 for R, 2 for G, 3 for B
        if p ==1
            Image_Stack=Image_Stack_R;
        elseif p ==2
            Image_Stack=Image_Stack_G;
        elseif p ==3
            Image_Stack=Image_Stack_B;
        end
        Number_of_Filter_Pixel_one_side=size(Image_Stack,3)*Filtering_Ratio/2;   %�`�@�����k���side

        Mean_Array_for_norm=mean(mean(Image_Stack,2),1); % �T�w�Υ��i�v����, ���n���a, �u�I: ���|�]��sample�I�ƤӤ֦Ӥ���, ���I: �Y�v�����������M�h����
        Mean_Stack_for_norm=repmat(Mean_Array_for_norm,[size(Image_Stack,1) size(Image_Stack,2) 1]);

        Averaged_Length=floor(size(Image_Stack,3)/Averaging_Factor);
        Temp=0;
        for q=1:Averaging_Factor
           Temp=Temp+Image_Stack(:,:,(Averaging_Factor-(q-1)):Averaging_Factor:(Averaging_Factor*Averaged_Length)-(q-1));
        end
        Image_Stack_Ave=Temp/Averaging_Factor;
        if Noise_Removal_Method == 1
            Mean_Array_Map_Temp=mean(Image_Stack(:,:,:),3);
            STD_Array_Map_Temp=std(Image_Stack(:,:,:),0,3);
        elseif Noise_Removal_Method == 2
            All_Mean=mean(Mean_Array_for_norm);
            Image_Stack_norm=Image_Stack./Mean_Stack_for_norm*All_Mean;
            Mean_Array_Map_Temp=mean(Image_Stack_norm(:,:,:),3);
            STD_Array_Map_Temp=std(Image_Stack_norm(:,:,:),0,3);    
        elseif Noise_Removal_Method == 3
            Image_Stack_FFT=fft((Image_Stack-Mean_Stack_for_norm),[],3);
            Image_Stack_FFT_Filter=Image_Stack_FFT;
            Image_Stack_FFT_Filter(:,:,1:Number_of_Filter_Pixel_one_side)=0;
            Image_Stack_FFT_Filter(:,:,(size(Image_Stack_FFT_Filter,1)-Number_of_Filter_Pixel_one_side+1):size(Image_Stack_FFT_Filter,1))=0;
            Image_Stack_Filter=real(ifft(Image_Stack_FFT_Filter,[],3))+Mean_Stack_for_norm;

            Mean_Array_Map_Temp=mean(Image_Stack_Filter(:,:,:),3);
            STD_Array_Map_Temp=std(Image_Stack_Filter(:,:,:),0,3);        
        elseif Noise_Removal_Method == 4

            Image_Stack_Npoint=zeros(size(Image_Stack_Ave,1),size(Image_Stack_Ave,2),floor(size(Image_Stack_Ave,3)/N));
            for q=1:size(Image_Stack_Npoint,3)
                Temp=Image_Stack_Ave(:,:,((q-1)*N+1):(q*N));
                Image_Stack_Npoint(:,:,q)=((N*sum(Temp.^2,3)-sum(Temp,3).^2).^0.5)*(2^0.5)/N;
            end
            Mean_Array_Map_Temp=mean(Image_Stack(:,:,:),3);
            STD_Array_Map_Temp=mean(Image_Stack_Npoint(:,:,3),3)*Unbias_Estimator;        
        end
        
        if p ==1
            Mean_Array_Map_R(:,:,QQQ)=Mean_Array_Map_Temp;
            STD_Array_Map_R(:,:,QQQ)=STD_Array_Map_Temp;
        elseif p ==2
            Mean_Array_Map_G(:,:,QQQ)=Mean_Array_Map_Temp;
            STD_Array_Map_G(:,:,QQQ)=STD_Array_Map_Temp;
        elseif p ==3
            Mean_Array_Map_B(:,:,QQQ)=Mean_Array_Map_Temp;
            STD_Array_Map_B(:,:,QQQ)=STD_Array_Map_Temp;
        end
        disp(p);
    end
end
 %%
 imagesc(STD_Array_Map_B(:,:,1));

%%
fclose('all');

%%
Size=25;
Downsampling_Ratio=100;
%It_Array=[1 2 3 4 6 8 10 12 14 16 18 20];
It_Array=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
%It_Array=[1 1.5 2 2.5 3 3.5 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

%It_Array=[5 10 20 40 50 60 70 80 100 150 200 250 300];
%It_Array=[1 2 3 4 5 6 8 10 12 15 20 25 30 35];
%It_Array=[5 10 20 40 50 60 70 80 100 150 200 250 300];
% %%
% NN=6;
% 
% plot(real(Data_Image_FFT(:,NN)));
% 

%

for p=1:3
    if p ==1
        Mean_Array_Map=Mean_Array_Map_R;
        STD_Array_Map=STD_Array_Map_R;
    elseif p ==2
        Mean_Array_Map=Mean_Array_Map_G;
        STD_Array_Map=STD_Array_Map_G;
    elseif p ==3
        Mean_Array_Map=Mean_Array_Map_B;
        STD_Array_Map=STD_Array_Map_B;
    end
    Mean_Array_Map_1D=Mean_Array_Map(:);
    Mean_Array=squeeze(mean(mean(Mean_Array_Map,1),2));
    STD_Array_Map_1D=STD_Array_Map(:);
    STD_Array=squeeze(mean(mean(STD_Array_Map,1),2));
    Mean_Array_Map_1D(isnan(Mean_Array_Map_1D))=[];
    STD_Array_Map_1D(isnan(STD_Array_Map_1D))=[];
    Mean_Array_Map_1D=downsample(Mean_Array_Map_1D,Downsampling_Ratio);
    STD_Array_Map_1D=downsample(STD_Array_Map_1D,Downsampling_Ratio);

    VAR_Array_Map_1D=STD_Array_Map_1D.^2;
    if p==1
        scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'r','filled');
    elseif p==2
        scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'g','filled');
    elseif p==3
        scatter(Mean_Array_Map_1D,VAR_Array_Map_1D,Size,'b','filled');
    end
    xlabel('Signal (DN)','fontsize',15);
    ylabel('Variance (DN^2)','fontsize',15);
    set(gca,'fontsize',15)
    %ylim([0 200]);
    xlim([0 256]);
    if p ==1
        saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g_R',Noise_Removal_Method) '.png']);
    elseif p==2
        saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g_G',Noise_Removal_Method) '.png']);
    elseif p==3
        saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g_B',Noise_Removal_Method) '.png']);
    end
    if p==1
        plot(It_Array,Mean_Array,'-or','LineWidth',2,'MarkerFaceColor','r');
    elseif p==2
        plot(It_Array,Mean_Array,'-og','LineWidth',2,'MarkerFaceColor','g');
    elseif p==3
        plot(It_Array,Mean_Array,'-ob','LineWidth',2,'MarkerFaceColor','b');

    end
    xlabel('Exposure Time (ms)','fontsize',15);
    ylabel('Signal (DN)','fontsize',15);
    set(gca,'fontsize',15)
    ylim([0 256]);
    %xlim([0 256]);
    Slope=(Mean_Array(5)-Mean_Array(1))/(It_Array(5)-It_Array(1));
    if p ==1
        saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g_Linearity_R_Slope%g',Noise_Removal_Method,Slope) '.png']);
        dlmwrite([Data_Save_Folder last_folder_name sprintf('_Method_%g_Linearity_R_Slope%g',Noise_Removal_Method,Slope) '.txt'],[It_Array' Mean_Array] ,'delimiter','\t','newline','pc','precision', '%.6f');
    elseif p==2
        saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g_Linearity_G_Slope%g',Noise_Removal_Method,Slope) '.png']);
        dlmwrite([Data_Save_Folder last_folder_name sprintf('_Method_%g_Linearity_G_Slope%g',Noise_Removal_Method,Slope) '.txt'],[It_Array' Mean_Array] ,'delimiter','\t','newline','pc','precision', '%.6f');

    elseif p==3
        saveas(gcf,[Data_Save_Folder last_folder_name sprintf('_Method_%g_Linearity_B_Slope%g',Noise_Removal_Method,Slope) '.png']);
        dlmwrite([Data_Save_Folder last_folder_name sprintf('_Method_%g_Linearity_B_Slope%g',Noise_Removal_Method,Slope) '.txt'],[It_Array' Mean_Array] ,'delimiter','\t','newline','pc','precision', '%.6f');
    end
end