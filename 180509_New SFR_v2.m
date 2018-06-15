clear all
clf
addpath('C:\TuanShu\MATLAB\TSLib');
%%
Data_Save_Folder='C:\TuanShu\180509_SFR related Plot\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%%
Image_Path='C:\TuanShu\180509_SFR related Plot\Simulated Slant-Edge_A95_SR0.45_R1_XFOV20_YFOV20.tif';

Input_Image=imread(Image_Path);
imagesc(Input_Image);
colormap(gray);
axis off
TSTightMargin

%% 
Error_Acceptance_micron=2E-2;
Known_Angle=5;  %degree
Known_SR=0.45;  %micron
Error_Acceptance=1E-5;

Interline_Ratio_Dfference=tan(Known_Angle/180*pi);

Required_Line_to_Comp_1Cycle=round(1/Interline_Ratio_Dfference);

Error_Ratio=[];
Error_micron=[];
N=1;
while N==1 || (Error_micron(N-1)>Error_Acceptance_micron)
    Required_Line_to_Comp_NCycle=round(N/Interline_Ratio_Dfference);
    Error_Ratio(N)=abs(Required_Line_to_Comp_NCycle*Interline_Ratio_Dfference-N)/N;
    Error_micron(N)=Error_Ratio(N)*N;
    if (Error_micron(N)<Error_Acceptance_micron)
        break;
    end
    N=N+1;
end
Required_Line_to_Comp_NCycle
clf
plot(Error_Ratio)
%% Not doing interpolation
Y_ROI_Start=1;
Y_ROI=[Y_ROI_Start:(Y_ROI_Start+Required_Line_to_Comp_NCycle-1)];
Image_ROI=Input_Image(Y_ROI,:);
%% Offset and average in case the N>1
ESF_Cell=cell(N,1);
X_Cell=cell(N,1);
ESF_ori_Cell=cell(N,1);
X_ori_Cell=cell(N,1);
for p=1:N
    Start_X=p;
    SuperResolution=Known_SR*N/Required_Line_to_Comp_NCycle;
    Number_of_Line=floor((size(Image_ROI,2)-Start_X)/N);
    Image_ROI_Interlace=Image_ROI(:,Start_X:N:((Number_of_Line-1)*N+Start_X));
    ESF_Cell{p}=reshape(Image_ROI_Interlace,[length(Image_ROI_Interlace(:)) 1]);
    X_Cell{p}=0:SuperResolution:SuperResolution*(length(ESF_Cell{p})-1);
    
    ESF_ori_Cell{p}=Image_ROI((p-1)*Required_Line_to_Comp_1Cycle+1,:);  %只是大約, 不同line之間的error會比上面的error大
    X_ori_Cell{p}=0:Known_SR:Known_SR*(size(Image_ROI,2)-1);    %const
end
%%
SetNumber=1;
plot(X_Cell{SetNumber},ESF_Cell{SetNumber},'-*');
hold on
plot(X_ori_Cell{SetNumber},ESF_ori_Cell{SetNumber},'*');    %概念是在圖中找一條線對應
hold off
%% LSF
LSF_Cell=cell(N,1);
LSF_ori_Cell=cell(N,1);
for p=1:N
    LSF_Cell{p}=[0; diff(ESF_Cell{p})];
    LSF_ori_Cell{p}=[0 diff(ESF_ori_Cell{p})]*N/Required_Line_to_Comp_NCycle;
end
%%
SetNumber=1;
plot(X_Cell{SetNumber},LSF_Cell{SetNumber},'-*');
hold on
plot(X_ori_Cell{SetNumber},LSF_ori_Cell{SetNumber},'*');    %概念是在圖中找一條線對應
hold off
%% Generating F for MTF
F_Cell=cell(N,1);
for p=1:N
    Spatial_Frequency_SR=1000/(SuperResolution*length(LSF_Cell{p}));
    Max_Spatial_Frequency=1000/SuperResolution-Spatial_Frequency_SR;
    F_Cell{p}=0:Spatial_Frequency_SR:Max_Spatial_Frequency;     %注意必須從0開始!!! 原本的程式可能也得改
end
%% MTF
MTF_Cell=cell(N,1);
for p=1:N
    MTF_notnorm=abs(fft(LSF_Cell{p}));
    MTF_Cell{p}=MTF_notnorm/MTF_notnorm(1);
end
%%
SetNumber=1;
plot(F_Cell{SetNumber},MTF_Cell{SetNumber}*100,'LineWidth',4);
xlim([0 1000]);
ylim([0 100]);
set(gca, 'FontSize', 14)
grid on
xlabel('Spatial Frequency (lp/mm)');
ylabel('MTF (%)');