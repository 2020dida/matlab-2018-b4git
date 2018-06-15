clear all
clf
addpath('C:\TuanShu\MATLAB\TSLib');
%%
Data_Save_Folder='C:\TuanShu\180510_SFR ROI for ML data\Processed Data\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%%
Image_Path='C:\TuanShu\180510_SFR for standard image\Simulated Slant-Edge_A100_SR0.45_R1_XFOV20_YFOV20.tif';

Input_Image=double(imread(Image_Path));
imagesc(Input_Image);
colormap(gray);
axis equal off
TSTightMargin

%%
Spatial_Frequency_SR_Limit=100;      %lp/mm
Error_Acceptance_micron=2E-2;
Known_Angle=10;  %degree
Known_SR=0.45;  %micron
Error_Acceptance=1E-5;

Interline_Ratio_Dfference=abs(tan(Known_Angle/180*pi));

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
%%
SuperResolutionRatio=round(Required_Line_to_Comp_NCycle/N);
SuperResolutionRatio=4;
%%
SuperResolution=Known_SR/SuperResolutionRatio;

%% Projection
ESF=TSprojection(Input_Image,Known_Angle,SuperResolutionRatio);
X=0:SuperResolution:SuperResolution*(length(ESF)-1);
ESF_ori=Input_Image(1,:);  %只是大約, 不同line之間的error會比上面的error大
X_ori=[0:Known_SR:Known_SR*(size(Input_Image,2)-1)]+0.5;    %const
%
plot(X,ESF,'-*');
hold on
plot(X_ori,ESF_ori,'*');    %概念是在圖中找一條線對應
hold off
%% LSF
Maximum_Window_Size=floor(1000/Spatial_Frequency_SR_Limit/SuperResolution);
LSF_full=[0; diff(ESF)];
[LSF_peakvalue LSF_peakindex]=max(LSF_full);
LSF=LSF_full((LSF_peakindex-round(Maximum_Window_Size/2)):(LSF_peakindex+round(Maximum_Window_Size/2)-1));
LSF_ori=[0 diff(ESF_ori)]/SuperResolutionRatio;
%
plot(X,LSF_full,'-*');
% hold on
% plot(X_ori,LSF_ori,'*');    %概念是在圖中找一條線對應
% hold off
%% Generating F for MTF
Spatial_Frequency_SR=1000/(SuperResolution*length(LSF));
Max_Spatial_Frequency=1000/SuperResolution-Spatial_Frequency_SR;
F=0:Spatial_Frequency_SR:Max_Spatial_Frequency;     %注意必須從0開始!!! 原本的程式可能也得改

%% MTF
MTF_notnorm=abs(fft(LSF));
MTF=MTF_notnorm/MTF_notnorm(1);

%%
plot(F,MTF*100,'LineWidth',4);
xlim([0 1000]);
ylim([0 100]);
set(gca, 'FontSize', 14)
grid on
xlabel('Spatial Frequency (lp/mm)');
ylabel('MTF (%)');