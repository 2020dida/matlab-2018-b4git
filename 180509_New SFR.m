clear all
clf
addpath('C:\TuanShu\MATLAB\TSLib');
%%
Data_Save_Folder='C:\TuanShu\180509_SFR related Plot\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%%
Image_Path='C:\TuanShu\180509_SFR related Plot\Simulated Slant-Edge_A95_SR0.45_R0.5_XFOV20_YFOV20.tif';

Input_Image=imread(Image_Path);
imagesc(Input_Image);
colormap(gray);
axis off
TSTightMargin

%% 
Error_Acceptance_micron=3E-3;
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
Required_Line_to_Comp_NCycle_Used=Required_Line_to_Comp_NCycle;
Y_ROI_Start=4;
Ori_YLine=1;
Y_ROI=[Y_ROI_Start:(Y_ROI_Start+Required_Line_to_Comp_NCycle_Used-1)];
Image_ROI=Input_Image(Y_ROI,:);
imagesc(Image_ROI);
colormap(gray);
grid on 
grid minor
ax=gca;
ax.MinorGridColor = [0.1, 0.1, 0.1]*8;  % [R, G, B]
ax.GridColor = [0.1, 0.1, 0.1]*8;  % [R, G, B]


%% Offset and average in case the N>1
Start_X=1;
SuperResolution=Known_SR*N/Required_Line_to_Comp_NCycle;
Number_of_Line=floor((size(Image_ROI,2)-Start_X)/N);
Image_ROI_Interlace=Image_ROI(:,Start_X:N:((Number_of_Line-1)*N+Start_X));
Binned_Array=reshape(Image_ROI_Interlace,[length(Image_ROI_Interlace(:)) 1]);
X=0:SuperResolution:SuperResolution*(length(Binned_Array)-1);
X_ori=0:Known_SR:Known_SR*(size(Image_ROI,2)-1);
LSF=[0;diff(Binned_Array)];
plot(X,Binned_Array,'-*');
hold on
plot(X_ori,Image_ROI(Start_X,:),'*');
hold off
