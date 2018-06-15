clear all
fclose all
addpath('C:\TuanShu\MATLAB\TSLib');
%%
Size_X=1024;
Size_Y=1024;
Target_YFOV=11;
If_Sub_Stray=1;
SR=10.6/200*9;
%%
Folder_Number=[1];  
File_Number=[351];      %160¬°peak
Case=0;
switch Case
    case 4 
        Reference_Peak_Value=2568.7;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Selected Data\Alignment 4\'; %4(asym) 5 6 7 8 ok
    case 5
        Reference_Peak_Value=2308.1;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Selected Data\Alignment 5\'; %4(asym) 5 6 7 8 ok
    case 6
        Reference_Peak_Value=2332.0;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Selected Data\Alignment 6\'; %4(asym) 5 6 7 8 ok
    case 7
        Reference_Peak_Value=2341.0;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Selected Data\Alignment 7\'; %4(asym) 5 6 7 8 ok
    case 8
        Reference_Peak_Value=2313.2;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Selected Data\Alignment 8\'; %4(asym) 5 6 7 8 ok
    case -1 %for Change CY case
        Reference_Peak_Value=2416;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Change CY\Alignment 1\'; %4(asym) 5 6 7 8 ok
    case -3 %for Change CY case
        Reference_Peak_Value=2486;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Change CY\Alignment 3\'; %4(asym) 5 6 7 8 ok
    case 0 %for Change CY case
        Reference_Peak_Value=100;
        Root_Folder_Path='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Change MMF\'; %4(asym) 5 6 7 8 ok
end

% Alignment 4: 2568.7
% Alignment 5: 2308.1
% Alignment 6: 2332.0
% Alignment 7: 2341.0
% Alignment 8: 2313.2

Data_Save_Folder='C:\TuanShu\180604_Imaging Guiding Test\d10.8cm\Selected Data\';
Folder_List=dir(Root_Folder_Path);
Folder_List(1:2)=[];
SubFolder=[Root_Folder_Path Folder_List(Folder_Number).name]
File_List=dir(SubFolder);
File_List(1:2)=[];


fin=fopen([SubFolder '\' File_List(File_Number).name]);
Image=fread(fin,[Size_X,Size_Y],'uint16')';


imagesc(Image);
colormap(gray);
caxis([0 4095]);

%
X_Index=1:size(Image,1);
Y_Index=1:size(Image,2);
Y_Profile=mean(Image,2);
Y_Diff_Profile=abs(diff(smooth(Y_Profile,8)));
Y_Diff_Profile=[Y_Diff_Profile(1); Y_Diff_Profile];
plot(Y_Index,Y_Profile/max(Y_Profile),Y_Index,Y_Diff_Profile/max(Y_Diff_Profile));
%
[value Y_Peak_Index]=max(Y_Profile);
[value Y_Lower]=max(Y_Diff_Profile(1:Y_Peak_Index));
[value Y_Upper]=max(Y_Diff_Profile(Y_Peak_Index:end));
Y_Upper=Y_Upper+Y_Peak_Index-1;
Height=(Y_Upper-Y_Lower)*SR;
%
imagesc(Image)
hold on
rectangle('Position',[1 Y_Lower Size_X Y_Upper-Y_Lower+1],'EdgeColor','y');
hold off
%
X_Profile_SigStray=mean(Image(Y_Lower:Y_Upper,:),1);
plot(X_Profile_SigStray);

X_Profile_StrayLower=mean(Image(1:Y_Lower,:),1);
X_Profile_StrayUpper=mean(Image(Y_Upper:end,:),1);
plot(X_Index,X_Profile_SigStray,X_Index,X_Profile_StrayLower,X_Index,X_Profile_StrayUpper);
% Find Max Stray Light X Profile
Image_wo_Signal=Image;
Buffer=0.5*(Y_Upper-Y_Lower+1);
Image_wo_Signal((Y_Lower-Buffer):(Y_Upper+Buffer),:)=0;
Image_wo_Signal=imgaussfilt(Image_wo_Signal,2);
imagesc(Image_wo_Signal);
Max_Stray_X_Profile=max(Image_wo_Signal,[],1);
plot(Max_Stray_X_Profile);


% Find Best FOV
Best_Mean_Value=0;
Best_Edge_Value=0;
Best_Mean_Index=1;
Best_Variance_Index=1;
Stray=repmat(Max_Stray_X_Profile,[Target_YFOV 1]);
for p=1:(Y_Upper-Y_Lower+1-Target_YFOV+1)
    Current_Y_Start=Y_Lower+p-1;
    if If_Sub_Stray == 1
            Current_FOV=Image(Current_Y_Start:(Current_Y_Start+Target_YFOV-1),:)-Stray;
    else
            Current_FOV=Image(Current_Y_Start:(Current_Y_Start+Target_YFOV-1),:);
    end
    Current_Mean=mean(Current_FOV(:));
    Edge=[Current_FOV(:,1:200) Current_FOV(:,(end-200+1):end)];
    Current_Edge_Value=mean(Edge(:));
    if Current_Mean>Best_Mean_Value
        Best_Mean_Value=Current_Mean;
        Best_Mean_Index=p;
    end
    if Current_Edge_Value>Best_Edge_Value
        Best_Edge_Value=Current_Edge_Value;
        Best_Edge_Index=p;
    end
end

Best_Mean_Y_Range=[Y_Lower+Best_Mean_Index-1 Y_Lower+Best_Mean_Index-1+Target_YFOV-1];
Best_Edge_Y_Range=[Y_Lower+Best_Edge_Index-1 Y_Lower+Best_Edge_Index-1+Target_YFOV-1];

Best_Mean_Profile=mean(Image(Best_Mean_Y_Range(1):Best_Mean_Y_Range(2),:),1);
Best_Edge_Profile=mean(Image(Best_Edge_Y_Range(1):Best_Edge_Y_Range(2),:),1);


subplot(3,3,[1:6])
imagesc(Image)
hold on
rectangle('Position',[1 Y_Lower Size_X Y_Upper-Y_Lower+1],'EdgeColor','m');
rectangle('Position',[1 Best_Mean_Y_Range(1) Size_X Best_Mean_Y_Range(2)-Best_Mean_Y_Range(1)+1],'EdgeColor','b');
%rectangle('Position',[1 Best_Edge_Y_Range(1) Size_X Best_Edge_Y_Range(2)-Best_Edge_Y_Range(1)+1],'EdgeColor','r');
hold off
axis equal off
hsp1 = get(gca, 'Position');                     % Get 'Position' for (2,1,1)
caxis([0 4095])

subplot(3,3,[7 8 9])
plot(X_Index,Best_Mean_Profile,'LineWidth',3)
PtP_Var=(max(Best_Mean_Profile)-min(Best_Mean_Profile))/max(Best_Mean_Profile);
ylim([500 3000])
ylim([0 2000])

xlim([1 Size_X])
legend(sprintf('Peak Value=%g%% \nPeak-to-Peak Variance=%g%% \nFOI Width=%gmicron',round(max(Best_Mean_Profile)/Reference_Peak_Value*100),round(PtP_Var*100),round(Height)),'Location', 'Best')
grid on
hsp2 = get(gca, 'Position');                     % Get 'Position' for (2,1,2)
set(gca, 'Position', [hsp2(1:2)  hsp1(3) 0.5*hsp1(4)]);    % Use 2*(2,1,1) Height for (2,1,2)
TSsaveas([Data_Save_Folder '\' sprintf('%s_XProfile_Peak%g_Var%g_Height%g.tif',Folder_List(Folder_Number).name,round(max(Best_Mean_Profile)/Reference_Peak_Value*100),round(PtP_Var*100),round(Height))]);
%TSTightMargin
%% 
