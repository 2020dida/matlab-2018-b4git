clear all
addpath('C:\TuanShu\MATLAB\sfrmat3_post\');
addpath('C:\TuanShu\MATLAB\TsLib\');
%% Parameters;
SR=0.954;  %micron
%% Flags
If_Last=1;
If_Predefined_ROI=1;
%%  Paths
Image_Name='';
FolderPath='C:\TuanShu\180521_Image Folder for Test\';
%% predefined ROI
X_Start=434;
Y_Start=500;
X_Size=100;
Y_Size=100;
%%
LastImageName=TSLastModified(FolderPath);
if If_Last == 1
    Original_Image=double(imread([FolderPath '\' LastImageName]));
else
    Original_Image=double(imread([FolderPath '\' Image_Name]));
end
imagesc(Original_Image);
axis equal off
colormap(gray)
if If_Predefined_ROI == 1
hold on
rectangle('Position',[X_Start Y_Start X_Size Y_Size],'EdgeColor','r')
hold off
end
%%
if If_Predefined_ROI == 1
    Image=Original_Image(Y_Start:(Y_Start+Y_Size-1),X_Start:(X_Start+X_Size-1));
    imwrite(uint16(Image/max(Image(:))*65535),[FolderPath '\' 'Test_ROI.tif'],'tiff');

else
    Image=Original_Image;
end
%%
HalfSR = 1000/(2*SR);    % half-sampling frequency
[OutputData, ESF, Super_Ratio, SR2, Y_Size_Used]=TSsfrmat3(SR,Image);
%% Plot
clf
subplot(2,1,1)
plot(OutputData(:,1),OutputData(:,2)*100,'LineWidth',4)
set(gca, 'FontSize', 14)
xlim([0 HalfSR]);
ylim([0 100]);
grid on
xlabel('Spatial Frequency (cy/mm)');
ylabel('SFR-MTF (%)');

subplot(2,1,2)
plot(ESF,'LineWidth',4)
set(gca, 'FontSize', 14)
grid on
