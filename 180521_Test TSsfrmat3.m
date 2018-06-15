clear all
%%
addpath('C:\TuanShu\MATLAB\sfrmat3_post\');
ImagePath='C:\TuanShu\180517_Image Normalization for Comparison\MySystem\Plan N 20X\Processed Data\ROI_Image_X100_Y100_A5.39457_SR0.954.tif';
Image=imread(ImagePath);
SR=0.1908;  %micron
HalfSR = 1000/(2*SR);    % half-sampling frequency
[OutputData, ESF, Super_Ratio, SR2, Y_Size_Used]=TSsfrmat3(SR,Image);
%% Plot
clf
subplot(1,1,1)
plot(OutputData(:,1),OutputData(:,2)*100,'LineWidth',4)
set(gca, 'FontSize', 14)
xlim([0 HalfSR]);
ylim([0 100]);
grid on
xlabel('Spatial Frequency (cy/mm)');
ylabel('SFR-MTF (%)');