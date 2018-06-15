clear all
clf
cd('C:\TuanShu\MATLAB\TSLib');
%%
If_Plot_Signal=0;
If_Plot_LSF=1;
%%
Data_Save_Folder='C:\TuanShu\180507_MTF related Plot\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%% LSF related
Resolution_LSF=1;     %micron
LSF_Type='Airy';    %'Gaussian' 'Airy'
Airy_1stZero_toFWHM_Coef=0.8420;
%% Input Signal Related
Spatial_Period=1;                       %micron
Spatial_Phase=0;                     %radian
Input_Contrast=0.5;
%%
Spatial_Frequency=1000/Spatial_Period;  %lp/mm
%% Domain related
Spatial_SR=0.01;                         %micron
Spatial_Width=100;                      %micron
%%
X=(Spatial_SR:Spatial_SR:floor(Spatial_Width/Spatial_SR)*Spatial_SR)-Spatial_Width/2;
%%
Spatial_Signal=(1+Input_Contrast*sin(2*pi*X/Spatial_Period+Spatial_Phase))/2;
%%
switch LSF_Type
    case 'Gaussian'
        LSF=gaussmf(X,[Resolution_LSF/(8*log(2))^0.5 0]);
    case 'Airy'
        LSF=(besselj(1,X/Resolution_LSF*3.8317*Airy_1stZero_toFWHM_Coef)./(X/Resolution_LSF*3.8317*Airy_1stZero_toFWHM_Coef)).^2;
        LSF=LSF./max(LSF(:));
        [value NAN_Index]=max(isnan(LSF));
        LSF(NAN_Index)=(LSF(NAN_Index-1)+LSF(NAN_Index+1))/2;
end
%% Plot Input Signal
if If_Plot_Signal ==1
    plot(X,Spatial_Signal,'LineWidth',4);
    daspect([2.5 1 1])
    set(gca, 'FontSize', 14)
    set(gca,'XTick',[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[0:0.1:1]);
    ylim([0 1])
    grid on
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('Input Spatial Signal_P%g_C%g.png',Spatial_Period,Input_Contrast)]);
    %% Show Input Signal as Bar Chart
    imagesc(Spatial_Signal)
    colormap(gray)
    caxis([0 1])
    axis off
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('Input Bar Chart_P%g_C%g.png',Spatial_Period,Input_Contrast)]);
end
%% Plot LSF
if If_Plot_LSF ==1
    clf
    plot(X,LSF,'LineWidth',4);
    daspect([5 1 1])
    set(gca, 'FontSize', 14)
    %set(gca,'XTick',[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]);
    %set(gca,'XTickLabel',[]);
    ylim([0 1])
    xlim([-2*Resolution_LSF 2*Resolution_LSF])
    %xlim([-10*Resolution_LSF 10*Resolution_LSF])

    grid on
    xlabel('X Position (\mum))');
    ylabel('Normalized LSF');
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('LSF_FWHM%g.png',Resolution_LSF)]);
    LSF_FWHM=abs(find(LSF>0.5,1,'last')-find(LSF>0.5,1,'first'))*Spatial_SR
end
%%
clf
Max_Spatial_Frequency=1000/Spatial_SR;
Spatial_Frequency_SR=1000/Spatial_Width;
F=Spatial_Frequency_SR:Spatial_Frequency_SR:Max_Spatial_Frequency;
MTF_notnorm=abs(fft(LSF));
MTF=MTF_notnorm/MTF_notnorm(1);
plot(F,MTF*100,'LineWidth',4);
xlim([0 1000/Resolution_LSF]);
set(gca, 'FontSize', 14)
grid on
xlabel('Spatial Frequency (lp/mm)');
ylabel('MTF (%)');