clear all
clf
cd('C:\TuanShu\MATLAB\TSLib');
%%
MTF_Type='LSFbased';    %'LSFbased': ifo & Zemax notation; 'PSFbased': test, maybe SFR notation
If_Plot_Signal=0;
If_Plot_PSF=0;
If_Plot_LSF=0;
If_Plot_MTF=1;
If_Plot_Output=0;
%%
Data_Save_Folder='C:\TuanShu\180507_MTF related Plot\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%% LSF related
Resolution_LSF=1;     %micron
LSF_Type='Airy';    %'Gaussian' 'Airy'
%% Input Signal Related
Spatial_Period=2;                       %micron
Spatial_Phase=0;                     %radian
Input_Contrast=1;
%%
Spatial_Frequency=1000/Spatial_Period;  %lp/mm
%% Domain related
Spatial_SR=0.01;                         %micron
Spatial_Width=100;                      %micron
%%
X=(Spatial_SR:Spatial_SR:floor(Spatial_Width/Spatial_SR)*Spatial_SR)-Spatial_Width/2;
Central_Index=round(length(X)/2);
X_ROI=(Central_Index-round(Spatial_Period/Spatial_SR)):(Central_Index+round(Spatial_Period/Spatial_SR));
%%
Spatial_Signal=(1+Input_Contrast*sin(2*pi*X/Spatial_Period+Spatial_Phase))/2;
%%
switch LSF_Type
    case 'Gaussian'
        LSF=gaussmf(X,[Resolution_LSF/(8*log(2))^0.5 0]);
    case 'Airy'
        Y=X;
        R=(repmat(Y',[1 length(X)]).^2+repmat(X,[length(Y) 1]).^2).^0.5;
        PSF=((besselj(1,R/Resolution_LSF*3.8317)./(R/Resolution_LSF*3.8317)).^2)';

        PSF(Central_Index,Central_Index)=(PSF(Central_Index-1,Central_Index)+PSF(Central_Index+1,Central_Index)+PSF(Central_Index,Central_Index-1)+PSF(Central_Index,Central_Index+1))/4;  
        PSF_XSlice=PSF(:,5000);
        PSF_XSlice=PSF_XSlice./max(PSF_XSlice(:));
        LSF_X=sum(PSF,2);
        LSF_X=(LSF_X-min(LSF_X(:)))./(max(LSF_X(:))-min(LSF_X(:)));
        LSF=LSF_X;
        LSF_FWHM=abs(find(LSF>0.5,1,'last')-find(LSF>0.5,1,'first'))*Spatial_SR
end
%% Plot Input Signal
if If_Plot_Signal ==1
    clf
    plot(X,Spatial_Signal,'LineWidth',4);
    daspect([5 1 1])
    set(gca, 'FontSize', 14)
    set(gca,'XTick',[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1]*Spatial_Period);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[0:0.1:1]);
    xlim([-1 1]*Spatial_Period)
    ylim([0 1])
    grid on
    xlabel('X Position (\mum))');
    ylabel('Input Signal');
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('Input Spatial Signal_P%g_C%g.png',Spatial_Period,Input_Contrast)]);
    %% Show Input Signal as Bar Chart
    imagesc(X,Y,Spatial_Signal)
    colormap(gray)
    caxis([0 1])
    xlim([-1*Spatial_Period 1*Spatial_Period]);
    axis off
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('Input Bar Chart_P%g_C%g.png',Spatial_Period,Input_Contrast)]);
end
%% Plot PSF
if If_Plot_PSF ==1
    clf
    xlim([X_ROI(1) X_ROI(end)]);
    ylim([X_ROI(1) X_ROI(end)]);
    axis equal off
    colormap(gray)
    TSTightMargin
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
    xlim([-1 1]*Spatial_Period)
   
    %xlim([-10*Resolution_LSF 10*Resolution_LSF])

    grid on
    xlabel('X Position (\mum))');
    ylabel('Normalized LSF');
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('LSF_FWHM%g_R%g.png',LSF_FWHM,Resolution_LSF)]);
end
%%

Max_Spatial_Frequency=1000/Spatial_SR;
Spatial_Frequency_SR=1000/Spatial_Width;
F=Spatial_Frequency_SR:Spatial_Frequency_SR:Max_Spatial_Frequency;
switch MTF_Type
    case 'LSFbased'
        MTF_notnorm=abs(fft(LSF));
    case 'PSFbased'
        MTF_notnorm=abs(fft(PSF_XSlice));
end
MTF=MTF_notnorm/MTF_notnorm(1);
%%
if If_Plot_MTF == 1
    clf
    plot(F,MTF*100,'LineWidth',4);
    xlim([0 1000/Spatial_Period*2]);
    set(gca, 'FontSize', 14)
    grid on
    xlabel('Spatial Frequency (lp/mm)');
    ylabel('MTF (%)');

    TSsaveas([Data_Save_Folder sprintf('MTF_FWHM%g_R%g_%s.png',LSF_FWHM,Resolution_LSF,MTF_Type)]);

end
%% Calculate Output Contrast Through Convolution
if If_Plot_Output == 1
    clf
    LSF_uni=LSF/sum(LSF);
    Output_Signal=conv(Spatial_Signal,LSF_uni,'same')
    plot(X,Output_Signal,'LineWidth',4);
    daspect([2.5 1 1])
    set(gca, 'FontSize', 14)
    set(gca,'XTick',[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1]*Spatial_Period);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[0:0.1:1]);
    xlim([-1 1])
    ylim([0 1])
    grid on
    xlabel('X Position (\mum))');
    ylabel('Output Signal');
    TSTightMargin
    MTFatPeriod=MTF(find(F>Spatial_Frequency,1,'first'))
    Output_Contrast=max(Output_Signal(X_ROI))-min(Output_Signal(X_ROI))
    TSsaveas([Data_Save_Folder sprintf('Output Spatial Signal_P%g_IC%g_MTF%g_OC%g.png',Spatial_Period,Input_Contrast,MTFatPeriod,Output_Contrast)]);
    %% Show Output Signal as Bar Chart
    imagesc(X,Y,Output_Signal)
    colormap(gray)
    caxis([0 1])
    xlim([-1*Spatial_Period 1*Spatial_Period]);
    axis off
    TSTightMargin
    TSsaveas([Data_Save_Folder sprintf('Output Bar Chart_P%g_IC%g_MTF%g_OC%g.png',Spatial_Period,Input_Contrast,MTFatPeriod,Output_Contrast)]);
end