clear all
clf
%%
Data_Save_Folder='C:\TuanShu\180507_MTF related Plot\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%%
If_Plot_PSF=0;
Convolution_Method='normal';  %'normal', 'fast'
%% Source Related
Optical_SR=0.01;            %micron
Optical_Resolution_LSF=1;   %micron
Optical_Window_Size=10;     %micron
%% Sample Plane Related
Sample_SR=Optical_SR;       %micron, 兩者不同在convolution時會比較麻煩
FOV=Optical_Window_Size;    %micron, should be the same for fast convolute
Edge_Angle=90+5;            %degree
%% Camera Related
Magnification=200/10;       %micron
Camera_Pixel_Size=10.6;     %micron
%% Source Related Calculation
X_Optical=Optical_SR:Optical_SR:floor(Optical_Window_Size/Optical_SR)*Optical_SR;
Central_Index_Source=round(length(X_Optical)/2);
X_Optical=X_Optical-X_Optical(Central_Index_Source);
Y_Optical=X_Optical;
R_Optical=(repmat(Y_Optical',[1 length(X_Optical)]).^2+repmat(X_Optical,[length(Y_Optical) 1]).^2).^0.5;
PSF=((besselj(1,R_Optical/Optical_Resolution_LSF*3.8317)./(R_Optical/Optical_Resolution_LSF*3.8317)).^2)';
PSF(Central_Index_Source,Central_Index_Source)=(PSF(Central_Index_Source-1,Central_Index_Source)+PSF(Central_Index_Source+1,Central_Index_Source)+PSF(Central_Index_Source,Central_Index_Source-1)+PSF(Central_Index_Source,Central_Index_Source+1))/4;  

%PSF=PSF/sum(PSF(:));
%LSF_Theo=sum(PSF,2);
%LSF_Theo=(LSF_Theo-min(LSF_Theo(:)))/(max(LSF_Theo(:))-min(LSF_Theo(:)));
%LSF_Theo_FWHM=abs(find(LSF_Theo>0.5,1,'last')-find(LSF_Theo>0.5,1,'first'))*Optical_SR

%% Plot PSF
if If_Plot_PSF == 1
    clf
    imagesc(X_Optical,Y_Optical,PSF);
    colormap(gray);
    axis equal
    xlim([X_Optical(1) X_Optical(end)])
    ylim([Y_Optical(1) Y_Optical(end)])
    TSTightMargin
end
%% Sample Related Calculation
switch Convolution_Method
    case 'normal'
        XY_Temp=Optical_SR:Optical_SR:(floor((FOV+Optical_Window_Size)/Optical_SR)-1)*Optical_SR;
    case 'fast'
        XY_Temp=Optical_SR:Optical_SR:(floor((FOV)/Optical_SR))*Optical_SR;
end

Central_Index_Temp=round(length(XY_Temp)/2);
Edge=double((repmat(XY_Temp',[1 length(XY_Temp)])-XY_Temp(Central_Index_Temp))>(repmat(XY_Temp,[length(XY_Temp) 1])-XY_Temp(Central_Index_Temp))*tan(Edge_Angle/180*pi));
imagesc(Edge);
%% Blurred Sample
switch Convolution_Method
    case 'normal'
        Edge_Blur=conv2(Edge,PSF,'valid');
    case 'fast'
        Edge_Blur=abs(ifft2(fft2(Edge).*fft2(PSF)));
end
X_Sample=Optical_SR:Optical_SR:floor(FOV/Optical_SR)*Optical_SR;
X_Sample=X_Sample-X_Sample(round(length(X_Sample)/2));
Y_Sample=X_Sample;
X_Sample_Grid=repmat(X_Sample,[length(Y_Sample) 1]);
Y_Sample_Grid=repmat(Y_Sample',[1 length(X_Sample)]);
imagesc(X_Sample,Y_Sample,Edge_Blur);
imwrite(Edge_Blur./max(Edge_Blur(:)),[Data_Save_Folder sprintf('Simulated SlantEdge_A%g_SR%g_R%g_FOV%g_Raw_%s.tif',Edge_Angle,Optical_SR,Optical_Resolution_LSF,FOV,Convolution_Method)]);
%% Calculate LSF and MTF
Max_Spatial_Frequency=1000/Optical_SR;
Spatial_Frequency_SR=1000/FOV;
F=Spatial_Frequency_SR:Spatial_Frequency_SR:Max_Spatial_Frequency;
%%
ESF=Edge_Blur(:,round(length(Y_Sample)/2));
plot(ESF)
LSF=[0;diff(ESF)];
plot(X_Sample,ESF)
%plot(X_Sample,LSF/max(LSF(:)),X_Optical,LSF_Theo/max(LSF_Theo(:)));
plot(X_Sample,LSF/max(LSF(:)));
%%
%LSF_Theo_Long=zeros(length(F),1);
%LSF_Theo_Long((length(LSF_Theo_Long)/2-length(LSF_Theo)/2):(length(LSF_Theo_Long)/2+length(LSF_Theo)/2-1))=LSF_Theo;
MTF_notnorm=abs(fft(LSF));
MTF=MTF_notnorm/MTF_notnorm(1);
plot(F,MTF*100,'LineWidth',4);
xlim([0 1000]);
set(gca, 'FontSize', 14)
grid on
xlabel('Spatial Frequency (lp/mm)');
ylabel('MTF (%)');
%% Camera Related Calculation
Camera_SR=Camera_Pixel_Size/Magnification;
X_Pixel=floor(FOV/Camera_SR);
Y_Pixel=X_Pixel;
X_Camera=Camera_SR:Camera_SR:floor(FOV/Camera_SR)*Camera_SR;
X_Camera=X_Camera-X_Camera(round(length(X_Camera)/2));
Y_Camera=X_Camera;
X_Camera_Grid=repmat(X_Camera,[length(Y_Camera) 1]);
Y_Camera_Grid=repmat(Y_Camera',[1 length(X_Camera)]);
Edge_Camera=interp2(X_Sample_Grid,Y_Sample_Grid,Edge_Blur,X_Camera_Grid,Y_Camera_Grid);
imagesc(X_Camera,Y_Camera,Edge_Camera);
imwrite(Edge_Camera,[Data_Save_Folder sprintf('Simulated Slant-Edge_A%g_SR%g_R%g.tif',Edge_Angle,Camera_SR,Optical_Resolution_LSF)]);

%%
%%
