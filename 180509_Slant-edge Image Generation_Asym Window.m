clear all
clf
addpath('C:\TuanShu\MATLAB\TSLib');
%%
Data_Save_Folder='C:\TuanShu\180507_MTF related Plot\';
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%%
If_Plot_PSF=0;
%% Source Related
Optical_SR=0.01;            %micron
Optical_Resolution_LSF=0.5;   %micron
Optical_Window_Size_X=20;     %micron
Optical_Window_Size_Y=20;     %micron

%% Sample Plane Related
Sample_SR=Optical_SR;       %micron, 兩者不同在convolution時會比較麻煩
FOV_X=Optical_Window_Size_X;    %micron, should be the same for fast convolute
FOV_Y=Optical_Window_Size_Y;    %micron, should be the same for fast convolute
Edge_Angle=90+7.125;            %degree
%% Camera Related
Camera_SR=0.45;     %micron
%% Source Related Calculation
X_Optical=Optical_SR:Optical_SR:floor(Optical_Window_Size_X/Optical_SR)*Optical_SR;
X_Optical=X_Optical-X_Optical(round(length(X_Optical)/2));
Y_Optical=Optical_SR:Optical_SR:floor(Optical_Window_Size_Y/Optical_SR)*Optical_SR;
Y_Optical=Y_Optical-Y_Optical(round(length(Y_Optical)/2));
R_Optical=(repmat(Y_Optical',[1 length(X_Optical)]).^2+repmat(X_Optical,[length(Y_Optical) 1]).^2).^0.5;
PSF=((besselj(1,R_Optical/Optical_Resolution_LSF*3.8317)./(R_Optical/Optical_Resolution_LSF*3.8317)).^2);
PSF(round(length(Y_Optical)/2),round(length(X_Optical)/2))=(PSF(round(length(Y_Optical)/2)-1,round(length(X_Optical)/2))+PSF(round(length(Y_Optical)/2)+1,round(length(X_Optical)/2))+PSF(round(length(Y_Optical)/2),round(length(X_Optical)/2)-1)+PSF(round(length(Y_Optical)/2),round(length(X_Optical)/2)+1))/4;  

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
clf
X_Edge=Optical_SR:Optical_SR:floor((FOV_X+Optical_Window_Size_X)/Optical_SR-1)*Optical_SR;
X_Edge=X_Edge-X_Edge(round(length(X_Edge)/2));
Y_Edge=Optical_SR:Optical_SR:floor((FOV_Y+Optical_Window_Size_Y)/Optical_SR-1)*Optical_SR;
Y_Edge=Y_Edge-Y_Edge(round(length(Y_Edge)/2));

Edge=double(repmat(Y_Edge',[1 length(X_Edge)])>repmat(X_Edge,[length(Y_Edge) 1])*tan(Edge_Angle/180*pi));
imagesc(X_Edge,Y_Edge,Edge);
%% Blurred Sample
Edge_Blur=conv2(Edge,PSF,'valid');
X_Blur=Optical_SR:Optical_SR:floor(FOV_X/Optical_SR)*Optical_SR;
X_Blur=X_Blur-X_Blur(round(length(X_Blur)/2));
Y_Blur=Optical_SR:Optical_SR:floor(FOV_Y/Optical_SR)*Optical_SR;
Y_Blur=Y_Blur-Y_Blur(round(length(Y_Blur)/2));

X_Blur_Grid=repmat(X_Blur,[length(Y_Blur) 1]);
Y_Blur_Grid=repmat(Y_Blur',[1 length(X_Blur)]);
imagesc(X_Blur,Y_Blur,Edge_Blur);
axis equal
imwrite(Edge_Blur./max(Edge_Blur(:)),[Data_Save_Folder sprintf('BlurEdge_A%g_SR%g_R%g_FOVX%g_FOVY%g_Raw.tif',Edge_Angle,Optical_SR,Optical_Resolution_LSF,FOV_X,FOV_Y)],'tiff');
%% Calculate LSF and MTF (X)
Spatial_Frequency_SR=1000/FOV_X;
Max_Spatial_Frequency=1000/Optical_SR-Spatial_Frequency_SR;
F=0:Spatial_Frequency_SR:Max_Spatial_Frequency;     %注意必須從0開始!!! 原本的程式可能也得改
%%
ESF=Edge_Blur(round(length(Y_Blur)/2),:);
plot(ESF)
if ESF(end)>ESF(1)
    LSF=[0 diff(ESF)];
else
    LSF=[0 diff(ESF(end:-1:1))];
end
plot(X_Blur,ESF)
%plot(X_Sample,LSF/max(LSF(:)),X_Optical,LSF_Theo/max(LSF_Theo(:)));
plot(X_Blur,LSF/max(LSF(:)));

%%
%LSF_Theo_Long=zeros(length(F),1);
%LSF_Theo_Long((length(LSF_Theo_Long)/2-length(LSF_Theo)/2):(length(LSF_Theo_Long)/2+length(LSF_Theo)/2-1))=LSF_Theo;
clf
MTF_notnorm=abs(fft(LSF));
MTF=MTF_notnorm/MTF_notnorm(1);
plot(F,MTF*100,'LineWidth',4);
xlim([0 1000]);
ylim([0 100]);
set(gca, 'FontSize', 14)
grid on
xlabel('Spatial Frequency (lp/mm)');
ylabel('MTF (%)');

TSTightMargin
TSsaveas([Data_Save_Folder sprintf('MTF_FOVX%02g_FOVY%02g_WX%02g_WY%02g_Angle%g.png',FOV_X,FOV_Y,Optical_Window_Size_X,Optical_Window_Size_Y,Edge_Angle)]);
%% Camera Related Calculation
X_Pixel=floor(FOV_X/Camera_SR);
X_Pixel=floor(FOV_Y/Camera_SR);
X_Camera=Camera_SR:Camera_SR:floor(FOV_X/Camera_SR)*Camera_SR;
X_Camera=X_Camera-X_Camera(round(length(X_Camera)/2));
Y_Camera=Camera_SR:Camera_SR:floor(FOV_Y/Camera_SR)*Camera_SR;
Y_Camera=Y_Camera-Y_Camera(round(length(Y_Camera)/2));
X_Camera_Grid=repmat(X_Camera,[length(Y_Camera) 1]);
Y_Camera_Grid=repmat(Y_Camera',[1 length(X_Camera)]);
Edge_Camera=interp2(X_Blur_Grid,Y_Blur_Grid,Edge_Blur,X_Camera_Grid,Y_Camera_Grid);
imagesc(X_Camera,Y_Camera,Edge_Camera);
imwrite(uint16(Edge_Camera/max(Edge_Camera(:))*65535),[Data_Save_Folder sprintf('Simulated Slant-Edge_A%g_SR%g_R%g_XFOV%g_YFOV%g.tif',Edge_Angle,Camera_SR,Optical_Resolution_LSF,FOV_X,FOV_Y)],'tiff');