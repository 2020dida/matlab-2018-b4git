clear all
addpath('C:\TuanShu\MATLAB\TSLib');
clf
%% Target Image Format
% L: black
% R: white
% Edge: RT to LB
%% Default Size
If_Use_Default_Size=1;
Default_Xwidth=100; %200 100 for vertical, 100 200 for horizontal
Default_Ywidth=100;
%% Data Parameter
folder_path='C:\TuanShu\180517_Image Normalization for Comparison\MySystem\Plan N 20X\';
file_name='170518_172322_Exp26';
%% System Parameter
System='My4x';
if strcmp(System,'MLdata180430')
    Camera_Pixel_Size=5.5; %micron
    Magnification=12.21;
elseif  strcmp(System,'My4x')
    Camera_Pixel_Size=10.6; %micron
    Magnification=500/45;
end
%% Estimation Method Parameter
Sampling_Resolution=Camera_Pixel_Size/Magnification;
%% AOI (Automatic Object Inspection) Parameter
ContrastTH=0.7;
%%
Data_Save_Folder=[folder_path 'Processed Data'];
if ~exist(Data_Save_Folder,'dir')
    mkdir(Data_Save_Folder);
end
%%
Data_Image=double(imread([folder_path '\' file_name],'tiff'));%.*FF;
%% Selection
go=1;
while go
    %% Main
    AxMain=subplot(1,2,1);
    imagesc(Data_Image);
    axis equal
    axis off
    colormap(gray);
    %% Sub
    h = imrect(AxMain);
    Position=h.getPosition ;
    if (Position(1)+Position(3))>size(Data_Image,2) || (Position(2)+Position(4))>size(Data_Image,1) || Position(1)<1 || Position(2)<1 
        break;
    end
    X_min=round(Position(1));
    Y_min=round(Position(2));
    if If_Use_Default_Size == 0
        X_max=round(Position(1)+Position(3));
        Y_max=round(Position(2)+Position(4));
        Last_Position=Position;

    else
        X_max=round(Position(1)+Default_Xwidth);
        Y_max=round(Position(2)+Default_Ywidth);
        Last_Position=[Position(1) Position(2) Default_Xwidth Default_Ywidth];
    end
    ROI_Image=TSnormalizeSFRimage(Data_Image(Y_min:Y_max,X_min:X_max));
    Angle=TSSimpleEdgeAngle(ROI_Image,ContrastTH);
    if ~isnan(Angle)
        disp(sprintf('The Edge angle is %g,XFOV=%g, YFOV=%g',Angle,Position(3),Position(4)));
    else
       disp('Image is Not Simple');
        %Angle=TSAngleDetection(ROI_Image,Type);
    end
    AxSub=subplot(1,2,2);
    imagesc(ROI_Image);
    axis equal
    axis off
    if If_Use_Default_Size == 0
        Last_Width=Position(3);
        Last_Height=Position(4);
    else
        Last_Width=Default_Xwidth;
        Last_Height=Default_Ywidth;       
    end
    pause(1)
end
%%
imwrite(uint16(ROI_Image/max(ROI_Image(:))*65535),[Data_Save_Folder '\' sprintf('ROI_Image_X%g_Y%g_A%g_SR%g.tif',Last_Width,Last_Height,Angle,Sampling_Resolution)],'tiff');
%% Save ROI
subplot(1,1,1)
imagesc(Data_Image);
axis equal off
hold on
rectangle('Position',Last_Position,'EdgeColor','r')
hold off

TSTightMargin
TSsaveas([Data_Save_Folder '\' sprintf('Marked Image_X%g_Y%g_A%g_SR%g.tif',Last_Width,Last_Height,Angle,Sampling_Resolution)]);
