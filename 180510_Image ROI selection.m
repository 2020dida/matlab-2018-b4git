clear all
addpath('C:\TuanShu\MATLAB\TSLib');
%% Data Parameter
folder_path='C:\TuanShu\180510_SFR ROI for ML data\';
file_name='#3-RB';
%% System Parameter
System='MLdata180430';
if strcmp(System,'MLdata180430')
    Camera_Pixel_Size=5.5; %micron
    Magnification=12.21;
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
    AxMain=subplot(2,1,1);
    imagesc(Data_Image);
    axis equal
    axis off
    colormap(gray);
    %% Sub
    h = imrect(AxMain);
    Position=h.getPosition ;
    if abs(Position(3))>size(Data_Image,1) || abs(Position(4))>size(Data_Image,2) || Position(1)<1 || Position(2)<1 
        break;
    end
    X_min=round(Position(1));
    X_max=round(Position(1)+Position(3));
    Y_min=round(Position(2));
    Y_max=round(Position(2)+Position(4));
    ROI_Image=Data_Image(Y_min:Y_max,X_min:X_max);
    Angle=TSSimpleEdgeAngle(ROI_Image,ContrastTH);
    if ~isnan(Angle)
        disp(sprintf('The Edge angle is %g,XFOV=%g, YFOV=%g',Angle,Position(3),Position(4)));
    else
       disp('Image is Not Simple');
        %Angle=TSAngleDetection(ROI_Image,Type);
    end
    AxSub=subplot(2,1,2);
    imagesc(ROI_Image);
    axis equal
    axis off
    Last_Width=Position(3);
    Last_Height=Position(4);
    pause(1)
end
%%
imwrite(uint16(ROI_Image/max(ROI_Image(:))*65535),[Data_Save_Folder '\' sprintf('ROI_Image_X%g_Y%g_A%g_SR%g.tif',Last_Width,Last_Height,Angle,Sampling_Resolution)],'tiff');