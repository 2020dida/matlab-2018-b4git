clear all
clf
cd('C:\TuanShu\180530_MTF Data Analysis\');
addpath('C:\TuanShu\MATLAB\sfrmat3_post\');
addpath('C:\TuanShu\MATLAB\TsLib\');
%% Parameters;
SR=0.0495;%5.5/12.21*1.12;%0.0495;%0.099;%0.0495;%1.59;%0.2385;%1.1925;%0.954;  %micron
FF_Blur=2;
%% Flags
If_Normalized_MTF=1;
FF_Measurement=0;
If_Apply_FF=0;
Fie_Select_Mode='Bulk';
ROI_Mode='Predefined';
If_Use_Default_Size=1;
%%  Paths
Predefined_Image_Name='';
Predefined_FolderPath='C:\TuanShu\180531_MTF\WIO 40x\';
FF_Folder=[Predefined_FolderPath '\FF\'];
if exist(FF_Folder)==0
    mkdir(FF_Folder);
end
LastFFName=TSLastModified(FF_Folder);
if isnan(LastFFName) ~=1
    LastFFNumber=sscanf(LastFFName,'FF%g.tif');
else
    LastFFNumber=0;
end
%% predefined ROI
X_Start=450;%434
Y_Start=425;
Default_Xwidth=150;
Default_Ywidth=150;

%% File Selection Modes
switch Fie_Select_Mode
    case 'Last'
        FolderPath=Predefined_FolderPath;
        Image_Name=TSLastModified(FolderPath);
        Original_Image=double(imread([FolderPath '\' Image_Name]));
        Exe_Times=1;
    case 'Predefined'
        FolderPath=Predefined_FolderPath;
        Image_Name=Predefined_Image_Name;
        Original_Image=double(imread([FolderPath '\' Image_Name]));
        Exe_Times=1;
    case 'Manual'
        [Image_Name,FolderPath] = uigetfile('*.*');
        Manual_Image_Name=[FolderPath '\' Image_Name];
        Original_Image=double(imread([Manual_Image_Name]));
        Exe_Times=1;
    case 'Bulk'
        FolderPath=Predefined_FolderPath;
        Image_List=dir([Predefined_FolderPath '\*.tif']);
        Exe_Times=length(Image_List);
end

for p=1:Exe_Times
    if strcmp(Fie_Select_Mode,'Bulk')
        Image_Name=Image_List(p).name;
        Original_Image=double(imread([Predefined_FolderPath '\' Image_Name]));
    end
    %% ROI Selection
    switch ROI_Mode
        case 'Predefined'
        Image=Original_Image(Y_Start:(Y_Start+Default_Ywidth-1),X_Start:(X_Start+Default_Xwidth-1));
        %imwrite(uint16(Image/max(Image(:))*65535),[FolderPath '\' 'Test_ROI.tif'],'tiff');
        case 'Manual'
        go=1;
        while go
            %% Main
            AxMain=subplot(2,1,1);
            imagesc(Original_Image);
            axis equal
            axis off
            colormap(gray);
            if go ==2
                hold on
                rectangle('Position',Last_Position,'EdgeColor','r');
                hold off
            end
            %% Sub
            h = imrect(AxMain);
            Position=h.getPosition ;
            if (Position(1)+Position(3))>size(Original_Image,2) || (Position(2)+Position(4))>size(Original_Image,1) || Position(1)<1 || Position(2)<1 
                break;
            end
            go=2;
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
            ROI_Image=TSnormalizeSFRimage(Original_Image(Y_min:Y_max,X_min:X_max));
            AxSub=subplot(2,1,2);
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
        Image=ROI_Image;
        case 'None'
        Image=Original_Image;
    end
    %% MTF
    [OutputData, ESF, Super_Ratio, SR2, Y_Size_Used,Angle]=TSsfrmat3(SR,Image);
    ESF=TSnorm(ESF);
    LSF=(smooth(diff(ESF),10));
    LSF=LSF/max(LSF);
    LSF=[LSF(1); LSF];
    %% Interpolation for Uniform sampling
    Required_MTF_Sampling=1;
    F_array=0:Required_MTF_Sampling:max(OutputData(:,1));
    MTF=interp1(OutputData(:,1),OutputData(:,2),F_array);
    if If_Normalized_MTF == 1
        Norm_MTF=MTF/max(MTF(1:find(MTF<0.5,1,'first')));
    end
    %% Plot
    clf
    HalfSR = 1000/(2*SR);    % half-sampling frequency
    subplot(2,2,[1 2])
    plot(F_array,MTF*100,'LineWidth',4)
    if If_Normalized_MTF == 1 && max(MTF)>1
        hold on
        plot(F_array,Norm_MTF*100,'LineWidth',4)
        hold off
    end
    title(sprintf('MTF50=%glp/mm',F_array(find(MTF<0.5,1,'first'))),'Color','b')
    set(gca, 'FontSize', 14)
    %xlim([0 HalfSR*2]);
    xlim([0 2500]);
    ylim([0 100]);
    grid on
    xlabel('Spatial Frequency (lp/mm)');
    ylabel('SFR-MTF (%)');

    subplot(2,2,4)
    plot(ESF,'LineWidth',4)
    hold on
    plot(LSF,'LineWidth',2)
    hold off
    xlim([0 size(ESF,1)]);
    ylim([-0.1 1]);

    set(gca, 'FontSize', 14)
    ylabel('ESF/LSF');
    grid minor
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);

    subplot(2,2,3)
    imagesc(Original_Image);
    axis equal off
    colormap(gray)
    hold on
    switch ROI_Mode
        case 'Manual'
            rectangle('Position',Last_Position,'EdgeColor','r')
        otherwise
            rectangle('Position',[X_Start Y_Start Default_Xwidth Default_Ywidth],'EdgeColor','r')
    end
    hold off
    TSTightMargin
    TSsaveas([FolderPath '\' sprintf('SFR-MTF_SR%g_%s.tif',SR,Image_Name)]);
    dlmwrite([FolderPath sprintf('SFR-MTF_SR%g_%s.txt',SR,Image_Name)],[F_array;MTF*100]','delimiter','\t','newline','pc','precision', '%.6f');
    disp(sprintf('%g/%g',p,Exe_Times));
    %%
    fclose all
end
