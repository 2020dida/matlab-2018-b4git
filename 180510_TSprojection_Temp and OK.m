function OutputArray=TSprojection(InputImage,Angle,SuperRatio)
%%
    Image_Path='C:\TuanShu\180510_SFR for standard image\Simulated Slant-Edge_A97.125_SR0.45_R0.5_XFOV20_YFOV20.tif';
    SR=0.45;
    InputImage=double(imread(Image_Path));
    Angle=7.125;
    SuperRatio=8;
    YROI=44;
    InputImage=InputImage(1:YROI,:);
    Max_Error=0.01;
    Max_Count=6;
    Cos=cos(Angle/180*pi);
    Sin=sin(Angle/180*pi);
    
    ArrayLength_Spixel=floor(size(InputImage,2)/Cos*SuperRatio);
    ArrayLength_Opixel=ArrayLength_Spixel/SuperRatio;
    Array_Spacing=(1/Cos/SuperRatio);   %unit: original pixel
    New_SR=SR*Array_Spacing;
    L=0:New_SR:New_SR*(ArrayLength_Spixel-1);
    Array_Coor_UnitOPixel=0:Array_Spacing:Array_Spacing*(ArrayLength_Spixel-1);
    Array_Coordinate_ThisPixel_Exact=zeros(size(InputImage,1),size(InputImage,2));%unit: original pixel
    Array_Coordinate_ThisPixel_Nearest=zeros(size(InputImage,1),size(InputImage,2));%unit: original pixel
    ErrorMap=zeros(size(InputImage,1),size(InputImage,2));           %unit: original pixel
    Array_Value=zeros(ArrayLength_Spixel,1);
    Array_Count=zeros(ArrayLength_Spixel,1);
    Array_MaxError=zeros(ArrayLength_Spixel,1);
    for p=1:size(InputImage,1)
        for q=1:size(InputImage,2)
            Array_Coordinate_ThisPixel_Exact(p,q)=(p)*Sin+(q)*Cos;      % 1st: Y, 2nd: X
            if Array_Coordinate_ThisPixel_Exact(p,q)>ArrayLength_Opixel
                continue;
            end
            [minvalue minindex]=min(abs(Array_Coor_UnitOPixel-Array_Coordinate_ThisPixel_Exact(p,q)));
            Array_Coordinate_ThisPixel_Nearest(p,q)=Array_Coor_UnitOPixel(minindex);%unit: original pixel
            ErrorMap(p,q)=abs(Array_Coordinate_ThisPixel_Exact(p,q)-Array_Coordinate_ThisPixel_Nearest(p,q))./Array_Coordinate_ThisPixel_Nearest(p,q);
            ErrorMap(isnan(ErrorMap))=0;  
            if ErrorMap(p,q)>Max_Error || Array_Count(minindex)>=Max_Count
                continue;
            end
            Array_Value(minindex)=Array_Value(minindex)+InputImage(p,q);
            Array_Count(minindex)=Array_Count(minindex)+1;
            Array_MaxError(minindex)=max(Array_MaxError(minindex),ErrorMap(p,q));            
        end
    end
    ErrorMap(isnan(ErrorMap))=0;    
    %% Solve for Zero Count
    
    Array_Mean=Array_Value./Array_Count;
    Array_Mean(isnan(Array_Mean))=0;
    
    for p=1:length(Array_Count)
        if Array_Count(p) ==0
            switch p
                case 1
                    Array_Mean(1)=Array_Mean(2);
                case length(Array_Count)
                    Array_Mean(length(Array_Count))=Array_Mean(length(Array_Count)-1);
                otherwise
                    Array_Mean(p)=mean(Array_Mean(p-1)+Array_Mean(p+1));
            end
        end
    end
    %%
    subplot(2,2,1)
    imagesc(InputImage)
    axis equal
    subplot(2,2,2)
    imagesc(ErrorMap)
    caxis([0 0.02])
    axis equal
    subplot(2,2,3)
    plot(L,Array_Count)
    subplot(2,2,4)
    plot(L,Array_Mean/mean(Array_Mean)*0.4)
    %
    ESF=Array_Mean/mean(Array_Mean);
    LSF=[0;diff(ESF)];
    subplot(1,1,1)    
    plot(L,LSF,'-*')
    %%
    Spatial_Frequency_SR=1000/(New_SR*length(LSF));
    Max_Spatial_Frequency=1000/New_SR-Spatial_Frequency_SR;
    F=0:Spatial_Frequency_SR:Max_Spatial_Frequency;
    MTF_notnorm=abs(fft(LSF));
    MTF=MTF_notnorm/MTF_notnorm(1);
    subplot(1,1,1)    
    plot(F,MTF)
    xlim([0 1000])
    ylim([0 1])
