clear all
clf
%%
Folder_Path_Oly='C:\TuanShu\180504_MLOBJ stray\ML_min distance\';
Folder_Path_ML='C:\TuanShu\180503_MLOBJ_Stray\Run_ML_Worst\';

Name_List={'025-050-Worst' '025-100-Worst' '025-150-Worst' '025-200-Worst' '025-250-Worst' '025-300-Worst'};
Number=1:6;
if_Oly=1;
if_ML=1;
Smooth_Size=5;
Max_Matrix=zeros(length(Number),2);
for p=1:length(Number)
    Name=Name_List{p};
    if if_Oly == 1
        Oly=imread([Folder_Path_Oly '\' Name '_Oly.tif'])-imread([Folder_Path_Oly '\' Name '_Ref.tif']);
        Oly_Smooth=imgaussfilt(Oly,Smooth_Size);
        Max_Matrix(p,1)=max(Oly_Smooth(:));
    end
    if if_ML == 1
        ML=imread([Folder_Path_ML '\' Name '_ML.tif'])-imread([Folder_Path_ML '\' Name '_Ref.tif']);
        ML_Smooth=imgaussfilt(ML,Smooth_Size);
        Max_Matrix(p,2)=max(ML_Smooth(:));
    end
    if if_Oly && if_ML == 1
        %subplot(length(Number),2,1+(p-1)*2)
        subplot('Position',[(p-1)*1/length(Number) 0.375 1/length(Number) 0.2])

        imagesc(ML);
        axis equal off
        colormap(gray)
        caxis([0 400])
        %subplot(length(Number),2,2+(p-1)*2)
        subplot('Position',[(p-1)*1/length(Number) 0.625 1/length(Number) 0.2])
        imagesc(Oly);
    else
        subplot(1,1,1)
        if if_Oly ==1
            imagesc(Oly);

        else
            imagesc(ML);
        end
    end
    axis equal off
    colormap(gray)
    caxis([0 400])
end
plot(Max_Matrix)