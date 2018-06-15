clear all
%%
Folder_Path='C:\TuanShu\180503_MLOBJ_Stray\Run_Oly_Worst\';

Name='025-050-Worst';
if_Oly=1;
if_ML=0;
if if_Oly == 1
    Oly=imread([Folder_Path '\' Name '_Oly.tif'])-imread([Folder_Path '\' Name '_Ref.tif']);
end
if if_ML == 1
ML=imread([Folder_Path '\' Name '_ML.tif'])-imread([Folder_Path '\' Name '_Ref.tif']);
end
if if_Oly && if_ML == 1
    subplot(2,1,1)
    imagesc(ML);
    axis equal off
    colormap(gray)
    caxis([0 400])
    subplot(2,1,2)
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