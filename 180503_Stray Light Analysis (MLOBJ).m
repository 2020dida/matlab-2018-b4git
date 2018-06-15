clear all
%%
Folder_Path='C:\TuanShu\180503_MLOBJ_Stray\';

REF='Ref';
DUT='ML';

Number=1;

Stray=imread([Folder_Path '\' DUT sprintf('%g',Number) '.tif'])-imread([Folder_Path '\' REF sprintf('%g',Number) '.tif']);

imagesc(Stray);
axis equal off
colormap(gray)
caxis([0 500])