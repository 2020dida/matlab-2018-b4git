clear all
clf
%%
Folder_List={'C:\TuanShu\180504_MLOBJ stray\ML worst_8.4mW_0.16it_2' 'C:\TuanShu\180504_MLOBJ stray\ML_min distance\' 'C:\TuanShu\180504_MLOBJ stray\Oly worst_8.4mW_0.16it'};
Ratio=11/8.4;
Name_List={'025-050' '025-100' '025-150' '025-200' '025-250' '025-300'};
Smooth_Size=5;
Max_Matrix=zeros(length(Name_List),length(Folder_List));
Max_DN=3800;
for q=1:length(Folder_List)
    for p=1:length(Name_List)
        Name=Name_List{p};
        Image=Ratio*(imread([Folder_List{q} '\' Name_List{p} '_Sam.tif'])-imread([Folder_List{q} '\' Name_List{p} '_Ref.tif']));
        Image_Smooth=imgaussfilt(Image,Smooth_Size);
        Max_Matrix(p,q)=max(Image_Smooth(:));
        subplot('Position',[(p-1)*1/length(Name_List) (length(Folder_List)-q)*1/length(Folder_List) 1/length(Name_List) 1/length(Folder_List)])
        imagesc(Image);
        axis equal off
        colormap(gray)
        caxis([0 600])
    end
end
subplot(1,1,1)
plot(Max_Matrix,'-*','LineWidth',4)

XTickArray=1:length(Name_List);
XTickLabelArray={50 100 150 200 250 300};
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)

set(gca, 'FontSize', 14)
ylabel('Stray Light (DN)');
xlabel('Cylindrical Lens Focal Length (mm)');
legend('MLOBJ, +5cm spacing','MLOBJ, min. spacing','Olympus, min. spacing');
%%
Max_Ratio=1-(Max_Matrix/Max_DN);
Max_dB=-20*log10(Max_Ratio);


subplot(1,1,1)
plot(Max_dB,'-*','LineWidth',4)

XTickArray=1:length(Name_List);
XTickLabelArray={50 100 150 200 250 300};
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)

set(gca, 'FontSize', 14)
ylabel('SNR Loss (dB)');
xlabel('Cylindrical Lens Focal Length (mm)');
legend('MLOBJ, +5cm spacing','MLOBJ, min. spacing','Olympus, min. spacing');

%%
Speed_Loss_Ratio=(100-(10.^(-1*Max_dB/10))*100);

subplot(1,1,1)
plot(Speed_Loss_Ratio,'-*','LineWidth',4)

XTickArray=1:length(Name_List);
XTickLabelArray={50 100 150 200 250 300};
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)

set(gca, 'FontSize', 14)
ylabel('Speed Loss (%)');
xlabel('Cylindrical Lens Focal Length (mm)');
legend('MLOBJ, +5cm spacing','MLOBJ, min. spacing','Olympus, min. spacing');
