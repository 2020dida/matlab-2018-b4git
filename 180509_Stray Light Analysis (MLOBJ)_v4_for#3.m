clear all
clf
%%
Folder_List={'C:\TuanShu\180509_#3MLOBJ stray\MLOBJ#2_New Worst+1_016it_with silicon' 'C:\TuanShu\180509_#3MLOBJ stray\MLOBJ#3_New Worst+1_016it_with silicon\' 'C:\TuanShu\180509_#3MLOBJ stray\Oly_New worst'};
Ratio=11/8.4;
Name_List={'025-050' '025-100' '025-150' '025-200' '025-250' '025-300'};
Smooth_Size=5;
Max_Matrix=zeros(length(Name_List),length(Folder_List));
Sum_Matrix=zeros(length(Name_List),length(Folder_List));

Max_DN=3800;
for q=1:length(Folder_List)
    for p=1:length(Name_List)
        Name=Name_List{p};
        Image=Ratio*(imread([Folder_List{q} '\' Name_List{p} '_Sam.tif'])-imread([Folder_List{q} '\' Name_List{p} '_Ref.tif']));
        Image_Smooth=imgaussfilt(Image,Smooth_Size);
        Max_Matrix(p,q)=max(Image_Smooth(:));
        Sum_Matrix(p,q)=sum(Image_Smooth(:));        
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
ylabel('Maximum Stray Light (DN)');
xlabel('Achromatic Lens Focal Length (mm)');
legend('OBJ#2','OBJ#3','Olympus');

%%
subplot(1,1,1)
plot(Sum_Matrix,'-*','LineWidth',4)

XTickArray=1:length(Name_List);
XTickLabelArray={50 100 150 200 250 300};
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)

set(gca, 'FontSize', 14)
ylabel('Sum Stray Light (DN)');
xlabel('Achromatic Lens Focal Length (mm)');
legend('OBJ#2','OBJ#3','Olympus');
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
grid on
ylabel('SNR Loss (dB)');
xlabel('Achromatic Lens Focal Length (mm)');
legend('OBJ#2','OBJ#3','Olympus');

%%
Speed_Loss_Ratio=(100-(10.^(-1*Max_dB/10))*100);

subplot(1,1,1)
plot(Speed_Loss_Ratio,'-*','LineWidth',4)

XTickArray=1:length(Name_List);
XTickLabelArray={50 100 150 200 250 300};
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)
grid on

set(gca, 'FontSize', 14)
ylabel('Speed Loss (%)');
xlabel('Achromatic Lens Focal Length (mm)');
legend('OBJ#2','OBJ#3','Olympus');
