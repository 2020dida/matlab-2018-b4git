clear all
%%
Main_Folder=('C:\TuanShu\180301_Oil and Elast Comp (txt files)');

Sets{1}=['Oil'];

Sets{2}=['Elastomer'];

SubSets{1}='R';
SubSets{2}='G';
SubSets{3}='B';

ArrayLabel{1}='Exposure Time (ms)';
ArrayLabel{2}='Signal (DN)';
Column1_Value_Considered=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];

Mean_Cell=cell(length(Sets),length(SubSets));
STD_Cell=cell(length(Sets),length(SubSets));
for p=1:length(Sets)
    for q=1:3
        Folder=[Main_Folder '\' Sets{p} '\' SubSets{q} '\'];
        File_List=dir(Folder);
        Original_Lnegth_folder_list=length(File_List);
        for r=1:Original_Lnegth_folder_list
            if File_List(Original_Lnegth_folder_list+1-r).isdir == 1
                File_List(Original_Lnegth_folder_list+1-r)=[];
            elseif strcmp(File_List(Original_Lnegth_folder_list+1-r).name,'.') == 1
                File_List(Original_Lnegth_folder_list+1-r)=[];
            elseif strcmp(File_List(Original_Lnegth_folder_list+1-r).name,'..') == 1
                File_List(Original_Lnegth_folder_list+1-r)=[];
            end
        end
        Temp_Matrix=zeros(length(Column1_Value_Considered),length(ArrayLabel),length(File_List));
        for r=1:length(File_List)
            Data_Temp=dlmread([Folder '\' File_List(r).name]);
            Original_Length=size(Data_Temp,1);
            for w=1:Original_Length
                if max(Data_Temp(Original_Length-w+1,1)==Column1_Value_Considered)~=1
                   Data_Temp(Original_Length-w+1,:)=[];
                end
            end
        Temp_Matrix(:,:,r)=Data_Temp;
        end
        Mean_Cell{p,q}=mean(Temp_Matrix,3);
        STD_Cell{p,q}=std(Temp_Matrix,[],3);

    end
    
end

%%
Set_Number=2;
SubSet_Number=1;
If_Error=0;
Color_Array_Subset{1}='red';
Color_Array_Subset{2}='green';
Color_Array_Subset{3}='blue';

subplot(1,1,1);
%plot(Mean_Cell{Set_Number,SubSet_Number}(:,1),Mean_Cell{Set_Number,SubSet_Number}(:,2),['-o' Color_Array_Subset{SubSet_Number}],'LineWidth',2,'MarkerFaceColor',Color_Array_Subset{SubSet_Number});
%errorbar(Mean_Cell{Set_Number,SubSet_Number}(:,1),Mean_Cell{Set_Number,SubSet_Number}(:,2),STD_Cell{Set_Number,SubSet_Number}(:,2),['-s' Color_Array_Subset{SubSet_Number}],'LineWidth',2,'MarkerFaceColor',Color_Array_Subset{SubSet_Number});
for p=1:3
    if If_Error ==1
        errorbar(Mean_Cell{Set_Number,p}(:,1),Mean_Cell{Set_Number,p}(:,2),STD_Cell{Set_Number,p}(:,2),['-s' Color_Array_Subset{p}],'LineWidth',2,'MarkerFaceColor',Color_Array_Subset{p});
    else
        plot(Mean_Cell{Set_Number,p}(:,1),Mean_Cell{Set_Number,p}(:,2),['-s' Color_Array_Subset{p}],'LineWidth',2,'MarkerFaceColor',Color_Array_Subset{p});
    end

    hold on
end
hold off
ylim([0 260])
xlim([0 20])
title(['Exposure vs Signal (' Sets{Set_Number} ')'],'FontSize',16);
xlabel('Exposure Time (ms)','FontSize',12);
ylabel('Signal (DN)','FontSize',12);
set(gca,'fontsize',12);

pbaspect auto

saveas(gcf,[Main_Folder '\' sprintf('%s_EvS',Sets{Set_Number}) '.png']);
%% Ratio Calc
NNN=2;
Temp=Mean_Cell{2,NNN}./Mean_Cell{1,NNN}
Temp(1,2)