clear all
%%
Number_of_Line=18;
File_Path='C:\TuanShu\180411_MiiS\OPL.txt';
D = dlmread(File_Path);
fclose('all');
%%
OPL=sum(abs(D(:,1)).*D(:,2))
