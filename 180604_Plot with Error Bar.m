clear all
%%
Number_of_Line=24;
File_Path='C:\TuanShu\180604_Imaging Guiding Test\Summarizing Data.txt';
fin=fopen(File_Path);
C = textscan(fin,'%f %f %f %f %f %f %f %f',Number_of_Line,'Delimiter','\t');
fclose('all');
%
%1 f: #
%2 f: Alignment #
%3 f: AC1
%4 f: CY
%5 f: Peak Ratio
%6 f: PtP Variance
%7 f: FOI Width
%8 f: Expected FOI Width
%% Get Name Cell & Data Arrays
Number=ones(Number_of_Line,1);
Alignment=ones(Number_of_Line,1);
AC1=ones(Number_of_Line,1);        %4
CY=ones(Number_of_Line,1);               %4
PeakRatio=ones(Number_of_Line,1);               %5
PtPVariance=ones(Number_of_Line,1);   %6
FOIwidth=ones(Number_of_Line,1);       %7
FOIwidthExpected=ones(Number_of_Line,1);       %7


for p=1:Number_of_Line
    Number(p)=C{1}(p);
    Alignment(p)=C{2}(p);               %4
    AC1(p)=C{3}(p);               %5
    CY(p)=C{4}(p);   %6
    PeakRatio(p)=C{5}(p);;       %7
    PtPVariance(p)=C{6}(p);           %8
    FOIwidth(p)=C{7}(p);            %10
    FOIwidthExpected(p)=C{8}(p);;          %9
end

