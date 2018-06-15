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
%% Plot per Alignment
Plot_Alignment=5:8;
for p=1:length(Plot_Alignment)
    plot(AC1(Alignment==Plot_Alignment(p)),100*PeakRatio(Alignment==Plot_Alignment(p)),AC1(Alignment==Plot_Alignment(p)),100*PtPVariance(Alignment==Plot_Alignment(p)),'LineWidth',3);
    hold on
end
ylim([0 100]);
xlabel('AC1 Focal Length (mm)');
ylabel('Percentage (%)');
set(gca, 'FontSize', 14);
hold off
grid on
%% Averaging by Alignment
Plot_Alignment=5:8;
Peak_Ratio_Temp=zeros(length(PeakRatio(Alignment==Plot_Alignment(1))),length(Plot_Alignment));
PtPVariance_Temp=zeros(length(PtPVariance(Alignment==Plot_Alignment(1))),length(Plot_Alignment));
for p=1:length(Plot_Alignment)
    Peak_Ratio_Temp(:,p)=PeakRatio(Alignment==Plot_Alignment(p));
    PtPVariance_Temp(:,p)=PtPVariance(Alignment==Plot_Alignment(p));
end
Peak_Ratio_Mean=mean(Peak_Ratio_Temp,2);
PtPVariance_Mean=mean(PtPVariance_Temp,2);
Peak_Ratio_STD=std(Peak_Ratio_Temp,0,2);
PtPVariance_STD=std(PtPVariance_Temp,0,2);
errorbar(AC1(Alignment==Plot_Alignment(p)),100*Peak_Ratio_Mean,100*Peak_Ratio_STD,'LineWidth',3);
hold on
errorbar(AC1(Alignment==Plot_Alignment(p)),100*PtPVariance_Mean,100*PtPVariance_STD,'LineWidth',3);
hold off
ylim([0 110]);
xlim([25 40]);
xlabel('AC1 Focal Length (mm)');
ylabel('Percentage (%)');
set(gca, 'FontSize', 14);
grid on
title(sprintf('Averaged from %g Sets of Data',length(Plot_Alignment)),'Color','b')
legend('Peak Intensity Ratio (%)','Peak-to-Valley Ratio (%)')



