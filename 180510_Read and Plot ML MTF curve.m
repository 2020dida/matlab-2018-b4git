clear all
%%
T=dlmread('C:\TuanShu\180430_ML data\central-S.txt');
%%
plot(T(:,1)/0.45*1000,T(:,2),'LineWidth',4);
set(gca, 'FontSize', 14)
grid on
xlabel('Spatial Frequency (lp/mm)');
ylabel('MTF (%)')
xlim([0 1700])