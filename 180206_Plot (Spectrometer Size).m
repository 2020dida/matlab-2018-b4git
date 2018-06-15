clear all
%%

Size=[27.5 31 19 21 18.3 12.5 11.5];% 29.5];

Spectrl_Res=[0.024 0.04 0.06 0.1 0.1 0.5 0.068];% 0.068];
Imaging_Depth=840*840/4./Spectrl_Res/1000/1000/1.35;
Model={{'Tornado OCT-Prime',sprintf('%gnm/pix',Spectrl_Res(1))}, 
       {'Wasatch CS800-840-80',sprintf('%gnm/pix',Spectrl_Res(2))}, 
       {'Wasatch CS800-840-114',sprintf('%gnm/pix',Spectrl_Res(3))}, 
       {'BaySpec DeepView',sprintf('%gnm/pix',Spectrl_Res(4))}, 
       {'G&H'},%{,sprintf('%gnm/pix',Spectrl_Res(5))}, 
       {'Ibsen EAGLE',sprintf('%gnm/pix',Spectrl_Res(6))}, 
       {'Tornado OCTANE', sprintf('%gnm/pix',Spectrl_Res(7))}};
       %{'New Span SP-850', sprintf('%gnm/pix',Spectrl_Res(8))}};


scatter(Size(1:(end)),Spectrl_Res(1:(end)),'filled');
set(gca,'yscale','log')
xlabel('Largest Dimension (cm)')
ylabel('Wavelength Dispersion (nm/pix)')
dx = 0.3; 
dy = 0;
text(Size+dx, Spectrl_Res+dy, Model,'FontSize',12);
xlim([10 40]);
ylim([0.015 0.7]);
set(gca,'fontsize',14)

Model={{'Tornado OCT-Prime',sprintf('%.01fmm',Imaging_Depth(1))}, 
       {'Wasatch CS800-840-80',sprintf('%.01fmm',Imaging_Depth(2))}, 
       {'Wasatch CS800-840-114',sprintf('%.01fmm',Imaging_Depth(3))}, 
       {'BaySpec DeepView',sprintf('%.01fmm',Imaging_Depth(4))}, 
       {'G&H'},%{,sprintf('%gnm/pix',Spectrl_Res(5))}, 
       {'Ibsen EAGLE',sprintf('%.01fmm',Imaging_Depth(6))}, 
       {'Tornado', 'OCTANE (On-Chip)', sprintf('%.01fmm',Imaging_Depth(7))}};
       %{'New Span SP-850', sprintf('%.01fnm/pix',Imaging_Depth(8))}};


scatter(Size(1:(end)),Imaging_Depth(1:(end)),'filled');
%set(gca,'yscale','log')
xlabel('Largest Dimension (cm)')
ylabel('Maximum Imaging Depth (mm)')
dx = 0.3; 
dy = 0.1;
text(Size+dx, Imaging_Depth+dy, Model,'FontSize',12);
xlim([10 40]);
%ylim([0.015 0.7]);
set(gca,'fontsize',14)
