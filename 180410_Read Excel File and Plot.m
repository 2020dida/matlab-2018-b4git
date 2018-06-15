clear all
%%
NotAvaliable_Value=-1;
Tol=0.01;
%%
Number_of_Line=18;
File_Path='C:\TuanShu\180410_Light Source Related\Light Sources.txt';
fin=fopen(File_Path);
C = textscan(fin,'%s %s %s %f %f %f %f %f %f %f %f %f',Number_of_Line,'Delimiter','\t');
fclose('all');
%
%1 s 1: Type
%2 s 2: Manufacterer
%3 s 3: Model
%4 f 1: Price (NTD)
%5 f 2: Power (mW)
%6 f 3: Wavelength (nm)
%7 f 4: BW
%8 f 5: Fiber dia
%9 f 6: Fiber area (calculated)
%10 f 7: NA
%11 f 8: Theta (calculated)
%12 f 9: Solid Angle  (calculated)
%% Get Name Cell & Data Arrays
Manufacturer=cell(Number_of_Line,1);
Model=cell(Number_of_Line,1);

% Necessary
Index=ones(Number_of_Line,1);               %4
Price=ones(Number_of_Line,1);               %4
Power=ones(Number_of_Line,1);               %5
Wavelength_Center=ones(Number_of_Line,1);   %6
Wavelength_BW=ones(Number_of_Line,1);       %7
% Optional
Fiber_Dia=ones(Number_of_Line,1);           %8
Fiber_NA=ones(Number_of_Line,1);            %10
% Calculated
Fiber_Area_Read=ones(Number_of_Line,1);          %9
Fiber_Theta_Read=ones(Number_of_Line,1);         %11
Fiber_SR_Read=ones(Number_of_Line,1);            %12


for p=1:Number_of_Line
    Manufacturer{p}=char(C{2}(p));
    Model{p}=char(C{3}(p));; %char(10)
    Price(p)=C{4}(p);
    % Necessary
    Index(p)=p;
    Price(p)=C{4}(p);               %4
    Power(p)=C{5}(p);               %5
    Wavelength_Center(p)=C{6}(p);   %6
    Wavelength_BW(p)=C{7}(p);;       %7
    % Optional
    Fiber_Dia(p)=C{8}(p);;           %8
    Fiber_NA(p)=C{10}(p);;            %10
    % Calculated
    Fiber_Area_Read(p)=C{9}(p);;          %9
    Fiber_Theta_Read(p)=C{11}(p);;         %11
    Fiber_SR_Read(p)=C{12}(p);;            %12

end

%% Handling Fiber_Area
Mask=(Fiber_Dia~=NotAvaliable_Value);
Fiber_Area_Calculated=pi*(Fiber_Dia/2).^2;
Check=((abs((Fiber_Area_Calculated-Fiber_Area_Read)./Fiber_Area_Calculated)).*Mask)<Tol;
if min(Check) == 0
    Index=find(Check==0,1);
    warning(sprintf('Fiber Area Error at #%g etc.',Index));
end
Fiber_Area=Fiber_Area_Calculated.*Mask+Fiber_Area_Read.*(1-Mask);

%% Handling Fiber_Theta
Mask=(Fiber_NA~=NotAvaliable_Value);
Fiber_Theta_Calculated=asin(Fiber_NA);
Check=((abs((Fiber_Theta_Calculated-Fiber_Theta_Read)./Fiber_Theta_Calculated)).*Mask)<Tol;
if min(Check) == 0
    Index=find(Check==0,1);
    warning(sprintf('Fiber Theta Error at #%g etc.',Index));
end
Fiber_Theta=Fiber_Theta_Calculated.*Mask+Fiber_Theta_Read.*(1-Mask);

%% Handling Fiber_SR
Mask=(Fiber_NA~=NotAvaliable_Value);
Fiber_SR_Calculated=2*pi()*(1-cos(Fiber_Theta));
Check=((abs((Fiber_SR_Calculated-Fiber_SR_Read)./Fiber_SR_Calculated)).*Mask)<Tol;
if min(Check) == 0
    Index=find(Check==0,1);
    warning(sprintf('Fiber SR Error at #%g etc.',Index));
end
Fiber_SR=Fiber_SR_Calculated.*Mask+Fiber_SR_Read.*(1-Mask);

%% Cacalcuating Other Parameters
% Radiance
Radiance=(Power*0.001)./(Fiber_Area.*Fiber_SR*1E-12);
Radiance_Norm=Radiance/Radiance(1);
Radiance_Norm_dB=10*log10(Radiance_Norm);   %10x, ~shot-limited SNR
% Bandwidth %假設在freq-based時為對稱, 在wavelength-based不對稱
Frequency_BW=3E8./(Wavelength_Center.^2).*Wavelength_BW;
Related_BW=(Frequency_BW-Frequency_BW(1))/Frequency_BW(1)*100;
%% Generating Scatterer Label
Lebal=cell(Number_of_Line,1);
for p=1:Number_of_Line
    %Lebal{p}=[sprintf('%g.',Index(p)) '\fontsize{8}' Model{p} '\fontsize{5}' Manufacturer{p}];
    Lebal{p}=['\fontsize{8}' Model{p} '\fontsize{5}' Manufacturer{p}];
end
%% Data Exclusion
Index_Excluded=[];
Mask=Index;               %4
Lebal_Masked=Lebal;

for p=1:length(Index_Excluded)
    Mask(Index_Excluded(length(Index_Excluded)-p+1))=[];
    Lebal_Masked(Index_Excluded(length(Index_Excluded)-p+1),:)=[];
end
%% Data Grouping
%% Scatterer Plot
%Line_Excluded
scatter(Radiance_Norm_dB(Mask),Related_BW(Mask),'filled');
%set(gca,'yscale','log')
grid on
grid minor
xlabel('Related SNR (dB)')
ylabel('Related Bandwidth (%)')
ax=gca;
ax.YTick =[-100 -75 -50 -25 0 25 50 75 100];
ax.XTick =[-40 -30 -20 -10 0 10 20 30 40];
dx = 0; 
dy = 0;
textfit(Radiance_Norm_dB(Mask)+dx, Related_BW(Mask)+dy, Lebal_Masked,0.5,0,'Color',[0.5 0.5 0.5],'FontSize',12,'HorizontalAlignment','center');
xlim([-45 45]);
ylim([-120 120]);
set(gca,'fontsize',14)