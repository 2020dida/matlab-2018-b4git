clear all
%%
N=3;
New_NA=0.22;
Min_New_Core=2;
Core_Cladding_Ratio=0;
Pixelation=0:0.1:2;
for p=1:length(Core_Cladding_Ratio)
Number_of_Core=sum(N:-1:((N-1)/2+1))*2-N;
Original_Core=106;
Original_Area=pi*(Original_Core/2)^2;
Original_NA=0.22;
Original_SR=2*pi*(1-cos(asin(Original_NA)));
New_Core=max(Original_Core/N,Min_New_Core);
Pixelation=New_Core*Core_Cladding_Ratio(p);
New_Area=pi*(New_Core/2)^2;
New_SR=2*pi*(1-cos(asin(Original_NA)));
Packing_Density(p)=New_Area*Number_of_Core/(pi*((New_Core*N+Pixelation*(N-1))/2)^2);

Resulting_Dia=Number_of_Core*New_Core;

Radiance_Gain=Packing_Density*(Original_Core*asin(Original_NA))/(New_Core*asin(New_NA));
Radiance_Gain_Expected=Radiance_Gain*0.7;
end

plot(Core_Cladding_Ratio*100,Packing_Density*100,'LineWidth',2);

grid on
grid minor
xlabel('Cladding/Core Ratio (%)')
ylabel('Packing Density (%)')
set(gca,'fontsize',14)