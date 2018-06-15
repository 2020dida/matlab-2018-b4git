clear all
%%

Material='Water';
Wavelength_Range=[0.7 1];    %micron

%%
Wavelength_Array=Wavelength_Range(1):0.01:Wavelength_Range(2);

if strcmp(Material,'Water')
    A1=0.5672526103;
    A2=0.1736581125;
    A3=0.02121531502;
    A4=0.1138493213;
    B1=0.005085550461;
    B2=0.01814938654;
    B3=0.02617260739;
    B4=10.73888649;
end

RI=(1+(A1*Wavelength_Array.^2)./(Wavelength_Array.^2-B1)+(A2*Wavelength_Array.^2)./(Wavelength_Array.^2-B2)+(A3*Wavelength_Array.^2)./(Wavelength_Array.^2-B3)+(A4*Wavelength_Array.^2)./(Wavelength_Array.^2-B4)).^0.5;
plot(Wavelength_Array,RI);
% 
% Pole1=(1+(A1*Wavelength_Array.^2)./(Wavelength_Array.^2-B1)).^0.5;
% Pole2=(1+(A2*Wavelength_Array.^2)./(Wavelength_Array.^2-B2)).^0.5;
% Pole3=(1+(A3*Wavelength_Array.^2)./(Wavelength_Array.^2-B3)).^0.5;
% Pole4=(1+(A4*Wavelength_Array.^2)./(Wavelength_Array.^2-B4)).^0.5;
% plot(Wavelength_Array,Pole1,Wavelength_Array,Pole2,Wavelength_Array,Pole3,Wavelength_Array,Pole4);
