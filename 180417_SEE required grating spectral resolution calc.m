clear all
%%
Wavelength_Center=0.57;  %1.3   0.57
Wavelength_Short=0.48;  %1.22   0.48
Wavelength_Long=0.66;   %1.36   0.66
Corresponding_FOV=550;  %micron 223 550
%%
Order=1;
Grating_Pitch=1000; %1144   1000
Grating_Diameter=1.35;   %%2.2 1.35
Total_Pitch=Grating_Pitch*Grating_Diameter;
Wavelength_Resolution=Wavelength_Center/Total_Pitch/Order;
Spectral_Spatial_Resolution=Wavelength_Resolution/(Wavelength_Long-Wavelength_Short)*Corresponding_FOV;
Number_of_Resolvable_Point=(Wavelength_Long-Wavelength_Short)/Wavelength_Resolution;
%%
