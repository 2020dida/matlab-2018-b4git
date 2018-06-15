clear all
%%
Focal_Length=1.25;  %micron Zemax notation, physical focal length/image space RI
Order=-1;
Wavelength=[1.22:0.01:1.36]+0;    
Theta_In=31.24;
Grating_Pitch=1144; %1144   1000 950
Aperture_Diameter=2.2;   %%2.2 1.35
Total_Pitch=Grating_Pitch*Aperture_Diameter/cos(Theta_In/180*pi);
Wavelength_Resolution=abs(mean(Wavelength)/Total_Pitch/Order);
RI_Inc=1.85;
RI_Out=1;
d=1000/Grating_Pitch;
FOV_Simulated=550;  %micron 223 550
Number_of_Resolvable_Point=abs(Wavelength(end)-Wavelength(1))/Wavelength_Resolution;

%
Theta_Out=zeros(length(Wavelength),1);
for p=1:length(Wavelength)
    Theta_Out(p)=asin((RI_Inc*sin(Theta_In/180*pi)+Order*Wavelength(p)/d)/RI_Out)/pi*180;
end
Divergence_Angle=abs(Theta_Out(end)-Theta_Out(1));
FOV_Estimated=abs(Focal_Length*tan(Divergence_Angle/2/180*pi))*2;
%Forward View
plot(Wavelength,Theta_Out-Theta_In);
Spectral_Spatial_Resolution=abs(Wavelength_Resolution/(Wavelength(end)-Wavelength(1))*FOV_Estimated);
Theta_Out_Mean=mean(Theta_Out);
%Littrow
plot(Wavelength,Theta_Out+Theta_In);
%Theta_D=abs(max(Theta_Out-Theta_In))-abs(min(Theta_Out-Theta_In))
Theta_D2=abs(max(Theta_Out+Theta_In))-abs(min(Theta_Out+Theta_In))
%% Image Space NA
RI=1;
R=0.52;
f=2;
NAi=RI*sin(atan(R/f));
