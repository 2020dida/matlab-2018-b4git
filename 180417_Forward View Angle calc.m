clear all
%%
Focal_Length=2000;  %micron
Order=-1;
Wavelength=[0.5:0.01:0.6]+0.2;

Theta_In=62.5;
%Theta_Out
Cypermm=1000;
RI_Inc=1.84;
RI_Out=1;
d=1000/Cypermm;

Theta_Out=zeros(length(Wavelength),1);
for p=1:length(Wavelength)
    Theta_Out(p)=asin((RI_Inc*sin(Theta_In/180*pi)+Order*Wavelength(p)/d)/RI_Out)/pi*180;
end
%Forward View
plot(Wavelength,Theta_Out-Theta_In);

%Littrow
%plot(Wavelength,Theta_Out+Theta_In);
FOV=Focal_Length*(tan(min(Theta_Out)/180*pi)-tan(max(Theta_Out)/180*pi));

%% Image Space NA
RI=1;
R=0.52;
f=2;
NAi=RI*sin(atan(R/f));