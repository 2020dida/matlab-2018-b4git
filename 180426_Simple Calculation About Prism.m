clear all
%%
RI_IN=1;

RI_Prism=1.85;

Theta_1=30;

Theta_2=asin((RI_IN.*sin(Theta_1/180*pi)/RI_Prism))*180/pi;


Tilt_Angle=Theta_1-Theta_2;
