clear all
%%
Order=1;
Wavelength=0.8;
Cypermm=1144;
RI_Inc=1.5;
RI_Out=1;

d=1000/Cypermm;
Theta_L=asin(Order*Wavelength/d/(RI_Inc+RI_Out))/pi*180; %�J�g��=�Ϯg�� + �ʶq�u��