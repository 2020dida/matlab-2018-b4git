clear all

Slow_Ratio=128;



Old_a0=10;

Old_a1=1.68E-6;

Old_a2=-4.4;

Old_aex=1.199;


New_a0=Old_a0

New_a1=Old_a1*((1/Slow_Ratio)^Old_aex)

New_a2=Old_a2*(1/Slow_Ratio)


New_aex=Old_aex