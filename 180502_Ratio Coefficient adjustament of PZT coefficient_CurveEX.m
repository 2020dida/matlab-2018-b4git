clear all

Slow_Ratio=0.5;



Old_a0=0.00;

Old_a1=3.78E-03;

Old_a2=1.72E-05;

Old_aex=0.500;


New_a0=Old_a0

New_a1=Old_a1*((1/Slow_Ratio)^Old_aex)

New_a2=Old_a2*(1/Slow_Ratio)


New_aex=Old_aex