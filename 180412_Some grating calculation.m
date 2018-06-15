clear all
%%
Length=10;
Grating_Order=1;
Center_Wavelength=780;
Grating_Density=1300; %cycle/mm
Spectral_FWHM_nm=200;
Ang_Inc=asin((Grating_Order*Center_Wavelength/(1000000/Grating_Density))/2)/pi*180;

Grating_Pitch=1000/Grating_Density; %micron
Diverging_Agnle_FWHM=asin(Grating_Order*(Center_Wavelength/1000+Spectral_FWHM_nm/2000)/Grating_Pitch-sin(Ang_Inc/180*pi))/pi*180-asin(Grating_Order*(Center_Wavelength/1000-Spectral_FWHM_nm/2000)/Grating_Pitch-sin(Ang_Inc/180*pi))/pi*180;
Diverging_Agnle_FWHM_Rad=Diverging_Agnle_FWHM/180*pi;
FOV=sin(Diverging_Agnle_FWHM/2/180*pi)*Length*2;