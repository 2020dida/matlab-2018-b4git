clear all;

%%
RI=1.33;

%% Simple relation between spectral resolution and roll-off, assume k-space
Center_Wavelength=840;      %nm
Grating_Density=1800;   %line per mm

if Center_Wavelength==840
    Fiber_NA=0.13;
    Fiber_Core_Size=5;  %micron
    Spectral_Sampling_Resolution_nm=0.04;   %nm, FWHM, optical, not working, only for determination of Sampling Res
    if Grating_Density==1200
        Collimator_f=100; %mm
    elseif Grating_Density==1800
        Collimator_f=60; %mm
    end
elseif Center_Wavelength==1050
    Fiber_NA=0.14;
    Fiber_Core_Size=6.2;  %micron
    Spectral_Sampling_Resolution_nm=0.06;   %nm, FWHM, optical, not working, only for determination of Sampling Res
    if Grating_Density==1200
        Collimator_f=75; %mm
    elseif Grating_Density==1800
        Collimator_f=35; %mm
    end
end
Padding_Style=2;    %0 for no padding, 1 for pad to zero, 2 for pad to high freq
Padding_Ratio=8;
Center_Frequency=3E8/(Center_Wavelength*1E-9);
Target_Axial_Resolution_in_Tissue=10/1.35;
Spectral_FWHM_nm=0.44*((Center_Wavelength)^2)/RI/(Target_Axial_Resolution_in_Tissue*1E3);        %nm
Spectral_FWHM_Hz=3E8/((Center_Wavelength-Spectral_FWHM_nm/2)*1E-9)-3E8/((Center_Wavelength+Spectral_FWHM_nm/2)*1E-9);
Spectral_Sampling_Resolution_Hz=3E8/((Center_Wavelength-Spectral_Sampling_Resolution_nm/2)*1E-9)-3E8/((Center_Wavelength+Spectral_Sampling_Resolution_nm/2)*1E-9);
Ratio_CameravsSpectrum=2.5;
Spectrometer_Coverrange_Hz=Spectral_FWHM_Hz*Ratio_CameravsSpectrum;
Spectral_FWHM_Pixel=Spectral_FWHM_Hz/Spectral_Sampling_Resolution_Hz;
Number_of_Pixel_required=round(Spectral_FWHM_Hz/Spectral_Sampling_Resolution_Hz*Ratio_CameravsSpectrum);
Spectral_Frequency_Array=[Spectral_Sampling_Resolution_Hz:Spectral_Sampling_Resolution_Hz:Spectral_Sampling_Resolution_Hz*Number_of_Pixel_required]-Spectral_Sampling_Resolution_Hz*round(Number_of_Pixel_required/2)+round(Center_Frequency/Spectral_Sampling_Resolution_Hz)*Spectral_Sampling_Resolution_Hz;
Spectral_Frequency_Array_PadtoZero=Spectral_Sampling_Resolution_Hz:Spectral_Sampling_Resolution_Hz:Spectral_Frequency_Array(end);
Spectral_Frequency_Array_PadtoHighFreq=Spectral_Sampling_Resolution_Hz:Spectral_Sampling_Resolution_Hz:Spectral_Frequency_Array(end)*Padding_Ratio;

Spectral_Sine_period_nm=10; %nm
%Free spectral range
Spectral_Sine_period_Hz=3E8/((Center_Wavelength-Spectral_Sine_period_nm/2)*1E-9)-3E8/((Center_Wavelength+Spectral_Sine_period_nm/2)*1E-9); %nm
L=3E8/RI/Spectral_Sine_period_Hz;
l=L/2;
l_micron=l*1E6;
Spectral_Sine_period_Pixel=Spectral_Sine_period_Hz/Spectral_Sampling_Resolution_Hz;

Dummy_Spectrum=gaussmf(Spectral_Frequency_Array,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array/Spectral_Sine_period_Hz*2*pi);
Dummy_Spectrum_PadtoZero=gaussmf(Spectral_Frequency_Array_PadtoZero,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array_PadtoZero/Spectral_Sine_period_Hz*2*pi);
Dummy_Spectrum_PadtoHighFreq=gaussmf(Spectral_Frequency_Array_PadtoHighFreq,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array_PadtoHighFreq/Spectral_Sine_period_Hz*2*pi);

if Padding_Style == 0
    Ascan_Carrier=real(fft(Dummy_Spectrum));
    Ascan_Envelope=abs(fft(Dummy_Spectrum));
elseif Padding_Style == 1
    Ascan_Carrier=real(fft(Dummy_Spectrum_PadtoZero));
    Ascan_Envelope=abs(fft(Dummy_Spectrum_PadtoZero));
elseif Padding_Style == 2
    Ascan_Carrier=real(fft(Dummy_Spectrum_PadtoHighFreq));
    Ascan_Envelope=abs(fft(Dummy_Spectrum_PadtoHighFreq));
end
Number_of_Pixel_in_Position_Array=length(Ascan_Carrier);
Maximum_Position=3E8/Spectral_Sampling_Resolution_Hz/2/RI; %/2 for roundtrip
Position_SamplingResolution_micron=Maximum_Position/Number_of_Pixel_in_Position_Array*1E6;
Position_Array_micron=0:Position_SamplingResolution_micron:(Number_of_Pixel_in_Position_Array-1)*Position_SamplingResolution_micron;

plot(Position_Array_micron,Ascan_Envelope,Position_Array_micron,Ascan_Carrier);

%% Rolloff Calculation
OPD=1000; %micron
Max_Signal_Array=zeros([length(OPD) 1]);
for p=1:length(OPD)
    Spectral_Sine_period_Hz_Loop=3E8/RI/(OPD(p)/1E6*2);
    Spectral_Sine_period_Pixel=Spectral_Sine_period_Hz/Spectral_Sampling_Resolution_Hz;

    Dummy_Spectrum=gaussmf(Spectral_Frequency_Array,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array/Spectral_Sine_period_Hz_Loop*2*pi);
    Dummy_Spectrum_PadtoZero=gaussmf(Spectral_Frequency_Array_PadtoZero,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array_PadtoZero/Spectral_Sine_period_Hz_Loop*2*pi);
    Dummy_Spectrum_PadtoHighFreq=gaussmf(Spectral_Frequency_Array_PadtoHighFreq,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array_PadtoHighFreq/Spectral_Sine_period_Hz_Loop*2*pi);

    if Padding_Style == 0
        Ascan_Carrier=real(fft(Dummy_Spectrum));
        Ascan_Envelope=abs(fft(Dummy_Spectrum));
    elseif Padding_Style == 1
        Ascan_Carrier=real(fft(Dummy_Spectrum_PadtoZero));
        Ascan_Envelope=abs(fft(Dummy_Spectrum_PadtoZero));
    elseif Padding_Style == 2
        Ascan_Carrier=real(fft(Dummy_Spectrum_PadtoHighFreq));
        Ascan_Envelope=abs(fft(Dummy_Spectrum_PadtoHighFreq));
    end
    Number_of_Pixel_in_Position_Array=length(Ascan_Carrier);
    Maximum_Position=3E8/Spectral_Sampling_Resolution_Hz/2/RI; %/2 for roundtrip
    Position_SamplingResolution_micron=Maximum_Position/Number_of_Pixel_in_Position_Array*1E6;
    Position_Array_micron=0:Position_SamplingResolution_micron:(Number_of_Pixel_in_Position_Array-1)*Position_SamplingResolution_micron;
    Max_Signal_Array(p)=max(Ascan_Envelope);
    disp(p);
end
if length(OPD)>1
    plot(OPD,Max_Signal_Array);
else
    plot(Position_Array_micron,Ascan_Envelope,Position_Array_micron,Ascan_Carrier);
end
%% Some Basic Calculation
Camera_Pixel_Size=10;   %micron
Collimated_Beam_Size=Collimator_f*Fiber_NA*2;
Grating_Order=1;
Ang_Inc=asin((Grating_Order*Center_Wavelength/(1000000/Grating_Density))/2)/pi*180;

Grating_Pitch=1000/Grating_Density; %micron
Grating_Resolving_Power_Hz=Center_Frequency/(Collimated_Beam_Size*Grating_Density);
Diverging_Agnle_FWHM=asin(Grating_Order*(Center_Wavelength/1000+Spectral_FWHM_nm/2000)/Grating_Pitch-sin(Ang_Inc/180*pi))/pi*180-asin(Grating_Order*(Center_Wavelength/1000-Spectral_FWHM_nm/2000)/Grating_Pitch-sin(Ang_Inc/180*pi))/pi*180;
Diverging_Agnle_FWHM_Rad=Diverging_Agnle_FWHM/180*pi;
Required_Spectrum_Length=Spectral_FWHM_Hz/Spectral_Sampling_Resolution_Hz*Camera_Pixel_Size/1000;

Required_2nd_Focal_for_Spectrum_Span=Required_Spectrum_Length/2/tan(Diverging_Agnle_FWHM_Rad/2);

Optical_Resolving_Power_Hz=(Fiber_Core_Size/Collimator_f*Required_2nd_Focal_for_Spectrum_Span)/Camera_Pixel_Size*Spectral_Sampling_Resolution_Hz;
Total_Resolving_Power_Hz=(Grating_Resolving_Power_Hz^2+Optical_Resolving_Power_Hz^2)^0.5;
Total_Resolving_Power_Hz_nm=(3E8/((Center_Frequency-Total_Resolving_Power_Hz/2))-3E8/((Center_Frequency+Total_Resolving_Power_Hz/2)))*1E9;

%% Rolloff with finite optical resolution
OPD=10:10:3000; %micron
Max_Signal_Array=zeros([length(OPD) 1]);
for p=1:length(OPD)
    Spectral_Sine_period_Hz_Loop=3E8/RI/(OPD(p)/1E6*2);
    Spectral_Sine_period_Pixel=Spectral_Sine_period_Hz/Spectral_Sampling_Resolution_Hz;

    Dummy_Spectrum=gaussmf(Spectral_Frequency_Array,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array/Spectral_Sine_period_Hz_Loop*2*pi);
    Dummy_Spectrum_PadtoZero=gaussmf(Spectral_Frequency_Array_PadtoZero,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array_PadtoZero/Spectral_Sine_period_Hz_Loop*2*pi);
    Dummy_Spectrum_PadtoHighFreq=gaussmf(Spectral_Frequency_Array_PadtoHighFreq,[Spectral_FWHM_Hz/(2*sqrt(2*log(2))) Center_Frequency]).*sin(Spectral_Frequency_Array_PadtoHighFreq/Spectral_Sine_period_Hz_Loop*2*pi);
    Spectral_PSF=gaussmf(1:(round(Total_Resolving_Power_Hz/Spectral_Sampling_Resolution_Hz*3/2)*2+1),[Total_Resolving_Power_Hz/Spectral_Sampling_Resolution_Hz/(2*sqrt(2*log(2))) round(Total_Resolving_Power_Hz/Spectral_Sampling_Resolution_Hz*3/2)+1]);

    if Padding_Style == 0
        Ascan_Carrier=real(fft(conv(Dummy_Spectrum,Spectral_PSF,'same')));
        Ascan_Envelope=abs(fft(conv(Dummy_Spectrum,Spectral_PSF,'same')));
    elseif Padding_Style == 1
        Ascan_Carrier=real(fft(conv(Dummy_Spectrum_PadtoZero,Spectral_PSF,'same')));
        Ascan_Envelope=abs(fft(conv(Dummy_Spectrum_PadtoZero,Spectral_PSF,'same')));
    elseif Padding_Style == 2
        Ascan_Carrier=real(fft(conv(Dummy_Spectrum_PadtoHighFreq,Spectral_PSF,'same')));
        Ascan_Envelope=abs(fft(conv(Dummy_Spectrum_PadtoHighFreq,Spectral_PSF,'same')));
    end
    Number_of_Pixel_in_Position_Array=length(Ascan_Carrier);
    Maximum_Position=3E8/Spectral_Sampling_Resolution_Hz/2/RI; %/2 for roundtrip
    Position_SamplingResolution_micron=Maximum_Position/Number_of_Pixel_in_Position_Array*1E6;
    Position_Array_micron=0:Position_SamplingResolution_micron:(Number_of_Pixel_in_Position_Array-1)*Position_SamplingResolution_micron;
    Max_Signal_Array(p)=max(Ascan_Envelope);
    disp(p);
end
if length(OPD)>1
    plot(OPD,Max_Signal_Array);
    plot(OPD,20*log10(Max_Signal_Array)-20*log10(Max_Signal_Array(1)));
else
    plot(Position_Array_micron,Ascan_Envelope,Position_Array_micron,Ascan_Carrier);
end

%% Sensitivty Rolloff in Literature
% J.A. Izatt
zmax=((Center_Wavelength^2)/4/Spectral_Sampling_Resolution_nm/RI)*1E-6;  %mm
z6dB=log(2)/pi/Total_Resolving_Power_Hz_nm*(Center_Wavelength^2)*1E-6;

Number_of_Pixel_required
zmax
z6dB
Spectral_Sampling_Resolution_nm
Total_Resolving_Power_Hz_nm
Collimated_Beam_Size
Ang_Inc
Spectral_FWHM_nm
Required_2nd_Focal_for_Spectrum_Span
%% 6dB