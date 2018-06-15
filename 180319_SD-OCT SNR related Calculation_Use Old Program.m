clear all
%% Mode
Mode='SD';  %SD or TD
%% System Parameter
FWC_e=200000;
Bit=12;
DN_max=2^Bit;
Dark_Noise_dn=0.5*2^4; %DN @12bit
Gain=DN_max/FWC_e;
X_number=1024;
Y_number=1;
Sturation_Level=0.9;
DN_mean=DN_max*Sturation_Level;
PreAVE=2;
Npoint_Corr=0.5;   %experiemental unbiased estimator for N=4
Shot_e=((PreAVE*FWC_e*Sturation_Level)^0.5)/PreAVE;
Shot_dn=Shot_e*Gain;
Shot_PostN_dn=Shot_dn*Npoint_Corr;
IE=0.02;
OPD=5;  %micron
Fringe_Period_frequency=3E8/((OPD*1E-6)*2);

%% Simulation Setting
Repeat_Times=2000;  %for Noise Calculation
%% Spectrum Parameter
X_spectrum_number=1024;
Center_Wavelengh=0.84;  %micron
Wavelength_FWHM=0.05;   %micron
Spectrum_Sampling_Ratio=1;    %?Spectrum?FWHM????????
ZP_Ratio=8;
%% SD-OCT Part
%% Simulating Spectrum (Linear with k, simplified) without Noise
Frequency_FWHM=3E8*Wavelength_FWHM/(Center_Wavelengh^2)/(1E-6);
Calculated_Axial_Resolution=3E8/Frequency_FWHM/2*1E6;   %micron
Center_Frequency=3E8/(Center_Wavelengh*1E-6);
Frequency_Sampling_Resolution=Frequency_FWHM/Spectrum_Sampling_Ratio/X_spectrum_number;

Frequency_Array=Frequency_Sampling_Resolution*((-1*X_spectrum_number/2+1):X_spectrum_number/2)+round(Center_Frequency/Frequency_Sampling_Resolution)*Frequency_Sampling_Resolution;
% DC
Spectrum_e=gaussmf(Frequency_Array,[Frequency_FWHM/2.355 Center_Frequency]).*FWC_e.*Sturation_Level;
Spectrum_dn=Spectrum_e*Gain;    %?????????Ave spectrum?, ????round
% Inter
Inter_Spectrum_e=Spectrum_e*IE.*cos((Frequency_Array-Center_Frequency)/Fringe_Period_frequency*2*pi);
Spectrum_with_Inter_e=Spectrum_e+Inter_Spectrum_e;
% Full Frequency Range
Full_Frequency_Array=[0:Frequency_Sampling_Resolution:ZP_Ratio*max(Frequency_Array)]';
Min_Frequency_Index=find(Full_Frequency_Array>min(Frequency_Array),1,'first');
%% Generating Position Array
Time_max=1/Frequency_Sampling_Resolution;    %sec
Temporal_Resolution=Time_max/length(Full_Frequency_Array);
Time=[Temporal_Resolution:Temporal_Resolution:Time_max];
Position_Full=3E8*Time/2;  %/2 for round-trip
Position=Position_Full(1:round(size(Full_Frequency_Array)/2));
Position_micron=Position*1E6;
Position_Sampling_Resolution_micron=Position_micron(2)-Position_micron(1);
Index_Theoritical_Peak_Signal=round(OPD/Position_Sampling_Resolution_micron);

%% Simulating Spectrum (Linear with k, simplified) with Noise
Signal_SDOCT_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    % Noise
    Shot_Noise_Spectrum_e=randn(X_spectrum_number,1)'.*Spectrum_with_Inter_e.^0.5;

    Spectrum_e_with_Noise=Spectrum_with_Inter_e+Shot_Noise_Spectrum_e;
    Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
    % plot(Frequency_Array,Spectrum_dn_with_Noise);
    %% SD-OCT Calculation (Tranditional)
    Full_Spectrum=zeros(size(Full_Frequency_Array,1),size(Full_Frequency_Array,2));
    Full_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Spectrum_dn_with_Noise)-1))=Spectrum_dn_with_Noise-Spectrum_dn;
    % plot(Full_Frequency_Array,Full_Spectrum);

    S_Full=fft(Full_Spectrum);
    S=S_Full(1:round(size(S_Full)/2));
    % 
    % plot(abs(S))
    % xlim([50 300]/8*ZP_Ratio)
    % ylim([-4E4 4E4])
    
    % Record Signal
    Signal_SDOCT_Array(p)=abs(S(Index_Theoritical_Peak_Signal));
    disp(p);
end
%% Plot Signal

plot(Position_micron,abs(S),Position_micron,real(S))
xlim([OPD-5 OPD+10])
ylim([-4E4 4E4])

%% SNR Calculation for SD-OCT
plot(Signal_SDOCT_Array);

Signal_SDOCT=mean(Signal_SDOCT_Array);
Noise_SDOCT=std(Signal_SDOCT_Array);
SNR_SDOCT=20*log10(Signal_SDOCT/Noise_SDOCT);
