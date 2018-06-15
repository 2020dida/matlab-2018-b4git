clear all
cd('C:\TuanShu\MATLAB\TSLib');
%% Hardware-related Parameters
Bit=12;         
FWC_e=200000;
X_Spectral_Pixel_Number=2048*8;

Center_Wavelengh=0.84;  %micron
Wavelength_FWHM=0.055;   %micron
%% Dirived Parameter (Hardware)
% Camera
DN_max=2^Bit;
Gain=DN_max/FWC_e;
% Light source
Frequency_FWHM=3E8*Wavelength_FWHM/(Center_Wavelengh^2)/(1E-6);
Calculated_Axial_Resolution=3E8/Frequency_FWHM/2*1E6;   %micron
Center_Frequency=3E8/(Center_Wavelengh*1E-6);
Pixel_Number_Array=1:X_Spectral_Pixel_Number;

%% System-related Parameters
IE=0.15;
Spectrum_Sampling_Ratio=0.5;    %繵眯FWHMcameraぇゑㄒ
%% Dirived Parameter (System)
Frequency_Sampling_Resolution=Frequency_FWHM/Spectrum_Sampling_Ratio/X_Spectral_Pixel_Number;
Frequency_Array=Frequency_Sampling_Resolution*((-1*X_Spectral_Pixel_Number/2+1):X_Spectral_Pixel_Number/2)+round(Center_Frequency/Frequency_Sampling_Resolution)*Frequency_Sampling_Resolution;
Wavelength_Array=3E8./Frequency_Array*1E6;
Central_Wavelength_Index=find(Wavelength_Array<Center_Wavelengh,1,'first');
%% Measurement-related Parameters
Sturation_Level=0.9;
OPD=2500;  %micron
%% Dirived Parameter (Measurement)
Fringe_Period_frequency=3E8/((OPD*1E-6)*2);

%% Calculation-related Parameter
ZP_Ratio=16;
%% Dirived Parameter (Calculation)
% Spectral
Max_Frequency_ZP=round((ZP_Ratio*max(Frequency_Array))/Frequency_Sampling_Resolution/2)*2*Frequency_Sampling_Resolution+Frequency_Sampling_Resolution; % 絋﹚计
Full_Frequency_Array=[0:Frequency_Sampling_Resolution:Max_Frequency_ZP]'; %絋﹚案计bin, 安砞ノ场spectrum, ぃノ疭よ猭
Min_Frequency_Index=find(Full_Frequency_Array>min(Frequency_Array),1,'first');

% Temporal
Time_max=1/Frequency_Sampling_Resolution;    %sec
Temporal_Resolution=Time_max/length(Full_Frequency_Array);
Time=[Temporal_Resolution:Temporal_Resolution:Time_max];
Position_Full=3E8*Time/2;  %/2 for round-trip
Position=Position_Full(1:round(size(Full_Frequency_Array)/2));
Position_micron=Position*1E6;
Position_Sampling_Resolution_micron=Position_micron(2)-Position_micron(1);
Index_Theoritical_Peak_Signal=round(OPD/Position_Sampling_Resolution_micron);
%% Simulating Spectrum
% DC
DC_Spectrum_e=gaussmf(Frequency_Array,[Frequency_FWHM/2.355 Center_Frequency]).*FWC_e.*Sturation_Level;
DC_Spectrum_dn=DC_Spectrum_e*Gain;    %?????????Ave spectrum?, ????round
subplot(2,1,1)
plot(Frequency_Array/1E12,DC_Spectrum_e);
xlabel('Optical Frequency (THz)');
ylabel('Number of Electron per Pixel');
xlim([min(Frequency_Array/1E12) max(Frequency_Array/1E12)]);
subplot(2,1,2)
plot(Pixel_Number_Array,DC_Spectrum_dn);
xlabel('Pixel Number');
ylabel('Digital Number (DN)');
xlim([min(Pixel_Number_Array) max(Pixel_Number_Array)]);
% Inter
Inter_Spectrum_e=DC_Spectrum_e*IE.*cos((Frequency_Array-Center_Frequency)/Fringe_Period_frequency*2*pi);
Complete_Spectrum_e=DC_Spectrum_e+Inter_Spectrum_e;
Complete_Spectrum_dn=Complete_Spectrum_e*Gain;
subplot(2,1,1)
plot(Frequency_Array/1E12,Complete_Spectrum_e);
xlabel('Optical Frequency (THz)');
ylabel('Number of Electron per Pixel');
xlim([min(Frequency_Array/1E12) max(Frequency_Array/1E12)]);
subplot(2,1,2)
plot(Pixel_Number_Array,Complete_Spectrum_e);
xlabel('Pixel Number');
ylabel('Digital Number (DN)');
xlim([min(Pixel_Number_Array) max(Pixel_Number_Array)]);
%% Sub DC
Complete_Spectrum_DC_removed_dn=Complete_Spectrum_dn-DC_Spectrum_dn;
%%
S_Full=TSZPFFT(Complete_Spectrum_DC_removed_dn,Full_Frequency_Array,Min_Frequency_Index);
S=S_Full(1:round(length(S_Full)/2));

plot(Position_micron,abs(S),Position_micron,real(S));
xlabel('Depth (micron)');
ylabel('Digital Number (DN)');
xlim([450 550]);
%% Generate Complex Spectrum
S_Temp=S_Full;
S_Temp(((length(S_Temp)/2+1):end))=0;
Full_Complex_Spectrum_dn=ifft(S_Temp)*2;
Complex_Spectrum_dn=Full_Complex_Spectrum_dn(Min_Frequency_Index:(Min_Frequency_Index+length(Complete_Spectrum_DC_removed_dn)-1));

%% Change to Linear-to-Wavelength Sensor Array
Sensor_Spectral_Pixel_Number=2048;
Sensor_Spectrum_Sampling_Ratio=0.5;
Sensor_Wavelength_Sampling_Resolution=Wavelength_FWHM/Sensor_Spectrum_Sampling_Ratio/Sensor_Spectral_Pixel_Number;
Sensor_Wavelength_Array=Sensor_Wavelength_Sampling_Resolution*((-1*Sensor_Spectral_Pixel_Number/2+1):Sensor_Spectral_Pixel_Number/2)+round(Center_Wavelengh/Sensor_Wavelength_Sampling_Resolution)*Sensor_Wavelength_Sampling_Resolution;

Sensor_Central_Wavelength_Index=find(Sensor_Wavelength_Array<Center_Wavelengh,1,'first');

plot(Sensor_Wavelength_Array,'.');

%% Calculate the "Measured" Spectrum
Optical_Wavelength_Resolution=0.02*1E-3; %micron
Optical_Frequency_Resolution=Optical_Wavelength_Resolution/Center_Wavelengh*Center_Frequency;
Sigma=Optical_Frequency_Resolution/2/(2*log(2))^0.5;    %Used for Filtering, ぃノ
Input_Spectrum=real(Complex_Spectrum_dn);    %Varied_Comepensated_Spectrum_dn, Comepensated_Spectrum_dn, Dispersed_Spectrum_dn, Complex_Spectrum_dn
Sensor_Spectrum_AVE=zeros(Sensor_Spectral_Pixel_Number,1);
Sensor_Spectrum_noAVE=zeros(Sensor_Spectral_Pixel_Number,1);

for p=1:Sensor_Spectral_Pixel_Number
    Index_min=find(Wavelength_Array<(Sensor_Wavelength_Array(p)-Sensor_Wavelength_Sampling_Resolution/2),1,'first');
    Index_max=find(Wavelength_Array<(Sensor_Wavelength_Array(p)+Sensor_Wavelength_Sampling_Resolution/2),1,'first');
    Index_center=find(Wavelength_Array<(Sensor_Wavelength_Array(p)),1,'first');
    if ~isempty(Index_min) & ~isempty(Index_max)
        Sensor_Spectrum_AVE(p)=mean(Input_Spectrum(min(Index_min,Index_max):max(Index_min,Index_max)));
    else
        Sensor_Spectrum_AVE(p)=0;
    end
    if ~isempty(Index_center)
        Sensor_Spectrum_noAVE(p)=Input_Spectrum(Index_center);
    else
        Sensor_Spectrum_noAVE(p)=0;
    end
    disp(p);
end
subplot(1,1,1)
plot(Sensor_Wavelength_Array,Sensor_Spectrum_AVE,Sensor_Wavelength_Array,Sensor_Spectrum_noAVE,Wavelength_Array,Input_Spectrum);
%% SNR Calculation for SD-OCT
plot(Signal_SDOCT_Array);

Signal_SDOCT=mean(Signal_SDOCT_Array);
Noise_SDOCT=std(Signal_SDOCT_Array);
SNR_SDOCT=20*log10(Signal_SDOCT/Noise_SDOCT);


%% Speed Related Calculation
% Calculate the number of electrons required for above TD and SD schme
Total_Electron_SDOCT=sum(DC_Spectrum_e)*X_number;
Total_Electron_TDOCT=FWC_e.*Sturation_Level*X_number;

FrameRate_TDOCT=4000;   %Partial Frame, 1024*11, fps
FrameRate_SDOCT=FrameRate_TDOCT/Total_Electron_SDOCT*Total_Electron_TDOCT;    %Full Frame, 1024*1024, fps

Depth_per_Frame_SDOCT=Center_Wavelengh/2;    %micron
Depth_per_Frame_TDOCT=Center_Wavelengh/2;   %micron

Target_Depth=400;   %micron

Frame_Required_for_Target_Depth_SDOCT=Target_Depth/Depth_per_Frame_SDOCT;
Frame_Required_for_Target_Depth_TDOCT=Target_Depth/Depth_per_Frame_TDOCT;

Time_Required_for_Target_Depth_SDOCT=Frame_Required_for_Target_Depth_SDOCT/FrameRate_SDOCT;
Time_Required_for_Target_Depth_TDOCT=Frame_Required_for_Target_Depth_TDOCT/FrameRate_TDOCT;

Avaliable_Averaging_TDOCT=round(Time_Required_for_Target_Depth_SDOCT/Time_Required_for_Target_Depth_TDOCT);
% Ave Factor = 8; ?????Include 4-point and PreAVE???(????PostAVE, ???Depth
% per Frame??????PostAVE)


%% SD-OCT (Roll-off optimization, reducing Spectral Resolution)
Spectral_AVE_Factor=4;
Frequency_Sampling_Resolution_SRreduced=Frequency_Sampling_Resolution*Spectral_AVE_Factor;
Full_Frequency_SRreduced_Array=[0:Frequency_Sampling_Resolution_SRreduced:ZP_Ratio*max(Frequency_Array)]';
Min_Frequency_SRreduced_Index=find(Full_Frequency_SRreduced_Array>min(Frequency_Array),1,'first');

%% Generating Position Array (Roll-offed)
Time_max_SRreduced=1/Frequency_Sampling_Resolution_SRreduced;    %sec
Temporal_Resolution_SRreduced=Time_max_SRreduced/length(Full_Frequency_SRreduced_Array);%same
Time_SRreduced=[Temporal_Resolution_SRreduced:Temporal_Resolution_SRreduced:Time_max_SRreduced];
Position_Full_SRreduced=3E8*Time_SRreduced/2;  %/2 for round-trip
Position_SRreduced=Position_Full_SRreduced(1:round(size(Full_Frequency_SRreduced_Array)/2));
Position_micron_SRreduced=Position_SRreduced*1E6;
Position_Sampling_Resolution_micron_SRreduced=Position_micron_SRreduced(2)-Position_micron_SRreduced(1);
Index_Theoritical_Peak_Signal_SRreduced=round(OPD/Position_Sampling_Resolution_micron_SRreduced);

Signal_SDOCT_SRreduced_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    % Add Noise
    Shot_Noise_Spectrum_e=randn(X_Spectral_Pixel_Number,1)'.*Complete_Spectrum_e.^0.5;
    Spectrum_e_with_Noise=Complete_Spectrum_e+Shot_Noise_Spectrum_e;
    Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
    plot(Frequency_Array,Spectrum_dn_with_Noise);
    % Reducing Spectral Resolution
    Reduced_Spectrum_Width=floor(X_Spectral_Pixel_Number/Spectral_AVE_Factor);
    Temp_Spectrum=0;
    Temp_Spectrum_DC=0;
    for q=1:Spectral_AVE_Factor
        Temp_Spectrum=Temp_Spectrum+Spectrum_dn_with_Noise(q:Spectral_AVE_Factor:(Reduced_Spectrum_Width*Spectral_AVE_Factor-Spectral_AVE_Factor+q));
        Temp_Spectrum_DC=Temp_Spectrum_DC+DC_Spectrum_dn(q:Spectral_AVE_Factor:(Reduced_Spectrum_Width*Spectral_AVE_Factor-Spectral_AVE_Factor+q));
    
    end
    Spectrum_dn_with_Noise_SRreduced=Temp_Spectrum/Spectral_AVE_Factor;
    Spectrum_dn_SRreduced=Temp_Spectrum_DC/Spectral_AVE_Factor;

    %% SD-OCT Calculation (Tranditional)
    Full_Spectrum_SRreduced=zeros(size(Full_Frequency_SRreduced_Array,1),size(Full_Frequency_SRreduced_Array,2));
    Full_Spectrum_SRreduced(Min_Frequency_SRreduced_Index:(Min_Frequency_SRreduced_Index+length(Spectrum_dn_with_Noise_SRreduced)-1))=Spectrum_dn_with_Noise_SRreduced-Spectrum_dn_SRreduced;
    % plot(Full_Frequency_Array,Full_Spectrum);

    S_Full_SRreduced=fft(Full_Spectrum_SRreduced);
    S_SRreduced=S_Full_SRreduced(1:round(size(S_Full_SRreduced)/2));
    % 
    % plot(abs(S))
    % xlim([50 300]/8*ZP_Ratio)
    % ylim([-4E4 4E4])
    
    % Record Signal
    Signal_SDOCT_SRreduced_Array(p)=abs(S_SRreduced(Index_Theoritical_Peak_Signal_SRreduced));
    disp(p);
end
%% Plot Signal

plot(Position_micron,abs(S),Position_micron,real(S),Position_micron_SRreduced,abs(S_SRreduced),Position_micron_SRreduced,real(S_SRreduced))
xlim([OPD-5 OPD+10])
ylim([-4E4 4E4])

%% SNR Calculation for SD-OCT (SRreduced)
plot(Signal_SDOCT_SRreduced_Array);

Signal_SDOCT_SRreduced=mean(Signal_SDOCT_SRreduced_Array);
Noise_SDOCT_SRreduced=std(Signal_SDOCT_SRreduced_Array);
SNR_SDOCT_SRreduced=20*log10(Signal_SDOCT_SRreduced/Noise_SDOCT_SRreduced);


%% SD-OCT (Roll-off optimization, Pixel smoothing)

Signal_SDOCT_Smooth_Array=zeros(Repeat_Times,1);
for p=1:Repeat_Times
    % Add Noise
    Shot_Noise_Spectrum_e=randn(X_Spectral_Pixel_Number,1)'.*Complete_Spectrum_e.^0.5;
    Spectrum_e_with_Noise=Complete_Spectrum_e+Shot_Noise_Spectrum_e;
    Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
    plot(Frequency_Array,Spectrum_dn_with_Noise);
    % Smoothing
    Spectrum_dn_with_Noise_Smooth=smooth(Spectrum_dn_with_Noise,Spectral_AVE_Factor);
    % SD-OCT Calculation (Tranditional)
    Full_Spectrum_Smooth=zeros(size(Full_Frequency_Array,1),size(Full_Frequency_Array,2));
    Full_Spectrum_Smooth(Min_Frequency_Index:(Min_Frequency_Index+length(Spectrum_dn_with_Noise)-1))=Spectrum_dn_with_Noise-DC_Spectrum_dn;
    % plot(Full_Frequency_Array,Full_Spectrum);

    S_Full_Smooth=fft(Full_Spectrum_Smooth);
    S_Smooth=S_Full_Smooth(1:round(size(S_Full_Smooth)/2));
    % 
    % plot(abs(S))
    % xlim([50 300]/8*ZP_Ratio)
    % ylim([-4E4 4E4])
    
    % Record Signal
    Signal_SDOCT_Smooth_Array(p)=abs(S_Smooth(Index_Theoritical_Peak_Signal));
    disp(p);
end
%%
plot(Position_micron,abs(S_Smooth),Position_micron,real(S_Smooth))
xlim([OPD-5 OPD+10])
ylim([-4E4 4E4])

%% SNR Calculation for SD-OCT (SRreduced)

Signal_SDOCT_Smooth=mean(Signal_SDOCT_Smooth_Array);
Noise_SDOCT_Smooth=std(Signal_SDOCT_Smooth_Array);
SNR_SDOCT_Smooth=20*log10(Signal_SDOCT_Smooth/Noise_SDOCT_Smooth);

%% Final SNR Comparision
SNR_diff_final=SNR_SDOCT_SRreduced-SNR_TDOCT_AVE_BF
Ratio_final=10^(SNR_diff_final/20)
