clear all
cd('C:\TuanShu\MATLAB\TSLib');
%% Hardware-related Parameters
Bit=12;         
FWC_e=200000;
X_Spectral_Pixel_Number=2048;

Center_Wavelengh=0.84;  %micron
Wavelength_FWHM=0.050;   %micron
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
Spectrum_Sampling_Ratio=0.5;    %只頻譜FWHM占camera之比例
%% Dirived Parameter (System)
Frequency_Sampling_Resolution=Frequency_FWHM/Spectrum_Sampling_Ratio/X_Spectral_Pixel_Number;
Frequency_Array=Frequency_Sampling_Resolution*((-1*X_Spectral_Pixel_Number/2+1):X_Spectral_Pixel_Number/2)+round(Center_Frequency/Frequency_Sampling_Resolution)*Frequency_Sampling_Resolution;
Wavelength_Array=3E8./Frequency_Array*1E6;
Central_Wavelength_Index=find(Wavelength_Array<Center_Wavelengh,1,'first');
%% Measurement-related Parameters
Sturation_Level=0.9;
OPD=500;  %micron
%% Dirived Parameter (Measurement)
Fringe_Period_frequency=3E8/((OPD*1E-6)*2);

%% Calculation-related Parameter
ZP_Ratio=256;
%% Dirived Parameter (Calculation)
% Spectral
Max_Frequency_ZP=round((ZP_Ratio*max(Frequency_Array))/Frequency_Sampling_Resolution/2)*2*Frequency_Sampling_Resolution+Frequency_Sampling_Resolution; % 確定為奇數倍
Full_Frequency_Array=[0:Frequency_Sampling_Resolution:Max_Frequency_ZP]'; %確定為偶數bin, 先假設用全部的spectrum, 不用特殊方法
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

%% Zero Padding
Full_Spectrum_dn=zeros(length(Full_Frequency_Array),1);
Full_Spectrum_dn(Min_Frequency_Index:(Min_Frequency_Index+length(Complete_Spectrum_dn)-1))=Complete_Spectrum_dn-DC_Spectrum_dn;

%
subplot(1,1,1)
plot(Full_Frequency_Array/1E12,Full_Spectrum_dn);
xlabel('Optical Frequency (THz)');
ylabel('Digital Number (DN)');
xlim([min(Full_Frequency_Array/1E12) max(Full_Frequency_Array/1E12)]);

%% Fourier Transform to TD
S_Full=fft(Full_Spectrum_dn);
S=S_Full(1:length(Position_micron));
plot(Position_micron,abs(S),Position_micron,real(S))
xlabel('Depth (micron)');
ylabel('Digital Number (DN)');
xlim([450 550]);
%% Generate Complex Spectrum
S_Temp=S_Full;
S_Temp(((length(S_Temp)/2+1):end))=0;
Complex_Spectrum_dn=ifft(S_Temp)*2;

subplot(2,1,1)
plot(Full_Frequency_Array/1E12,Full_Spectrum_dn);
xlabel('Optical Frequency (THz)');
ylabel('Digital Number (DN)');
xlim([min(Full_Frequency_Array/1E12) max(Full_Frequency_Array/1E12)]);

subplot(2,1,2)
plot(Full_Frequency_Array/1E12,Complex_Spectrum_dn);
xlabel('Optical Frequency (THz)');
ylabel('Digital Number (DN)');
xlim([min(Full_Frequency_Array/1E12) max(Full_Frequency_Array/1E12)]);

%% Plot Material Dispersion
%RI=TSRI(Wavelength_Array,'S-TIH53')
subplot(1,1,1)
plot(Wavelength_Array,TSRI(Wavelength_Array,'Silica'));
xlim([Wavelength_Array(end) Wavelength_Array(1)]);
%% Plot Total Material Dispersion Difference (以Center_Wavelengh為基準)
%Varying_Thickness_Array=([28.03 22.615 17.2 14.28 11.36]-17.2)*1E3; %micron
Material_List={'S-TIH53' 'L-LAL13' 'Water'};% 'S-NPH2' 'Silica'}; %Cell
Thickness_List=[5 14.93 0.6+3+4+17.2];% 0 0];%[0.6+3+4+17.2 14.93 5 0 0];%[0.6+3+4+17.2 14.93 5 3.42 0];    %array, mm, single-trip thickness
%Material_List={'L-LAL13'}; %Cell
%Thickness_List=[14.93];%[0.6+3+4+17.2 14.93 5 0 0];%[0.6+3+4+17.2 14.93 5 3.42 0];    %array, mm, single-trip thickness

OPD_Spectrum=zeros(1,length(Wavelength_Array)); %Unit: micron
for p=1:length(Material_List)
    OPD_Spectrum=OPD_Spectrum+TSRI(Wavelength_Array,Material_List{p})*Thickness_List(p)*1E3*2;  %micron, roundtrip
end
Phase_Spectrum=OPD_Spectrum./Wavelength_Array*2*pi;

plot(Frequency_Array,Phase_Spectrum); 
%Q=[Frequency_Array' ones(length(Frequency_Array),1)]
%Q2=Frequency_Array';
%beta=Q\Phase_Spectrum';
% To find the suitable OPD(air) for reference arm to compensate
FitResult = fit(Frequency_Array',Phase_Spectrum','poly1');
Phase_Spectrum_Modi=Phase_Spectrum-Frequency_Array.*FitResult.p1+FitResult.p2;      %這步減掉了2個東西: 1. constant phase (本來就沒差) 2. 前述List中之原件對干涉訊號位置之影響 (假設可用reference arm調回來)
Phase_Spectrum_Modi=Phase_Spectrum_Modi-Phase_Spectrum_Modi(Central_Wavelength_Index);  

plot(Wavelength_Array,Phase_Spectrum_Modi);
xlabel('Optical Wavelength (\mum)');
ylabel('Residual Phase (Radian)');
xlim([min(Wavelength_Array) max(Wavelength_Array)]);
xlim([min(Wavelength_Array) max(Wavelength_Array)]);

%% Zero Padding Phase Spectrum

Full_Phase_Spectrum=zeros(length(Full_Frequency_Array),1);
Full_Phase_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Complete_Spectrum_dn)-1))=Phase_Spectrum_Modi;

%% Generate Dispersed Spectrum
Dispersed_Spectrum_dn=Complex_Spectrum_dn.*exp(i*Full_Phase_Spectrum);

subplot(1,1,1)
plot(Full_Frequency_Array/1E12,Dispersed_Spectrum_dn);
xlabel('Optical Frequency (THz)');
ylabel('Digital Number (DN)');
xlim([min(Full_Frequency_Array/1E12) max(Full_Frequency_Array/1E12)]);


%% Diserpsed PSF

S_Dispersed_Full=fft(Dispersed_Spectrum_dn);
S_Dispersed=S_Dispersed_Full(1:length(Position_micron));

subplot(2,1,1)
plot(Position_micron-OPD,abs(S),Position_micron-OPD,real(S))
xlabel('Depth (micron)');
ylabel('Digital Number (DN)');
xlim([-100 100]);

subplot(2,1,2)
plot(Position_micron-OPD,abs(S_Dispersed),Position_micron-OPD,real(S_Dispersed))
xlabel('Depth (micron)');
ylabel('Digital Number (DN)');
xlim([-100 100]);

%% Dispersion Compensation with Single Material
DC_Material={'Silica' 'L-LAL13' 'S-TIH53' 'S-NPH2'};
Thickness_Tol=250;   %micron, +/-  250
delta_Thickness=10;%10
DC_Thickness_Array=[-1*Thickness_Tol:delta_Thickness:Thickness_Tol];
PSF_FWHM_Matrix=zeros(length(DC_Thickness_Array),length(DC_Material));
DC_Thickness_Matrix=zeros(length([-1*Thickness_Tol:delta_Thickness:Thickness_Tol]),length(DC_Material));
for q=1:length(DC_Material)
    % 先用Phase spectrum的2nd order derivative找接近的thickness, 然後在附近搜尋
    Phase_Spectrum_Modi_2d=diff(Phase_Spectrum_Modi,2);
    %Phase_Spectrum_Modi_2d=[Phase_Spectrum_Modi_2d(1) Phase_Spectrum_Modi_2d(1) Phase_Spectrum_Modi_2d];
    DCPhase_UnityLength_2d=diff(TSRI(Wavelength_Array,DC_Material{q})./Wavelength_Array*2*pi*2,2); %Unity for ROUNDTRIP!!

%     subplot(2,1,1)
%     plot(Wavelength_Array(3:end),Phase_Spectrum_Modi_2d);
%     xlabel('Optical Wavelength (\mum)');
%     ylabel('2nd Derivative of Phase (vs Frequency)');
%     xlim([min(Wavelength_Array) max(Wavelength_Array)]);
%     xlim([min(Wavelength_Array) max(Wavelength_Array)]);
% 
%     subplot(2,1,2)
%     plot(Wavelength_Array(3:end),DCPhase_UnityLength_2d);
%     xlabel('Optical Wavelength (\mum)');
%     ylabel('2nd Derivative of Phase (vs Frequency)');
%     xlim([min(Wavelength_Array) max(Wavelength_Array)]);
%     xlim([min(Wavelength_Array) max(Wavelength_Array)]);
% 
     DC_Thickness_Est(q)=DCPhase_UnityLength_2d'\Phase_Spectrum_Modi_2d';
% 
%     %%
%     Coef=1;
%     subplot(1,1,1)
%     plot(Frequency_Array(3:end),Phase_Spectrum_Modi_2d,Frequency_Array(3:end),DCPhase_UnityLength_2d*DC_Thickness_Est(q)*Coef);
%     xlabel('Optical Frequency (THz)');
%     ylabel('2nd Derivative of Phase (vs Frequency)');
%     xlim([min(Frequency_Array) max(Frequency_Array)]);
%     xlim([min(Frequency_Array) max(Frequency_Array)]);

    %% Dispersion Compensation with Selected DC
    DC_Thickness_Matrix(:,q)=DC_Thickness_Array;
    for p=1:length(DC_Thickness_Array)
        DC_Phase_Spectrum=TSRI(Wavelength_Array,DC_Material{q}).*(DC_Thickness_Est(q)+DC_Thickness_Array(p))./Wavelength_Array*2*pi*2;  %roundtrip
        FitResult = fit(Frequency_Array',DC_Phase_Spectrum','poly1');
        DC_Phase_Spectrum_Modi=DC_Phase_Spectrum-Frequency_Array.*FitResult.p1+FitResult.p2;      %這步減掉了2個東西: 1. constant phase (本來就沒差) 2. 前述List中之原件對干涉訊號位置之影響 (假設可用reference arm調回來)
        DC_Phase_Spectrum_Modi=DC_Phase_Spectrum_Modi-DC_Phase_Spectrum_Modi(Central_Wavelength_Index);  

        Full_DC_Phase_Spectrum=zeros(length(Full_Frequency_Array),1);
        Full_DC_Phase_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Complete_Spectrum_dn)-1))=DC_Phase_Spectrum_Modi;

        Comepensated_Spectrum_dn=Dispersed_Spectrum_dn.*exp(-i.*Full_DC_Phase_Spectrum);

        S_Compensated_Full=fft(Comepensated_Spectrum_dn);
        S_Compensated=S_Compensated_Full(1:length(Position_micron));
        PSF_FWHM_Matrix(p,q)=Position_micron(abs(find(abs(S_Compensated)>0.5*max(abs(S_Compensated)),1,'first')-find(abs(S_Compensated)>0.5*max(abs(S_Compensated)),1,'last')));
        plot(Position_micron-OPD,abs(S_Compensated));%,Position_micron,real(S_Compensated))
%         hold on
%         xlabel('Depth (micron)');
%         ylabel('Digital Number (DN)');
%         xlim([-60 60]);
%         disp(p);
    end
%     hold off
end
subplot(1,1,1)
plot(DC_Thickness_Array/1E3,PSF_FWHM_Matrix/TSRI(Center_Wavelengh,'Water'),'LineWidth',4);
%xlabel(sprintf('Thickness of %s Block (mm)',DC_Material));
xlabel('Relative Thickness of DC Block (mm)');
ylabel('Axial Resolution in Water (\mum)');
set(gca, 'FontSize', 14)

PSF_FWHM_Matrix_Folded=zeros((size(PSF_FWHM_Matrix,1)-1)/2+1,size(PSF_FWHM_Matrix,2));
DC_Tol_Array=DC_Thickness_Array((size(PSF_FWHM_Matrix,1)-1)/2+1:end);
for p=1:size(PSF_FWHM_Matrix,2)
    Right_Half=PSF_FWHM_Matrix(((size(PSF_FWHM_Matrix,1)-1)/2+1:end),p);
    Left_Half=PSF_FWHM_Matrix((((size(PSF_FWHM_Matrix,1)-1)/2+1):-1:1),p);
    PSF_FWHM_Matrix_Folded(:,p)=max([Right_Half Left_Half],[],2);
end

subplot(1,1,1)
plot(DC_Tol_Array/1E3,PSF_FWHM_Matrix_Folded/TSRI(Center_Wavelengh,'Water'),'LineWidth',4);
%xlabel(sprintf('Thickness of %s Block (mm)',DC_Material));
xlabel('DC Thickness Tolerance (mm)');
ylabel('Axial Resolution (\mum)');
xlim([0 0.15]);
ylim([4.85 5]);

set(gca, 'FontSize', 14)
grid on
DC_Thickness_Est/1000
legend(DC_Material,'Location','northwest');
% xlim([0 max(DC_Thickness_Array/1E3)]);
% ylim([min(PSF_FWHM_Matrix/TSRI(Center_Wavelengh,'Water')) max(PSF_FWHM_Matrix/TSRI(Center_Wavelengh,'Water'))]);

%% Change Material Length (i.e. Water of Eye), in SAMPLE ARM
Varying_Material='Water';
Varying_Thickness_Array=([28.03 22.615 17.2 14.28 11.36]-17.2)*1E3; %micron
PSF_FWHM_Matrix_2=zeros(length(Varying_Thickness_Array),1);

for p=1:length(Varying_Thickness_Array)
    Varying_Phase_Spectrum=TSRI(Wavelength_Array,Varying_Material).*Varying_Thickness_Array(p)./Wavelength_Array*2*pi*2;  %roundtrip
    FitResult = fit(Frequency_Array',Varying_Phase_Spectrum','poly1');
    Varying_Phase_Spectrum_Modi=Varying_Phase_Spectrum-Frequency_Array.*FitResult.p1+FitResult.p2;      %這步減掉了2個東西: 1. constant phase (本來就沒差) 2. 前述List中之原件對干涉訊號位置之影響 (假設可用reference arm調回來)
    Varying_Phase_Spectrum_Modi=Varying_Phase_Spectrum_Modi-Varying_Phase_Spectrum_Modi(Central_Wavelength_Index);  

    Full_Varying_Phase_Spectrum=zeros(length(Full_Frequency_Array),1);
    Full_Varying_Phase_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Complete_Spectrum_dn)-1))=Varying_Phase_Spectrum_Modi;

    Varied_Comepensated_Spectrum_dn=Comepensated_Spectrum_dn.*exp(i.*Full_Varying_Phase_Spectrum);

    S_Varied_Full=fft(Varied_Comepensated_Spectrum_dn);
    S_Varied=S_Varied_Full(1:length(Position_micron));
    PSF_FWHM_Matrix_2(p)=Position_micron(abs(find(abs(S_Varied)>0.5*max(abs(S_Varied)),1,'first')-find(abs(S_Varied)>0.5*max(abs(S_Varied)),1,'last')));
    plot(Position_micron-OPD,abs(S_Varied));%,Position_micron,real(S_Compensated))
    hold on
    xlabel('Depth (micron)');
    ylabel('Digital Number (DN)');
    xlim([-60 60]);
    disp(p);
end
hold off

plot(Varying_Thickness_Array/1E3+17.2,PSF_FWHM_Matrix_2/TSRI(Center_Wavelengh,'Water'),'LineWidth',4);
xlabel(sprintf('Thickness of %s Block (mm)',Varying_Material));
ylabel('Axial Resolution in Water (\mum)');
xlim([min(Varying_Thickness_Array/1E3+17.2) max(Varying_Thickness_Array/1E3+17.2)]);
ylim([min(PSF_FWHM_Matrix_2/TSRI(Center_Wavelengh,'Water')) max(PSF_FWHM_Matrix_2/TSRI(Center_Wavelengh,'Water'))]);

Diopter_Array=[-20 -10 0 10 20];

plot(Diopter_Array,PSF_FWHM_Matrix_2/TSRI(Center_Wavelengh,'Water'),'LineWidth',4);
xlabel('Diopter');
ylabel('Axial Resolution (\mum)');
xlim([min(Diopter_Array) max(Diopter_Array)]);
ylim([min(PSF_FWHM_Matrix_2/TSRI(Center_Wavelengh,'Water')) max(PSF_FWHM_Matrix_2/TSRI(Center_Wavelengh,'Water'))]);

set(gca, 'FontSize', 16)
grid on
%% Change to Linear-to-Wavelength Sensor Array
Sensor_Spectral_Pixel_Number=2048;
Sensor_Spectrum_Sampling_Ratio=0.5;
Sensor_Wavelength_Sampling_Resolution=Wavelength_FWHM/Sensor_Spectrum_Sampling_Ratio/Sensor_Spectral_Pixel_Number;
Sensor_Wavelength_Array=Sensor_Wavelength_Sampling_Resolution*((-1*Sensor_Spectral_Pixel_Number/2+1):Sensor_Spectral_Pixel_Number/2)+round(Center_Wavelengh/Sensor_Wavelength_Sampling_Resolution)*Sensor_Wavelength_Sampling_Resolution;
Sensor_Central_Wavelength_Index=find(Sensor_Wavelength_Array<Center_Wavelengh,1,'first');

plot(Sensor_Wavelength_Array,'.');

%% Calculate the "Measured" Spectrum
Optical_Wavelength_Resolution=0.02;

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
