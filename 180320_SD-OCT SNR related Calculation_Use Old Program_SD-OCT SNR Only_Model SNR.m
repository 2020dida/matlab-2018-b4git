clear all

If_Model_SNR=1;
If_AxialBinning=0;
% Power Related Parameter
SLD_Power=7.5;       %mW
Target_Power_on_Sample=0.75;    %mw
Line_Rate=80000;    %Hz
Sample_SinglePass_CE=0.7;
Sample_R=0.005;
Reference_R=0.001:0.0002:0.05;
Reference_SinglePass_CE=0.7;
Spectrometer_CE=0.6*0.6;
Sturation_Level_Max=1;
IE_Max=1;
Target_Axial_Sampling_Resolution=2.5;   %micron
RI=1.33;
%% Calculating BS ratio
BS_Ratio=(1-Target_Power_on_Sample/Sample_SinglePass_CE/SLD_Power);
%%
Plank=6.626069934*1E-34; %J*s
% System Parameter
FWC_e=140000;
Bit=12;
DN_max=2^Bit;
Dark_Noise_dn=0.5*2^4; %DN @12bit
Gain=DN_max/FWC_e;
OPD=50;  %micron, free space
Fringe_Period_frequency=3E8/((OPD*1E-6)*2);

% Simulation Setting
Repeat_Times=10000;  %for Noise Calculation
% Spectrum Parameter
X_spectrum_number=2048;
Center_Wavelengh=0.84;  %micron
Wavelength_FWHM=0.05;   %micron
Spectrum_Sampling_Ratio=0.7;    %?Spectrum?FWHM????????
ZP_Ratio=8;
% SD-OCT Part
% Simulating Spectrum (Linear with k, simplified) Normalzed 
Frequency_FWHM=3E8*Wavelength_FWHM/(Center_Wavelengh^2)/(1E-6);
Calculated_Axial_Resolution=3E8/Frequency_FWHM/2*1E6;   %micron
Center_Frequency=3E8/(Center_Wavelengh*1E-6);
Frequency_Sampling_Resolution=Frequency_FWHM/Spectrum_Sampling_Ratio/X_spectrum_number;

Frequency_Array=Frequency_Sampling_Resolution*((-1*X_spectrum_number/2+1):X_spectrum_number/2)+round(Center_Frequency/Frequency_Sampling_Resolution)*Frequency_Sampling_Resolution;
% DC
Spectrum_Norm_e=gaussmf(Frequency_Array,[Frequency_FWHM/2.355 Center_Frequency]);
%% Interference Efficiency and Satunration Level Estimation
Exp_Time=1/Line_Rate;
for q=1:length(Reference_R)
    Reference_Power(q)=SLD_Power*BS_Ratio*(1-BS_Ratio)*(Reference_SinglePass_CE^2)*Reference_R(q);
    Ratio_R=BS_Ratio*(1-BS_Ratio)*(Reference_SinglePass_CE^2)*Reference_R(q);

    Sample_Power(q)=SLD_Power*(1-BS_Ratio)*BS_Ratio*(Sample_SinglePass_CE^2)*Sample_R;
    Ratio_S=(1-BS_Ratio)*BS_Ratio*(Sample_SinglePass_CE^2)*Sample_R;

    Sample_Power_at_Sample=SLD_Power.*(1-BS_Ratio).*Sample_SinglePass_CE;

    %%
    Number_of_Electron_per_Line=Spectrometer_CE*(Reference_Power(q)+Sample_Power(q))*1E-3.*Exp_Time/(Plank*3E8/(Center_Wavelengh.*1E-6));
    Sturation_Level_Power=Number_of_Electron_per_Line./sum(Spectrum_Norm_e)/FWC_e;
    Sturation_Level(q)=min(Sturation_Level_Max,Sturation_Level_Power);
    IE(q)=min(IE_Max,2*((Sample_Power(q)*Reference_Power(q))^0.5)/(Sample_Power(q)+Reference_Power(q)));

    DN_mean=DN_max*Sturation_Level(q);
    Spectrum_Norm_e=gaussmf(Frequency_Array,[Frequency_FWHM/2.355 Center_Frequency]);
    
    % Use V100 to Estimate
    %((11*1E-3)*(0.16E-3)*0.25*0.01*0.17)/(Plank*3E8/(Center_Wavelengh.*1E-6))
    %%
    Spectrum_e=gaussmf(Frequency_Array,[Frequency_FWHM/2.355 Center_Frequency]).*FWC_e.*Sturation_Level(q);
    Spectrum_dn=Spectrum_e*Gain; 
    % Inter
    Inter_Spectrum_e=Spectrum_e*IE(q).*cos((Frequency_Array-Center_Frequency)/Fringe_Period_frequency*2*pi);
    Spectrum_with_Inter_e=Spectrum_e+Inter_Spectrum_e;
    % Full Frequency Range
    Full_Frequency_Array=[0:Frequency_Sampling_Resolution:ZP_Ratio*max(Frequency_Array)]';
    Min_Frequency_Index=find(Full_Frequency_Array>min(Frequency_Array),1,'first');
    %% Generating Position Array
    Time_max=1/Frequency_Sampling_Resolution;    %sec
    Temporal_Resolution=Time_max/length(Full_Frequency_Array);
    Time=[Temporal_Resolution:Temporal_Resolution:Time_max];
    Position_Full=3E8*Time/2/RI;  %/2 for round-trip
    Position=Position_Full(1:round(size(Full_Frequency_Array)/2));
    Position_micron=Position*1E6;
    Position_Sampling_Resolution_micron=Position_micron(2)-Position_micron(1);
    Index_Theoritical_Peak_Signal=round(OPD/Position_Sampling_Resolution_micron/RI);

    %% Simulating Spectrum (Linear with k, simplified) with Noise
    Full_Spectrum=zeros(size(Full_Frequency_Array,1),size(Full_Frequency_Array,2));
    if If_Model_SNR == 1
            Spectrum_with_Inter_dn=round(Spectrum_with_Inter_e.*Gain);
            Full_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Spectrum_with_Inter_dn)-1))=Spectrum_with_Inter_dn-Spectrum_dn;
            % plot(Full_Frequency_Array,Full_Spectrum);
            % Signal
            S_Full=fft(Full_Spectrum);
            S=S_Full(1:round(size(S_Full)/2));
            if If_AxialBinning == 1
                Axial_Binning_Factor=floor(Target_Axial_Sampling_Resolution/Position_Sampling_Resolution_micron);
                S_Binned=TSBinning(abs(S),1,Axial_Binning_Factor);
                Index_Theoritical_Peak_Signal_Binned=round(Index_Theoritical_Peak_Signal/Axial_Binning_Factor);
                Signal_SDOCT=S_Binned(Index_Theoritical_Peak_Signal_Binned);
            else
                Signal_SDOCT=abs(S(Index_Theoritical_Peak_Signal));
            end
            Spectrum_noise_dn=round((Spectrum_with_Inter_e.^0.5)*Gain);
            Noise_SDOCT=mean(Spectrum_noise_dn)*(X_spectrum_number/1.8)^0.5;    %real:2, abs: 4.68, sim: 1.8
    else

        Signal_SDOCT_Array=zeros(Repeat_Times,1);
        for p=1:Repeat_Times
            % Noise
            Shot_Noise_Spectrum_e=randn(X_spectrum_number,1)'.*Spectrum_with_Inter_e.^0.5;

            Spectrum_e_with_Noise=Spectrum_with_Inter_e+Shot_Noise_Spectrum_e;
            Spectrum_dn_with_Noise=round(Spectrum_e_with_Noise.*Gain);
            % plot(Frequency_Array,Spectrum_dn_with_Noise);
            %% SD-OCT Calculation (Tranditional)
            Full_Spectrum(Min_Frequency_Index:(Min_Frequency_Index+length(Spectrum_dn_with_Noise)-1))=Spectrum_dn_with_Noise-Spectrum_dn;
            % plot(Full_Frequency_Array,Full_Spectrum);

            S_Full=fft(Full_Spectrum);
            S=S_Full(1:round(size(S_Full)/2));
            %plot(Position_micron,S);
            %% Axial Binning
            if If_AxialBinning == 1
                Axial_Binning_Factor=floor(Target_Axial_Sampling_Resolution/Position_Sampling_Resolution_micron);
                S_Binned=TSBinning(abs(S),1,Axial_Binning_Factor);
                Index_Theoritical_Peak_Signal_Binned=round(Index_Theoritical_Peak_Signal/Axial_Binning_Factor);
            end

            if If_AxialBinning == 1
                Signal_SDOCT_Array(p)=S_Binned(Index_Theoritical_Peak_Signal_Binned);
            else
                Signal_SDOCT_Array(p)=abs(S(Index_Theoritical_Peak_Signal));
            end
            disp(p);
        end

        Signal_SDOCT=mean(Signal_SDOCT_Array);
        Noise_SDOCT=std(Signal_SDOCT_Array);
    end

    %% SNR Calculation for SD-OCT
    %plot(Signal_SDOCT_Array);

    SNR_SDOCT(q)=20*log10(Signal_SDOCT/Noise_SDOCT);
end
%% Simple Est
SNR_SD=10*log10(140000*2048/8)
10*log(4*(0.7E-3)*1/50000*0.6/(6.6E-34)/(3E8/(0.84E-6)))
%%
Reference_Sample_Power_Ratio=Reference_Power./Sample_Power;
Data_Save_Folder=['C:\TuanShu\180316_Setup MiiS Probe for OPD Calc\'];

subplot(2,2,1)
plot(Reference_R*100,SNR_SDOCT,'r','LineWidth',3);
xlabel('Reference arm Reflectance (%)')
ylabel('Signal-to-Noise Ratio (dB, shot)')
%dlmwrite([Data_Save_Folder 'SNRvsRr.txt'],[Reference_R'*100 SNR_SDOCT'],'delimiter','\t','newline','pc','precision', '%.6f','-append');

%%

    subplot(2,2,2)
    plot(Reference_R*100,Reference_Sample_Power_Ratio,'b','LineWidth',3);
    xlabel('Reference Arm Reflectance (%)')
    ylabel('Reference/Sample Power Ratio')    
    dlmwrite([Data_Save_Folder 'RSratio.txt'],[Reference_R'*100 Reference_Sample_Power_Ratio'],'delimiter','\t','newline','pc','precision', '%.6f','-append');

    subplot(2,2,3)
    plot(Reference_R*100,Sturation_Level*100,'b','LineWidth',3);
    xlabel('Reference Arm Reflectance (%)')
    ylabel('Saturation Level (%)')    
    dlmwrite([Data_Save_Folder 'SaturationLevel.txt'],[Reference_R'*100 Sturation_Level'*100],'delimiter','\t','newline','pc','precision', '%.6f','-append');

    subplot(2,2,4)
    plot(Reference_R*100,IE*100,'b','LineWidth',3);
    xlabel('Reference Arm Reflectance (%)')
    ylabel('Interference Efficiency (%)')    
    dlmwrite([Data_Save_Folder 'IE.txt'],[Reference_R'*100 IE'*100],'delimiter','\t','newline','pc','precision', '%.6f','-append');
%% 
subplot(1,1,1)
[hAx,hLine1,hLine2] = plotyy(Reference_R*100,SNR_SDOCT,Reference_R*100,IE*100);
hLine1.LineWidth = 4;
hLine2.LineWidth = 2.5;
hLine1.LineStyle = '-';
hLine2.LineStyle = ':';

hAx(1).FontSize=12;
hAx(2).FontSize=12;
xlabel(hAx(1),'Reference Arm Reflectance (%)');
ylabel(hAx(1),'Signal-to-Noise Ratio (dB)');
ylabel(hAx(2),'Reference/Sample Power Ratio');

%% 

subplot(1,1,1)
plot(Reference_R*100,SNR_SDOCT,'LineWidth',4,'LineStyle','-');
Ax=gca;
Ax.FontSize=12;
xlabel(Ax,'Reference Arm Reflectance (%)');
ylabel(Ax,'Signal-to-Noise Ratio (dB)');
