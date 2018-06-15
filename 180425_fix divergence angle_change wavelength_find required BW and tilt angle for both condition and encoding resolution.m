clear all
%%
Central_Wavelength=[0.48:0.01:0.80];
Focal_Length=1.25;  %micron Zemax notation, physical focal length/image space RI
Grating_Pitch=500:100:1000; %1144   1000 950
Order=-1;
RI_Inc=1.85;
RI_Out=1;
Aperture_Diameter=2;   %%2.2 1.35
Target_Divergence_Angle=12; %degree

BW_Increment=0.00001;  %micron
%%
Theta_FW=zeros(length(Central_Wavelength),length(Grating_Pitch));
Theta_LI=zeros(length(Central_Wavelength),length(Grating_Pitch));
Required_BW_FW=zeros(length(Central_Wavelength),length(Grating_Pitch));
Required_BW_LI=zeros(length(Central_Wavelength),length(Grating_Pitch));
Wavelength_Resolution_FW=zeros(length(Central_Wavelength),length(Grating_Pitch));
Wavelength_Resolution_LI=zeros(length(Central_Wavelength),length(Grating_Pitch));
Number_of_Resolvable_Point_FW=zeros(length(Central_Wavelength),length(Grating_Pitch));
Number_of_Resolvable_Point_LI=zeros(length(Central_Wavelength),length(Grating_Pitch));
Required_Grating_Diameter_FW=zeros(length(Central_Wavelength),length(Grating_Pitch));
Required_Grating_Diameter_LI=zeros(length(Central_Wavelength),length(Grating_Pitch));
Total_Pitch_FW=zeros(length(Central_Wavelength),length(Grating_Pitch));
Total_Pitch_LI=zeros(length(Central_Wavelength),length(Grating_Pitch));
%%
for q=1:length(Grating_Pitch)
    for p=1:length(Central_Wavelength)
        %% Forward-view Condition (Theta_IN = Theta_OUT)
        Theta_FW(p,q)=asin(((RI_Inc-RI_Out)^-1)*Central_Wavelength(p)/1000*Grating_Pitch(q))/pi*180;
        Current_Wavelength_BW=0;
        Current_Divergence_Angle=0;
        while (Current_Divergence_Angle < Target_Divergence_Angle)
            S_Wavelength=Central_Wavelength(p)-Current_Wavelength_BW/2;
            L_Wavelength=Central_Wavelength(p)+Current_Wavelength_BW/2;
            S_Angle=asin((RI_Inc*sin(Theta_FW(p,q)/180*pi)+Order*S_Wavelength/1000*Grating_Pitch(q))/RI_Out)/pi*180;
            L_Angle=asin((RI_Inc*sin(Theta_FW(p,q)/180*pi)+Order*L_Wavelength/1000*Grating_Pitch(q))/RI_Out)/pi*180;
            Current_Divergence_Angle=abs(S_Angle-L_Angle);
            Current_Wavelength_BW=Current_Wavelength_BW+BW_Increment;
        end
        Required_BW_FW(p,q)=Current_Wavelength_BW;
        Required_Grating_Diameter_FW(p,q)=Aperture_Diameter/cos(Theta_FW(p,q)/180*pi);
        Total_Pitch_FW(p,q)=Grating_Pitch(q)*Required_Grating_Diameter_FW(p,q);
        Wavelength_Resolution_FW(p,q)=abs(Central_Wavelength(p)/Total_Pitch_FW(p,q)/Order);
        Number_of_Resolvable_Point_FW(p,q)=Required_BW_FW(p)/Wavelength_Resolution_FW(p,q);

        %% Littrow Condition (Theta_IN = -Theta_OUT)
        Theta_LI(p,q)=asin(((RI_Inc+RI_Out)^-1)*Central_Wavelength(p)/1000*Grating_Pitch(q))/pi*180;
        Current_Wavelength_BW=0;
        Current_Divergence_Angle=0;
        while (Current_Divergence_Angle < Target_Divergence_Angle)
            S_Wavelength=Central_Wavelength(p)-Current_Wavelength_BW/2;
            L_Wavelength=Central_Wavelength(p)+Current_Wavelength_BW/2;
            S_Angle=asin((RI_Inc*sin(Theta_LI(p,q)/180*pi)+Order*S_Wavelength/1000*Grating_Pitch(q))/RI_Out)/pi*180;
            L_Angle=asin((RI_Inc*sin(Theta_LI(p,q)/180*pi)+Order*L_Wavelength/1000*Grating_Pitch(q))/RI_Out)/pi*180;
            Current_Divergence_Angle=abs(S_Angle-L_Angle);
            Current_Wavelength_BW=Current_Wavelength_BW+BW_Increment;
        end
        Required_BW_LI(p,q)=Current_Wavelength_BW;
        Required_Grating_Diameter_LI(p,q)=Aperture_Diameter/cos(Theta_LI(p,q)/180*pi);
        Total_Pitch_LI(p,q)=Grating_Pitch(q)*Required_Grating_Diameter_LI(p,q);
        Wavelength_Resolution_LI(p,q)=abs(Central_Wavelength(p)/Total_Pitch_LI(p,q)/Order);
        Number_of_Resolvable_Point_LI(p,q)=Required_BW_LI(p,q)/Wavelength_Resolution_LI(p,q);
    end
    disp(q);
end
% plot(Central_Wavelength,Theta_FW,Central_Wavelength,Theta_LI);
% plot(Central_Wavelength,Required_BW_FW,Central_Wavelength,Required_BW_LI);
% plot(Central_Wavelength,Number_of_Resolvable_Point_FW,Central_Wavelength,Number_of_Resolvable_Point_LI);
% plot(Required_BW_LI,Number_of_Resolvable_Point_LI);

imagesc(abs(Number_of_Resolvable_Point_LI./Required_Grating_Diameter_LI))
TickSampling_Ratio_X=1;
TickSampling_Ratio_Y=5;

XTickArray=TickSampling_Ratio_X:TickSampling_Ratio_X:floor(length(Grating_Pitch)/TickSampling_Ratio_X)*TickSampling_Ratio_X;
YTickArray=TickSampling_Ratio_Y:TickSampling_Ratio_Y:floor(length(Central_Wavelength)/TickSampling_Ratio_Y)*TickSampling_Ratio_Y;

XTickLabelArray=round(Grating_Pitch(XTickArray));
YTickLabelArray=round(Central_Wavelength(YTickArray)*1000);
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)
set(gca,'YTick',YTickArray)
set(gca,'YTickLabel',YTickLabelArray)
set(gca, 'FontSize', 14)
xlabel('Grating Density (cycles/mm)');
ylabel('Central Wavelength (nm)');
