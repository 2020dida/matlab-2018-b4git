clear all
%%
Number_of_Set=3;
Central_Wavelength=[0.78 0.82 1.03];
Wavelength_BW=[0.14 0.062 0.12];      %micron
Focal_Length=1.25;  %micron Zemax notation, physical focal length/image space RI
Order=-1;   %-1
RI_Inc=1.92;    %1.85 1.92
RI_Out=1;
Aperture_Diameter=2;   %%2.2 1.35
Target_Divergence_Angle=5:0.1:25; %degree
Target_Number_of_Resolvable_Point=400;
Grating_Density_Increment=5;  %cycles/mm
%%
Theta_FW=zeros(length(Target_Divergence_Angle),Number_of_Set);
Theta_LI=zeros(length(Target_Divergence_Angle),Number_of_Set);
Required_Grating_Density_FW=zeros(length(Target_Divergence_Angle),Number_of_Set);
Required_Grating_Density_LI=zeros(length(Target_Divergence_Angle),Number_of_Set);
Required_Grating_Diameter_FW=zeros(length(Target_Divergence_Angle),Number_of_Set);
Required_Grating_Diameter_LI=zeros(length(Target_Divergence_Angle),Number_of_Set);
Required_Aperture_Diameter_FW=zeros(length(Target_Divergence_Angle),Number_of_Set);
Required_Aperture_Diameter_LI=zeros(length(Target_Divergence_Angle),Number_of_Set);
%%
for Set=1:Number_of_Set
    p=Set;
    q=Set;
    for r=1:length(Target_Divergence_Angle)

    S_Wavelength=Central_Wavelength(p)-Wavelength_BW(q)/2;
    L_Wavelength=Central_Wavelength(p)+Wavelength_BW(q)/2;
    %% Forward-view Condition (Theta_IN = Theta_OUT)

    % Detect TIR
    Current_Grating_Density=0;
    Current_Divergence_Angle=0;
    Current_Theta_FW=0;
    while (Current_Divergence_Angle < Target_Divergence_Angle(r)) && isreal(Current_Theta_FW)
        Current_Theta_FW=asin(((RI_Inc-RI_Out)^-1)*Central_Wavelength(p)/1000*Current_Grating_Density)/pi*180;

        S_Angle=asin((RI_Inc*sin(Current_Theta_FW/180*pi)+Order*S_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
        L_Angle=asin((RI_Inc*sin(Current_Theta_FW/180*pi)+Order*L_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
        Current_Divergence_Angle=abs(S_Angle-L_Angle);
        Current_Grating_Density=Current_Grating_Density+Grating_Density_Increment;
    end
    Theta_FW(r,Set)=Current_Theta_FW;
    Required_Grating_Density_FW(r,Set)=Current_Grating_Density;
    Required_Grating_Diameter_FW(r,Set)=abs(Central_Wavelength(p)/Wavelength_BW(q)/Order/Required_Grating_Density_FW(r,Set)*Target_Number_of_Resolvable_Point);
    Required_Aperture_Diameter_FW(r,Set)=Required_Grating_Diameter_FW(r,Set)*cos(Theta_FW(r,Set)/180*pi);
    %% Littrow Condition (Theta_IN = -Theta_OUT)
    Current_Grating_Density=0;
    Current_Divergence_Angle=0;
    Current_Theta_LI=0;
    while (Current_Divergence_Angle < Target_Divergence_Angle(r)) && isreal(Current_Theta_LI)
        Current_Theta_LI=asin(((RI_Inc+RI_Out)^-1)*Central_Wavelength(p)/1000*Current_Grating_Density)/pi*180;
        S_Angle=asin((RI_Inc*sin(Current_Theta_LI/180*pi)+Order*S_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
        L_Angle=asin((RI_Inc*sin(Current_Theta_LI/180*pi)+Order*L_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
        Current_Divergence_Angle=abs(S_Angle-L_Angle);
        Current_Grating_Density=Current_Grating_Density+Grating_Density_Increment;
    end
    Theta_LI(r,Set)=Current_Theta_LI;
    Required_Grating_Density_LI(r,Set)=Current_Grating_Density;
    Required_Grating_Diameter_LI(r,Set)=abs(Central_Wavelength(p)/Wavelength_BW(q)/Order/Required_Grating_Density_LI(r,Set)*Target_Number_of_Resolvable_Point);
    Required_Aperture_Diameter_LI(r,Set)=Required_Grating_Diameter_LI(r,Set)*cos(Theta_LI(r,Set)/180*pi);
    end

end

plot(tan(Target_Divergence_Angle/2/180*pi)*2*Focal_Length*1000,real(Theta_LI)*2,'LineWidth',4);
%xlabel(sprintf('Thickness of %s Block (mm)',DC_Material));
xlabel('Field of View Diameter (\mum)');
ylabel('Required Y-Tilt Angle (degree)');
xlim([min(tan(Target_Divergence_Angle/2/180*pi)*2*Focal_Length*1000) max(tan(Target_Divergence_Angle/2/180*pi)*2*Focal_Length*1000)])
ylim([0 150])
set(gca,'YTick',[0 30 60 90 120 150])

set(gca, 'FontSize', 14)
grid on
title('Light Source Comparison','Color','b')
legend('Ti:sapphire','QSDM-820-10D','QSDMI-1030-15');

%% Old Relation
%         Required_Grating_Diameter_LI(p,q)=Aperture_Diameter/cos(Theta_LI(p,q)/180*pi);
%         Total_Pitch_LI(p,q)=Grating_Pitch(q)*Required_Grating_Diameter_LI(p,q);
%         Wavelength_Resolution_LI(p,q)=abs(Central_Wavelength(p)/Total_Pitch_LI(p,q)/Order);
%         Number_of_Resolvable_Point_LI(p,q)=Required_BW_LI(p,q)/Wavelength_Resolution_LI(p,q);
