clear all
%%
Central_Wavelength=[0.65 0.78 0.84 1.03];
Wavelength_BW=[0.04:0.01:0.2];      %micron
Focal_Length=1.25;  %micron Zemax notation, physical focal length/image space RI
Order=-1;   %-1
RI_Inc=1.92;    %1.85 1.92
RI_Out=1;
Aperture_Diameter=2;   %%2.2 1.35
Target_Divergence_Angle=20; %degree
Target_Number_of_Resolvable_Point=400;
Grating_Density_Increment=5;  %cycles/mm
%%
Theta_FW=zeros(length(Central_Wavelength),length(Wavelength_BW));
Theta_LI=zeros(length(Central_Wavelength),length(Wavelength_BW));
Required_Grating_Density_FW=zeros(length(Central_Wavelength),length(Wavelength_BW));
Required_Grating_Density_LI=zeros(length(Central_Wavelength),length(Wavelength_BW));
Required_Grating_Diameter_FW=zeros(length(Central_Wavelength),length(Wavelength_BW));
Required_Grating_Diameter_LI=zeros(length(Central_Wavelength),length(Wavelength_BW));
Required_Aperture_Diameter_FW=zeros(length(Central_Wavelength),length(Wavelength_BW));
Required_Aperture_Diameter_LI=zeros(length(Central_Wavelength),length(Wavelength_BW));
%%
for q=1:length(Wavelength_BW)
    for p=1:length(Central_Wavelength)
        S_Wavelength=Central_Wavelength(p)-Wavelength_BW(q)/2;
        L_Wavelength=Central_Wavelength(p)+Wavelength_BW(q)/2;
        %% Forward-view Condition (Theta_IN = Theta_OUT)

        % Detect TIR
        Current_Grating_Density=0;
        Current_Divergence_Angle=0;
        Current_Theta_FW=0;
        while (Current_Divergence_Angle < Target_Divergence_Angle) && isreal(Current_Theta_FW)
            Current_Theta_FW=asin(((RI_Inc-RI_Out)^-1)*Central_Wavelength(p)/1000*Current_Grating_Density)/pi*180;
            
            S_Angle=asin((RI_Inc*sin(Current_Theta_FW/180*pi)+Order*S_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
            L_Angle=asin((RI_Inc*sin(Current_Theta_FW/180*pi)+Order*L_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
            Current_Divergence_Angle=abs(S_Angle-L_Angle);
            Current_Grating_Density=Current_Grating_Density+Grating_Density_Increment;
        end
        Theta_FW(p,q)=Current_Theta_FW;
        Required_Grating_Density_FW(p,q)=Current_Grating_Density;
        Required_Grating_Diameter_FW(p,q)=abs(Central_Wavelength(p)/Wavelength_BW(q)/Order/Required_Grating_Density_FW(p,q)*Target_Number_of_Resolvable_Point);
        Required_Aperture_Diameter_FW(p,q)=Required_Grating_Diameter_FW(p,q)*cos(Theta_FW(p,q)/180*pi);
        %% Littrow Condition (Theta_IN = -Theta_OUT)
        Current_Grating_Density=0;
        Current_Divergence_Angle=0;
        Current_Theta_LI=0;
        while (Current_Divergence_Angle < Target_Divergence_Angle) && isreal(Current_Theta_LI)
            Current_Theta_LI=asin(((RI_Inc+RI_Out)^-1)*Central_Wavelength(p)/1000*Current_Grating_Density)/pi*180;
            S_Angle=asin((RI_Inc*sin(Current_Theta_LI/180*pi)+Order*S_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
            L_Angle=asin((RI_Inc*sin(Current_Theta_LI/180*pi)+Order*L_Wavelength/1000*Current_Grating_Density)/RI_Out)/pi*180;
            Current_Divergence_Angle=abs(S_Angle-L_Angle);
            Current_Grating_Density=Current_Grating_Density+Grating_Density_Increment;
        end
        Theta_LI(p,q)=Current_Theta_LI;
        Required_Grating_Density_LI(p,q)=Current_Grating_Density;
        Required_Grating_Diameter_LI(p,q)=abs(Central_Wavelength(p)/Wavelength_BW(q)/Order/Required_Grating_Density_LI(p,q)*Target_Number_of_Resolvable_Point);
        Required_Aperture_Diameter_LI(p,q)=Required_Grating_Diameter_LI(p,q)*cos(Theta_LI(p,q)/180*pi);
    end
    disp(q);
end
% plot(Central_Wavelength,Theta_FW,Central_Wavelength,Theta_LI);
% plot(Central_Wavelength,Required_BW_FW,Central_Wavelength,Required_BW_LI);
% plot(Central_Wavelength,Number_of_Resolvable_Point_FW,Central_Wavelength,Number_of_Resolvable_Point_LI);
% plot(Required_BW_LI,Number_of_Resolvable_Point_LI);

imagesc(real(Theta_LI))
TickSampling_Ratio_X=1;
TickSampling_Ratio_Y=1;

XTickArray=TickSampling_Ratio_X:TickSampling_Ratio_X:floor(length(Wavelength_BW)/TickSampling_Ratio_X)*TickSampling_Ratio_X;
YTickArray=TickSampling_Ratio_Y:TickSampling_Ratio_Y:floor(length(Central_Wavelength)/TickSampling_Ratio_Y)*TickSampling_Ratio_Y;

XTickLabelArray=round(Wavelength_BW(XTickArray)*1000);
YTickLabelArray=round(Central_Wavelength(YTickArray)*1000);
set(gca,'XTick',XTickArray)

set(gca,'XTickLabel',XTickLabelArray)
set(gca,'YTick',YTickArray)
set(gca,'YTickLabel',YTickLabelArray)
set(gca, 'FontSize', 14)
xlabel('Bandwidth (nm)');
ylabel('Central Wavelength (nm)');
Required_Grating_Density_LI
Required_Aperture_Diameter_LI
Theta_LI


%%
plot(Wavelength_BW*1000,real(Theta_LI)*2,'LineWidth',4);
%xlabel(sprintf('Thickness of %s Block (mm)',DC_Material));
xlabel('Avaliable Optical Bandwidth (nm)');
ylabel('Y-Tilt Angle (degree)');
xlim([min(Wavelength_BW*1000) max(Wavelength_BW*1000)])
ylim([0 150])

set(gca, 'FontSize', 14)
grid on
title(sprintf('Objective Lens Field of View: %g degree',Target_Divergence_Angle),'Color','b')
legend(sprintf('Central Wavelength: %g nm',Central_Wavelength(1)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(2)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(3)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(4)*1000));
%%


%%

plot(Wavelength_BW*1000,Required_Grating_Density_LI,'LineWidth',4);
%xlabel(sprintf('Thickness of %s Block (mm)',DC_Material));
xlabel('Avaliable Optical Bandwidth (nm)');
ylabel('Required Grating Density (cycle/mm)');
xlim([min(Wavelength_BW*1000) max(Wavelength_BW*1000)])
ylim([0 4000])

set(gca, 'FontSize', 14)
grid on
title(sprintf('Objective Lens Field of View: %g degree',Target_Divergence_Angle),'Color','b')
legend(sprintf('Central Wavelength: %g nm',Central_Wavelength(1)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(2)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(3)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(4)*1000));

%%
plot(Wavelength_BW*1000,Required_Aperture_Diameter_LI,'LineWidth',4);
%xlabel(sprintf('Thickness of %s Block (mm)',DC_Material));
xlabel('Avaliable Optical Bandwidth (nm)');
ylabel('Required System Aperture Diameter (mm)');
xlim([min(Wavelength_BW*1000) max(Wavelength_BW*1000)])
ylim([0 3])
set(gca, 'FontSize', 14)
grid on
title(sprintf('Objective Lens Field of View: %g degree \n  Number of Resolving Points: %g',Target_Divergence_Angle,Target_Number_of_Resolvable_Point),'Color','b')
legend(sprintf('Central Wavelength: %g nm',Central_Wavelength(1)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(2)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(3)*1000),sprintf('Central Wavelength: %g nm',Central_Wavelength(4)*1000));

%% Old Relation
%         Required_Grating_Diameter_LI(p,q)=Aperture_Diameter/cos(Theta_LI(p,q)/180*pi);
%         Total_Pitch_LI(p,q)=Grating_Pitch(q)*Required_Grating_Diameter_LI(p,q);
%         Wavelength_Resolution_LI(p,q)=abs(Central_Wavelength(p)/Total_Pitch_LI(p,q)/Order);
%         Number_of_Resolvable_Point_LI(p,q)=Required_BW_LI(p,q)/Wavelength_Resolution_LI(p,q);
