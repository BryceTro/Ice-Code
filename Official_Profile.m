%Pressure Profile

clc;
close all;
clear all;

addpath('C:\Users\Bryce Troncone\Desktop\Matlab Folder\SeaFreeze-master\SeaFreeze-master\Matlab')

global start_height
global resolution
global final_height

%Input the below values to run

T_surf=200; %Surface Temp of planet (in K)
Base_Flux=100/1000;  %Heat Flux (in W/m^2)
Timescale=10000000; %ex. 1,000,000 year timescale; 10,000,000 year timescale;100,000,000 timescale etc. 
start_height=0; %from 0 km onward
resolution=1;   %height of each step in km (can be decimals)
final_height=80; %up to 80 km
g=23.7; %gravity of planet (in m/s^2)
Patmosphere=0.101325; %atmospheric pressure of planet (in MPa)

%------------------------------------------------------------------------
%------------------------------------------------------------------------


Patmosphere=0.101325*10^6;

dt=(86400*365.25*Timescale)/1000;

rho1=917; %density of ice I
rho2=1170; %density of ice II
rho5=1230; %density of ice V
rho6=1310; %density of ice VI
 
 
h1=13.804; %height of glacier up to Ice I and Ice II boundary
h2=21.0166; %height of glacier between Ice II and Ice V boundary
h3=41.5991; %height of glacier between Ice V and Ice VI boundary
h4=80.2502; %height of glacier between Ice VI and Ice VII boundary

Height_list=[];
Pressure_list=[];
Phase_list=[];
Density_list=[];
Specific_Heat_list=[];
Thermal_Conductivity_list=[];
Melting_Temperature_list=[];
 
Height_list=start_height:resolution:final_height;

for i=1:length(Height_list)
    
h=Height_list(i);
 
%F(h,g,dz) -> [P(z),n(z),c(z),k(z),rho(z),Tm(z)]
 
%For each depth (h) in km with a given resolution (dh) gives:
 
%Pressure (P) 
%Ice type (n) 
%Specific heat (c)
%Thermal conductivity (k) 
%Density (rho)
%Melting temperature (Tm) 
 
%P=rho*g*h
    
if h==0 %is 0 km
    P1=Patmosphere;
    
elseif (h<13.804);
    
    P1=Patmosphere +(rho1*g*h)*1000; %pressure of glacier up to Ice I and Ice II boundary
 
elseif (h>=13.804) && (h<21.0166);
    
    P1=Patmosphere +((rho1*g*h1) + (rho2*g*(h-h1)))*1000; %pressure of glacier between Ice II and Ice V boundary
    
elseif (h>=21.0166) && (h<41.5991);
    
    P1=Patmosphere +((rho1*g*h1) + (rho2*g*(h2-h1))+(rho5*g*(h-h2)))*1000; %pressure of glacier between Ice V and Ice VI boundary
    
else %(h>=41.5991) && (h<80.2502);
    
    P1=Patmosphere +((rho1*g*h1) + (rho2*g*(h2-h1))+(rho5*g*(h3-h2)+(rho1*g*(h-h3))))*1000; %pressure of glacier between Ice VI and Ice VII boundary
 
end

 
%Calculator to know ice phase
 
Pressure_list(i)=P1.*(10^-6); %give pressure in MPa
 
Pressure=Pressure_list(i);
 
T=T_surf;   %give temperature in Kelvin-assume constant temp for now (surface temp of LHS 1140 b)

Phase= SF_WhichPhase({Pressure,T}); 

Phase_list(i)=Phase;
 
 
%output gives:
% 0=liquid
% 1= 'Ih' for ice Ih 
% 2= 'II' for ice II 
% 3= 'III' for ice III 
% 5= 'V' for ice V 
% 6 = 'VI' for ice VI
 
%SeaFreeze Calculator to get the specific heat and density

if Phase==0;
    PT = {Pressure,T};
    out=SeaFreeze(PT,'water1');
 
elseif Phase==1;
    PT = {Pressure,T};
    out=SeaFreeze(PT,'Ih');
 
elseif Phase==2;
    PT = {Pressure,T};
    out=SeaFreeze(PT,'II');
 
elseif Phase==3;
    PT = {Pressure,T};
    out=SeaFreeze(PT,'III');
    
elseif Phase==5;
    PT = {Pressure,T};
    out=SeaFreeze(PT,'V');
 
else %Phase==6;
    PT = {Pressure,T};
    out=SeaFreeze(PT,'VI');
 
end
   

Densities=out.rho;
Specific_Heats=out.Cp;
 
 
Density_list(i)=Densities;
 
Specific_Heat_list(i)=Specific_Heats;
 
% Cp = gives specific heat J/kg K
% rho = gives density in kg/m^3
 
 
%Thermal Conductivity Profile (Anderson 1 Paper Figure 2)

if (Phase == 0);
    Thermal_Conductivity=.60; %(W/(m*K))
 
elseif (Phase == 1);
    Thermal_Conductivity= 2.4;
    
elseif (Phase == 2);
    Thermal_Conductivity= 1.8;
 
elseif (Phase == 3);
    Thermal_Conductivity= 1.1;
 
elseif (Phase == 5);
    Thermal_Conductivity=1.5;
 
elseif (Phase == 6);
    Thermal_Conductivity= 1.9;
 
else
    Thermal="Not Available";
    Thermal_Conductivity= str2double(Thermal);
end

 
Thermal_Conductivity_list(i)=Thermal_Conductivity;
 
 
%Water Melting Points Diagram

minP=0; %MPa
stepP=10;
maxP=2300;
minT=200; %K
maxT=400;
T1 = [minT:maxT]; %temperature range
P1 = [minP:maxP]; %pressure range
out = SF_WhichPhase({minP:stepP:maxP,minT:maxT});   

%imagesc(minT:maxT,minP:maxP,out);
%ylabel('Pressure (MPa)')
%xlabel('Temperature (K)')
%hcb=colorbar;
%title(hcb,'Ice Phase')
%set(gca,'FontSize',24)
%set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
%figure;
 
 
nT1=repmat(T1,231,1);
 
mT=nT1(find(out==0));
 
P1_short=minP:stepP:maxP;
 
nP1=reshape(repmat(P1_short',201,1)', 231,201);
 
mP=nP1(find(out==0));
 
% scatter(mT,mP);
% ylabel('Pressure (MPa)')
% xlabel('Temperature (K)')
% title('Melting Points for Water')
% hold on

%Get melting points of that pressure
 
%Psimple=round(P,-1)
%plot(mT,Psimple,'r*')

%plot(mT,Pressure,'r*');
 
tol = 5;
MeltingTemps=mT(abs(mP-Pressure) < tol);
Melting_Temperature=MeltingTemps(1);
 
Melting_Temperature_list(i)= Melting_Temperature;
 
end


Pressure_Profile = cat(1,Height_list,Pressure_list, Phase_list, Density_list,Specific_Heat_list,...
Thermal_Conductivity_list,Melting_Temperature_list);
 
 
Pressure_Profile_Final = array2table(Pressure_Profile,'RowNames',{'Height(km)','Pressure(MPa)',...
    'Phase','Density(kg/m^3)','Specific Heat((J/kg/K)','Thermal Conductivity(W/m/K)','Melting Temperature(K)'})

%Making a figure of all properties
figure('units','normalized','position',[.1 .1 .6 .6])
subplot(2,3,1)
plot(Pressure_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',12)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel("Pressure (MPa)")
ylabel("Depth Below Surface (km)")
subplot(2,3,2)
scatter(Phase_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',12)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel("Ice Phase")
ylabel("Depth Below Surface (km)")
xlim([0 7])
subplot(2,3,3)
plot(Density_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',12)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel("Density (kg/m^3)")
ylabel("Depth Below Surface (km)")
subplot(2,3,4)
plot(Specific_Heat_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',12)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel("Specific Heat (J/(kg K))")
ylabel("Depth Below Surface (km)")
subplot(2,3,5)
plot(Thermal_Conductivity_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',12)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel("Thermal Conductivity (W/(m K))")
ylabel("Depth Below Surface (km)")
xlim([1.4 2.6])
subplot(2,3,6)
plot(Melting_Temperature_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',12)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel("Melting Temperature (K)")
ylabel("Depth Below Surface (km)")
sgtitle('Steady State Temperature Plots')

%Heat Diffusion Equation

% vectors
T_Start=T_surf*ones(1,length(Height_list));      %Starting Temp.
Phi_Start=0*ones(1,length(Height_list));      %Porosity Start
T_Start_hold=T_Start;
Phi_Start_hold=Phi_Start;
k_i=Thermal_Conductivity_list;          %Thermal Conductivites of Ices
rho_i=Density_list;                 %Densities of Ices
c_i=Specific_Heat_list;             %Specific Heats of Ices
k_w=0.6*ones(1,length(Height_list));          %Thermal Conductivity of Water
rho_w=1000*ones(1,length(Height_list));       %Density of Water
c_w=4180*ones(1,length(Height_list));     %Specific Heat Water
Tm=Melting_Temperature_list;        %Melting Temp

dz=1000*resolution; %resolution in (m)
L=334778;   %latent heat of fusion
TTol=0.01;      %Temp. Tolerance
PhiTol=0.001;       %Porosity Tolerance

% plotting mats
hold_temps=[];
hold_phi=[];
hold_phase=[];
i=0;

% loop over time
    %Ice phase and Temperature plots showing progressive change with depth and time
   skip = 1;
   ss = 0;
for time=[dt:dt:1000*dt]  
    i=i+1;
   
  if rem(i,100)==0
      ss = ss+1;
           [k_i,rho_i,c_i,Tm,Evolve_Phase_list]=SeaFreezeCall(T_Start,Pressure_list);
                     hold_phase=[hold_phase Evolve_Phase_list];
                     
    hold_temps=[hold_temps T_new];
    hold_phi=[hold_phi Phi_new];  
    test2(:,ss)  = Evolve_Phase_list';
   end
   
       [T_new,Phi_new]=HP_Ice_Evolve_v3(T_Start,Phi_Start,k_i,rho_i,c_i,k_w,rho_w,...
    c_w,dt,dz,T_surf,Base_Flux,Tm,L,TTol,PhiTol);

    % redefine T_start & Phi_Start for next time loop  
    T_Start=T_new';
    T_Start(1)=T_surf;
      test(:,skip)  = T_Start';

      Phi_Start=Phi_new';
 skip = skip+1;
end

figure
contourf([dt:dt:1000*dt],[Height_list],([test]))
cb=colorbar
figure
imagesc([(dt:dt:1000*dt)],[Height_list],([test2]))
cb=colorbar

figure('units','normalized','position',[.1 .1 .6 .6])
scatter(Phase_list,-Height_list,'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
title("Phase vs. Depth Below Surface Final at Time=T")
xlabel("Ice Phase")
ylabel("Depth Below Surface (km)")
xlim([0 7])

% plots over time
subplot(2,2,1)
plot(T_Start,-[Height_list],'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
hold on
i=1;
while i<size(hold_temps,2)
    plot(hold_temps(:,i),-[Height_list],'LineWidth',1.9)
    i=i+1;
end
xlabel('Temp (K)')
ylabel('Height (km)')

subplot(2,2,2)
plot(Phi_Start,-[Height_list],'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
hold on
i=1;
while i<size(hold_phi,2)
    plot(hold_phi(:,i),-[Height_list],'LineWidth',1.9)
    i=i+1;
end
xlabel('Porosity')
ylabel('Height (km)')

% final profile
subplot(2,2,3)
plot(T_new,-[Height_list],'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel('Temp (K)')
ylabel('Height (km)')

subplot(2,2,4)
plot(Phi_new,-[Height_list],'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel('Porosity')
ylabel('Height (km)')

figure('units','normalized','position',[.1 .1 .6 .6])


yyaxis left

plot((-[Height_list]), [Melting_Temperature_list],'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlabel('Height( km)')
ylabel('Melting Temperature (K)')
ylim([200 500])
hold on 

yyaxis right 

plot((-[Height_list]), [T_new],'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]); 
xlabel('Height (km)')
ylabel('Temperature In Glacier (K)')
ylim([200 500])
