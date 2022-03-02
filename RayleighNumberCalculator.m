%Newtonian Convection
clc;
close all;
clear all;
addpath ('C:\Users\Bryce Troncone\Desktop\Matlab Folder\SeaFreeze-master\SeaFreeze-master\Matlab')
addpath 'C:\Users\Bryce Troncone\Desktop\Ice-Code-main\Ice-Code-main'

global start_height
global resolution
global final_height

%Input the below values to run
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
minP=0; %MPa
maxP=2300; %MPa
minT=200; %K
maxT=400; %K

T_surf=235; %Surface Temp of planet (in K)
Base_Flux=30/1000;  %Heat Flux (in W/m^2)
Timescale=10*1E6; %ex. 1,000,000 year timescale; 10,000,000 year timescale;100,000,000 timescale etc.
start_height=0; %from 0 km onward
resolution=0.1;   %height of each step in km (can be decimals)
final_height=80 ; %up to 80 km
Patmosphere=0.101325; %atmospheric pressure of planet (in MPa)
Mx = 6.38; %Mass factor of the planet 
Rx = 1.64; %Radius factor of the planet
g=24; % gravity on planet surface (m/s^2)

%Get an initial adiabatic profile for the ice sheet and depth dependent
%melting temperature using the two subroutines below

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

[T_Start,dtdz,Pressure,rho,alpha,Cp,K,phase] = adiabat_profile(Patmosphere,T_surf,(final_height-start_height+1)*1000,resolution,Mx,Rx);
step=(maxP-minP+1)/numel(K);
[MeltT] = findmeltT(minP,maxP,minT,maxT,step,Pressure);
MeltT = MeltT';
Height_list=linspace(0,80,numel(T_Start));
Height_list=Height_list';
Pressure_list=Pressure;
alpha;        %thermal expansion coefficient  
Cp;         %specific heat
k=K';       %thermal conductivity
rho;        %density
k_w=0.6*ones(1,length(Height_list))';          %Thermal Conductivity of Water
rho_w=1000*ones(1,length(Height_list))';       %Density of Water
c_w=4180*ones(1,length(Height_list))';     %Specific Heat Water
Tm=MeltT;        %Melting Temp

%Rayleigh Number Calculator
%Ra=Rayleigh Number

R= 8.3144598; %universal gas constant J/mol/K
r= 6563*Rx;  %radius (km) of Planet
Area= 4*pi*r^2;  %Surface Area of Planet km^2
T_Surface=235; %K
h=Height_list; %km

skip=1;
for i=1:length(h)
       
td(skip)= k(i)/(rho(i)*Cp(i)); %thermal diffusivity
 
tau(skip)= Pressure(i)/Area; %basal stress

shearv(skip)=(tau(i)/rho(i))^(1/2); %shear velocity

gamma(skip)= (tau(i)*h(i))/shearv(i); %dynamic viscosity %check this

kv(skip)=gamma(i)/rho(i); %kinematic viscosity


if T_Start(i)>=T_Surface
    
Ra(skip)= (rho(i)*g*alpha(i)*(T_Start(i)-T_Surface)*h(i)^3)/(td(i)*kv(i)*1e9); %Rayleigh

elseif T_Start(i)<T_Surface
    
Ra(skip)= (rho(i)*g*alpha(i)*(T_Surface-T_Start(i))*h(i)^3)/(td(i)*kv(i)*1e9); %Rayleigh
end

Pr(skip)=(kv(i))/td(i);    %Prandtl Number

% for free convection, the average Nusselt number is expressed as a 
%function of the Rayleigh number and the Prandtl number-Churchill and Chu
%Churchill–Bernstein equation
%Nu >> 1 , More effective Convection
%Nu = 1, For fluid layer represent heat transfer by pure conduction across layer

%Nusselt Calculator

Var1(skip) = 1 + (Ra(i)./282000).^(5/8);
Var2(skip) = Var1(i).^(4/5);

Var3(skip) = 0.62*Ra(i).^(1/2)*Pr(i).^(1/3);
Var4(skip) = (1 + (0.4/Pr(i)).^(2/3))^(1/4);

Nu(skip) = 0.3 + Var3(i)/Var4(i)*Var2(i);

skip=skip+1;
end

figure
plot(Ra,h)
xlabel('Rayleigh Number')
ylabel('Depth(km)')
xline(10^6,'--','Critical Rayleigh Threshold')

figure
plot(Pr,h)
xlabel('Prandtl Number')
ylabel('Depth(km)')

figure
plot(Nu,h)
xlabel('Nusselt Number')
ylabel('Depth(km)')

Resolution_Profile= [h Ra' Pr' Nu' T_Start];

%add ice evolve temperature and phases to it now?