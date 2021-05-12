%For each height and melting temperature, calculates the heat flow needed to reach that temp.

%Pressure Profile

clc;
close all;

tic
 
addpath('C:\Users\Bryce Troncone\Desktop\Matlab Folder\SeaFreeze-master\SeaFreeze-master\Matlab')

%Only need to input these three values and calculates everything else

start_height=1;
resolution=1;
final_height=80;

%------------------------------------------------------------------------

rho1=917; %density of ice I
rho2=1170; %density of ice II
rho5=1230; %density of ice V
rho6=1310; %density of ice VI
 
g=23.7; %m/s^2
 
h1=13.804; %height of glacier up to Ice I and Ice II boundary
h2=21.0166; %height of glacier between Ice II and Ice V boundary
h3=41.5991; %height of glacier between Ice V and Ice VI boundary
h4=80.2502; %height of glacier between Ice VI and Ice VII boundary

%Depth (in km) Resolution Goes Below:
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
    
if (h<13.804);
    
    P1=(rho1*g*h)*1000; %pressure of glacier up to Ice I and Ice II boundary
 
elseif (h>=13.804) && (h<21.0166);
    
    P1=((rho1*g*h1) + (rho2*g*(h-h1)))*1000; %pressure of glacier between Ice II and Ice V boundary
    
elseif (h>=21.0166) && (h<41.5991);
    
    P1=((rho1*g*h1) + (rho2*g*(h2-h1))+(rho5*g*(h-h2)))*1000; %pressure of glacier between Ice V and Ice VI boundary
    
else %(h>=41.5991) && (h<80.2502);
    
    P1=((rho1*g*h1) + (rho2*g*(h2-h1))+(rho5*g*(h3-h2)+(rho1*g*(h-h3))))*1000; %pressure of glacier between Ice VI and Ice VII boundary
 
end

 
%Calculator to know ice phase
 
Pressure_list(i)=P1.*(10^-6); %give pressure in MPa
 
Pressure=Pressure_list(i);
 
T=235;   %give temperature in Kelvin-assume constant temp for now (surface temp of LHS 1140 b)

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

% imagesc(minT:maxT,minP:maxP,out);
% ylabel('Pressure (MPa)')
% xlabel('Temperature (K)')
% hcb=colorbar;
% title(hcb,'Ice Phase')
% figure;
 
 
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

i=1
for z= [Height_list]*1000 % Heights (depths) in meters
    
for Tz= [Melting_Temperature_list]   % Melting Temperature (K) at depth z 
    
x1= 13800;    %thickness of Ice I (m)
k1= 2.4;      %thermal conductivity of Ice I (W/m K)
x2= 7300;     %thickness of Ice II
k2= 1.8;    %thermal conductivity of Ice II
x5= 20600;    %thickness of Ice V
k5= 1.5;      %thermal conductivity of Ice V
x6= z-41000;  %thickness of Ice VI
k6= 1.9;      %thermal conductivity of Ice VI
    
%r= 9110014  % Planetary radius, ex. LHS 1140 b     
    
%A= 4*pi*(r.^2) % Surface Area of Planet

Ts= 235; % Planetary Surface Temperature

if (z>0) & (z<13800)
    
    q(i) = ((Tz-Ts))/(x1/k1)
    
elseif (z>=13800) & (z<21000)
    
    q(i)= ((Tz-Ts))/((x1/k1)+(x2/k2))

elseif (z>=21000) & (z<41600)
    
    q(i)= ((Tz-Ts))/((x1/k1)+(x2/k2)+(x5/k5))
    
else %(z>=41600) && (z<80250)
    
q(i)= (Tz-Ts)/((x1/k1)+(x2/k2)+(x5/k5)+(x6/k6))

i=i+1
end
end
end

figure
qq= reshape(q,[],80);
qq=qq*1000;
imagesc([Height_list],[Melting_Temperature_list], qq)
 
xlabel('Depth of Ice (km)')
ylabel('Melting Temperature (K)')
title('Heat Flow (mW/m^2)')
colorbar

