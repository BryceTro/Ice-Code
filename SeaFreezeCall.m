
function [k_i,rho_i,c_i,Tm,Phase_list]=SeaFreezeCall(T,Pressure)

tic
 
addpath('C:\Users\Bryce Troncone\Desktop\Matlab Folder\SeaFreeze-master\SeaFreeze-master\Matlab')

%Only need to input these three values and calculates everything else

global start_height; %in km
global resolution; %in km
global final_height; %in km

%------------------------------------------------------------------------

%start_height=start_height+1;
%final_height=final_height+1;

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

Phase= SF_WhichPhase({Pressure(i),T(i)}); 

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
    PT = {Pressure(i),T(i)};
    out=SeaFreeze(PT,'water1');
 
elseif Phase==1;
    PT = {Pressure(i),T(i)};
    out=SeaFreeze(PT,'Ih');
 
elseif Phase==2;
    PT = {Pressure(i),T(i)};
    out=SeaFreeze(PT,'II');
 
elseif Phase==3;
    PT = {Pressure(i),T(i)};
    out=SeaFreeze(PT,'III');
    
elseif Phase==5;
    PT = {Pressure(i),T(i)};
    out=SeaFreeze(PT,'V');
 
else %Phase==6;
    PT = {Pressure(i),T(i)};
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
 
 
nT1=repmat(T1,231,1);
 
mT=nT1(find(out==0));
 
P1_short=minP:stepP:maxP;
 
nP1=reshape(repmat(P1_short',201,1)', 231,201);
 
mP=nP1(find(out==0));
 
 
tol = 5;
MeltingTemps=mT(abs(mP-Pressure(i)) < tol);
Melting_Temperature=MeltingTemps(1);
 
Melting_Temperature_list(i)= Melting_Temperature;
 
end

k_i=Thermal_Conductivity_list;
rho_i=Density_list;
c_i=Specific_Heat_list;
Tm=Melting_Temperature_list;

end

