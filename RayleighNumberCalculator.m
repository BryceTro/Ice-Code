%Newtonian Convection

function [Ra,Nu,td] = RayleighNumberCalculator(T_Start,T_surf,Base_Flux,g,Rx,k_i,rho_i,c_i,alpha,Pressure,Height_list);
g=g;
T_surf=T_surf; %K
Height_list=Height_list; %km
k_i=k_i;
rho_i=rho_i;
c_i=c_i;
Rx=Rx;
Pressure=Pressure;
Base_Flux=Base_Flux;
v=10e12;   %constant viscosity (Pa)

%^from Water-rich planets: How habitable is a water layer deeper than on Earth'
skip=1;
for i=1:length(Height_list)
       
td(skip)= k_i(i)/(rho_i(i)*c_i(i)); %thermal diffusivity (in m^2/s)
%rho (kg/m^3)

if T_Start(i)>=T_surf
    
Ra(skip)= (rho_i(i)*g*alpha(i)*(T_Start(i)-T_surf)*((Height_list(i)*1000)^3))/(td(i)*v); %Rayleigh

elseif T_Start(i)<T_surf
    
Ra(skip)= (rho_i(i)*g*alpha(i)*(T_surf-T_Start(i))*((Height_list(i)*1000)^3))/(td(i)*v); %Rayleigh

end

% for free convection, the average Nusselt number is expressed as a 
%function of the Rayleigh number and the Prandtl number-Churchill and Chu
%Nu >> 1 , More effective Convection
%Nu = 1, For fluid layer represent heat transfer by pure conduction across layer
%Nu ∼ 1, sluggish convection
%1-10 laminar flow
%#>100 turbulent flow

%Nusselt Calculator

 %h=heat transfer coefficient, L=characteristic length, k=thermal cond.
 
%Nusselt Calculator
if T_Start(i)>T_surf
flux(skip)= (Base_Flux)/(T_Start(i)-T_surf);
elseif T_Start(i)<T_surf
flux(skip)= (Base_Flux)/(T_surf-T_Start(i));
end

%1/3 power law rayleigh estimate, up to ~10^15
%Raest(skip)=Ra(i)^(1/3);

Nu(skip)=(Base_Flux*(Height_list(i)*1000))/k_i(i);
skip=skip+1;
end


Resolution_Profile= [Height_list' Ra' Nu' T_Start];

%add ice evolve temperature and phases to it now?