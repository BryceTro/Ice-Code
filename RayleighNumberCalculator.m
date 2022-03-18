%Newtonian Convection

function [Ra,Nu,td,tau,shearv,gamma,kv] = RayleighNumberCalculator(T_Start,T_surf,Base_Flux,g,Rx,k_i,rho_i,c_i,alpha,Pressure,Height_list);
g=g;
T_surf=T_surf; %K
Height_list=Height_list; %km
k_i=k_i;
rho_i=rho_i;
c_i=c_i;
Rx=Rx;
Pressure=Pressure;
Base_Flux=Base_Flux;
f=.8;   %shape factor

skip=1;
for i=1:length(Height_list)

theta(skip)=f*atand(Height_list(i)/1); %slope angle

%^From 'Estimating the volume of glaciers in the Himalayan–Karakoram region'
       
td(skip)= k_i(i)/(rho_i(i)*c_i(i)); %thermal diffusivity (in m^2/s)
%rho (kg/m^3)

tau(skip)= (Pressure(i)*1e6)*(sind(theta(i))); %basal shear (in (kg⋅m/s^2)/m^2)

%how can calculate shearing on bottom of glacier
%^think this is now resolved

shearv(skip)=(tau(i)/rho_i(i))^(1/2); %shear velocity (in m/s)

gamma(skip)=((tau(i)*Height_list(i))/shearv(i))*1e3; %dynamic viscosity (in kg/m*s)
%^fix this
kv(skip)=gamma(i)/rho_i(i); %kinematic viscosity (in m^2/s)

if T_Start(i)>=T_surf
    
Ra(skip)= (rho_i(i)*g*alpha(i)*(T_Start(i)-T_surf)*((Height_list(i)*1000)^3))/(td(i)*kv(i)); %Rayleigh

elseif T_Start(i)<T_surf
    
Ra(skip)= (rho_i(i)*g*alpha(i)*(T_surf-T_Start(i))*((Height_list(i)*1000)^3))/(td(i)*kv(i)); %Rayleigh
end

% for free convection, the average Nusselt number is expressed as a 
%function of the Rayleigh number and the Prandtl number-Churchill and Chu
%Nu >> 1 , More effective Convection
%Nu = 1, For fluid layer represent heat transfer by pure conduction across layer

%Nusselt Calculator

 %h=heat transfer coefficient, L=characteristic length, k=thermal cond.
 
%Nusselt Calculator
if T_Start(i)>T_surf
flux(skip)= (Base_Flux)/(T_Start(i)-T_surf);
elseif T_Start(i)<T_surf
flux(skip)= (Base_Flux)/(T_surf-T_Start(i));
end
%Var1(skip)=.663*(Ra(i).^25)
%Var2(skip)=(1+((.492/Pr(i)).^(9/16))).^(4/9)

%Nu(skip)=.68+(Var1(i)./Var2(i))

Nu(skip)=(Base_Flux*(Height_list(i)*1000))/k_i(i);
skip=skip+1;
end

Resolution_Profile= [Height_list' Ra' Pr' Nu' T_Start];

%add ice evolve temperature and phases to it now?