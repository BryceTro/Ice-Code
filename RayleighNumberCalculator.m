%Newtonian Convection

function [Ra,Pr,Nu] = RayleighNumberCalculator(T_s,Rx,k_i,rho_i,c_i,Pressure,radius,Height_list);

R= 8.3144598; %universal gas constant J/mol/K
radius=radius;  %radius (km) of Planet
Area= 4*pi*radius^2;  %Surface Area of Planet km^2
T_Surface=235; %K
Height_list=h; %km
k_i=k;
rho_i=rho;
c_i=Cp;
Rx=Rx;
Pressure=Pressure;

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
