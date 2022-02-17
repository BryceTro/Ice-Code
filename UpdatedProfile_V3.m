clc;
close all;
clear all;
addpath '/Users/luju/Dropbox/Mac/Documents/LHS1140B'

global start_height
global resolution
global final_height

%Input the below values to run
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
minP=0; %MPa
maxP=2300; %MPa
minT=200; %K
maxT=400;

T_surf=235; %Surface Temp of planet (in K)
Base_Flux=30/1000;  %Heat Flux (in W/m^2)
Timescale=10*1E6; %ex. 1,000,000 year timescale; 10,000,000 year timescale;100,000,000 timescale etc.
start_height=0; %from 0 km onward
resolution=0.1;   %height of each step in km (can be decimals)
final_height=80 ; %up to 80 km
% res = 100;
% g=24; %gravity of planet (in m/s^2)
Patmosphere=0.101325; %atmospheric pressure of planet (in MPa)
Mx = 6.38; %Mass factor of the planet 
Rx = 1.64; %Radius factor of the planet

%Get an initial adiabatic profile for the ice sheet and depth dependent
%melting temperature using the two subroutines below
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
[T_Start,dtdz,Pressure,rho,alpha,Cp,K,phase] = adiabat_profile(Patmosphere,T_surf,(final_height-start_height+1)*1000,resolution,Mx,Rx);
step=(maxP-minP+1)/numel(K);
[MeltT] = findmeltT(minP,maxP,minT,maxT,step,Pressure);
MeltT = MeltT';
Height_list=linspace(0,80,numel(T_Start));

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%Do you want to plot the intial profiles of parameters?%
%plotpaam = 1 plots the relevant initial profiles
plotparam = 1;
if plotparam == 1
figure('units','normalized','position',[.1 .1 .5 .6])
subplot(2,2,1)
plot(Height_list,T_Start,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Temperature (K)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 80])
subplot(2,2,2)
plot(Height_list,rho,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Density (kg m^{-3})')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 80])
subplot(2,2,3)
plot(Height_list,Cp,'k','LineWidth',1.9)
xlabel('Depth (km)')
ylabel('Cp (J kg^{-1} K^{-1})')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
xlim([0 80])
subplot(2,2,4)
scatter(Height_list,phase,'k','LineWidth',1.9)
xlabel('Depth (km)')
% ylabel('K (W m^{-1} K^{-1})')
ylabel('Phase')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
box on
xlim([0 80])
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

dt=(86400*365.25*Timescale)/1000;
Phi_Start=0*ones(1,length(Height_list))';      %Porosity Start
k_i=K';          %Thermal Conductivites of Ices
rho_i=rho;                 %Densities of Ices
c_i=Cp;             %Specific Heats of Ices
k_w=0.6*ones(1,length(Height_list))';          %Thermal Conductivity of Water
rho_w=1000*ones(1,length(Height_list))';       %Density of Water
c_w=4180*ones(1,length(Height_list))';     %Specific Heat Water
Tm=MeltT;        %Melting Temp

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
%Work on making an ice phase plot showing progressive change with depth and time
skip = 1;
ss = 0;
for time=[0:dt:1000*dt]
    i=i+1;
    if rem(i,100)==0
        ss = ss+1;
        PT = [Pressure T_Start];
        phasenew = SF_WhichPhase(PT);
        [out,oo] = compute_params(PT,phasenew);
        k_i =oo;
        rho_i = out.rho;
        c_i = out.Cp;
        hold_phase=[hold_phase phasenew];
        hold_temps=[hold_temps T_new];
        hold_phi=[hold_phi Phi_new];
        test2(:,ss)  = phasenew';
    end

    [T_new,Phi_new]=HP_Ice_Evolve_v3(T_Start,Phi_Start,k_i,rho_i,c_i,k_w,rho_w,...
        c_w,dt,dz,T_surf,Base_Flux,Tm,L,TTol,PhiTol);
    T_Start=T_new;
    T_Start(1)=T_surf;
    test(:,skip)  = T_Start';
    Phi_Start=Phi_new';
    [MeltT] = findmeltT(minP,maxP,minT,maxT,step,Pressure);
    Tm = MeltT';
    Tmtest(:,skip) = Tm;
    skip = skip+1;
end

time = 86400*365.25*Timescale./dt;
time = time./100000;
figure('units','normalized','position',[.1 .1 .3 .6])
ax(1) = subplot(2,1,1)
imagesc([time:time:1000*time],[Height_list],([test]))
hh = colorbar
ylabel(hh,'Temperature (K)','FontSize',24)
xlabel('Time (Myr)')
ylabel('Depth (km)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
colormap hot
% title("T_s = " + Ts + "K, " +"q_s = " + q + "mW m^{-2}") %---------_Change this value
ax(2) = subplot(2,1,2)
imagesc([time:time:1000*time],[Height_list],([test2]))
cb=colorbar
xlabel('Time (Myr)')
ylabel('Depth (km)')
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
ylabel(cb,'Ice Phase','FontSize',24)
colormap parula
% cmap = [0 0 0; lbmap(1,'RedBlue')];
% colormap(ax(2),flipud(cmap))
%  cb.Ticks = linspace(0, 1, 2) ; %Create 8 ticks from zero to 1
%  cb.TickLabels = num2cell(0:1) ;



