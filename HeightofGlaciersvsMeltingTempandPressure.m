close all
clear all
rawTable = readtable('meltpoints.xlsx','Sheet','Sheet1');
x = rawTable.Pressure; %excel column, Header1
%From pressure (corresponds to melting temp) get height...then
rho1=917; %density of ice I
rho2=1170; %density of ice II
rho5=1230; %density of ice V
rho6=1310; %density of ice VI

g=23.7; %m/s^2

P1=3e8; %pressure up to ice I and ice II boundary...300 MPa
P2=2e8; %added pressure up to ice II and ice V boundary...500 MPa
P3=6e8; %added pressure up to ice V and ice VI boundary...1.1 GPa
P4=1.2e9; %added pressure up to ice VI and ice VII boundary...2.3

%H= P/(rho*g)...Pressure=normal stress

for i=1:length(x)
P=x(i)
if (P<= 3e8)
    H=(P/(g*rho1))/1000 
    
elseif (P>3e8) & (P<=5e8)
        H= (P1/(g*rho1)/1000)+(((P-3e8)/(g*rho2))/1000)
        
elseif  (P>5e8) & (P<=1.1e9);
            H= (P1/(g*rho1)/1000)+(P2/(g*rho2)/1000)+(((P-5e8)/(g*rho5))/1000)
            
else % P>1.1e9 && P<=2.3e9
            H= + (P1/(g*rho1)/1000)+(P2/(g*rho2)/1000)+(P3/(g*rho5)/1000)+(((P-1.1e9)/(g*rho6))/1000)
                 
end

list(i)=H

end

Glacier_Heights=list.';

figure('units','normalized','position',[.1 .1 .6 .6])
yyaxis left

x = rawTable.Height; 
y = rawTable.Pressure./1e6;
plot(x,y,'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);
ylim([0 2300])  
xlim([0 85])
xlabel('Height of Glacier (km)')
ylabel('Pressure(MPa)')

hold on 

yyaxis right 
ylabel('Melting Temperature (K)')

x = rawTable.Height; %excel column, Header1
y = rawTable.MeltT; %excel column, Header2 
plot(x,y,'LineWidth',1.9)
set(gca,'FontSize',24)
set(gca,'LineWidth',2,'TickLength',[0.03 0.03]);









