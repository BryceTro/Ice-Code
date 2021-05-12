%Pressure in Pa
Pinitial=500000 %lowest pressure
P0=1.5e8 %mid pressure value in ice I for graph
P1=3e8 %added pressure up to ice I and ice II boundary...300 MPa
P2=2e8 %added pressure up to ice II and ice V boundary...500 MPa
P3=6e8 %added pressure up to ice V and ice VI boundary...1.1 GPa
P4=1.2e9 %added pressure up to ice VI and ice VII boundary...2.3 GPa
P5=59.7e9 %added pressure up to ice VII and X boundary...62 GPa

rho1=917 %density of ice I
rho2=1170 %density of ice II
rho5=1230 %density of ice V
rho6=1310 %density of ice VI
rho7=1650 %density of ice VII
g=23.7 %m/s^2

%P=rho*g*h...Pressure=normal stress
%Height given in km
Hstart=0
Pstart=0
Hinitial=(Pinitial/(g*rho1))/1000 %height at first height input
H0= (P0/(g*rho1))/1000 %height at ice 1 input 
H1= (P1/(g*rho1))/1000 %height of ice glacier up to ice I and II boundary
H2= H1+((P2/(g*rho2))/1000) %height of ice glacier up to ice II and ice V
H3= H2+((P3/(g*rho5))/1000) %height of ice glacier up to ice V and Ice VI boundary
H4= H3+((P4/(g*rho6))/1000) %height of ice glacier up to ice VI and Ice VII boundary
H5= H4+((P5/(g*rho7))/1000) %height of ice glacier up to ice VII and Ice X boundary

yyaxis left

plot(([Pstart Pinitial P0 P1 P1+P2 P1+P2+P3 P1+P2+P3+P4 P1+P2+P3+P4+P5]), [Hstart Hinitial H0 H1 H2 H3 H4 H5])
xlim([0 2.3e9])  
xlabel('Pressure of Glacier (Pa)')
ylabel('Height of Glacier (km)')

hold on 

yyaxis right 
ylabel('Melting Temperature (K)')
rawTable = readtable('meltpoints.xlsx','Sheet','Sheet1');
x = rawTable.Pressure; %excel column, Header1
y = rawTable.MeltT; %excel column, Header2 
plot(x,y)


