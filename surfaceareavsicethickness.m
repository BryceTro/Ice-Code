clc;
clear all;
close all;
Mearth = 5.972e24;  %kg
Mlhs = Mearth*6.98;
Rearth= 6371;   %km
Rlhs=1.727*Rearth*1000; %km to m
%density of ice~.917 g/cm^3 to 917 kg/m^3

i = 1
for w = [1e-4:1e-4:7e-4];      %X axis-Water Mass Fraction of .0001 to .00075 (up to 1.5 times water content of Earth)
    for s = [0:.01:1];           %Y axis-Surface Area Percentage of 0-100%
        d(i)= ((w.*((Mlhs))/917))/(s.*(4*pi*(Rlhs)^2))  %Depth=Volume/Surface Area
        i=i+1;
    end
end

dd= reshape(d,101,[]);
dd=dd/1000
imagesc([1e-4:1e-4:7e-4],[0:.01:1]*100, dd)
caxis([0 225]) %use line for narrowing results
colorbar
xlabel('Water Mass Fraction (%)')
ylabel('Surface Area Covered with Ice/Water (%)')
title('Uniform Ice/Water Thickness (km)')
