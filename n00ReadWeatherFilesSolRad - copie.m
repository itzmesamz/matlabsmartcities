clear all
clc %clear console
% Load weaher data
M = csvread('FRA_Lyon.074810_IWEC.csv');
from = 5*24;         %start time
period = [1:3*24]';  %simulation period

Time = 3600*period;
% Day = M(from+period,3); Hour = M(from+period,4);
Temp = M(from+period,7);      %Dry bulb temperarure [°C]
RadNDir =  M(from+period,15); %Direct normal solar radiation [Wh/m2]
RadHDif =  M(from+period,16); %Diffuse horizontal solar radiation [Wh/m2]
WDir = M(from+period,21);     %Wind direcrion: N=00; E=90°; S=180°; W=270°
WSpeed = M(from+period,22);   %wind speed [m/s]

month = M(from+period,2); %
day = M(from+period,3); %
hour = M(from+period,4); %
minute = M(from+period,5); %

clear M

B = 90; Z = 0; L = 45; albedo = 0.2;
[PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);
  
%Z = 90;
%[PhiDir1, PhiDif1, PhiRef1] = fSolRadTiltSurf(month, day, hour, minute, ...
% RadNDir, RadHDif, B, Z, L, albedo);
% 
%plot(Time/(24*3600), PhiDir,'b'), hold on %direct on surface
%plot(Time/(24*3600), PhiDir1,'r.'), hold on %direct on surface`

%plot(Time/(24*3600), RadNDir,'g') % direct on normal to sun
plot(Time/(24*3600), PhiDif,'r')  % diffusif on surface
%plot(Time/(24*3600), RadHDif,'k') % diffusif on horizontal surface
%plot(Time/(24*3600), PhiRef,'m')  % reflected on surface
legend('\Phi_d_i_r','\Phi_N_d_i_r','\Phi_D_i_f', '\Phi_N_D_i_f')
