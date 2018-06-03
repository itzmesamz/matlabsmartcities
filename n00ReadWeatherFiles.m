% clear all
clc %clear console
% Load weaher data
% M = csvread('FRA_Lyon.074810_IWEC.csv'); % if the header is removed

filename = 'FRABordeaux075100IWEC.txt';
delimiterIn = '\t';
hederlinesIn = 8;
M = importdata(filename,delimiterIn, hederlinesIn);
M = M.data;

%M = dlmread('FRA_Lyon.csv');

from = 25*24;             %start time
period = [1:31*24]';  %simulation period

Time = 3600*period;
% Day = M(from+period,3); Hour = M(from+period,4);
% Temp = M(from+period,7-1); %Dry bulb temperarure [째C]
% PhiDirN =  M(from+period,15-1); %Direct horizontal solar radiation [Wh/m2]
% PhiDiff =  M(from+period,16-1); %Diffuse solar radiation [Wh/m2]
% WDir = M(from+period,21-1);   %Wind direcrion: N=00; E=90째; S=180째; W=270째
% WSpeed = M(from+period,22-1); %wind speed [m/s]
%%  radition on walls
% south wall
%Albedo of aged concrete is 0.2
Temp = M(from+period,7-1);      %Dry bulb temperarure [캜]
RadNDir =  M(from+period,15-1); %Direct normal solar radiation [Wh/m2]
RadHDif =  M(from+period,16-1); %Diffuse horizontal solar radiation [Wh/m2]


month = M(from+period,2); %
day = M(from+period,3); %
hour = M(from+period,4); %
minute = M(from+period,5); %
%%
%south Wall 
B = 90; Z = 0; L = 44.83; albedo = 0.2;
[PhiDir0, PhiDiff0,PhiRef0] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);

%west Wall
B = 90; Z = 180; L = 44.83; albedo = 0.2;
[PhiDir1, PhiDiff1,PhiRef1] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);

%east wall
B = 90; Z = -180; L = 44.83; albedo = 0.2;
[PhiDir2, PhiDiff2,PhiRef2] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);


%south window 
B = 90; Z = 0; L = 44.83; albedo = 0.05;
[PhiDir3, PhiDiff3,PhiRef3] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);

%west window
B = 90; Z = 180; L = 44.83; albedo = 0.05;
[PhiDir4, PhiDiff4,PhiRef4] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);

%east window
B = 90; Z = -180; L = 44.83; albedo = 0.2;
[PhiDir5, PhiDiff5,PhiRef5] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);
%% Total radiation
PhiDir= PhiDir0+PhiDir1+PhiDir2+PhiDir3+PhiDir4+PhiDir5;
PhiDiff = PhiDiff0+PhiDiff1+PhiDiff2+PhiDiff3+PhiDiff4+PhiDiff5;
PhiRef = PhiRef0+ PhiRef1+PhiRef2+ PhiRef3+ PhiRef4+PhiRef5;
WDir = M(from+period,21-1);     %Wind direcrion: N=00; E=90째; S=180째; W=270�
WSpeed = M(from+period,22-1);   %wind speed [m/s]

[PhiDir, PhiDif, PhiRef] = fSolRadTiltSurf(month, day, hour, minute, ...
  RadNDir, RadHDif, B, Z, L, albedo);

clear M