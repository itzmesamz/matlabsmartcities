% 5 walls: concrete (4 layesrs) + insulation (2 layers)
% glass (1 node)
% ventilation
% Inputs: outdoor temperature, solar radiation, HVAC heat flow rate
% Free running (FR) / without control

clc, clear all
% pkg load control
% Data
% ******
Sc = 5*3*3; Si = Sc; Sg = 3*3; %surface [m2]: concrete, insulation, glass
Va = 3*3*3; %air volume[m3]

rhoa = 1.2; ca = 1000;  %indoor air density; heat capacity
Vpa = 1*Va/3600;        %infiltration and ventilation air: volume/hour

% c: concrete; i: insulation;  g: glass
lamc = 2;       lami = 0.04;    lamg = 1.2;   %[W/m K]
rhoccc = 2.5e6; rhoici = 2e6;   rhogcg = 2e6; %[J/K m3]
wc = 0.2;       wi = 0.08;      wg = 0.01;    %[m]
epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.8;   %short wave wall emmisivity

epsgLW = 0.9;   %long wave glass emmisivity
taugSW = 0.8;   %short wave glass transmittance
alphagSW = 0.2; %short wave glass absortivity

% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]

% MODEL
% *****
nth = 19; nq = 22; % # of temperature node, # of flow nodes

Tm7_15=19+273; Tm8_9=19+273; Tm9_10=19+273; %mean temp for radiative exchange
sigma = 5.67e-8; %[W/m2 K4]
Fwg = 1/5;

% G conductance matrix
G = zeros(nq,nq);
G(1,1)=ho*Sc; G(2,2)=lamc/(wc/8)*Sc; for i=3:9; G(i,i)=G(2,2); end;
G(10,10)=lami/(wi/4)*Si; for i=11:13; G(i,i)=G(10,10); end; G(14,14)=hi*Si;
G(15,15)=epswLW/(1-epswLW)*Si*4*sigma*Tm7_15^3;
G(16,16)=Fwg*Si*4*sigma*Tm8_9^3;
G(17,17)=epsgLW/(1-epsgLW)*Sg*sigma*Tm9_10^3;
G(18,18)=ho*Sg; G(19,19)=lamg/(wg/2)*Sg; G(20,20)=G(19,19); G(21,21)=hi*Sg;
G(22,22)=Vpa*rhoa*ca;

% C capacity matrix
C = zeros(nth);
C(12,12)=1/4*Sc*wc*rhoccc; for i=13:15; C(i,i)=C(12,12); end;
C(16,16)=1/2*Si*wi*rhoici; C(17,17)=C(16,16);
C(18,18)=Sg*wg*rhogcg; C(19,19)=Va*rhoa*ca;

% A adjancy matrix
A = zeros(nq,nth);
A(1,1)=1;
A(2,1)=-1;  A(2,12)=1;
A(3,2)=1;   A(3,12)=-1;
A(4,2)=-1;  A(4,13)=1;
A(5,3)=1;   A(5,13)=-1;
A(6,3)=-1;  A(6,14)=1;
A(7,4)=1;   A(7,14)=-1;
A(8,4)=-1;  A(8,15)=1;
A(9,5)=1;   A(9,15)=-1;
A(10,5)=-1; A(10,16)=1;
A(11,6)=1;  A(11,16)=-1;
A(12,6)=-1; A(12,17)=1;
A(13,7)=1;  A(13,17)=-1;
A(14,7)=-1; A(14,19)=1;
A(15,7)=-1; A(15,8)=1;
A(16,8)=-1; A(16,9)=1;
A(17,9)=-1; A(17,10)=1;
A(18,11)=1;
A(19,11)=-1;A(19,18)=1;
A(20,10)=1; A(20,18)=-1;
A(21,10)=-1;A(21,19)=1;
A(22,19)=1;

%State-space representation
%State-space model
nnodes = size(C,1); %n° total nodes
nC = rank(C);       %n° nodes with capacity
n0 = nnodes - nC;   %n° of nodes with zero capacity

K = -A'*G*A;
K11 = K(1:n0,1:n0);
K12 = K(1:n0,n0+1:end);
K21 = K(n0+1:end,1:n0);
K22 = K(n0+1:end,n0+1:end);

Kb = A'*G;
Kb1 = Kb(1:n0,:);
Kb2 = Kb(n0+1:end,:);

CC = C(n0+1:end,n0+1:end);

As = inv(CC)*(-K21*inv(K11)*K12 + K22);
Bs = inv(CC)*[-K21*inv(K11)*Kb1+Kb2 -K21*inv(K11) eye(nnodes-n0,nnodes-n0)];

%Select relevant inputs and outputs
Bs = Bs(:,[[1 18 22] nq+[1 7 11 19]]); %inputs: [To To To Phiw Phii Phig Qh]
Cs = zeros(1,nC);Cs(nC)=1;  %output
Ds = 0;

% SIMULATION
% **********

n00ReadWeatherFiles;  %read weather files

n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);
u = [Temp Temp Temp ...
  epswSW*Sc*PhiDiff taugSW*epswSW*Sg*PhiDiff alphagSW*Sg*PhiDiff ...
  Qh];

x0 = zeros(nth,n);
sys = ss(As,Bs,Cs,Ds);
[y, t, x] = lsim(sys, u, Time);
plot(Time/3600,y,'g', Time/3600, Temp,'b'), xlabel('Time [h]')
title('Simulation lsim'), ylabel('T [C]')
