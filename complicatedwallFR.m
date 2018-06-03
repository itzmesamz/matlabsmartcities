close all;
clear all;
clc;
%% FREE RUN

% Physical properties
% *******************

Sc1 = (3.6576+0.19+15.24)*3; Si1 = Sc1; Sg = 13.24*1; %surface [m2]: concrete, insulation, glass
Va = 15.24*3*(3.6576+0.19)/2; %air volume[m3]
wc=0.20;
wi=0.08;
wg=0.05;

% concrete:c ; insulation:i; glass:g
lamc =2;   lami =0.039; lamg =0.96; %[W/m K] thermal conductivity
Rhoc= 2400; Rhoi= 19;    Rhog=2.6; %[kg/m3]
Cc=750; Ci=1700;    Cg=840; %[J/kg*K]

epswLW = 0.9;   %long wave wall emmisivity
epswSW = 0.8;   %short wave wall emmisivity

epsgLW = 0.9;   %long wave glass emmisivity
taugSW = 0.8;   %short wave glass transmittance
alphagSW = 0.2; %short wave glass absortivit

RhocCc = Rhoc*Cc; RhoiCi= Rhoi*Ci; RhogCg = Rhog*Cg;  %[J/K m3]
rho3c3 = 1.2e3; %rho air capacity air

% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]
% Thermal resistances
% concrete
Rc = wc/(lamc*Sc1);   
Cc = Sc1*wc*RhocCc;
% insulation
Ri = wi/(lami*Si1);   
Ci = Si1*wi*RhoiCi;
%glass
Rg = wg/(lamg*Sg);
% convection
Rvi = 1/(hi*Si1); Rvo = 1/(ho*Sc1);

dt = 3600/15;
ntemp =11;
nq = 14;
%% Matrix A
    A1 = zeros(nq,ntemp);

    A1(1,1)= 1;

    A1 (2,1)= -1;
    A1 (2,2)= 1;

    A1 (3,2)= -1;
    A1 (3,3)= 1;

    A1 (4,3)= -1;
    A1 (4,4)= 1;

    A1 (5,4)= -1;
    A1 (5,5)= 1;

    A1 (6,5)= -1;
    A1 (6,6)= 1;

    A1 (7,6)= -1;
    A1 (7,7)= 1;

    A1 (8,7)= -1;
    A1 (8,8)= 1;

    A1 (9,8)= -1;
    A1 (9,9)= 1;

    A1 (10,9)= -1;
    A1 (10,10)= 1;

    A1 (11,10)= -1;

    A1 (12,8)= 1;

    A1 (13,8)= 1;
    A1 (13,11)= -1;

    A1 (14,11)= 1;

    %end of matrix A
%% R matrix
% Thermal circuit
%****************
nth = 11; nq = 14;  % # of temperature node, # of flow nodes
% resistances
R = zeros(nq,nq);
R(1,1) = Rvo;
R(2,2) = Rc/8; 
R(3,3) = Rc/4; R(4,4)=R(3,3); R(5,5)=R(3,3);
R(6,6) = Rc/8 + Ri/4; 
R(7,7) = Ri/2; 
R(8,8) = Ri/4 + Rvi;
R(9,9) = Rg/2;
R(10,10) = Rg/2+Rvi;
R(11,11) = Va*rho3c3;
R(12,12) = Va*rho3c3;
R(13,13) = Rg/2;
R(14,14) = Rg/2+Rvo;

G = inv(R);
%% Capacitance
% capacitances
C1 = zeros(nth,nth);
C1(2,2)=Cc/4; 
C1(3,3)=C1(2,2); 
C1(4,4)=C1(2,2);
C1(5,5) = C1(2,2);
C1(6,6)=Ci/2;
C1(7,7) = C1(6,6);
C1(8,8)= Va*rho3c3;
C1(9,9)= Ci;
C1(10,10)= Va*rho3c3;
C1(11,11) = Ci;
%% % State-space representation
nnodes = size(C1,1); %nÂ° total nodes
nC = rank(C1);       %nÂ° nodes with capacity
n0 = nnodes - nC;   %nÂ° of nodes with zero capacity

K = -A1'*G*A1;
K11 = K(1:n0,1:n0);
K12 = K(1:n0,n0+1:end);
K21 = K(n0+1:end,1:n0);
K22 = K(n0+1:end,n0+1:end);


Kb = A1'*G;
Kb1 = Kb(1:n0,:);
Kb2 = Kb(n0+1:end,:);

CC = C1(n0+1:end,n0+1:end);

As = inv(CC)*(-K21*inv(K11)*K12 + K22);
Bs = inv(CC)*[-K21*inv(K11)*Kb1+Kb2 -K21*inv(K11) eye(nnodes-n0,nnodes-n0)];

%Select relevant inputs and outputs
Bs = Bs(:,[[1 10 11] nq+[1 7 10 8]]); %inputs: [To To To Phiw Phii Phig Qh]
Cs = zeros(1,nC);Cs(nC-3)=1;Cs(nC)=1;  %output
Ds = 0;

% STEP RESPONSE
% *************

dh = 200;
dt = 3600/dh; % simulation step 1h/dt = 3600s / dt
n = dh*24*30; % n° of time samples for 30 days
n = floor(n);
Time1 = 0:dt:(n-1)*dt; % time
th = zeros(nth-n0,n);
B = Bs (:, [1 2]);
u1 = [ones(1,n); zeros(1,n)];
for k = 1:n-1
 th(:,k+1) = (eye(nth-n0) + dt*As)*th(:,k) + dt*B*u1(:,k);
end
subplot(2,1,1)
plot(Time1/3600,th(7,1:n),'r'), xlabel('Time [h]')
xlim([0 200]);
title('Step response for T_o = 1 C'), ylabel('T [C]')

u1 = [zeros(1,n); ones(1,n)];
for k = 1:n-1
 th(:,k+1) = (eye(nth-n0) + dt*As)*th(:,k) + dt*B*u1(:,k);
end
subplot(2,1,2)
plot(Time1/3600,th(7,1:n),'r'), xlabel('Time [h]')
xlim([0 200]);
title('Step response for Q_h = 1 W'), ylabel('T [C]')

clear dh dt n Time1 th B u1 

%% SIMULATION
% 
% **********

n00ReadWeatherFiles;  %read weather files
n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);
TintSP = 20*ones(n,1);
% Inputs
u = [Temp Temp Temp ...
  epswSW*Sc1*PhiDiff taugSW*epswSW*Sg*PhiDiff alphagSW*Sg*PhiDiff ...
  Qh];

  % Integrate using lsim (linear systems)
x0 = zeros(nth,n);
sys = ss(As,Bs,Cs,Ds);
[y, t, x] = lsim(sys, u, Time);
figure (2)
subplot(2,1,1)
plot(Time/3600,y,'g', Time/3600, Temp,'b'), 
xlabel('Time [h]'), ylabel('T [C]')
