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

Kp = 10000; %P-controller gain: large for precision

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
nq = 15;
%% model
nth = 11; nq = 15; % # of temperature node, # of flow nodes

Tm = 20 + 273;   %mean temp for radiative exchange
sigma = 5.67e-8;  %[W/m2 K4]
Fwg = 1/5;        %view factor wall - glass
%% Matrix A
    A1 = zeros(nq,nth);

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

    A1 (9,9)= -1;
    A1 (9,8)= 1;

    A1 (10,10)= -1;
    A1 (10,9)= 1;

    A1 (11,11)= 1;

    A1 (12,8)= 1;

    A1 (13,8)= 1;
    A1 (13,11)= -1;

    A1 (14,11)= 1;
    A1 (15,8)=1;

    %end of matrix A

%% R matrix
% Thermal circuit
%****************

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
R(15,15) = 1/Kp;
G = inv(R);
%% Capacitance
% capacitances
C1 = zeros(nth);
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
%State-space model
nnodes = size(C1,1); %n° total nodes
nC = rank(C1);       %n° nodes with capacity
n0 = nnodes - nC;   %n° of nodes with zero capacity

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
Bs = Bs(:,[[1 11 12 14 15] nq+[1 8]]); %inputs: [To To To Phiw Phii Phig Qh]
Cs = zeros(1,nC);Cs(8)=1;  %output
Ds = zeros(1,7);

% SIMULATION
% **********
n00ReadWeatherFiles;
n = size(Time,1);
th = zeros(nth,n);
Qa = zeros(n,1);  %auxiliary sources (electrical, persons, etc.)
TintSP = 20*ones(n,1);
% Inputs
PhiTot = PhiDir + PhiDif + PhiRef;
u = [Temp Temp Temp Temp TintSP ...
  epswSW*Sc1*PhiTot ...
  Qa];

% Integrate using lsim (linear systems)
x0 = zeros(nth,n);
sys = ss(As,Bs,Cs,Ds);
[y, t, x] = lsim(sys, u, Time);
subplot(2,1,1)
plot(Time/3600,y,'g', Time/3600, Temp,'b'), 
xlabel('Time [h]'), ylabel('T [C]')
Qh = Kp*(TintSP - y); Qh(1)=0; % thermal load
subplot(212), plot(Time/3600,Qh,'r')
xlabel('Time [h]'), ylabel('Q_h_v_a_c [W]')

% Integrate at each time step 0:dt:dt
% needs to iterate
dt = 3600/1; % simulation step 1h/dt = 3600s / dt
th = zeros(size(As,2),n);
for k = 1:n-1
 x0 = th(:,k);
 % x = lsode (@(x, t) ssm(x, t, As, Bs, u(k,:)'), x0, 0:dt:dt);
 [y, t, x] = lsim(sys, u(k:k+1,:), 0:dt:dt,x0);
 th(:,k+1) = x(2,:)';
 th(:,k+1) = min(TintSP(k+1),th(:,k+1)); % limit th_int to set point
 Qhvac(k+1) = Kp*(TintSP(k+1) - th(4,k+1));
 
%  th(:,k) = x(2,:)';
%  th(:,k) = min(TintSP(k),th(:,k)); % limit th_int to set point
%  Qhvac(k+1) = Kp*(TintSP(k) - th(4,k));
end
figure(2)
subplot(211), hold on, plot(Time/3600, th(4,:),'r'), hold off
xlabel('Time [h]'), ylabel('T [C]')
subplot(212), hold on, plot(Time/3600,Qhvac), hold off
xlabel('Time [h]'), ylabel('Q_h_v_a_c [W]')
