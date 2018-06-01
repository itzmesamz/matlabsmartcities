% Simple wall with capacities in all temperature nodes
% Inputs: outdoor temperature, indoor convection heat flow rate (from HVAC))
clc, clear all
pkg load control

% Physical properties
% *******************
Sw = 3*3;   %wall surface [m2]
Va = 3*3*3; %air volume[m3]
% concrete:1 ; insulation:2
lam1 = 2;   lam2 = 0.04; %[W/m K]
rho1c1 = 2.5e6; rho2c2 = 2e6; rho3c3 = 1.2e3;  %[J/K m3]
w1 = 0.2;   w2 = 0.08;  %[m]
x1 = 0.05;  x2 = 0.04;  %[m]
% convection coefficents
hi = 4; ho = 10;  %[W/m2 K]
% Thermal resistances
% concrete
Rc = w1/(lam1*Sw);   Cc = Sw*w1*rho1c1;
% insulation
Ri = w2/(lam2*Sw);   Ci = Sw*w2*rho2c2;
% convection
Rvi = 1/(hi*Sw); Rvo = 1/(ho*Sw);

% Thermal circuit
%****************
nth = 7; nq = 7;  % # of temperature node, # of flow nodes
% resistances
R = zeros(nq,nq);
R(1,1) = Rvo + Rc/8; 
R(2,2) = Rc/4; R(3,3)=R(2,2); R(4,4)=R(2,2);
R(5,5) = Rc/8 + Ri/4; 
R(6,6) = Ri/2; 
R(7,7) = Ri/4 + Rvi;
G = inv(R);
% capacitances
C = zeros(nth,nth);
C(1,1) = 1/4*Sw*w1*rho1c1; C(2,2)=C(1,1); C(3,3)=C(1,1); C(4,4)=C(1,1);
C(5,5) = 1/2*Sw*w2*rho2c2; C(6,6)=C(5,5);
C(7,7) = Va*rho3c3;
% adjancecy matrix
A = eye(nq+1,nth);
A = -diff(A,1,1)';
% steady-state with
b = [1 0 0 0 0 0 0]'; f = [0 0 0 0 0 0 0]';
thsteadyTo = inv(A'*G*A)*(A'*G*b + f);
b = [0 0 0 0 0 0 0]'; f = [0 0 0 0 0 0 1]';
thsteadyQi = inv(A'*G*A)*(A'*G*b + f);
[thsteadyTo(7) thsteadyQi(7)]

% State-space representation
B = inv(C)*[A'*G eye(nth,nth)]; % inputs u = [b; f] size(b)=nq, size(f)=nth;
B = B(:,[1 14]);        % select the 2 relevant inputs: 1->To and 14->Qh
A = inv(C)*(-A'*G*A); 
C = zeros(1,7);C(7)=1;  % output: th(7)

% Time integration using Euler forward
% ************************************

% Step response
dh = 1.97;       %if dh = 1.98,stable if dh = 1.97, instable
dt = 3600/dh; % simulation step 1h/dt = 3600s / dt
n = floor(dh*24*30); % n° of time samples for 30 days
Time = 0:dt:(n-1)*dt; % time
th = zeros(nth,n); % Euler explicit
thi = zeros(nth,n); % Euler implicit
u = [ones(1,n); zeros(1,n)];
for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*A)*th(:,k) + dt*B*u(:,k);
 thi(:,k+1) = inv((eye(nth) - dt*A))*(thi(:,k) + dt*B*u(:,k));
end
subplot(2,2,1)
plot(Time/3600,th(7,1:n),'r', Time/3600, thi(7,1:n),'k'), xlabel('Time [h]')
title('Step response for T_o = 1 C'), ylabel('T [C]')

u = [zeros(1,n); ones(1,n)];
th = zeros(nth,n); % Euler explicit
thi = zeros(nth,n); % Euler implicit
for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*A)*th(:,k) + dt*B*u(:,k);
 thi(:,k+1) = inv((eye(nth) - dt*A))*(thi(:,k) + dt*B*u(:,k));
end
subplot(2,2,2)
plot(Time/3600,th(7,1:n),'r', Time/3600, thi(7,1:n),'k'), xlabel('Time [h]')
title('Step response for Q_h = 1 W'), ylabel('T [C]')

%Stability
lambda = eig(A);
min(min(lambda))*dt; %needs to be in [-2 0]
disp('Stability: min eigenvalue of A * dt in [-2 0]'), disp(min(min(lambda))*dt)
disp('lam/(rho*c)/w^2 :')
[(lam1/rho1c1)/(w1/4)^2*dt lam2/rho2c2/(w2/2)^2*dt] %needs to be <1/2

% Simulation with outoor temperature
n00ReadWeatherFiles;  %read weather files

Temp = interp1(Time, Temp, [Time(1):dt:Time(end)]'); %interpolate for dt
Time = [Time(1):dt:Time(end)]';

n = size(Time,1);
th = zeros(nth,n);
Qh = zeros(n,1);
u = [Temp'; Qh'];
for k = 1:n-1
 th(:,k+1) = (eye(nth) + dt*A)*th(:,k) + dt*B*u(:,k);
end
subplot(2,2,3)
plot(Time/3600,th(7,1:n),'r', Time/3600, Temp,'b'), xlabel('Time [h]')
title('Simulation Euler'), ylabel('T [C]')

disp('Mean temperatures In Out'), disp([mean(th(7,1:n)) mean(Temp)])

% Integration with library functions
x0 = zeros(nth,n);
sys = ss(A,B,C);
[y, t, x] = lsim(sys, u', Time);
hold on, plot(t/3600,y,'g.'), hold off
subplot(2,2,4)
plot(Time/3600,y,'g', Time/3600, Temp,'b'), xlabel('Time [h]')
title('Simulation lsim'), ylabel('T [C]')