clc, clear allnth = 7; nq = 7;  % # of temperature node, # of flow nodesdh = 2;       %if dh = 1.98,stable if dh = 1.97, instabledt = 3600/dh; % simulation step 1h/dt = 3600s / dtn = dh*24*30; % n° of time samples for 30 daysTime = 0:dt:(n-1)*dt; % timeth = zeros(nth,n);u = [ones(1,n); zeros(1,n)];% Integrate at each time step 0:dt:dt% needs to itteratefor k = 1:n-1 x0 = th(:,k); x = lsode (@(x, t) sw02(x, t, u(:,k)), x0, 0:dt:dt); th(:,k+1) = x(2,:)';endplot(Time, th(7,:))% Integrate for all periode Time% needs to pass all inputs and interpolate the inputs for the current time tx0 = zeros(nth,1);x = lsode (@(x, t) sw03(x, t, Time, u), x0, Time);hold on, plot(Time, x(:,7),'.r'), hold off