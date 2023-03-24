function HIV_Prototypecode 
x0 = [ 1000 0 0.001]; 
tspan = [0 1000]; 
options = odeset('Refine',10, 'RelTol',1e-4);
[t,x] = ode45(@HIV, tspan, x0, options);
  
%% Plot
figure(1)
hold on
   plot(t,x(:,1),'b','linewidth',2)
   plot(t,x(:,2),'r','linewidth',2)
   legend('T','I')
   title ('CD4^+ populations')
xlabel('time in days')
ylabel('Population size')

figure(2)
   plot(t,x(:,3),'--r','linewidth',2)
   legend('V')
   title ('Virus population')
xlabel('time in days')
ylabel('Population size')

figure(3)
   plot(t,x,'linewidth',2)
   legend('T','I','V')
   title ('All populations')
xlabel('time in days')
ylabel('Population size')
  end

%initial values for model parameters 
function dx = HIV(~, x)

lambda= 20;
mu = 0.02;
beta = 1.7714e-05;
delta = 0.1;
N = 500;
c=2.4;
betav = 0.00002;
%% Thresholds
R0 = N*beta*delta*lambda/((delta+mu)*(betav*lambda+c*mu));
betastar = ((c+350*betav)*(delta))/(350*(delta*N));
%% Model equations 
dx = zeros(3,1); % 
%% systems of differential equations
dx(1) = lambda - mu*x(1) - beta*x(1)*x(3);
dx(2) = beta*x(1)*x(3) - (mu+delta)*x(2);
dx(3) = delta*N*x(2) - c*x(3) - betav*x(1)*x(3);





end