% Created by Farai Chirove
% 3MC workshop Northwest University
% Contour plots and surface plots for R_0
clc
close all

% step size and time span
h = 0.0001;
tspan = 1:h:10;

%parameter values 
lambda= 20;
mu = 0.02;
%beta = 1.7714e-05;
%delta = 0.1;
N = 500;
c=2.4;
betav = 0.002;

U1 = 1.7714e-05:0.001:0.1;
U2 = 0.1:0.001:0.2;
U3 = 0:0.001:1;

for i = 1:length(U1)
    for j = 1:length(U2)
        
        u1 = U1(i); u2 = U2(j);
        R0(i,j) = N*u1*u2*lambda/((u2+mu)*(betav*lambda+c*mu));

         %R0(i,j) = ((1-u1)*(1-u2)*lambda*delta)/(mu*(mu+((1-u2)*gamma*phi*tau)/(ps*T)));
        %R0(i,j) = ((1-u1)*(1-u2)*lambda*beta)/(mu*(mu+d+(gamma*vartheta*tau)/(pi*T)));
        
    end
end

[X,Y] = meshgrid(U1,U2);

Ro = R0; P=X; Q=Y;

figure(1)
%subplot(1,2,1)
mesh(P,Q,Ro'), hold on
[C,h] = contour(P,Q,Ro','color','k'); axis tight

 text_handle1 = clabel(C,h);
 set(text_handle1,'BackgroundColor',[1 1 .6],...
     'Edgecolor',[.7 .7 .7])
% Xlabel('\beta_1'), Ylabel('\beta_2'), Zlabel('R_0')
% Xlabel('\phi_1'), Ylabel('\phi_2'), Zlabel('R_0')
title('Change in $\mathcal{R}_{0}$ with $\beta$ and $\delta$', 'Interpreter', 'Latex','FontSize',20)
xlabel('$\beta$', 'Interpreter', 'Latex','FontSize',20),
ylabel('$\delta$', 'Interpreter', 'Latex','FontSize',20)
zlabel('$\mathcal{R}_{0}$', 'Interpreter', 'Latex','FontSize',20)
set(gcf, 'PaperUnits', 'inches')
papersize = get(gcf, 'PaperSize');
width = 6;         % Initialize a variable for width.
height = 6;          % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);
opt = struct('FontSize',20);
set(get(gcf,'CurrentAxes'), opt)
 print -depsc FigureR0.eps

figure(2)
%subplot(1,2,2)
[C,h] = contourf(P,Q,Ro');
text_handle = clabel(C,h);
set(text_handle,'BackgroundColor',[1 1 .6],...
   'Edgecolor',[.7 .7 .7])
% Xlabel('\beta_1'), Ylabel('\beta_2')
% Xlabel('\phi_1'), Ylabel('\phi_2')
title('Level lines','FontSize',20)
xlabel('$\beta$', 'Interpreter', 'Latex','FontSize',20) 
ylabel('$\delta$', 'Interpreter', 'Latex','FontSize',20)
set(gcf, 'PaperUnits', 'inches')
papersize = get(gcf, 'PaperSize');
width = 6;         % Initialize a variable for width.
height = 6;          % Initialize a variable for height.
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);
opt = struct('FontSize',20);
set(get(gcf,'CurrentAxes'), opt)
title('b','FontSize',20)
 print -depsc FigureRV.eps
