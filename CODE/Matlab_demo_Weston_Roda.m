%This set of codes is prepared by Mr. Weston Roda, PhD student 
%at the University of Alberta, 2018

%Contents

%1. Vectors and matrices
%2. Plotting functions in 2D
%3. Newton's method
%4. Plotting functions in 3D
%5. Plotting a list of points
%6. Contour plot
%7. Basic matrix operations
%8. Eigenvalues and eigenvectors
%9. Solving nonlinear equations
%10. Solving ODEs
%11. Plotting solutions of ODEs in time plots
%12. 2D Phase Plots
%13. Solving PDE heat equation and plotting solution as a surface
%14. Solving PDE advection equation and plotting solutions

%%

%Part 1
%Vectors and matrices

x=[1 3 2]               %create a row vector
x=[1; 3; 2]             %create a column vector
x=0:2:10                %automatically generate vector

A=[1 3 2; 5 6 7; 8 3 1] %create a 3 by 3 matrix

A(2,3)                %the 2,3 entry of A
A(2,:)                % the second row of A
A(1:2,1:2)            % the first 2by2 block of A

%%

%Part 2
%Plotting functions in 2D

clear all % clear all variables 
close all % close all figures

x=-10:0.1:10; % define plotting range (x column)
y=sin(x); % define function to plot (y column)
plot(x,y) % plot the function (x and y columns)

%Alternatively, can create the same plot with "fplot" function

clear all % clear all variables 
close all % close all figures

fplot(@(x) sin(x), [-10 10])

%more plots

clear all % clear all variables 
close all % close all figures

x1=linspace(-10,10,64); % generates 64 evenly spaced points between -10 and 
                        %10
y1=x1.*sin(x1); % the operator '.*' multiplies the vectors element-wise
plot(x1,y1)

%Plot in the same frame, with different color

clear all % clear all variables 
close all % close all figures

x=-10:0.1:10; % define plotting range (x column)
y=sin(x);  %define function (y column)

x1=linspace(-10,10,64); 
y1=x1.*sin(x1);

plot(x,y,'g*',x1,y1,'mo:') %'g*' plots green *'s
                           %'mo:' plots magenta points connected by dashed
                           %lines


%Plot with labels

clear all % clear all variables 
close all % close all figures

x=-10:0.1:10; % define plotting range (x column)
y=sin(x);  %define function (y column)

x1=linspace(-10,10,64); 
y1=x1.*sin(x1);

plot(x,y,'g*',x1,y1,'mo:')
xlabel('x values');
ylabel('y values');
title('Example Graph');
%%

%Part3
%Newton's method

clear all   % clear all variables
close all % close all figures

x(1)=0.5;    % initial guess
for n=1:1000 % n is the iteration variable
x(n+1)=x(n)-((x(n)*sin(x(n)))/(sin(x(n))+x(n)*cos(x(n))));
fc=x(n+1)*sin(x(n+1));
display(fc)
if abs(fc)<10^(-8)
break       % quit the loop
end
end
x(n+1)      % print value of root, 5.7592e-05
fc          % print value of function (error), 3.3168e-09

%Alternatively, use "diff" function to differentiate the function f(x)
%in Newton's method

clear all % clear all variables 
close all % close all figures

syms x;
f(x) = x*sin(x);
df = diff(f, x);

F(x) = x - f(x)/df(x);

F(0.5)  %1/2 - sin(1/2)/(2*(cos(1/2)/2 + sin(1/2)))

%vpa function uses variable-precision floating point arithmetic (VPA)

vpa(F(0.5)) %0.2389366387126959524848716084988

k = 0.5;

for i = 1:1000
    
    y = vpa(F(k));  
    k = y;
    
    if abs(f(k))<10^(-8)
        break
    end
end

k     % print value of root, 5.7592e-05
f(k)  % print value of function (error), 3.3168e-09
%%

%Part 4
%Plotting functions in 3D

%Create an elliptic paraboloid
%z = 4 x^2 + y^2

clear all % clear all variables
close all % close all figures

[X, Y] = meshgrid(-20:0.5:20);

N = size(X);

Z = zeros(N(1), N(2));

F = @(x, y) 4*(x^2) + (y^2);

for i = 1:N(1)
    for j = 1:N(2)
        
        Z(i,j) = F(X(i,j), Y(i,j));
        
    end
end

surf(X, Y, Z)
xlabel('x values');
ylabel('y values');
zlabel('z values');
title('Example elliptic paraboloid');

%line of intersection of the planes
%x + 2 y + z - 1 = 0 
%2 x + 3 y - 2 z + 2 = 0

clear all % clear all variables
close all % close all figures

[X, Y] = meshgrid(-20:0.5:20, -20:0.5:20);

N = size(X);

G = zeros(N(1), N(2));
H = zeros(N(1), N(2));

g = @(x, y) 1 - x - 2*y;
h = @(x, y) x + (3/2)*y + 1;

for i = 1:N(1)
    for j = 1:N(2)
        
        G(i,j) = g(X(i,j), Y(i,j));
        H(i,j) = h(X(i,j), Y(i,j));
        
    end
end

hold on
surf(X, Y, G)
surf(X, Y, H)
hold off
%Note: You may need to use the Rotate 3D tool in the Figure window to see
%the three surfaces at a good perspective
%%

%Part 5
%Plotting a list of points

clear all % clear all variables
close all % close all figures

x = [1:10];
y = [22 37 65 100 111 181 190 150 185 180];

plot(x, y, 'r.') %'.' plots the points and 'r' plots the points as red
%%

%Part 6
%Contour plot

%%Creating a contour plot of an elliptic paraboloid
%z = 4 x^2 + y^2

clear all % clear all variables
close all % close all figures

[X, Y] = meshgrid(-20:0.5:20);

N = size(X);

Z = zeros(N(1), N(2));

F = @(x, y) 4*(x^2) + (y^2);

for i = 1:N(1)
    for j = 1:N(2)
        
        Z(i,j) = F(X(i,j), Y(i,j));
        
    end
end

[C, h] = contour(X, Y, Z, 10); %The '10' chooses 10 contour levels
clabel(C,h) %Labels the values of z choosen for the contour levels
xlabel('x values');
ylabel('y values');
title('Contour plot of example elliptic paraboloid');

v = [1; 10; 100; 1000];
[C, h] = contour(X, Y, Z, v); %Vector v contains the values of z selected
                              %for the contour plot
clabel(C,h)
xlabel('x values');
ylabel('y values');
title('Contour plot of example elliptic paraboloid');
%%

%Part 7
%Basic matrix operations

%Create 2x2 matrix A and B

clear all % clear all variables
close all % close all figures

A = [1 2; -3 4];
B = [10 2; 0 -8];

C = A*B;

%Solve a matrix equation of the form A x = b, where x and b are vectors

clear all % clear all variables
close all % close all figures

A = [1 1 1; 1 2 5; 2 5 -1];
b = [6; -4; 27];

x = mldivide(A, b)

%Create 3x3 identity matrix, I

I = eye(3)
%%

%Part 8
%Eigenvalues and eigenvectors

clear all % clear all variables
close all % close all figures

%Consider 2x2 matrix A
A = [-4 -2; -2 -2];

%Function 'eig' returns the eigenvalues and eigenvectors of a matrix
[V, D] = eig(A);

%Eigenvalues are returned on the diagonal entries of the matrix D
D

%The matrix V contains the right eigenvectors on the
%columns such that A * V = V * D
V

%Also, can obtain symbolic eigenvalues and eigenvectors if use a symbolic 
%matrix

clear all % clear all variables
close all % close all figures

%Create a 2x2 symbolic matrix J
A = sym([-4 -2; -2 -2]);

[V, D] = eig(A);

%Eigenvalues are returned on the diagonal entries of the matrix D
D

%The matrix V contains the right eigenvectors on the
%columns such that A * V = V * D
V
%%

%Part 9
%Solving nonlinear equations

clear all % clear all variables
close all % close all figures

syms x y z;

[solx, soly, solz] = solve([(x^2)+y+z == 1, x/y == 2, z == 3], x, y, z);

solx % (31^(1/2)*i)/4 - 1/4
     % -(31^(1/2)*i)/4 - 1/4
     
soly % (31^(1/2)*i)/8 - 1/8
     % -(31^(1/2)*i)/8 - 1/8
     
solz % 3
     % 3
%%

%Part 10
%Solving ODEs

%Numerically solve the following system of ODEs:
% x' = -x + x^3
% y' = -2y

clear all % clear all variables
close all % close all figures

f = @(t,x) [-1*x(1) + (x(1))^3; -2*x(2)];

[t, xa] = ode45(f, [0:1:10], [0.5, 0.5]);

%The point (0,0) is a stable node, and (-1, 0) and (1, 0) are saddles
%shown by the linearization of the system of ODEs.
%Since the initial condition was (0.5, 0.5), the solutions of x and y in
%time will tend to the stable node (0,0)

%Solution of x
solx = xa(:,1);
solx

%Solution of y
soly = xa(:,2);
soly
%%

%Part 11
%Plotting solutions of ODEs in time plots

clear all % clear all variables
close all % close all figures

% x' = -x + x^3
% y' = -2y

f = @(t,x) [-1*x(1) + (x(1))^3; -2*x(2)];

[t, xa] = ode45(f, [0:1:10], [0.5, 0.5]);

solx = xa(:,1);
soly = xa(:,2);

%Since the initial condition was (0.5, 0.5), the solutions of x and y in
%time will tend to the stable node (0,0).

%Plotting the solutions x and y with initial condition (0.5, 0.5)
hold on
plot(t, solx, 'b')
plot(t, soly, 'r')
hold off

xlabel('time');
ylabel('values');
title('ODE solutions for x and y over time with initial condition (0.5, 0.5)');
legend('x','y')
%%

%Part 12
%2D Phase Plots

clear all % clear all variables
close all % close all figures

%Consider Lotka-Volterra model of competition
% x' = x(3 - x - 2y)
% y' = y(2 - x - y)

f = @(t,x) [x(1)*(3 - x(1) - 2*x(2)); x(2)*(2 - x(1) - x(2))];

%Creating the vector field
[X, Y] = meshgrid(0:0.2:5);

u = zeros(size(X));
v = zeros(size(Y));

t = 0;

%Calculate the value of the derivatives at each point in the grid
for i = 1:length(X)
    
    for j = 1:length(Y)

        D = f(t, [X(i, j), Y(i, j)]);
        
        %Normalize vectors for plotting
        L = sqrt((D(1))^2 + (D(2))^2);
        
        u(i, j) = D(1)/L;
        v(i, j) = D(2)/L;
        
    end

end

%Plotting the solutions x and y with initial condition (0.5, 0.5)
[t, x1] = ode45(f, [0:1:10], [0.5, 0.5]);
solx1 = x1(:,1);
soly1 = x1(:,2);

%Plotting the solutions x and y with initial condition (0.5, 1)
[t, x2] = ode45(f, [0:1:10], [0.5, 1]);
solx2 = x2(:,1);
soly2 = x2(:,2);

%Plotting the solutions x and y with initial condition (3, 1)
[t, x3] = ode45(f, [0:1:10], [3, 1]);
solx3 = x3(:,1);
soly3 = x3(:,2);

%Plotting the solutions x and y with initial condition (2, 2)
[t, x4] = ode45(f, [0:1:10], [2, 2]);
solx4 = x4(:,1);
soly4 = x4(:,2);

%Plotting the solutions x and y with initial condition (0, 1)
[t, x5] = ode45(f, [0:1:10], [0, 1]);
solx5 = x5(:,1);
soly5 = x5(:,2);

%Plotting the solutions x and y with initial condition (1, 0)
[t, x6] = ode45(f, [0:1:10], [1, 0]);
solx6 = x6(:,1);
soly6 = x6(:,2);

hold on

quiver(X, Y, u, v, 0.5, 'r')
xlabel('x')
ylabel('y')
axis tight equal %This command keeps the axis strictly between the X and Y
                 %values

%Equilibrium points
a = [0; 0; 3; 1];
b = [0; 2; 0; 1];
scatter(a, b)

%Plotting the solutions x and y with four initial conditions: (0.5, 0.5), 
%(0.5, 1), (3, 1), (2, 2), (0, 1), (1, 0), (1, 1)

plot(solx1, soly1) %initial condition (0.5, 0.5)
plot(solx2, soly2) %initial condition (0.5, 1)
plot(solx3, soly3) %initial condition (3, 1)
plot(solx4, soly4) %initial condition (2, 2)
plot(solx5, soly5) %initial condition (0, 1)
plot(solx6, soly6) %initial condition (1, 0)

hold off

xlabel('x');
ylabel('y');
title('ODE solutions for x and y for Lotka-Volterra model with different initial conditions');
%%

%Part 13
%Solving PDE heat equation and plotting solution as a surface

%Consider PDE heat equation
%u_t = u_xx
%u(x, 0)= (2*x)/(1 + (x)^2), 0 <= x <= 1
%u(0, t) = 0, t > 0
%u(1, t) = 1, t > 0

%Use Matlab function 'pdepe'

clear all % clear all variables
close all % close all figures

m = 0;
x = linspace(0,1,100);
t = linspace(0,2,10);

sol = pdepe(m,@pde,@ic,@bc,x,t);
% Extract the first solution component as u.
u = sol(:,:,1);

%Surface plot
surf(x,t,u) 
title('Numerical solution computed with 100 mesh points in space and 10 mesh points in time')
xlabel('Distance x')
ylabel('Time t')

%Solution profile at t = 2
plot(x,u(10,:))
title('Solution at t = 2')
xlabel('Distance x')
ylabel('u(x,2)')

%Create movie of the evolution of the profile over time
h = animatedline(x,u(1,:));
for k = 2:length(t)
    clearpoints(h)
    addpoints(h,x,u(k,:))
    pause(2)
    drawnow
end
%%

%Part 14
%Solving PDE advection equation and plotting solutions

%Consider PDE advection equation
%u_t + u_x = 0, x > 0, t > 0
%u(0,x) = 0, x > 0
%u(t, 0) = 1
%u(t, 10) = 0

%Using the Lax-Wendroff scheme
%u_{i, n+1} = u_{i, n} - (C/2) (u_{i+1, n} - u_{i-1, n})
%             + ((C^2) / 2) (u_{i+1, n} - 2 u_{i, n} + u_{i-1, n})
%where C = k/h

clear all % clear all variables
close all % close all figures

%Choose a step size in space of h = 1/10 and a step size in time of k =
%1/20. 
%Then the Courant number C = (1/20)/(1/10) = 0.5 < 1.

%The solution u1 of the initial-boundary value problem with C = 0.5 < 1
%and with 1 <= n <= 100 and 1 <= j <= 99
u1 = u_solve(1/20, 1/10, 100, 99);

%The corresponding values of x_j for 1 <= j <= 99 with h = 1/10
x1 = xValues(99, 1/10);

%Now, choose a step size in space of h = 1/20 and a step size in time of 
%k = 1/20. %The solution u2 of the initial-boundary value problem with 
%C = 1 and with 1 <= n <= 100 and 1 <= j <= 199
u2 = u_solve(1/20, 1/20, 100, 199);

%The corresponding values of x_j for 1 <= j <= 99 with h = 1/20
x2 = xValues(199, 1/20);

%Lastly, choose a step size in space of h = 1/10 and a step size in time of 
%k = 1/5. %The solution u3 of the initial-boundary value problem with 
%C = 2 and with 1 <= n <= 25 and 1 <= j <= 99
u3 = u_solve(1/5, 1/10, 25, 99);

%The corresponding values of x_j for 1 <= j <= 99 with h = 1/10
x3 = xValues(99, 1/10);

%To plot solution u1
hold on

plot(x1, u1(100,:));
axis([0 10 0 1.5])
title('Approximate solution of PDE at t = 5 with C = 0.5 < 1')
xlabel('x')
ylabel('u_h(5, x)')

hold off

%To plot solution u2
hold on

plot(x2, u2(100,:));
axis([0 10 0 1.5])
title('Approximate solution of PDE at t = 5 with C = 1')
xlabel('x')
ylabel('u_h(5, x)')

hold off

%To plot solution u3
hold on

plot(x3, u3(25,:));
title('Approximate solution of PDE at t = 5 with C = 2 > 1')
xlabel('x')
ylabel('u_h(5, x)')

hold off

%Create movie of the evolution of the profile over time
h = animatedline(x2,u2(1,:));
for k = 2:100
    clearpoints(h)
    addpoints(h,x2,u2(k,:))
    drawnow
end