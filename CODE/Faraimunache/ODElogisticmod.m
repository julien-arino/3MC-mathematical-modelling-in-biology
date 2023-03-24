%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code solves the Logistic growth model
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ODElogisticmod
close all;
clc
%%%%%%%%%%%%%%%%%% Request parameter values from user %%%%%%%%%%%%%%%%%%%%%
  prompt={'Enter the value of r ',...
          'Enter the value of K ',...
          };
name='Parameter Values for Logistic Model';
numlines=1;%specifies the number of lines for each answer in numlines
   
defaultanswer={'0.02','1000'};
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';
% create a dialog box that returns user input for multiple prompts
%in a cell array
% inputdlg suspends execution until the user responds.
parameters=inputdlg(prompt,name,numlines,defaultanswer,options);
%%%%%%%%%%%%%%%%%%% Request Initial data from user %%%%%%%%%%%%%%%%%
 prompt1={'Enter the starting time, t_1',...
          'Enter the final time t_f',...
          'Enter the initial Population N_1','Enter the initial Population N_2','Enter the initial Population N_3'};
name1='Initial data';
numlines=2;
defaultanswer1={'0','500','100','600','1500'};

options1.Resize='on';
options1.WindowStyle='normal';
options1.Interpreter='tex';
initialDat=inputdlg(prompt1,name1,numlines,defaultanswer1,options1);      
%%%%%%%%%%%%%%%%%%%%%% Convert cell arrays to numbers %%%%%%%%%%%%%%%%%%%%%
r     = str2double(parameters{1}); 
K     = str2double(parameters{2});

t_1    = str2double(initialDat{1});
t_f    = str2double(initialDat{2});
N_0    = str2double(initialDat(3:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = [t_1 t_f];
for i=1:length(N_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,x] = ode45(@Logistic,tspan,N_0(i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold on
plot(t,x,'linewidth',2);
xlabel('time');
ylabel('N');
legend
%legend('Population N_1','Population N_2','Population N_3');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xdot = Logistic(~,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N= x(1);


dN = r*N*(1-(N/K));
  
xdot = dN;
end
end
