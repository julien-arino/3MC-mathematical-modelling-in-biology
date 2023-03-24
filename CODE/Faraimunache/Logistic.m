%% Get user input for logistic growth parameters and initial population values
prompt = {'Enter the value of r:', 'enter the value of K:', 'Initial time t_1:', 'Final time t_f:', 'Initial population N_1:', 'Initial population N_2:', 'Initial population N_3:', 'Initial population N_4:','Initial population N_5:'};
name = 'Logistic Growth Parameters and Initial Population Values';
tspan = [1 50];
defaultinput = {'0.02', '1000', '0', '500', '100', '600', '900','1100','1600'};
parameters = inputdlg(prompt, name, tspan, defaultinput);

%% Convert input to numeric values
r = str2double(parameters{1});
K = str2double(parameters{2});
t_1 = str2double(parameters{3});
t_f = str2double(parameters{4});
N0 = str2double(parameters(5:end));% using multiple initial conditions.

%% Set time range
t = t_1:0.1:t_f;

% Solve differential equation for each initial population value
for i = 1:length(N0)
    [t,N] = ode45(@(t,N) r*N*(1-N/K), t, N0(i));
    % Plot results for each initial population value on the same graph
    hold on
    plot(t,N,'linewidth',2)
end

% Set plot properties
%title('Logistic Growth')
xlabel('Time')
ylabel('Population Size, N(t)')
box on
legend('Population 1', 'Population 2', 'Population 3','Population 4', 'Population 5')
