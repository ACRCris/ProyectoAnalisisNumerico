clc

global g v0 pi_val a gamma z_crit m0

% Initialize the constants
g = 981;  % gravitational acceleration in cm/s^2
a = 0.1;  % radius
gamma = 0.01;  % damping coefficient
m0 = 2.7;  % initial mass
z0 = 0;  % initial height
z_crit = 5e-1;  % critical height

% Define the range of V0 values
V0_values = 1:0.0001:2;  % Adjust the range and step as needed

% Initialize arrays to store V0 and T_n
V0_list = [];
Tn_list = [];

for V0 = V0_values
    % Set v0 (global variable)
    v0 = V0;

    disp(V0)
    
    % Initialize initial_conditions
    initial_conditions = [z0; 0; m0];
    
    % Initialize periods array
    periods = [];
    
    N_per = 100;  % Total number of periods
    N_start = 30; % Start collecting data after this period to avoid transients
    
    for i = 1:N_per
        % Time span for the solution
        t_span = [0, 1000];
        
        % Solve the ODE system using ode45 with event handling
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Events', @event_breakoff);
        [t_solution, y_solution] = ode45(@system, t_span, initial_conditions, options);
        
        % Extract the solution components
        z = y_solution(:,1);
        v = y_solution(:,2);
        m = y_solution(:,3);
    
        periods(i) = round(t_solution(end),8);
        
        initial_conditions = [0; 0; m(end) * 0.6];
    end
    
    % Collect T_n after transients
    for n = N_start:N_per
        V0_list = [V0_list; V0];
        Tn_list = [Tn_list; periods(n)];
    end
end

% Plot the bifurcation diagram
figure;
scatter(V0_list, Tn_list, 10, '.');
xlabel('V_0');
ylabel('Period T_n');
title('Bifurcation Diagram: T_n vs V_0');
grid on;

% Define k(m) as a function
function k_val = k_m(m)
    if m < 4.61
        k_val = -11.4 * m + 52.5;
    else
        k_val = 0;
    end
end

% Define the ODE system
function dydt = system(t, y)
    global g v0 pi_val a gamma z_crit

    z = y(1);
    dzdt = y(2);
    m = y(3);
    k = k_m(m);
    flow_rate = pi * a^2 * v0;
    dmdt = flow_rate;
    dz2dt2 = (m * g - k * z - gamma * dzdt - flow_rate * (dzdt - v0)) / m;
    dydt = [dzdt; dz2dt2; dmdt];
end

% Event function to stop integration when z reaches z_crit
function [value, isterminal, direction] = event_breakoff(t, y)
    global z_crit
    z = y(1);
    value = z - z_crit;  % Event triggers when z = z_crit
    isterminal = 1;  % Stop the integration
    direction = 0;   % Trigger in both directions
end