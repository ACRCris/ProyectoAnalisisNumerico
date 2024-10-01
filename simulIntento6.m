clc

global g v0 pi_val a gamma z_crit m0,

% Initialize the constants
g = 981;  % gravitational acceleration in cm/s^2
v0 =0.04;  % constant velocity
a = 0.25;  % 7radiu0
gamma =0.001;  % damping coefficient
m0 = 0.5;  % initial mass
z0 = 0;  % initial height
z_crit = 5e-1;  % critical height where dro

initial_conditions = [z0; 0; m0];
periods = [];
for i=1:500
    % Initial conditions: z(0), dz/dt(0), and m(0)
    
    % Time span for the solution
    t_span = [0, 1000];
    
    % Solve the ODE system using ode45 with event handling
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Events', @event_breakoff);
    [t_solution, y_solution] = ode45(@system, t_span, initial_conditions, options);
    
    % Extract the solution components
    z=y_solution(:,1);
    v=y_solution(:,2);
    m=y_solution(:,3);

    %disp("t:")
    %disp(t_solution(length(t_solution)));
    periods(i) = round(t_solution(length(t_solution)),8);
    % disp("z")
    %disp(z(length(z)));
    % disp(v(length(v)));
    % disp("masa");
    %disp(m(length(m)));
    % 
    initial_conditions = [0,0,m(length(m))*.6];
    % disp("init")
    % disp(initial_conditions)

    % Extract the solution components
    z_solution = y_solution(:, 1);
    dzdt_solution = y_solution(:, 2);
    m_solution = y_solution(:, 3);
    % 
    % % Plot the results
    % figure;
    % subplot(3, 1, 1);
    % plot(t_solution, z_solution, 'DisplayName', 'Height (z)');
    % xlabel('Time [s]');
    % ylabel('z(t) [m]');
    % title('Height over Time');
    % yline(z_crit, '--r', 'Critical height z_{crit}');
    % legend();
    % grid on;

    % subplot(3, 1, 2);
    % plot(t_solution, dzdt_solution, 'DisplayName', 'Velocity (dz/dt)', 'Color', 'r');
    % xlabel('Time [s]');
    % ylabel('Velocity [m/s]');
    % title('Velocity over Time');
    % yline(0, '--r', 'Zero velocity at z_{crit}');
    % legend();
    % grid on;
    % 
    % subplot(3, 1, 3);
    % plot(t_solution, m_solution, 'DisplayName', 'Mass (m)', 'Color', 'g');
    % xlabel('Time [s]');
    % ylabel('Mass [kg]');
    % title('Mass over Time');
    % grid on;



end

% Crear la carpeta 'imagenes' si no existe
output_folder = 'imagenesDefinitivas';
% output_folder = 'interesantes';

if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Construir una cadena base para los nombres de archivo
% Reemplazar '.' por 'p' para evitar conflictos en nombres de archivo
filename_base = sprintf('g_%d_v0_%.8f_a_%.3f_gamma_%.8f_m0_%.2f_zcrit_%.1f', ...
    g, v0, a, gamma, m0, z_crit);
filename_base = strrep(filename_base, '.', 'p');  % Reemplazar '.' por 'p'

titulo = {
    sprintf('g = %.0f, v_0 = %.4f, a = %.2f', g, v0, a), ...
    sprintf('gamma = %.2f, m_0 = %.2f, z_{crit} = %.1f', gamma, m0, z_crit)
};
disp(titulo)
figure1 =figure;
plot(periods(41:end),'.');
xlabel('Drop Number');
ylabel('Period (s)');
title(titulo);
grid on;


% Guardar la primera gráfica como PDF
filename1 = fullfile(output_folder, [filename_base, '_PeriodosPorDrop.pdf']);
saveas(figure1, filename1);

disp(periods)


figure2 = figure;
plot(periods(41:end-1), periods(42:end), 'o', 'MarkerSize', 8);
xlabel('T_n (Period n)');
ylabel('T_{n+1} (Period n+1)');
title(titulo);
grid on;
hold on;

% Determinar los límites para la línea y = x basada en los datos
x_min = min([periods(10:end-1); periods(11:end)]);
x_max = max([periods(10:end-1); periods(11:end)]);
plot([x_min, x_max], [x_min, x_max], 'r--', 'LineWidth', 1.5, 'DisplayName', 'y = x');

% Guardar la segunda gráfica como PDF
filename2 = fullfile(output_folder, [filename_base, '_DiagramaPoincare.pdf']);
saveas(figure2, filename2);



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