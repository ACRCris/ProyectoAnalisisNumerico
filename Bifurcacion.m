function chaotic_dripping_faucet()
    % Parámetros físicos
    g = 9.81;            % aceleración gravitatoria (m/s^2)
    gamma = 0.1;         % coeficiente de amortiguación (Ns/m)
    a = 0.005;           % diámetro de la boquilla (m)
    
    % Función por partes para la constante de resorte k(m)
    k_func = @(m) (-11.4 * m + 52.5) .* (m < 4.61);  % k = -11.4m + 52.5 para m < 4.61, k = 0 en otro caso
    
    % Parámetros de tiempo
    tspan = [0 10000];      % intervalo de tiempo para la simulación (segundos)

    % Condiciones iniciales
    m0 = 0.1;            % masa inicial de la gota (kg)
    v0_droplet = 0;      % velocidad inicial de la gota (m/s)
    z0 = 0;              % altura inicial de la gota (m)

    % Tasas de flujo para las diferentes simulaciones
    flow_rates = [0.319, 0.35, 0.374, 0.456, 0.579]; % en mL/s
    
    %% FIGURA 4: Periodo-1 con tasa de flujo 0.319 mL/s
    v0 = flow_rates(1) / (pi * a^2 * 1000);  % convertir mL/s a m/s
    flow_rate = pi * a^2 * v0;  % tasa de flujo en kg/s

    % Simular para la tasa de flujo 0.319 mL/s (Periodo-1)
    simulate_and_plot(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, 'Figura 4: Periodo-1');

    %% FIGURA 5: Periodo-2 con tasa de flujo 0.35 mL/s
    v0 = flow_rates(2) / (pi * a^2 * 1000);  % convertir mL/s a m/s
    flow_rate = pi * a^2 * v0;  % tasa de flujo en kg/s

    % Simular para la tasa de flujo 0.35 mL/s (Periodo-2)
    simulate_and_plot(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, 'Figura 5: Periodo-2');

    %% FIGURA 6: Atractores caóticos con diferentes tasas de flujo
    figure('Name', 'Figura 6: Atractores Caóticos');
    
    % (a) Atractor extraño en tasa de flujo de 0.374 mL/s
    v0 = flow_rates(3) / (pi * a^2 * 1000);  % convertir mL/s a m/s
    flow_rate = pi * a^2 * v0;  % tasa de flujo en kg/s
    subplot(3,1,1);
    simulate_attractor(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, 'Atractor extraño a 0.374 mL/s');

    % (b) Atractor extraño en tasa de flujo de 0.456 mL/s
    v0 = flow_rates(4) / (pi * a^2 * 1000);  % convertir mL/s a m/s
    flow_rate = pi * a^2 * v0;  % tasa de flujo en kg/s
    subplot(3,1,2);
    simulate_attractor(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, 'Atractor extraño a 0.456 mL/s');

    % (c) Atractor extraño en tasa de flujo de 0.579 mL/s
    v0 = flow_rates(5) / (pi * a^2 * 1000);  % convertir mL/s a m/s
    flow_rate = pi * a^2 * v0;  % tasa de flujo en kg/s
    subplot(3,1,3);
    simulate_attractor(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, 'Atractor extraño a 0.579 mL/s');

    %% FIGURA 7: Periodo-3 con tasa de flujo 0.374 mL/s
    v0 = flow_rates(3) / (pi * a^2 * 1000);  % convertir mL/s a m/s
    flow_rate = pi * a^2 * v0;  % tasa de flujo en kg/s

    % Simular para la tasa de flujo 0.374 mL/s (Periodo-3)
    simulate_and_plot(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, 'Figura 7: Periodo-3');

end

function simulate_and_plot(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, fig_title)

    % Simulación del sistema para una tasa de flujo específica

    % Function to compute equations of motion for the droplet
    equations_of_motion = @(t, state) [
        state(2);                                                    % dz/dt = velocity
        (g - k_func(state(3)) * state(1)/state(3) - gamma * state(2)/state(3)) - flow_rate * (state(2) - v0) / state(3); % d^2z/dt^2 = acceleration with mass change
        flow_rate                                                   % dm/dt = flow_rate
    ];
    % Vector de estado inicial: [altura, velocidad, masa]
    initial_state = [z0, v0_droplet, m0];
    % Resolver el sistema de ecuaciones diferenciales con ODE45
    [t, states] = ode45(equations_of_motion, tspan, initial_state);
    % Extraer los datos de altura (z)
    z = states(:,1);
    v_droplet = states(:,2);   % droplet velocity over time
    m = states(:,3);   % droplet mass over time

    % Detectar los tiempos de caída (cuando la altura z cruza 0.02 desde positivo hacia negativo)

    drop_times = [];
    for i = 2:length(z)
        if z(i-1) > 0.02 && z(i) <= 0.02  % Umbral de ruptura de la gota
            drop_times = [drop_times, t(i)];
        end
    end

    % Calcular los periodos Tn como la diferencia entre tiempos de caída consecutivos

    Tn = diff(drop_times);
    % Si no hay suficientes tiempos de goteo, no se puede construir el mapa de Poincaré
    if isempty(Tn) || length(Tn) < 2
        error('No se han detectado suficientes tiempos de caída de gotas. Ajuste los parámetros.');
    end

    % Calcular los datos para el mapa de Poincaré (Tn+1 vs Tn)

    Tn_next = Tn(2:end);   % Tn+1
    Tn_current = Tn(1:end-1);  % Tn

    % Crear figura
    figure('Name', fig_title);
    
    % (a) Mapa de Poincaré
    subplot(3,2,1);
    plot(Tn_current, Tn_next, 'b.');
    xlabel('Tn');
    ylabel('Tn+1');
    title('Mapa de Poincaré: ');
    grid on;

    % (b) Serie temporal de los periodos
    subplot(3,2,2);
    plot(1:length(Tn), Tn, 'b.-');
    xlabel('Número de Gota');
    ylabel('Periodo Tn (s)');
    title('Periodo vs. Número de Gota: ');
    grid on;

    % (c) Altura de la gota a lo largo del tiempo
    subplot(3,2,3);
    plot(t, z);
    xlabel('Time (s)');
    ylabel('Height (m)');
    title('Height of Droplet Over Time ');
    grid on;
    
    % (d) Velocidad de la gota a lo largo del tiempo
    subplot(3,2,4);
    plot(t, v_droplet);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Velocity of Droplet Over Time ');
    grid on;
    
    % (e) Masa de la gota a lo largo del tiempo
    ax_mass = subplot(3,1,3); % Cambiado a la última fila
    plot(t, m);
    xlabel('Time (s)');
    ylabel('Mass (kg)');
    title('Mass of Droplet Over Time ');
    grid on;

    % Crear un nuevo eje para superponer el mapa de Poincaré
    ax_poincare = axes('Position', ax_mass.Position); % Usar la misma posición que el gráfico de masa
    hold(ax_poincare, 'on'); % Mantener el gráfico
    set(ax_poincare, 'Color', 'none'); % Hacer el fondo transparente

    % Graficar el mapa de Poincaré en el nuevo eje
    plot(ax_poincare, Tn_current, Tn_next, 'r.'); % Graficar el mapa de Poincaré
    ylabel('Tn+1 (s)');
    xlabel('Tn (s)');
    title('Poincaré Map on Mass Graph');
    hold(ax_poincare, 'off'); % Liberar el gráfico

    % Ajustar el espacio vertical entre subgráficas
    set(gcf, 'Position', [100, 100, 800, 800]); % Aumentar el tamaño de la figura
    sgtitle(fig_title); % Título general para la figura

    % Ajustar etiquetas para evitar cruces
    ax_mass.XLabel.Position(2) = ax_mass.XLabel.Position(2) - 0.5; % Mover etiqueta x hacia abajo
    ax_mass.YLabel.Position(1) = ax_mass.YLabel.Position(1) + 0.5; % Mover etiqueta y hacia arriba
    ax_poincare.XLabel.Position(2) = ax_poincare.XLabel.Position(2) - 0.5; % Mover etiqueta x hacia abajo
    ax_poincare.YLabel.Position(1) = ax_poincare.YLabel.Position(1) + 0.5; % Mover etiqueta y hacia arriba
end

function simulate_attractor(tspan, m0, v0_droplet, z0, g, gamma, flow_rate, k_func, v0, attractor_title)
    % Simulación y visualización del atractor caótico
    % Function to compute equations of motion for the droplet
    equations_of_motion = @(t, state) [
        state(2);                                                    % dz/dt = velocity
        (g - k_func(state(3)) * state(1)/state(3) - gamma * state(2)/state(3)) - flow_rate * (state(2) - v0) / state(3); % d^2z/dt^2 = acceleration with mass change
        flow_rate                                                   % dm/dt = flow_rate
    ];
    % Vector de estado inicial: [altura, velocidad, masa]
    initial_state = [z0, v0_droplet, m0];

    % Resolver el sistema de ecuaciones diferenciales con ODE45
    [t, states] = ode45(equations_of_motion, tspan, initial_state);

    % Extraer los datos de altura (z)
    z = states(:,1);

    % Detectar los tiempos de caída (cuando la altura z cruza 0.02 desde positivo hacia negativo)

    drop_times = [];

    for i = 2:length(z)
        if z(i-1) > 0.02 && z(i) <= 0.02  % Umbral de ruptura de la gota
            drop_times = [drop_times, t(i)];
        end
    end
    % Calcular los periodos Tn como la diferencia entre tiempos de caída consecutivos
    Tn = diff(drop_times);

    % Si no hay suficientes tiempos de goteo, no se puede construir el mapa de Poincaré
    if isempty(Tn) || length(Tn) < 2
        return;
    end

    % Calcular los datos para el mapa de Poincaré (Tn+1 vs Tn)
    Tn_next = Tn(2:end);   % Tn+1
    Tn_current = Tn(1:end-1);  % Tn

    % Graficar el atractor
    plot(Tn_current, Tn_next, 'b.');
    xlabel('Tn');
    ylabel('Tn+1');
    title(attractor_title);
end
