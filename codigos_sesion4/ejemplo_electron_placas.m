q = -1.6e-19;    % Carga del electrón en COULOMBS
m = 9.11e-31;    % Masa del electrón KILOGRAMOS

% Placas
d = 0.05;        % Separación entre placas (m)
v0 = 10;        % Velocidad inicial horizontal (m/s)

% Tiempo de simulación
dt = 1e-9;
t_max = 2e-7;
t = 0:dt:t_max;

% Posición inicial
x = zeros(size(t));
y = zeros(size(t));
x(1) = 0;
y(1) = d/2;

% Campo eléctrico no uniforme
E0 = 10;  % Magnitd del campo (V/m)

% Inicializar velocidad
vx = v0;
vy = 0;

% Simular movimiento
for i = 2:length(t)
    % Campo eléctrico
    E = E0 * (1 - y(i-1)/d); 
    
    % Fuerza eléctrica
    Fy = q * E;
    
   
    ay = Fy / m;
    
    % Actualizar posiciones y velocides
    vy = vy + ay * dt;
    x(i) = x(i-1) + vx * dt;
    y(i) = y(i-1) + vy * dt;
    
    % Condición de borde. Detenerse si choca con las placas
    if y(i) <= 0
        y(i) = 0;
        vy = 0;
        vx = 0;
    elseif y(i) >= d
        y(i) = d;
        vy = 0;
        vx = 0;
    end
end

% Animación
figure;
for i = 1:10:length(t)
    clf;
    
    % placas
    fill([0 x(end) x(end) 0], [0 0 0.001 0.001], [0.7 0.7 0.7]); % placa inferior
    hold on;
    fill([0 x(end) x(end) 0], [d-0.001 d-0.001 d d], [0.7 0.7 0.7]); % placa superior
    
    % Partícula
    plot(x(i), y(i), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    
    xlim([0 x(end)]);
    ylim([-0.002 d+0.002]);
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('Tiempo = %.1f ns', t(i)*1e9));
    grid on;
    drawnow;  %Actualizar en tiempo real la posición
end
