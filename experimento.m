x = linspace(-10, 10, 1000);
y = linspace(-10, 10, 1000);

% Crear la malla de puntos (grid)
[X, Y] = meshgrid(x, y);

% Definir la función z
f = @(x, y) sin(x) + sin(y);

% Calcular Z para cada punto de la malla
Z = f(X, Y);

% Graficar la superficie 3D
figure;
mesh(X, Y, Z);  % Usar surf para superficies
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Gráfico de la función f(x, y)');
