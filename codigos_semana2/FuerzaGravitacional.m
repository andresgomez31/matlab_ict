% Gráfica del "campo" gravitacional (ojo, en realidad es un campo de
% fuerzas, no un CAMPO como formalmente se define en Física).
% Por: Víctor M. Rico, 2025-Ene-26

N = 5;
fig = figure(1);
fig.Color = 'w';

% Definición de la malla
[x,y,z] = meshgrid(linspace(-1,1,N));
r = sqrt(x.^2 + y.^2 + z.^2);

% Constantes
G = 6.67e-11;   % N*m^2/kg^2
M = 5.972e24;   % kg, masa de la tierra
m = 60;         % kg, masa promedio de una persona adulta

F = -G * M*m/r.^2;
Fx = F.*x./r;
Fy = F.*y./r;
Fz = F.*z./r;

quiver3(x,y,z,Fx,Fy,Fz,'-b')
axis equal
xlabel('x (radios del planeta)')
ylabel('y (radios del planeta)')
zlabel('z (radios del planeta)')