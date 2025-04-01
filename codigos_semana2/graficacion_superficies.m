%Dibujar superficie de la forma z = f(x,y)
%{
%Generamos una malla
[X,Y] = meshgrid(1:0.2:20,0.2:20);
Z = sin(X) + cos(Y);
C = X.*Y;
%mesh,meshc, surf, surfc, plot3
mesh(X,Y,Z) %Grafica la función

%El parámetro C, determina qué valores corresponden al código de color, se
%toma Z por defecto.

%surf(X,Y,Z,C) %Rellena la malla con colores, surfc()
%colorbar

contour(X,Y,Z,3) %la última entrada controla No. de lineas
pcolor(X,Y,Z)
colorbar
shading interp  %Quitar la cuadrícula
%}

% Superficie paramétrica sigma(u,v)

%Ejemplo de Toroide
u = linspace(0,2*pi,50);
v = u;
[U,V] = meshgrid(u,v);

X = cos(U).*(2 + cos(V));
Y = sin(U).*(2 + cos(V));
Z = sin(V);
surf(X,Y,Z);
%shading interp; %Elimina las lineas y suaviza
axis equal;
