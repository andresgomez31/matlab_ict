%%%%%%%%%%%%%%%%%%%%%%
%Paso3_Sesion4_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%POR HACER%

% Cálculo de la Fuerza Eléctrica%

% (1) Dentro del bucle for interno, el que recorre las cargas en las placas (electrodos),
% calcula la fuerza total sobre el eritrocito.

% (2) Dentro del bucle while, después de abrir y cerrar el bucle for que usaste para calcular la 
% fuerza total sobre el eritrocito (y por encima de la parte de código delete(head) del paso anterior), 
% asume una masa de 1, y calcula (usando dinámica de aceleración constante) la aceleración (a), 
% velocidad en x (Vx), posición en x (xe) y posición en y (ye). Usa aquí el paso de tiempo "dt" predefinido.

clear;
clc;
clf;

%-----------------?------------------%
Lp=3.5;                      % 
Ln=2.5;                      % 
t=0.02;                      % 
d=0.4;                       % 
p=0.01;                      % 

% Definir las características de los electrodos
ke=1/(4*pi*8.85*10^-12);     % 
Q=1e-3;                      % ¿Qué es esto?
Nq=28;                       % ¿Y esto?

%-----------------?------------------%

xmin_original=-d/2-3*t;  xmax_original=-xmin_original;  % 
xmin=-d/2-3*t;  xmax=-xmin;                             % 
ymin=2*(-Lp/2);   ymax=-ymin;                           % 

%-----------------?------------------%
if ymin <= -1  
    if xmin >= -0.5 && xmax <= 0.5
            xmin = -1.5;
            xmax = -xmin;
    end
end

%-----------------?------------------%
Ny=30;  Nx=Ny;
x=linspace(xmin, xmax, Nx); y=linspace(ymin, ymax, Ny);

%-----------------?------------------%

vertices2d=[[-d/2-t,Lp
