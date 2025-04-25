%%%%%%%%%%%%%%%%%%%%%%
%Step0_Session3_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%{
Este código es muy similar al que realizaste la sesión de reto pasada, con
algunas lineas adicionales y modificaciones que resuleven el problema de la
visualización del campo eléctrico total independientemente de la forma en que lo
definamos; con el principio de superposición o usando el gradiente de el
potencial eléctrico.
%}

%POR HACER%
%{
Investiga, corre y comenta cada una de las lineas, averigua qué hacen los
bloques de código señalados con un  -----------------?------------------.
%}



clear;
clc;
clf;

%-----------------?------------------%
Lp=3.5;                      % 
Ln=2.5;                      % 
t=0.02;                      % 
d=0.4;                       % 
p=0.01;                      % 

% Define las características de los electrodos
ke=1/(4*pi*8.85*10^-12);     % 
Q=1e-3;                      % ¿Qué hace esto?
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

vertices2d=[[-d/2-t,Lp/2]    %1
    [-d/2,Lp/2]              %2
    [-d/2,-Lp/2]             %3  
    [-d/2-t,-Lp/2]           %4
    [d/2,Ln/2]               %5
    [d/2+t,Ln/2]             %6  
    [d/2+t,-Ln/2]            %7  
    [d/2,-Ln/2]];            %8  

%  
facesP=[1 2 3 4 1];
facesN=[5 6 7 8 5];

% 
colorP=[0.95,0,0];           % 
colorN=[0,0,0.7];            % 


% 

hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (No gradient)'
grid on


% 
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);


