%%%%%%%%%%%%%%%%%%%%%%
%Step3_Session3_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%POR HACER%
%{
(1) Usando el concepto de potencial eléctrico V(i,j) y su relación con el
campo eléctrico, calcula (usando los tres ciclos for anidados,
preferentemente), el campo eléctrico usando la función de gradiente de
MATLAB.

(2) Usando el comando streamslice nuevamente, crea un nuevo plot del campo
eléctrico antes mencionado, usando el gradiente.

(3) Prueba la robustez, cambiando Nq, Lp y Ln y observa qué sucede con las
líneas de campo.
%}



%(3) Test the robustness of you code via changing the Nq, Lp and Ln
%values (what happens with the field lines?) - Expect questions of this
%behavior in your Oral Exam.

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

%---------------------Inicia posicionado las cargas------------------%
% Define un diferencial de carga lineal
%dq=?/?;                     % Mangitud de diferencial de carga

% Define las posiciones de las cargas
yp=linspace(-(1-p)*Lp/2,(1-p)*Lp/2,Nq); % Posiciones Y de las cargas positivas
xp(1:Nq)=-d/2-t/2;                      % Posiciones X de las cargas positivas
yn=                                     % Posiciones Y de las cargas negativas
xn(?:?)=                                % Posiciones X de las cargas negativas

%{
Descomenta y completa la siguiente plantilla de código para averigur si
estás posicionando las cargas correctamente. De forma adicional, juega con
algunos de los parámetros, para asegurar de que las cargas discretas están
distribuidas en tus placas. Finalmente, cambia la longitud de las placas y
observa si el código responde adecuadamente.
%}

%plot(<vector de las X para las cargas positivas>, <vector de las Y para las cargas positivas>,'*')
%hold on
%plot(?,?,'*') Lo mismo pero para las cargas negativas


%-------Cáclulo del Campo eléctrico para cada punto XY (Sin el gradiente)-------%

% % Incializa el potencial aquí
 V(?,?)=0;

% Inicializar componentes del campo eléctrico
 Ex = zeros(?, ?);  %¿En dónde quieres calcular esta componente?
 Ey = zeros(?, ?);
 
% % Calcula las componentes del campo eléctrico
%Aquí empiezan los ciclos anidados...

 for 




%Aquí comienzan a calcular las componentes de campos eléctricos.
%Calcula también V(i,j), el potencial para la i-esima y j-ésima posición.



 end

%Los tres ciclos anidados terminan aquí...

% Calcula las componentes del campo eléctrico usando el gradiente del
% potencial
% aquí

[?, ?] = gradient(V'); %¿Por qué necesitamos V' en lugar de V? 

% 
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectrophoresis (No gradient)'
grid on

%Parámetros estéticos (para el campo eléctrico) usando los valores del potencial (sólo descomenta).
%pcolor(x,y,V')                % Mapa de color del potencial
%colormap bone                 % Color



%Usa aquí el comando streamslice
streamslice(?,?,?,?,?)  % Líneas de campo eléctrico sin el gradiente 
                        % ¿Por qué necesitamos Ex' y Ey' en lugar de Ex y Ey?

%Descomenta y añade comentarios ¿Para qué son estas nuevas líneas?
%shading interp;
%colorbar

% 
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);

%%%-- Grafica el campo eléctrico usando el gradiente --%%%
%Usa las líneas previas de código como plantilla, usando los mismos parámetros de graficación


