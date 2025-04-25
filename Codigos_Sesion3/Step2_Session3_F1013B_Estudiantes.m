%%%%%%%%%%%%%%%%%%%%%%
%Step2_Session3_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%POR HACER%
%{
(1) Usando tres ciclos for anidades, uno desde 1=1:Nx, otro desde j=1:Ny
para cubrir el espacio, o otro desde k=1:Nq para cubrir todas las cargas en
las placas. Calcula el campo eléctrico total dentro y alrededor de las
placas paralelas. Usa como pista el código dado y recuerda que un sólo
punto del espacio tiene contribuciones de todas y cada una de las cargas
eléctricas en las placas (i.e. positivas y negativas).
%}

%{
(2) Usando el comando de MATLAB llamado "streamslice", grafica el campo
eléctrico mencionado previamente, ya sabes cómo se debe ver.
%} .

%{
(3) Prueba la robustez del código cambiando Nq, Lp y Ln (¿Qué ocurre con las líneas de campo?)
este tipo de preguntas vendrán en tu Examen Oral.
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
% Inicializar componentes del campo eléctrico
 Ex = zeros(?, ?);  %¿En dónde quieres calcular esta componente?
 Ey = zeros(?, ?);
 
% % Calcula las componentes del campo
%Los tres ciclos for anidados comienzan aquí...

 for 




%Aqui debes calcular las componentes de E.




 end

%Los ciclos for terminan aquí...

% 
hold on
axis ([xmin xmax ymin ymax])
xlabel 'x position, mm'
ylabel 'y position, mm'
title 'Dielectroforesis (Sin gradiente)'
grid on

%Usa aquí el comando streamslice
streamslice(?,?,?,?,?)  % Líneas de campo eléctrico sin el gradiente 
                        % ¿Por qué necesitamos Ex' y Ey' en lugar de Ex y Ey?

% 
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);


