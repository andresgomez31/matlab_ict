%%%%%%%%%%%%%%%%%%%%%%
%Paso4_Sesion4_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%POR HACER%

% (1) Después de cerrar (finalizar) el bucle while, y antes de finalizar el bucle for de glóbulos rojos, introduce el siguiente contador. 

% ¿Qué hace esto? -describe-
     %if xe<0.07
     %    healthy=healthy+1;
     %else
     %    infected=infected+1;
     %end

% (2) Finalmente, después de finalizar el bucle for de glóbulos rojos, agrega los cálculos estadísticos finales,
% y describe cada línea.

% ¿Qué calcula esto?
%infected_probability = 100*infected/(healthy+infected);

% ¿Qué crea esto?
%str = ['Probabilidad de infección de  ', num2str(infected_probability), ' %'];
%annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', str, 'FitBoxToText', 'on');

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

%---------------------Inicio del posicionamiento de cargas------------------%
% Definir una carga diferencial lineal
%dq=?/?;                     % Magnitud del diferencial de carga

% Definir las posiciones de las cargas
yp=linspace(-(1-p)*Lp/2,(1-p)*Lp/2,Nq); % Posiciones Y de las cargas positivas
xp(1:Nq)=-d/2-t/2;                      % Posiciones X de las cargas positivas
yn=                                     % Posiciones Y de las cargas negativas
xn(?:?)=                                % Posiciones X de las cargas negativas

%Descomenta y completa el siguiente código plantilla para verificar si
%estás colocando correctamente tus cargas. Además, juega con algunos de los
%parámetros anteriores para asegurarte que las cargas discretas estén distribuidas
%entre tus placas. Finalmente, cambia la longitud de tus placas y verifica si
%tu código responde correctamente.

%plot(<vector de posiciones X de cargas positivas>, <vector de posiciones Y de cargas positivas>,'*')
%hold on
%plot(?,?,'*') Lo mismo pero para las cargas negativas


%-------Cálculo del campo eléctrico para cada punto XY (Sin gradiente)-------%

% % Inicializar el potencial aquí
 V(?,?)=0;

% % Inicializar componentes del campo eléctrico
 Ex = zeros(?, ?);  %¿Dónde quieres calcular esta componente?
 Ey = zeros(?, ?);
 
% % Calcular componentes del campo eléctrico
%Tres bucles for anidados comienzan aquí...

 for 

%Aquí deberías calcular las componentes del campo eléctrico.
%Y también calcular el potencial V(i,j) en la posición (i,j).


 end

%Terminan los tres bucles for anidados...

% Calcular las componentes del campo eléctrico usando el gradiente del potencial
% aquí

[?, ?] = gradient(V'); %¿Por qué necesitamos V' en lugar de V? 

% 
hold on
axis ([xmin xmax ymin ymax])
xlabel 'posición x, mm'
ylabel 'posición y, mm'
title 'Dielectroforesis (Sin gradiente)'
grid on

%Ajustar estética (para el gráfico del campo E) usando los valores del potencial recién añadidos (solo descomentar).
%pcolor(x,y,V')                % Mapa de color del voltaje
%colormap bone                 % Colores


%Usa aquí el comando streamslice para graficar el campo eléctrico calculado...
streamslice(?,?,?,?,?)  % Líneas del campo eléctrico sin gradiente 
                        % ¿Por qué necesitamos Ex' y Ey' en lugar de Ex y Ey?

%Descomenta y comenta (¿para qué sirven estas nuevas líneas?)
%shading interp;
%colorbar

% 
patch('Faces',facesP,'Vertices',vertices2d,'FaceColor',colorP);
patch('Faces',facesN,'Vertices',vertices2d,'FaceColor',colorN);

%%%-- Grafica aquí ahora el campo E usando el gradiente --%%%
%Usa las líneas anteriores como plantilla. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Paso 1 Sesión4 comienza aquí!%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------Configuración de los eritrocitos----------------------%
% Definir las características de los eritrocitos

%;                                 % Define un parámetro ajustable usando la letra "h". Considera que 
                                   % h<0.03 estará asociado a sangre sana y h>0.2 a sangre infectada 
                                   % (para un conjunto final de parámetros dados).

%rad=(xmax_original-xmin_original)/40;  % Define el radio del eritrocito como una fracción de la escala espacial en x.

%;                             % Crea una variable para definir el número de eritrocitos a modelar.

%;                             % Crea una variable de paso de tiempo "dt" en unidades arbitrarias (usa un valor de 0.2)

%;                             % Crea una variable "qe" para definir la carga eléctrica del eritrocito 
                               % (igualmente positivos y negativos - ¿por qué?),
                               % usa un valor de 1 micro Coulomb.
                         
                            


%-------------------------Desplazamiento de los eritrocitos-----------------------%
% Inicializar el conteo de eritrocitos sanos e infectados
%healthy=0; infected=0;

% Abre un bucle for que se ejecutará sobre el número total de eritrocitos (Ne) definido por el usuario, 
% y dentro de él inicializa la posición del eritrocito (xe, ye) y las velocidades (Vx , Vy).

% Luego, dentro de este bucle y después de la inicialización mencionada anteriormente, usa (define) un while que se 
% ejecute desde la coordenada Y superior del eritrocito, digamos "ye", hasta la inferior  
% (ye>ymin). 

% Finalmente, dentro de este while abre otro bucle for que se ejecutará desde la carga 1 hasta la última 
% carga en las placas (electrodos).


%-------------------------Configuración de los eritrocitos----------------------%
% Definir las características de los eritrocitos

%<faltan cosas aquí>


%-------------------------Desplazamiento de los eritrocitos----------------------%

%<faltan cosas aquí>


% Bucle sobre cada eritrocito
for ery=1:?
    

    % --Inicializar la posición y velocidades del eritrocito--
     xe=?; ye=?;    % Posición inicial del eritrocito
     Vy=?; Vx=?;    % Velocidades iniciales
      
     % Bucle hasta que el eritrocito alcance el fondo del área de trabajo
     while ? > ?

 
         % Bucle sobre las cargas en las placas (electrodos), para calcular la fuerza total sobre el eritrocito.

         for ?=?:?
                
         
         

         end

         % Calcular la nueva posición del eritrocito asumiendo masa=1
         
         %código delete(head) debe ir aquí

     end               % Fin del bucle while

 
end                 % Fin del bucle for
