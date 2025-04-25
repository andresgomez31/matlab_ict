%%%%%%%%%%%%%%%%%%%%%%
%Paso2_Sesion4_F1013B
%%%%%%%%%%%%%%%%%%%%%%

%POR HACER%

% (1) Dentro del bucle for que se ejecutará sobre el número total de
% eritrocitos (Ne), define una variable llamada path, usando el
% comando de Matlab "animatedline".

 %path=animatedline(...);
 
% Dentro de este mismo bucle, copia (antes de la inicialización de la posición y velocidad del eritrocito)
% y descubre para qué sirve el siguiente código:

 %dx=h*rand()*rad;       %¡Piensa profundamente en términos de dielectroforesis!

 %¿Qué está pasando aquí? ¡Completa el texto que falta!
 %fprintf('Eritrocito: %d tiene un ??????? asociado de : %f\n', ery, dx);
 
 
% (2) Dentro del while que se ejecutará desde la coordenada Y más alta del eritrocito, digamos "ye", hasta la más baja  
% (ye>ymin), agrega el siguiente código y comenta para qué sirve. Además, inicializa una variable "Fx",
% para calcular la componente en x de la fuerza eléctrica de coulomb
% experimentada por el glóbulo rojo.

 % ¿Para qué sirve esto?
         %addpoints(path,xe,ye); 
         %head=scatter(xe,ye, 100,'filled','o','red'); 
         %drawnow  
   
         %--------------Cálculo de la Fuerza-------------------%
         % Inicializar fuerza eléctrica
         % ;


% (3) Finalmente, también dentro del bucle while, pero después de abrir y cerrar el bucle for que se ejecutará desde la carga 1 hasta la última 
% carga en las placas (electrodos), agrega el siguiente código y como
% antes, explica.

% ¿Para qué sirve esto?
         %if ye>ymin+dx
         %    delete(head);    
         %end

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

%---------------------Iniciar posicionamiento de cargas------------------%
% Definir una carga diferencial lineal
%dq=?/?;                     % Magnitud del diferencial de carga

% Definir las posiciones de las cargas
yp=linspace(-(1-p)*Lp/2,(1-p)*Lp/2,Nq); % Posiciones en Y de cargas positivas
xp(1:Nq)=-d/2-t/2;                      % Posiciones en X de cargas positivas
yn=                                     % Posiciones en Y de cargas negativas
xn(?:?)=                                % Posiciones en X de cargas negativas

%Descomenta y completa el siguiente código para comprobar si
%estás colocando correctamente tus cargas. Además, juega con algunos de los
%parámetros anteriores, para asegurarte de que las cargas discretas están distribuidas
%entre tus placas. Finalmente, cambia la longitud de tus placas y verifica si
%tu código responde adecuadamente.

%plot(<vector de posiciones X de cargas positivas>, <vector de posiciones Y de cargas positivas>,'*')
%hold on
%plot(?,?,'*') Lo mismo pero para cargas negativas


%-------Cálculo del campo eléctrico para cada punto XY (sin gradiente)-------%

% % Inicializar potencial aquí
 V(?,?)=0;

% % Inicializar componentes del campo eléctrico
 Ex = zeros(?, ?);  %¿Dónde quieres calcular esta componente?
 Ey = zeros(?, ?);
 
% % Calcular componentes del campo eléctrico
%Comienzan tres bucles for anidados...

 for 




%En algún lugar aquí debes calcular las componentes del campo eléctrico.
%En algún lugar debes calcular el potencial V(i,j) para la posición i,j.



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

%Ajustar estética (para la gráfica del campo E) usando los valores de potencial recién añadidos (solo descomenta).
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
%Usa líneas anteriores como plantilla
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Paso 1 Sesión4 comienza aquí! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------Configuración de los eritrocitos----------------------%
% Define las características de los eritrocitos

%;                                 % Define un parámetro ajustable  usando la letra "h". Considera que 
                                   % h<0.03 estará asociado con sangre sana y h>0.2 con sangre infectada 
                                   % (para un conjunto final de parámetros dados).

%rad=(xmax_original-xmin_original)/40;  % Define el radio del eritrocito como una fracción de la escala del espacio en x.

%;                             % Crea una variable para definir el número de eritrocitos a modelar.

%;                             % Crea una variable de paso de tiempo "dt" en unidades arbitrarias (usa un valor de 0.2)

%;                             % Crea una variable "qe" para definir la carga eléctrica del eritrocito 
                               % (positivos y negativos por igual - ¿por qué?),
                               % usa un valor de 1 micro Coulomb.
                         
                            


%-------------------------Desplazamiento de los eritrocitos-----------------------%
% Inicializar el conteo de eritrocitos sanos e infectados
%healthy=0; infected=0;


% Abre un bucle for que se ejecutará sobre el número total de eritrocitos (Ne) definido por el usuario, 
% y dentro de él inicializa la posición (xe, ye) y velocidades (Vx , Vy) del eritrocito... :)

% Luego, dentro de este bucle y después de la inicialización mencionada anteriormente, usa (define) un while que se 
% ejecute desde la coordenada Y más alta del eritrocito, "ye", hasta la más baja  
% (ye>ymin). 

% Finalmente, dentro de este while abre otro bucle for que se ejecutará desde la carga 1 hasta la última 
% carga en las placas (electrodos).


%-------------------------Configuración de los eritrocitos----------------------%
% Define las características de los eritrocitos

%<falta contenido aquí>


%-------------------------Desplazamiento de los eritrocitos----------------------%

%<falta contenido aquí>


% Bucle sobre cada eritrocito
for ery=1:?
    

    % --Inicializar la posición y velocidades del eritrocito--
     xe=?; ye=?;    % Posición inicial del eritrocito
     Vy=?; Vx=?;    % Velocidades iniciales
      
     % Bucle hasta que el eritrocito llegue al fondo del área de trabajo
     while ? > ?

 
         % Bucle sobre las cargas en las placas (electrodos), para calcular la fuerza total sobre el eritrocito.

         for ?=?:?
                
            

            
         end

 

     end               % Fin del bucle while

 
end                 % Fin del bucle for
