
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Step 1 Session4 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%-------------------------Configuración de los eritrocitos----------------------%
% Define las características de los eritrocitos

%;                                 % Define el parámetro ajustable "h". Considera que 
                                   % h<0.03 implica un glóbulo rojo saludable y h>0.2 uno infectado 
                                   %

%rad=(xmax_original-xmin_original)/40;  % Define el radio del eritrocito como una fracción de la escala espacial x

%;                             % Crea una variable para definir el número
%                               de eritrocitos

%;                             % Crea una variable de paso en el tiempo "dt" con unidades aribitrarias (usa el valor 0.2)

%;                             % Crea una variable "qe" para definir la carga eléctrica del eritrocito 
                               % (Igualmente positivas y negativas (por qué?)),
                               % usa un valor de 1 micro Coulombs.
                         
                            


%-------------------------Desplazamiento de los eritrocitos-----------------------%
% Inicializa la cuenta de eritrocitos infectados y sanos
%healthy=0; infected=0;



% (3) Abre un ciclo for sobre el número total de eritrocitos (Ne) definido por el usuario, 
% y dentro de éste inicializa las posiciones de los eritrocitos y las
% velocidades.
% (4) Dentro del ciclo, y después de la inicialización anterior define un
% ciclo while que corre desde la coordenada más alta Y del eritrocito,
% digamos "ye" hasta la más abajo "ye>min".
% 
% (5)Finalmenta, dentro de este while abre otro ciclo for que corra sobre
% todas las cargas en las placas.


