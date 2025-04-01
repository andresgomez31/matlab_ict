%Graficar campos vectoriales con quiver
%{
x = linspace(-2,2,20);
y = x;
[X,Y] = meshgrid(x,y);

U = (1-2*X.^2).*exp(-X.^2-Y.^2);

V = (-2*Y.^2).*exp(-X.^2-Y.^2);

Z = X.*exp(-X.^2-Y.^2);

quiver(X,Y,U,V)
%}

%Campo vectorial sobre una esfera

[X, Y, Z] = sphere(40);  
figure;
surf(X, Y, Z, 'FaceAlpha', 0.5); %FaceAlpha modifica la opacidad de las caras de la malla
shading interp;
axis equal;
hold on;

[Xq, Yq, Zq] = sphere(30);

% Definir el campo vectorial en cada punto: F = (-y, x, 0)
U = Xq-Yq;
V = Yq-Xq;
W = zeros(size(Zq));

% Graficar el campo vectorial con quiver3
quiver3(Xq, Yq, Zq, U, V, W, 'k');

hold off;

