function [xev,K]=sliding_window(x,tamano,solap,ventana)
%
%   [xev,K]=sliding_window(x,tamano,solap,ventana)
%       xev:        señal "enventanada": cada columna un tramo
%       K:          numeros de muestra que comienza cada tramo
%       x:          señal original
%       tamano:     tamaño de ventana
%       solap:      solapamiento en tanto por uno
%       ventana:    ponderacion de cada elemento dentro de la ventana
%                   cuando se omite es 1 para todos (v. rectangular)

n=length(x);
delta=tamano-round(solap*tamano);
numtramos=floor((n-tamano)/delta)+1;
xev=zeros(tamano,numtramos);
K=zeros(1,numtramos);
x=x(:);
if(nargin>3)
    ventana=ventana(:);
    for k=1:numtramos
        inic=(k-1)*delta;
        K(k)=inic;
        xev(:,k)=x(inic+1:inic+tamano).*ventana;
    end
else
    for k=1:numtramos
        inic=(k-1)*delta;
        K(k)=inic;
        xev(:,k)=x(inic+1:inic+tamano);
    end
end