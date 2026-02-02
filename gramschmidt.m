function U=gramschmidt(V,Upart)
%
% U=gramschmidt(V,Upart)
%

%
% Author:
%	Abel Alberto Cuadrado Vega
%
% 	Grupo de Supervisión, Diagnóstico y Descubrimiento de Conocimiento en Procesos de Ingeniería
% 	(c) Universidad de Oviedo, 2023-2026
%

n=size(V,2);
if nargin==2
	U=Upart;
	ini=size(U,2)+1;
else
	ini=2;
	U=V(:,1)/norm(V(:,1));
end
for i=ini:n
    u_i=V(:,i);
    for j=1:i-1
        u_i=u_i-dot(V(:,i),U(:,j))*U(:,j);
    end
    U=[ U u_i/norm(u_i) ];
end
