function [f,fsig]=randaddfault(n,a,normalized)
% Random additive fault vector
%
%	[f,fsig]=randaddfault(n,a,normalized)
%
%		n: space dimension, number of components
%		a: number of nonzero-valued components [optional, default n]
%		normalized (bool): normalize fault [optional, default 0]
%		f: fault vector
%		fsig: fault signature [only if "a" is specified]
%

%
% Author:
%	Abel Alberto Cuadrado Vega
%
% 	Grupo de Supervisión, Diagnóstico y Descubrimiento de Conocimiento en Procesos de Ingeniería
% 	(c) Universidad de Oviedo, 2023-2026
%

if nargin<3
	normalized=0;
end
f=randn(n,1);				% Obtain random additive fault vector
if nargin>=2
	comp_idx=randperm(n);
	faultcomp_idx=comp_idx(1:a);
	fsig=zeros(n,1);
	fsig(faultcomp_idx)=1;		% Fault signature: nonzero fault components
	f=fsig.*f;			% Set zero components
end
if normalized
	f=f/norm(f);			% Normalize fault
end

