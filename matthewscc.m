function MCC = matthewscc(TP,TN,FP,FN)
%
%	MCC = matthewscc(TP,TN,FP,FN)
%

%
% Author:
%	Abel Alberto Cuadrado Vega
%
% 	Grupo de Supervisión, Diagnóstico y Descubrimiento de Conocimiento en Procesos de Ingeniería
% 	(c) Universidad de Oviedo, 2023-2026
%

% N = TN+TP+FN+FP;
% S = (TP+FN)/N;
% P = (TP+FP)/N;
% MCC = (TP/N-S*P)/sqrt(P*S*(1-S)*(1-P));

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
