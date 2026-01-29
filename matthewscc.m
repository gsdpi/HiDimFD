function MCC = matthewscc(TP,TN,FP,FN)
%
%	MCC = matthewscc(TP,TN,FP,FN)
%

% N = TN+TP+FN+FP;
% S = (TP+FN)/N;
% P = (TP+FP)/N;
% MCC = (TP/N-S*P)/sqrt(P*S*(1-S)*(1-P));

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
