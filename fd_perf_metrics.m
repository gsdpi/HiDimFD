function [f1s,fh,ppv,MCC] = fd_perf_metrics(f,fe,fthres)
%
%   function [f1s,fh,ppv,MMC] = fd_perf_mtrics(f,fe)
%

fs = (f~=0);        % Assume that the fault has actual zero components
fes = (abs(fe)>fthres);
TP = sum(fs & fes);
FP = sum(fs < fes);
FN = sum(fs > fes);
n = length(f);
TN = n-TP-FP-FN;
%TNc = sum(~fs & ~fes);
%assert(TN==TNc);
f1s = 2*sum(TP)/(2*sum(TP)+sum(FP)+sum(FN));
fh = all(fs==fes);
ppv = TP/(TP+FP);
MCC = matthewscc(TP,TN,FP,FN);