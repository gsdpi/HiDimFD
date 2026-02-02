% sporadic_and_random_faults:
%	Set the value of the variable "fault_type" below
%	Run this script to generate the results of the experiments in Section IV
%
% Title:
%	Generation of Interpretable Residuals for Fault Diagnosis based on
%	Projection Techniques: Leveraging Variable Redundancy
% Authors:
%	Abel Alberto Cuadrado Vega, Ignacio Díaz Blanco, José María Enguita González,
%	Diego García Pérez and Ana González Muñiz
%
% 	Grupo de Supervisión, Diagnóstico y Descubrimiento de Conocimiento en Procesos de Ingeniería
% 	(c) Universidad de Oviedo, 2023-2026
%

% Sporadic and random faults (dataicann)
printaspdf=0;

% Load data: http://hdl.handle.net/10651/53461
% Alternative: https://digibuo.uniovi.es/dspace/handle/10651/53461
load dataicann.mat

y=z{3};			% Only the third condition (normal operation)
saxn=y(:,2);		% horizontal vibration (acceleration)
sayn=y(:,3);		% vertical vibration (acceleration)
sirn=y(:,4);		% r-phase current
sisn=y(:,5);		% s-phase current
% Normalize data:
%	remove mean
saxn=saxn-mean(saxn);
sayn=sayn-mean(sayn);
sirn=sirn-mean(sirn);
sisn=sisn-mean(sisn);
%	scaling using standard deviation
saxn=saxn/std(saxn);
sayn=sayn/std(sayn);
sirn=sirn/std(sirn);
sisn=sisn/std(sisn);

M=size(y,1);		% Total number of samples
ftrain=0.8;		% Fraction of number of samples used for training
Mtrain=round(ftrain*M);
Nv=300;			% Window size (samples)
solap=0.998;		% Overlapping (fraction of window size)
% ... results in a stride of 1 sample for the sliding window
ax=sliding_window(saxn(1:Mtrain),Nv,solap);
ay=sliding_window(sayn(1:Mtrain),Nv,solap);
ir=sliding_window(sirn(1:Mtrain),Nv,solap);
is=sliding_window(sisn(1:Mtrain),Nv,solap);
% Training data matrix 15675 vector samples (rows) of 300 x 4 = 1200 components (columns)
X=[ax' ay' ir' is'];
axv=sliding_window(saxn(Mtrain+1:end),Nv,solap);
ayv=sliding_window(sayn(Mtrain+1:end),Nv,solap);
irv=sliding_window(sirn(Mtrain+1:end),Nv,solap);
isv=sliding_window(sisn(Mtrain+1:end),Nv,solap);
% Test data matrix: 3695 vector samples (rows) of 300 x 4 = 1200 components (columns)
Xv=[axv' ayv' irv' isv'];
m=size(X,1);		% Number of training samples
mv=size(Xv,1);		% Number of test samples
fprintf('Training samples: %d. Validation samples: %d\n',m,mv);
% PCA model
[V,S]=eig(X'*X/(m-1));
%figure;plot(diag(S))
l=200;			% Number of retained principal components
% Descending sorting of eigenvalues, retain "l" highest ones and their corresponding eigenvectors
lambda=diag(S);
[lambda,idxsort]=sort(lambda,'descend');
V=V(:,idxsort);
% Percentage of the variance explained by the model
var_exp=cumsum(lambda)/sum(lambda)*100;
fprintf('Model of dimension %d explains %.2f%% of the variance\n',l,var_exp(l));
P=V(:,1:l);		% Retained principal components (column vectors)
n=1200;			% Number of variables (ambient space dimension)
I=eye(n,n);		% Identity matrix n x n
F=I-P*P';		% Residual subspace projection matrix

% Significance level
alpha=0.0001;   % VALUE IN ORIGINAL SUBMISSION
%alpha=0.001;
fprintf('Confidence level: %.3f\n',(1-alpha)*100);

force_ct_comp = 0;
if ~exist('rth','var') || force_ct_comp   % Component thresholds
    % Residuals
    R=F*X';
    %R=F*X(1:100,:)';
    rth=prctile(abs(R),(1-alpha)*100,2);
    % LASSO
    lam = 0.002;
    LDX=[];
    for k=1:size(X,1)
	    ldx = lasso(F,F*X(k,:)',Lambda=lam);
        LDX = [ LDX ldx ];
        if ~mod(k,1000)
            fprintf('*');
        end
    end
    fprintf('\n');
    lth=prctile(abs(LDX),(1-alpha)*100,2);
    % SR3
    lam0 = 0.004;	% lambda recommended by SR3 authors
    SR3DX=[];
    for k=1:size(X,1)
	    [sr3dx, w0] = sr3(F,F*X(k,:)', 'mode', '0', 'lam',lam0,'ptf',0);
        SR3DX = [ SR3DX sr3dx ];
        if ~mod(k,1000)
            fprintf('*');
        end
    end
    fprintf('\n');
    sr3th=prctile(abs(SR3DX),(1-alpha)*100,2);
    %figure;plot(rth);hold on;
    %plot(lth);
    %plot(sr3th);
    %legend('Res thres.','LASSO thres.','SR3 thres.');
end

s=rng(44); % seed used for the results of submission to IEEE TCST (Feb 2025 and Aug 2025 on GNU Octave; Jan 2026 on Matlab, different results)

% Types of faults
SPORADIC_FAULT = 1;
RANDOM_FAULT = 2;
% Selected type
fault_type = SPORADIC_FAULT;
%fault_type = RANDOM_FAULT;

Nf=20;			% Number of fault samples
nplots=2;		% Show only first "nplots" faults
%nplots=Nf;
Xnewn=[X' Xv'];		% Total matrix of nonfaulty data (samples as columns). Faulty data added later.
LDX1=[];		% LASSO fault estimations (samples as columns)
DX1=[];			% SR3 fault estimations (samples as columns)
DXreal=[];		% Actual fault vector (samples as columns)
R1=[];			% Residual vectors (samples as columns)
tempo_r=[];		% Residual computation timings
tempo_sr3=[];		% SR3 estimation timings
tempo_lasso=[];		% LASSO estimation timings
idxtestfault=[];
f1s_sr3=[];
f1s_lasso=[];
f1s_res=[];

if fault_type == SPORADIC_FAULT
	% For sporadic faults
	fc=50;			% Number of consecutive components set to zero to simulate a fault
	iif=randi(n-fc+1,1,Nf);	% index of first component set to zero in each of the faults
	liif=length(iif);
	ilimvib=2*Nv-fc+1;
	ilimcor=2*Nv+1;
	fallosvib=find((iif<=ilimvib));
	falloscor=find((iif>=ilimcor));
end

for i=1:Nf
	ixv=randi(mv);	% Random index of test sample to generate fault with
	idxtestfault=[ idxtestfault ixv ];
	xo=Xv(ixv,:)';	% ... and corresponding test sample
	x=xo;
	if fault_type == SPORADIC_FAULT
		x(iif(i):iif(i)+fc-1)=zeros(1,fc);	% set components to zero
	elseif fault_type == RANDOM_FAULT
		frn = randn(n,1);   % random fault
            frn = frn.*(abs(frn)>sr3th);
    		frn = frn.*(abs(frn)>lth);
    		frn = frn.*(abs(frn)>rth);
    		nzc_rat=sum(frn~=0)/n*100;
    		fprintf('Fault %d # nonzero components ratio: %.2f%%\n',i,nzc_rat);
    		x = x + frn;
	else
		error('Unknown fault type');
	end
	Xnewn=[Xnewn x];			% add faulty sample to total matrix of data
	%figure;plot(x)
    % LASSO estimation
	lam = 0.002;
    tic;
	ldx = lasso(F,F*x,Lambda=lam);
    ttt=toc;
	tempo_lasso=[tempo_lasso ttt];	% SR3 estimation timings
    LDX1=[LDX1 ldx];			% LASSO fault estimations (samples as columns)
	% SR3 estimation: 0-norm mode, kappa = 1, C = I
	lam0 = 0.004;	% lambda recommended by SR3 authors
    tic;
	[dx, w0] = sr3(F, F*x, 'mode', '0', 'lam',lam0,'ptf',0);
    ttt=toc;
	tempo_sr3=[tempo_sr3 ttt];	% SR3 estimation timings
	DX1=[DX1 dx];			% SR3 fault estimations (samples as columns)
    dxreal=x-xo;
	DXreal=[DXreal dxreal];		% Actual fault vector (samples as columns)
	tic;
	r=F*x;				% Residual vector computation (projection on residual subspace)
	ttt=toc;
	tempo_r=[tempo_r ttt];		% Residual computation timings
	R1=[R1 r];			% Residual vectors (samples as columns)
    f1s=fd_perf_metrics(dxreal,dx,sr3th);
    f1s_sr3=[f1s_sr3 f1s*100];
    f1s=fd_perf_metrics(dxreal,ldx,lth);
    f1s_lasso=[f1s_lasso f1s*100];
    f1s=fd_perf_metrics(dxreal,r,rth);
    f1s_res=[f1s_res f1s*100];
	if(i<=nplots)
		figure;
		subplot(2,1,2);plot(xo)		% Plot normal signal
		hold on;plot(x-dx,'r')		% Plot SR3 estimation of normal signal
        plot(x-ldx,'m')		% Plot LASSO estimation of normal signal
		plot(x,'g')			% Plot faulty signal
		% show boundaries between original signals (ax,ay, ir, is) within the 1200-component vectors
		yr=ylim;
		plot([Nv Nv],yr,'k')
		plot(2*[Nv Nv],yr,'k')
		plot(3*[Nv Nv],yr,'k')
		% Plot labeling
		legend('Normal signal','normal-signal SR3 estimation','normal-signal LASSO estimation','Faulty signal')
		title('Signals');
		xlabel('Components');
		% Plot additive fault and its estimations (SR3 and residual vector)
		subplot(2,1,1);plot(x-xo)	% Plot actual fault
		ylim(yr);
		hold on;plot(dx,'r');		% Plot fault SR3 estimation
        plot(ldx,'m');		% Plot fault LASSO estimation
		plot(r,'c');			% Plot residual vector
		% show boundaries between original signals (ax,ay, ir, is) within the 1200-component vectors
		yr=ylim;
		plot([Nv Nv],yr,'k')
		plot(2*[Nv Nv],yr,'k')
		plot(3*[Nv Nv],yr,'k')
		% Plot labeling
		legend('Actual fault','SR3 estimation','LASSO estimation','Residual');
		title('Faults');
		xlabel('Components');
	end
	fprintf('-');
end
fprintf('\n');

% F1 score
figure;
h1=plot(f1s_res,'c.-');
hold on;
h2=plot(f1s_lasso,'m.-');
h3=plot(f1s_sr3,'r.-');
xlabel('Fault samples');
ylabel('%');
title('F1-score of estimated fault components')
if fault_type == SPORADIC_FAULT
	plot(fallosvib,f1s_res(fallosvib),'co');
	plot(falloscor,f1s_res(falloscor),'c+');
	plot(fallosvib,f1s_lasso(fallosvib),'mo');
	plot(falloscor,f1s_lasso(falloscor),'m+');
    plot(fallosvib,f1s_sr3(fallosvib),'ro');
	plot(falloscor,f1s_sr3(falloscor),'r+');
end
legend([h1,h2,h3],'F1-s Res','F1-s LASSO','F1-s SR3')
ylim([0 100]);


if fault_type == SPORADIC_FAULT
	media_f1s_res_vibra = mean(f1s_res(fallosvib))
	media_f1s_res_corr = mean(f1s_res(falloscor))
	media_f1s_lasso_vibra = mean(f1s_lasso(fallosvib))
	media_f1s_lasso_corr = mean(f1s_lasso(falloscor))
    media_f1s_sr3_vibra = mean(f1s_sr3(fallosvib))
	media_f1s_sr3_corr = mean(f1s_sr3(falloscor))
end

% Plot RMSE of fault estimations
figure;
% RMSE of SR3 estimation of the fault
RMSE_sr3=vecnorm(DX1-DXreal,2)/sqrt(n);
plot(RMSE_sr3,'r.-');
hold on;
RMSE_lasso=vecnorm(LDX1-DXreal,2)/sqrt(n);
plot(RMSE_lasso,'m.-');
% RMSE of residual-based estimation of the fault
RMSE_r=vecnorm(R1-DXreal,2)/sqrt(n);
plot(RMSE_r,'c.-');
title('Fault estimation RMSE');
legend('SR3','LASSO','res');
xlabel('Fault samples');
yy=ylim;
yy(1)=0;
ylim(yy);
% LaTeX table of RMSE values and % improvements from residuals to SR3/LASSO
subtabla=[100*RMSE_r;100*RMSE_sr3;100*RMSE_lasso;(RMSE_r-RMSE_sr3)./RMSE_r*100;(RMSE_r-RMSE_lasso)./RMSE_r*100];
TABLA=[ subtabla mean(subtabla,2)];
latex_table('RMSE.tex',TABLA,...
'1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 & 14 & 15 & 16 & 17 & 18 & 19 & 20 & AVG',...
'$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$',...
'c');

% Plot and show values of estimation timings
% figure;
% plot(tempo_sr3,'r.-');
% hold on;
% plot(tempo_lasso,'m.-');
% plot(tempo_r,'c.-');
% title('Estimation times');
% legend('SR3','LASSO','res');
% xlabel('Fault sample #');
mean_t_r=mean(tempo_r)
mean_t_sr3=mean(tempo_sr3)
mean_t_lasso=mean(tempo_lasso)
relacion_t_sr3=mean_t_sr3/mean_t_r
relacion_t_lasso=mean_t_lasso/mean_t_r

% PCA diagnosis

if length(alpha)==1
    alphat=alpha;
    alphas=alpha;
    alphaf=alpha;
elseif length(alpha)==3
    alphat=alpha(1);
    alphas=alpha(2);
    alphaf=alpha(3);
else 
    error('Parameter vector alpha must have 1 or 3 components');
end

Ptilde=V(:,l+1:n);			% Principal components discarded in PCA model
LAMBDA=diag(lambda(1:l));		% Diagonal matrix of retained eigenvalues
LAMBDAtilde=diag(lambda(l+1:n));	% Diagonal matrix of discarded eigenvalues
C=P*P';					% PCA model subspace projection matrix
Ctilde=Ptilde*Ptilde';			% Residual subspace projection matrix

% SPE control limit
theta1=sum(lambda(l+1:n));
theta2=sum(lambda(l+1:n).^2);
gSPE=theta2/theta1;
hSPE=theta1^2/theta2;
delta2=gSPE*chi2inv(1-alphas,hSPE);

% Hotelling's T2 control limit
tau2=chi2inv(1-alphat,l);

% Combined index phi control limit
gphi=(l/tau2^2+theta2/delta2^2)/(l/tau2+theta1/delta2);
hphi=(l/tau2+theta1/delta2)^2/(l/tau2^2+theta2/delta2^2);
dseta2=gphi*chi2inv(1-alphaf,hphi);

% Elements of PCA method
PCApar.P=P;
PCApar.LAMBDA=LAMBDA;
PCApar.D=P*inv(LAMBDA)*P';
PCApar.Ctilde=Ctilde;
PCApar.PHI=Ctilde/delta2+PCApar.D/tau2;
PCApar.l=l;
PCApar.alpha=alpha;
PCApar.delta2=delta2;
PCApar.tau2=tau2;
PCApar.dseta2=dseta2;

% Control indexes
figure;

% Custom (based on residuals)
r=F*Xnewn;
%rn=norm(r,2,'columns');
rn=vecnorm(r,2);
rntrain=rn(1:m);
pctl9x=prctile(rntrain,(1-alpha)*100);
subplot(2,2,1);
plot(rn(1:m),'.-');
hold on;
plot(m+1:m+mv,rn(m+1:m+mv),'g.-');
plot(m+mv+1:length(rn),rn(m+mv+1:end),'m-.');
ttl=sprintf('||r|| index (thres: %.2f)',pctl9x);
title(ttl);
xlabel('Samples');
FPRres=mean(rntrain>pctl9x);
fprintf('False positive rate with r index (train): %.3f\n',FPRres*100);
FPRtres=mean(rn(m+1:m+mv)>pctl9x);
fprintf('False positive rate with r index (test): %.3f\n',FPRtres*100);
hold on;
xr=xlim;
plot(xr,[pctl9x pctl9x],'r')
legend('training','test','fault','control limit','location','northwest');

% Standard T2
D=PCApar.D;
T2=diag(Xnewn'*D*Xnewn);	% Inefficient, probably
subplot(2,2,2);
plot(T2(1:m),'.-');
hold on;
plot(m+1:m+mv,T2(m+1:m+mv),'g.-');
plot(m+mv+1:length(T2),T2(m+mv+1:end),'m-.');
ttl=sprintf('T2 index (thres: %.2f)',tau2);
title(ttl);
xlabel('Samples');
hold on;
xr=xlim;
plot(xr,[tau2 tau2],'r')
legend('training','test','fault','control limit','location','northwest');
T2train=T2(1:m);
FPRT2=mean(T2train>tau2);
fprintf('False positive rate with T2 index (train): %.3f\n',FPRT2*100);
FPRtT2=mean(T2(m+1:m+mv)>tau2);
fprintf('False positive rate with T2 index (test): %.3f\n',FPRtT2*100);

% Standard SPE
Ctilde=PCApar.Ctilde;
SPE=diag(Xnewn'*Ctilde*Xnewn);	% Inefficient, probably
subplot(2,2,4);
plot(SPE(1:m),'.-');
hold on;
plot(m+1:m+mv,SPE(m+1:m+mv),'g.-');
plot(m+mv+1:length(SPE),SPE(m+mv+1:end),'m-.');
ttl=sprintf('SPE index (thres: %.2f)',delta2);
title(ttl);
xlabel('Samples');
hold on;
xr=xlim;
plot(xr,[delta2 delta2],'r')
legend('training','test','fault','control limit','location','northwest');
SPEtrain=SPE(1:m);
FPRSPE=mean(SPEtrain>delta2);
fprintf('False positive rate with SPE index (train): %.3f\n',FPRSPE*100);
FPRtSPE=mean(SPE(m+1:m+mv)>delta2);
fprintf('False positive rate with SPE index (test): %.3f\n',FPRtSPE*100);

FPRSPET2=mean((SPEtrain>delta2)|(T2train>tau2));
fprintf('False positive rate with SPE|T2 index: %.3f\n',FPRSPET2*100);
FPRtSPET2=mean((SPE(m+1:m+mv)>delta2)|(T2(m+1:m+mv)>tau2));
fprintf('False positive rate with SPE|T2 index: %.3f\n',FPRtSPET2*100);


% Standard phi
PHI=PCApar.PHI;
phi=diag(Xnewn'*PHI*Xnewn);	% Inefficient, probably
subplot(2,2,3);
plot(phi(1:m),'.-');
hold on;
plot(m+1:m+mv,phi(m+1:m+mv),'g.-');
plot(m+mv+1:length(phi),phi(m+mv+1:end),'m-.');
ttl=sprintf('\\phi index (thres: %.2f)',dseta2);
title(ttl);
xlabel('Samples');
hold on;
xr=xlim;
plot(xr,[dseta2 dseta2],'r')
legend('training','test','fault','control limit','location','northwest');
phitrain=phi(1:m);
FPRphi=mean(phitrain>dseta2);
fprintf('False positive rate with phi index (train): %.3f\n',FPRphi*100);
FPRtphi=mean(phi(m+1:m+mv)>dseta2);
fprintf('False positive rate with phi index (test): %.3f\n',FPRtphi*100);

if printaspdf
    if fault_type == SPORADIC_FAULT
	    figure(1); exportgraphics(gcf,'figures/sporadic-current-fault.pdf','ContentType','vector')
	    figure(2); exportgraphics(gcf,'figures/sporadic-vibration-fault.pdf','ContentType','vector')
        figure(3); exportgraphics(gcf,'figures/F1s-sporadic-fault.pdf','ContentType','vector')
	    figure(4); exportgraphics(gcf,'figures/RMSE-sporadic-fault.pdf','ContentType','vector')
	    figure(5); exportgraphics(gcf,'figures/indexes-sporadic-fault.pdf','ContentType','vector')
    end
    if fault_type == RANDOM_FAULT
        figure(1); exportgraphics(gcf,'figures/random-fault-1.pdf','ContentType','vector')
        figure(2); exportgraphics(gcf,'figures/random-fault-2.pdf','ContentType','vector')
        figure(3); exportgraphics(gcf,'figures/F1s-random-fault.pdf','ContentType','vector')
	    figure(4); exportgraphics(gcf,'figures/RMSE-random-fault.pdf','ContentType','vector')
	    figure(5); exportgraphics(gcf,'figures/indexes-random-fault.pdf','ContentType','vector')
    end
end
