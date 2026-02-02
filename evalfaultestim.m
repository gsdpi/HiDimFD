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

% Evaluation of additive-fault estimations, Monte Carlo method.
% Normal N(0,1) distribution of a few nonzero fault components
% Optional: nonzero initial state
% Optional: fault with norm 1

nonzeroinitialstate=1;	% It shouldn't matter
normalizefault=1;	% It shouldn't matter
TFPNrates=1;		% Replace TP,TN, FP, FN with their corresponding rates at the end
randprojrelation=1;	% Check factor sqrt((n-d)/n) in expected norm of projection residual vector
norm_est=2;		% Normalize fault estimation: 0=No, 1=magnitude, 2=components
			% norm_est>0 for threspcent to be actually a relative value

if ~exist('a','var')
	a=3;		% Number of non-zero fault vector components
end
if ~exist('threspcent','var')
	threspcent=12;	% Threshold for components in %
end
thres=threspcent/100;
if ~exist('estmode','var')
	estmode='res';
end
if ~exist('M','var')
	M=1000;
end

n=70;		% Ambient space dimension
I=eye(n,n);	% Identity matrix
fullhitrate=zeros(1,n-1);
meanTP=zeros(1,n-1);
meanFN=zeros(1,n-1);
meanFP=zeros(1,n-1);
meanTN=zeros(1,n-1);
F1score=zeros(1,n-1);
meanestnorm=zeros(1,n-1);
meandxnorm=zeros(1,n-1);
cosaver=zeros(1,n-1);
if strcmp(estmode,'res')
	modestr='projection residuals';
elseif strcmp(estmode,'sperbc')
	modestr='reconstruction based SPE contributions';
elseif strcmp(estmode,'lasso')
	lam = 0.002;
	modestr='LASSO';
elseif strcmp(estmode,'sr3')
	lam0 = 0.004;
	modestr='SR3 (sparse) solutions';
else
	error('Unknown method');
end
fprintf('Monte Carlo simulation: %s. Threshold: %.1f%%. # nonzero components: %d\n',modestr,threspcent,a);
for d=1:n-1	% Model space dimensions
	fprintf('.');
	[U,xp]=randaffinespace(n,d);	% Random affine subspace model
	if nonzeroinitialstate
		x0=U*randn(d,1)+xp;	% Nonzero normal initial state
	else
		x0=xp;			% Normal initial state: origin of affine subspace
	end
	fullhit=zeros(1,M);
	TP=zeros(1,M);
	FN=zeros(1,M);
	FP=zeros(1,M);
	TN=zeros(1,M);
	estnorm=zeros(1,M);
	dxnorm=zeros(1,M);
	cosa=zeros(1,M);
	F=(I-U*U');	% Fault space (model complementary) projection matrix
	for k=1:M
		[dx,faultsig]=randaddfault(n,a,normalizefault);	% Random additive fault
		dxnorm(k)=norm(dx);
		x=x0+dx;	% Obtain fault state
		if strcmp(estmode,'res')
			dx_est=F*(x-xp);			% Estimate fault as projection on residual space
		elseif strcmp(estmode,'sperbc')
			PCApar.Ctilde=F;
			spe_rbc=mspc_pca_spe_rbc(PCApar,x-xp,0);
			dx_est=sqrt(spe_rbc);
		elseif strcmp(estmode,'sr3')
			[dx_est, w0] = sr3(F, F*(x-xp), 'mode', '0', 'lam',lam0,'ptf',0);
		elseif strcmp(estmode,'lasso')
			%dx_est = LassoGaussSeidel(F,F*(x-xp),lam,'verbose',0);
            dx_est = lasso(F,F*(x-xp),Lambda=lam);
		else
			error('Unknown method');
		end
		estnorm(k)=norm(dx_est);
		cosa(k)=dot(dx,dx_est)/(dxnorm(k)*estnorm(k));
		if norm_est==1
			dx_est_orig=dx_est;
			dx_est=dx_est/estnorm(k);
		elseif norm_est==2
			dx_est_orig=dx_est;
			dx_est=dx_est/max(abs(dx_est));
		end
		[s,estcomp_idx]=sort(abs(dx_est),'descend'); % Sort estimated fault components
		% Fault-estimation components exceed threshold?
		thres_dx_est=(abs(dx_est)>=thres);
		% Which of the ones that exceed the threshold are involved in the actual fault?
		hits=(thres_dx_est==faultsig);
		% Is it a full hit?
		fullhit(k)=all(hits); % residual components over threshold match fault signature?
		% For F1-score computation
		TP(k)=sum(faultsig & thres_dx_est);
		FN(k)=sum(faultsig & ~thres_dx_est);
		FP(k)=sum(~faultsig & thres_dx_est);
		TN(k)=sum(~faultsig & ~thres_dx_est);
	end
	fullhitrate(d)=mean(fullhit);
	F1score(d)=2*sum(TP)/(2*sum(TP)+sum(FP)+sum(FN));
	meanestnorm(d)=mean(estnorm);
	meandxnorm(d)=mean(dxnorm);
	cosaver(d)=mean(cosa);
	% Extra info on F1-score results
	meanTP(d)=mean(TP);
	meanFN(d)=mean(FN);
	meanFP(d)=mean(FP);
	meanTN(d)=mean(TN);
	if TFPNrates
		meanTP(d)=meanTP(d)/a;
		meanFN(d)=meanFN(d)/a;
		meanFP(d)=meanFP(d)/(n-a);
		meanTN(d)=meanTN(d)/(n-a);
	end
end
fprintf('\n');
% Overall results based on areas under curves
AUCfullhitrate=sum(fullhitrate)/n*100;
AUC_TP=sum(meanTP)/n*100;
AUC_TN=sum(meanTN)/n*100;
AUC_FP=sum(meanFP)/n*100;
AUC_FN=sum(meanFN)/n*100;
AUC_F1score=sum(F1score)/n*100;
fprintf('AUC full-hite rate: %.2f%%\n',AUCfullhitrate);
fprintf('AUC True Positives: %.2f%%\n',AUC_TP);
fprintf('AUC True Negatives: %.2f%%\n',AUC_TN);
fprintf('AUC False Positives: %.2f%%\n',AUC_FP);
fprintf('AUC False Negatives: %.2f%%\n',AUC_FN);
fprintf('AUC F1-score: %.2f%%\n',AUC_F1score);
% Relative difference of dimensions (n-d)/n
reldiffdim=(n-(1:n-1))/n;

