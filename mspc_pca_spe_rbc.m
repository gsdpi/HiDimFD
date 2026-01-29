function spe_rbc=mspc_pca_spe_rbc(PCApar,X,normalize)
%
%   spe_rbc=mspc_pca_spe_rbc(PCApar,X,normalize)
%
%

% Alcala09

if nargin<3
  normalize=1;
end

if normalize
  meanX=PCApar.meanX;
  stdX=PCApar.stdX;

  norm_par.meanx=meanX;
  norm_par.stdx=stdX;
  norm_par.type='std';
  X=normalize_data('apply',X,norm_par);
end

Ct=PCApar.Ctilde;
n=size(Ct,1);
m=size(X,2);
I=eye(n,n);
spe_rbc=(Ct*X).^2./repmat(diag(Ct),1,m);

