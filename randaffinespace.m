function [U,x0]=randaffinespace(n,d)
% Random affine space
%	[U,x0]=randaffinespace(n,d)
%
%		n: ambient space dimension
%		d: affine subspace dimension
%
%		U (n x d matrix): subspace (orthonormalized) base vectors
%		x0: affine space origin
%

x0=randn(n,1);			% Affine space origin
V=[];
while(rank(V)~=d)
	V=randn(n,d);		% Generate random hyperplane base, dimension d
end
U=gramschmidt(V);		% Gram-Schmidt process (orthonormalizing)

