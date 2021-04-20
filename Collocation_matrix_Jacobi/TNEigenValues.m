% function a=TNEigenValues(B)
%
% Computes the eigenvalues of a TN matrix with bidiagonal decomposition
% stored in B
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written February 2003
%
% July 2015: Added the return of the eigenvector matrix V
%            Supported by SJSU's Woodward Fund



function [a,V]=TNEigenValues(B)


[B,V]=TNtridiag(B);
n=size(B,1);

% symmetrizing ... need to form the matrix explicitly first to compute the diagonal
%% symmetrizing factors
for i=1:n-1
   B(i+1,i)=sqrt(B(i+1,i)*B(i,i+1));
   B(i,i+1)=B(i+1,i);
end

d=sqrt(diag(B));
e=sqrt(diag(B,1).*diag(B,-1)).*d(1:length(d)-1);


a=mexdlasq1(d,e).^2;

  

[U,S,Z]=svd(diag(d)+diag(e,1));
V=V*Z;



