% Computes the singular values of a TN matrix A with bidiagonal 
% decomposition B=BD(A)
%
% Written February 2003
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
%
% function a=TNSingularValues(B);

function a=TNSingularValues(B);

B=TNbidiag(B);
d=diag(B);
[m,n]=size(B);
a=mexdlasq1(d,d(1:length(d)-1).*diag(B,(m>=n)-(m<n)));




