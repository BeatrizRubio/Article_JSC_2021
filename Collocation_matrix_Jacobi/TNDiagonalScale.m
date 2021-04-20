function B=TNDiagonalScale(f,B);

% function B=TNDiagonalScale(f,B);
%
% if B=BD(A), computes BD(diag(f)*A), works on rectangular matrices
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written September 28, 2004

[m,n]=size(B);
B(1,1)=B(1,1)*f(1);
for i=2:m
    if i<=n, B(i,i)=B(i,i)*f(i); end
    B(i,1:min(i-1,n))=B(i,1:min(i-1,n))*(f(i)/f(i-1));
end
