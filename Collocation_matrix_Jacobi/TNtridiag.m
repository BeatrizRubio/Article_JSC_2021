function [B,V]=TNtridiag(B)

% function [B,V]=TNtridiag(B)
%
% Using similarity transformation reduced a matrix A to tridiagonal form T
% Input: B=BD(A)
% Output: BD(T), tranformation matrix V, such that T=V^{-1}AV
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
%
% July 2015: Added the return of the transformation matrix V,
%            Supported by SJSU's Woodward Fund

n=size(B,1);
V=eye(n);

for i=1:n-2
   for j=n:-1:i+2
      x=B(j,i);
      B(j,i)=0;
      B=TNAddToPrevious(B,x,1,j);
      V(:,j-1)=V(:,j-1)+x*V(:,j);
      x=B(i,j);
      B(i,j)=0;
      B=(TNAddToPrevious(B',x,1,j))';
      V(:,j)=V(:,j)-x*V(:,j-1);
   end
end

% symmetrizing ... 

x=ones(n,1);
for i=1:n-1
   x(i+1)=x(i)*sqrt(B(i,i+1)/B(i+1,i));
   z=sqrt(B(i+1,i)*B(i,i+1));
   B(i+1,i)=z;
   B(i,i+1)=z;
   
   V(:,i+1)=V(:,i+1)/x(i+1);
end

