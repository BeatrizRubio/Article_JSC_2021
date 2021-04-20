function B=TNAddToPrevious(B,x,y,i) % i=2,...,n
% function B=TNAddToPrevious(B,x,y,i) % i=2,...,n
% An accurate, subtraction-free routine that computes BD(A*E), given 
% B=BD(A), where E=I_n, except E(i,i-1)=x, E(i-1,i-1)=y and E(i,i)=1/y
% In other words, given the bidiagonal decomposition B of a TN matrix A
% computes the bidiagonal decomposition of A*E. Multiplication by E means
% we add a multiple x of column i to i-1 and scale columns i and i-1 by y
% and 1/y respectively.
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written September 28, 2004

[m,n]=size(B);
if i<n, B(1,i+1)=B(1,i+1)*y; end
for j=1:min(i-1,m+1)
   if j>1
      B(j-1,i-1)=B(j-1,i-1)*y;
   end
   if (j<=m)
      y1=y+B(j,i)*x;
      B(j,i)=B(j,i)/y/y1;
      if (i<n) & (j<m)
         B(j+1,i+1)=B(j+1,i+1)*y1;
      end
      y=y1;
   end
end

% passing through the diagonal now ...
if i-1<=m
   B(i-1,i-1)=B(i-1,i-1)*y;
   if i<=m
      if B(i-1,i-1)~=0, x=B(i,i)*x/B(i-1,i-1); end;
      B(i,i)=B(i,i)/y; 
    
      % through the left now ...

      j=i;
      while (j<m) & (x>0)
         z=B(j,i-1);
         B(j,i-1)=z+x;
         y=B(j+1,i)/B(j,i-1);
         B(j+1,i)=y*z;
         x=y*x;
         j=j+1;
      end

      % fitting into the next lower triangular factor 
      B(m,i-1)=B(m,i-1)+x;
   end
end