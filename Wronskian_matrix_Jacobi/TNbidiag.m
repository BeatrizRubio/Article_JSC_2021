% function B=TNbidiag(B)
%
% If B=BD(A), reduces A to bidiagonal form C using Givens rotations 
% and returns BD(C).
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written October 19, 2004

function B=TNbidiag(B);

if size(B,1)<size(B,2), B=B'; end
[m,n]=size(B);
for i=1:n
   for j=m:-1:i+1
      x=B(j,i);
      B(j,i)=0;
      c=sqrt(1+x*x);
      B=(TNAddToPrevious(B',x/c,c,j))';
   end
   for j=n:-1:i+2
      x=B(i,j);
      B(i,j)=0;
      c=sqrt(1+x*x);
      B=TNAddToPrevious(B,x/c,c,j);
   end      
end
