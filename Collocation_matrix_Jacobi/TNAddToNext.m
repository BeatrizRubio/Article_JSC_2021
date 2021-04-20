% function B=TNAddToNext(B,x,i)
%
% Computes BD(E_i(x)*A) where B=BD(A); i>=2
% In other words, given the bidiagonal decomposition B of a matrix 
% A, computes the bidiagonal decomposition of C, where C is obtained from 
% A by adding a multiple x of row i-1 to i
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
%
% Written October 19, 2004

function B=TNAddToNext(B,x,i)


[m,n]=size(B);


z=0; 
while (z<min(i-1,n)) & ((z==0) | (B(i-1,z)==0))
   % put the appropriate diagonals of B in b and c

   b=zeros(m-i+1,1);
   c=zeros(m-i+1,1);
   for j=0:m-i
      if z+j>0
         if z+j<=n, 
            b(j+1)=B(j+i,z+j);
         end
      else b(j+1)=x;   
      end
      if z+j+1<=n
         c(j+1)=B(j+i,z+j+1);
      end
   end

   [b,c,q]=dqd2(b,c); % q=0 means no nonzeros in c were created
                      % if nonzeros were created in c, the new nonzeros may have to be chased

   % return the new b and c in the same diagonals of B
   for j=0:m-i
      if (z+j>0) & (z+j<=n), B(j+i,z+j) =b(j+1); end
      if (z+j+1<=n), B(j+i,z+j+1) =c(j+1); end
   end
      
   i=i+q-1;
   z=z+q;
end


function [b,c,i]=dqd2(b,c);

% takes 2 lower bidiagonal matrices with offdiagonal elements b and c
% and returns BD((I+diag(b,-1))*(I+diag(c,-1))
% if i>0 then c(i)=0 became positive and may need to be chased further

n=length(b);

t=c(1);
c(1)=b(1)+c(1);
d=b(1);
b(1)=0;

i=1;
while (i<n)&(b(i+1)>0) 
   e=b(i+1)/c(i);
   d=e*d;
   b(i+1)=e*t;
   t=c(i+1);
   c(i+1)=c(i+1)+d;
   i=i+1;
end
