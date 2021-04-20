% function x=TNSolve(B,b)
% Solves a TN linear system Ax=b, where B=BD(A).
%
% This yields Bjorck-Pereyra type methods for every TN matrix A for which
% we have B=BD(A). A (and thus B) must be square
%
% Plamen Koev, Massachusetts Institute of Technology, December 23, 2005
%
% Copyright Plamen Koev, see COPYRIGHT.TXT for details.

function x=TNSolve(B,b)
  x=b;
  [m,n]=size(B);
  if m~=n
      echo('The matrix must be square'); 
      exit; 
  end;
  for i=1:n-1
      for j=n-i+1:n
          x(j)=x(j)-x(j-1)*B(j,j-n+i);
      end
  end
  
  for i=1:n
      x(i)=x(i)/B(i,i);
  end
  
  for i=1:n-1
      for j=n-1:-1:i
          x(j)=x(j)-x(j+1)*B(j-i+1,j+1);
      end
  end
  