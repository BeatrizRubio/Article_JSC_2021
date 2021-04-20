% function A=TNInverseExpand(B)  
%
% Computes directly the inverse a square TN matrix whose bidiagonal
% bidiagonal decomposition is stored in B, using the 
% results on the factorization of A and its inverse presented in:
% Ana Marco, Jose-Javier Martinez:  Accurate computations with totally 
% positive Bernstein-Vandermonde matrices.
% Electronic Journal of Linear Algebra, Volume 26 (2013): 357--380.
%
% This is mathematically equivalent to TNExpand(TNInverse(B)), but
% is done directly without the computing the bidiagonal decomposition of
% the inverse first. As a result, this is an O(n^2) algorithm.
%
% Written by Ana Marco and Jose-Javier Martinez, 2018


function A=TNInverseExpand(B)
[m,n]=size(B);

if m~=n 
    error('The matrix is not square');
end
A=zeros(n);
if isequal(class(B(1,1)),'sym'), A=sym(A); end
for i=1:n
  A(i,i)=1;
end  

for i=1:n-1
   for j=n:-1:i+1
     A(:,j)=A(:,j)-B(i,j)*A(:,j-1);
   end
end
for i=1:n
   A(:,i)=A(:,i)/B(i,i);
end
for i=n-1:-1:1
   for j=i:n-1
     A(:,j)=A(:,j)-B(j+1,i)*A(:,j+1);
   end
end



