format long E


%Algorithm of the
%Bidiagonal decomposition of the collocation matrix of Jacobi polynomials presented in  
%Accurate computations with collocation and Wronskian matrices of Jacobi
%polynomials (2021), Scientific of Computing. To appear.
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Experimental results in Mathematica: CollocationMatrix_Jacobi_LinearSystem.nb

n=25

syms alpha
syms beta
syms t
f=((t-1)/2);

alpha=1 %(alpha, beta >-1) Corollary 1
beta=2


%Collocation matrix Jacobi polynomials
 
 for i=1: n
       t(i)= 1+i/(n+1);
 end 
 
 B = sym(zeros(n,n));
 
for col=1:size(B,2)
    n=col;
    v = sym(zeros(n,1)); 
    
    for i=1:n
        v(i)=((gamma(alpha+n)/(factorial(n-1)*gamma(alpha+beta+n))))*nchoosek(n-1,i-1)*gamma(alpha+beta+n-1+(i-1)+1)/gamma(alpha+(i-1)+1)*f^(i-1);
    end
    
    B(1,col)=simplify(sum(v));
end

X=B(1,:);

eval_t=t;

for i=1:n
    t=eval_t(i);
    A(i,:)=eval(X);    
end

     
 
%Bidiagonal factoriztion of de Collocation matrix of Jacobi polynomials.

%1. Bidiagonal factorization of  the matrix of Theorem 2 or Algorithm 1.  

n=25; %(same n as above)
BDA1=zeros(n,n);

%1.1 Computation of the multipliers m_{i,j}
BDA1= sym(zeros(n,n));
  
  for i=2: n 
	M=  (alpha+i-1)/(i-1); 
	BDA1( i,1)=M;
  
 	for j=2: i-1
  	     M=M* (alpha+beta+2*i-j)/(alpha+beta+2*i-j-2);   
 	     BDA1(i,j)= M ;
% 	   
      end 
  end
 
%1.2 Computation of the pivots p_{i,i}
 BDA1(1,1)=1;
 q=1;
  for  i=2:n
    BDA1(i,i)=(1/factorial(i-1));
     aux=1;
     for k=2:i
         aux=aux*(alpha+beta+i+(k-2));
     end
     BDA1(i,i)=BDA1(i,i)*aux;
  end 

    
%2. Bididagonal factorization Collocation matrix  polynomial basis (t-1)/2. See Algorithm 2

 for i=1: n
   t(i)= 1+i/(n+1); %Same t's as above.
 end 
 
 BDA2=zeros(n,n);

%2.1 Computation of the multipliers m_{i,j}
 
for i=2: n 
	M= 1; 
	BDA2( i,1)=M;
   
	for j=2: i-1
	     M=M*((t(i))-(t( i-j+1)))/((t(i-1))-(t(i-j)));   
	     BDA2(i,j)= M  ; 
	    
    end 
end


%2.2 Computation of the tilde multipliers m_{i,j}
 BDA2(1,1)=1;
 
  q=1;
   for  i=2:n
  	q=q;
   	aux=(1/2)^(i-1);
   	for k=1: i-1
   		aux=aux*(t(i)-t(k));
   	end 
        BDA2(i,i)=q*aux;     
   end

%2.3 Computation of pivots p_{i,j}

 for j=1:n-1	
    coef= (t(j)-1)/2;
    for i= j+1:n 
   	BDA2(j,i)=coef;
   	end
 end


%3. Computation of the bidiagonal factorization of Wronskian Jacobi.
%Remark 3 or Algorithm 3

BDA=TNProduct(BDA2,transpose(BDA1));

% function A=TNProduct(A,B)
%
% given BD(A) and BD(B), computes BD(AB)
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written September 29, 2004


    
%4. Inverse Matrix
IB=TNInverseExpand(BDA);
IM=inv(W);
dlmwrite('inverseCollocationMatrixJacobiB.csv',IB,'precision','%.45f');
dlmwrite('inverseCollocationMatrixJacobiM.csv',IM,'precision','%.45f');

 %function A=TNInverseExpand(B)  
%Computes directly the inverse a square TN matrix whose bidiagonal
% bidiagonal decomposition is stored in B, using the 
% results on the factorization of A and its inverse presented in:
% Ana Marco, Jose-Javier Martinez:  Accurate computations with totally 
% positive Bernstein-Vandermonde matrices.
% Electronic Journal of Linear Algebra, Volume 26 (2013): 357--380.


