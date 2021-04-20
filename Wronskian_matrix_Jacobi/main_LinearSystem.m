clear all
format long E


%Experiments results of the
%Bidiagonal decomposition of the Jacobi polynomials presented in  
%Accurate computations with collocation and Wronskian matrices of Jacobi
%polynomials (2021), Scientific of Computing. To appear.
%E. Mainar, J.M. Pe\~na, B. Rubio, 
%Experimental results in Mathematica: Wronskian_Jacobi_LinearSystem.nb

n=24

syms alpha
syms beta
syms t
f=((t-1)/2);

alpha=1 %(alpha, beta >-1) Corollary 1
beta=2


%Wronskian matrix Jacobi
 B = sym(zeros(n+1,n+1));
 
for col=1:size(B,2)
    n=col;
    v = sym(zeros(n,1)); 
    
    for i=1:n
        v(i)=((gamma(alpha+n)/(factorial(n-1)*gamma(alpha+beta+n))))*nchoosek(n-1,i-1)*gamma(alpha+beta+n-1+(i-1)+1)/gamma(alpha+(i-1)+1)*f^(i-1);
   
    end
    
    B(1,col)=simplify(sum(v));
end

for i=2:size(B,1)
    for j=1:size(B,2)
        B(i,j) = diff(B(1,j), i-1);
    end
end

 t=50; %(t>0)
 W=eval(B);  %Wronskian matrix Jacobi polynomials

 
%Bidiagonal factoriztion of de Wronskian matrix of Jacobi polynomials.

%1. Bidiagonal factorization  the matrix A. Theorem 2 or Algorithm 1.  

n=24; %(same n as above)
BDA1=zeros(n+1,n+1)

%1.1 Computation of the multipliers m_{i,j}
  for i=2: n+1
	M=  (alpha+i-1)/(i-1); 
	BDA1( i,1)=M;
  
 	for j=2: i-1
  	     M=M* (alpha+beta+2*i-j)/(alpha+beta+2*i-j-2);   
 	     BDA1(i,j)= M ;
 	   
      end 
  end
 
%1.2 Computation of the pivots p_{i,i}
 BDA1(1,1)=1;
 q=1;
  for  i=2:n+1
    BDA1(i,i)=(1/factorial(i-1));
     aux=1;
     for k=2:i
         aux=aux*(alpha+beta+i+(k-2));
     end
     BDA1(i,i)=BDA1(i,i)*aux;
  end 

    
%2. Bididagonal factorization Wronskian polynomial basis(t-1)/2. See Algorithm2
t=50 %(same t as above)

BDA2=zeros(n+1,n+1)

%2.1 Computation of the multipliers m_{i,j}
 
for i=2:n+1
     for j=1:i
         BDA2(i,j)= 0;
     end
end

%2.2 Computation of the tilde multipliers m_{i,j}
for i=1:n+1

    for j=i+1:n+1
        BDA2(i,j)= (t-1)/2;
    end
end

%2.3 Computation of pivots p_{i,j}

for i=0:n
    BDA2(i+1,i+1)=(1/2)^i*factorial(i);
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



%4. Solve linear system Wx=b


b=[17, -31, 77, -83, 27, -11, 96, -57, 70, -64, 29, -41,...
      46, -16, 74, -1, 2, -6, 7, -5, 1, -2, 6, -7, 5];

SolBD=transpose(TNSolve(BDA,b))
SolM = W\transpose(b)


dlmwrite('systemWronskianJacobiB.csv',SolBD,'precision','%.45f');
dlmwrite('systemWronskainJacobiM.csv',SolM,'precision','%.45f');

%function TNSolve(B,b)
%%Solves a TN linear system Ax=b, where B=BD(A). (see TNSolve of Plamen Koev https://math.mit.edu/~plamen/software/TNTool.html)

   
  


