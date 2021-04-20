function A=TNProduct(A,B)

% function A=TNProduct(A,B)
%
% given BD(A) and BD(B), computes BD(AB)
%
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for more details.
% Written September 29, 2004

[m,n]=size(A);
[n1,p]=size(B);
if (n~=n1)| ((m>n)&(n<p)), fprintf('matrix dimensions do not agree or product not TN\n'); 
   return; 
end
   
for i=1:n-1
    for j=n-i+1:min(n,n+p-i)
        A=TNAddToPrevious(A,B(j,i+j-n),1,j);
    end
end

if p<n, A=A(:,1:p); end
A=TNDiagonalScale(diag(B),A')';
A=[A zeros(m,p-n)];

for i=p-1:-1:1
    for j=p:-1:p-i+1
        A=TNAddToNext(A',B(i+j-p,j),j)';
    end
end
