function [X,D]=jeig(A,B);
L=chol(B,'lower');
G=inv(L);
C=G*A*G';
[Q,D]=schur(C);  %Q是一个正交矩阵  D是一个上三角矩阵
X=G'*Q; 