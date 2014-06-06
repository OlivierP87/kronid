A=rand(3,3);

kai=KronAI(A,3,2);
kai2=kron(A,eye(3,2));
fprintf(1,'Test A  kron I, error = %f\n',norm(kai-kai2))

kia=KronIA(A,3,2);
kia2=kron(eye(3,2),A);
fprintf(1,'Test I  kron A, error = %f\n',norm(kia-kia2))

A=rand(4,4);
B=A'*A; % symmetric

b=svecmex(B);
fprintf(1,'Test : trace(B*B) = %f,  svec(B)^T*svec(B) = %f\n',trace(B*B),b'*b)

B2=svecimex(b);
fprintf(1,'Test svecimex, error = %f\n',norm(B2-B))

% compute vec of B
b2=B(:);
P=vecPsvec(size(B,1));
fprintf(1,'Test vecPsvec, error = %f\n',norm(b2-P*b))

% compute symmetric kronecker product
% 1/2( CXB'+ BXC') = G --->  B kronsym C . svec(X) = svec(G)
X=rand(4,4);
X=X'*X;
G=1/2*(A*X+X*A');
g=svecmex(G);
mskron=skronAI(A);
g2=mskron*svecmex(X);
fprintf(1,'Test skronai, error = %f\n',norm(g2-g))
