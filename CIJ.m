function C_ij = CIJ(i,j,D,v,A,u,xi,yi,x)

% C_ij calculates the covariance matrix C_ij for a multioutput
% dependent gaussian.
% FORMAT
% DESC C_ij  calculates the covariance matrix C_ij for a single pair of
% outputs 
% DESC dependent gaussian.
% D : Input Dimention.
% v : Matrix Hyperparameters.
% A : Matrix Hyperparameters.
% u : Matrix Hyperparameters.
% phi : Hyperparameter.
% i: Position inside of covariance matrix
% j: Position inside of covariance matrix
% xi: Number of data X
% yi: Number of white noises
%
% 
% output : Matrix C_ij of covariance
%         /       D/2                                                                            \
%    M   ?    (2pi)  *V(m,i)*V(m,j)       /    1                                             \   ?
% SUM    ?  ----------------------   exp ?  - --- (ds-[u(m,i)-u(m,j))'*E*(ds-[u(m,i)-u(m,j)) ?  ?
%   m=1  ?                      ½         \    2                                             /   ?
%        ?    (?A(m,j)+A(m,i)?)                                                                 ?
%         \                                                                                       /
%
% COPYRIGHT: Wilmer Ariza, UTAS,2016
%

for si=1:1:xi
    sj=1;
    Sa=x(si,i);
    for sj=1:xi
        Sb=x(sj,j);
        d=Sa-Sb;
        SUM=0;
        for m=1:yi
            SIGMA=(A(i,m)*(A(i,m)+A(j,m))^(-1)*A(j,m));
            SUM=((2*pi)^(D/2)*v(i,m)*v(j,m))/(sqrt(abs(A(j,m)+A(i,m))))...
                *exp(-1/2*(d-(u(i,m)-u(j,m)))'*SIGMA*(d-(u(i,m)-u(j,m))))+SUM;
        end
        C_ij(si,sj)=SUM;
    end
end
