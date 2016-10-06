function C = covariance(D,v,A,u,phi)
% Covariance  calculates the covariance matrix u for a multioutput
% dependent gaussian.
% FORMAT
% DESC Covariance  calculates the covariance matrix u for a multioutput
% DESC dependent gaussian.
% D : Input Dimention.
% v : Matrix Hyperparameters.
% A : Matrix Hyperparameters.
% u : Matrix Hyperparameters.
% phi : Hyperparameter.
% 
% output : Matrix C of covariance
%                 /       D/2                                                                            \
%            M   ?    (2pi)  *V(m,i)*V(m,j)       /    1                                             \   ?
% C(1,1)= SUM    ?  ----------------------   exp ?  - --- (ds-[u(m,i)-u(m,j))'*E*(ds-[u(m,i)-u(m,j)) ?  ?
%           m=1  ?                      ½         \    2                                             /   ?
%                ?    (?A(m,j)+A(m,i)?)                                                                 ?
%                 \                                                                                       /
% 
%     ? C(1,2)  C(1,n)  ?
% C = ?                 ?
%     ? C(1,n)  C(n,n)  ?
%
% COPYRIGHT: Wilmer Ariza, UTAS,2016
%

global x
global y

xi=size(x,1);
yi=size(y,1);
C = cell(yi);
for i=1:yi
    for j=1:yi
        C_ij = CIJ(i,j,D,v,A,u,xi,yi,x);

        if i==j
            C{i,j}=C_ij+phi(i)^2*ones(xi);
        else
            C{i,j}=C_ij;
        end
    end
end

