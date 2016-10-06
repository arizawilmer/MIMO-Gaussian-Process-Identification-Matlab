function y_prime = prediction(D,x_prime,T,C_inv)
% DESC Prediction   calculates Y* based in the calculation of K(x*?x) .
% DESC giving X* and the inverse of the covariance.
% 
% FORMAT
%
% Input
% D : Input Dimention.(scalar)
% T : Vector of Hyperparameters(1,m)
% C_inv : inverse of covarriance matrix(nxn)
% Output
% y_prime=1,m
% COPYRIGHT: Wilmer Ariza, UTAS,06/10/2016
%
global x y ytotal compare
[n,~]=size(y);
[~,Tm]=size(T);
T=reshape(T,[n,Tm/n]);
T= mat2cell(T, n,[n n n 1]);
v=cell2mat(T(1));
A=cell2mat(T(2));
u=cell2mat(T(3));
phi=cell2mat(T(4));
[xi,yi]=size(x);
% Cx= @covariance;
% C=Cx(1,v,A,u,phi);
hi=size(x_prime,2);
%C = cell(2);
for i=1:yi
    for j=1:yi
        
        sj=1;
        Sa=x_prime(i);
        for sj=1:xi
            Sb=x(sj,i);
            d=Sa-Sb;
            SUM=0;
            for m=1:yi
            SIGMA=(A(i,m)*(A(i,m)+A(j,m))^(-1)*A(j,m));
            SUM=((2*pi)^(D/2)*v(i,m)*v(j,m))/(sqrt(abs(A(j,m)+A(i,m))))...
                *exp(-1/2*(d-(u(i,m)-u(j,m)))'*SIGMA*(d-(u(i,m)-u(j,m))))+SUM;
            end
            C_ij_prime(1,sj)=SUM;
            SUM=0;
        end
            if i==j
            C_prime{i,j}=C_ij_prime+phi(i)^2;
        else
            C_prime{i,j}=C_ij_prime;
        end
        
            

    end
end
C_prime= cell2mat(C_prime);


sum=0;

y_prime=C_prime*C_inv*ytotal;
