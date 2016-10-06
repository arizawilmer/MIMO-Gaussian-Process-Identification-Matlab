function Laux = likelihood(T)

% DESC LIkelihood   calculates the likelihood of a covariance matrix for 
% DESC multioutput gaussian processes.
% 
% FORMAT
%   log(?C?)   y'C^-1 y      N
%  - ------ - ---------  -  ---  log(2*pi)
%      2          2          2
% Input
% % T : Vector of Hyperparameters(1,m)
% % Output
% fun=likelihood
% COPYRIGHT: Wilmer Ariza, UTAS,06/10/2016
%
global x y 
[n,m]=size(y);
D=5;
ytotal=reshape(y,[1,m*n])';
[Tn,Tm]=size(T);
T=reshape(T,[n,Tm/n]);
T= mat2cell(T, n,[n n n 1]);
v=cell2mat(T(1));
A=cell2mat(T(2));
u=cell2mat(T(3));
phi=cell2mat(T(4));
Cx= @covariance;
[xi,~]=size(x);
C=Cx(D,v,A,u,phi);
C = cell2mat(C);
C = nearestSPD(C);
R=2 * sum(log(diag(chol(C)))); 
Laux=-1/2*ytotal'*(C^(-1)*ytotal)-1/2*R-1*xi/2*log(2*pi);

end