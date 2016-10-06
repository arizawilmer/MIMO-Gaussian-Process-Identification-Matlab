% Sample calculation of a MIMO gaussian process indetificacion  
% Based in:
% 
% Boyle, Phillip, and Marcus Frean. "Dependent gaussian processes." 
% Advances in neural information processing systems. 2004.
% 
%
clc;
clear;
%% Prepare figures to be plot
clf(figure(1))
clf(figure(2))
clf(figure(3))
%% Global variables
global y x xtotal ytotal 

%% Start data

X1=0:0.5:2;
X2=0:0.5:2;
X3=0:0.5:2;

D=5;
% Funtions who define output as example
Y1=2*sin(3*X1)-2*cos(X2);

Y2=2*sin(3*X2)+cos(Y1+X1);%+acos(Y1);
Y3=2*sin(3*X1)-2*cos(X2);%+acos(Y1);

x=[X1;X2;X3];
y=[Y1;Y2;Y3];
xtotal=[x(1,:) x(2,:) x(3,:)];
ytotal=[y(1,:) y(2,:) y(3,:)]';
x=[X1;X2;X3]';
%% Plot input data

figure(1)
scatter(X1,Y1,'*');hold on;
figure(2)
scatter(X2,Y2,'*');hold on;
figure(3)
scatter(X3,Y3,'*');hold on;

%%  Define an call funtion likelihood to find hyperparameters and X0 values
fun=@likelihood;

options = optimoptions(@lsqnonlin,...
     'MaxFunEvals',50000,'MaxIter',5000);

%  x0=[6.9962   -3.6025   -0.1201    3.5852   -0.3883   -9.5339    6.3237...
%      3.8123    3.7202   -1.2218   -2.6995    0.2454   -6.4137    2.2627...
%      2.2141    1.0104    5.8044   -1.3004   -1.2155    0.8543    1.1748...
%      3.5965    2.8772    6.1168    0.7694    6.4642    1.0811   -6.3159...
%      6.6458   -8.8765];

%options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',50000,'MaxIter',10000);
[n,m]=size(y);
x0=normrnd(0,2^2,[1,n*n*3+n]);


options1 = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxFunEvals',50000,'MaxIter',5000000);

%% call to optimization

param = lsqnonlin(fun,x0,[],[],options)
%param = fminsearch(fun,x0,options1) 

%% Output values
T=param;
[Tn,Tm]=size(T);
T=reshape(T,[n,Tm/n]);
T= mat2cell(T, n,[n n n 1]);
v=cell2mat(T(1))
A=cell2mat(T(2))
u=cell2mat(T(3))
phi=cell2mat(T(4))


%%  Call funtion covariance to calculate the covariance matrix 
% with the new hyperparameters

Cx= @covariance;
C_u=Cx(D,v,A,u,phi);
C=cell2mat(C_u);
C(1,:)*C^-1*ytotal
%C(43,:)*C^-1*ytotal
C_inv=(C)^-1;

%%  Prediction of values
X1=0:0.05:2;
X2=0:0.05:2;
X3=0:0.05:2;
x2=[X1;X2;X3]';
sample_s=size(X1,2);

for i=1:sample_s
y_prime(i,:)=prediction(D,x2(i,:) ,param,C_inv);
end


%% PLOT results
% Funtions who define output as example
Y1=2*sin(3*X1)-2*cos(X2);

Y2=2*sin(3*X2)+cos(Y1+X1);%+acos(Y1);
Y3=2*sin(3*X1)-2*cos(X2);%+acos(Y1);
x=[X1;X2;X3];
y=[Y1;Y2;Y3];


figure(1)
plot(X1,Y1);hold on;
plot(x2(:,1),y_prime(:,1));hold on;
figure(2)
plot(X2,Y2);hold on;
plot(x2(:,2),y_prime(:,2));hold on;
figure(3)
plot(X3,Y3);hold on;
plot(x2(:,3),y_prime(:,3));hold on;



