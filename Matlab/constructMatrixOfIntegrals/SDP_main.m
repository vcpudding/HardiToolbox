close all;
clear all;

UnitVectors;
g81=g([1:81],:);
g321=g([1:321],:);

b=1500;
order=4; delta=500;
sigma=0.018;
v1=[1/sqrt(2) 1/sqrt(2) 0];
orientation1=atan2(v1(2),v1(1));

s=SimulateData(b,g81,pi/6,orientation1);

%S=ricernd(S,sigma);

S=DWI_to_SH(s,g81,g321,8,0.006,sigma,false,0.05)'; % regularize and over-sample 81->321

C=constructMatrixOfIntegrals(g321,order,delta);

c=S*C;

G=constructMatrixOfMonomials(g321,order);

m=zeros(size(g321,1),1);

echo on
n=15;

cvx_begin
   variable x(n)
   dual variables y
   minimize(norm(C*x-S',2))
   subject to
     y : G * x > m;
     %x > 0;
cvx_end

echo off

TensorODF=x;

figure,plotTensors(TensorODF,1,[321 1 0]);
W=convert_DT4(TensorODF([15 5 1 12 3 10 11 8 7 14 13 9 4 6 2]));

[a,b,c]=cp_als(tensor(W),2,'tol',1e-6,'init','nvecs');










