
function W=estimate_ODF_with_SDP(S,GradientOrientations,order,delta)

% S=DWI_to_SH(s,g81,g321,8,0.006,sigma,false,0.05)'; % regularize and over-sample 81->321
n=15;


C=constructMatrixOfIntegrals_642(GradientOrientations,order,delta);

c=S*C;

G=constructMatrixOfMonomials(GradientOrientations,order);

m=zeros(size(GradientOrientations,1),1);

echo on

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

%figure,plotTensors(TensorODF,1,[321 1 0]);
W=convert_DT4(TensorODF([15 5 1 12 3 10 11 8 7 14 13 9 4 6 2]));

%[a,b,c]=cp_als(tensor(W),2,'tol',1e-6,'init','nvecs');










