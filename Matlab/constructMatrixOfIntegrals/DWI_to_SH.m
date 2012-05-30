function Snew=DWI_to_SH1(s,x,y,l,gamma,sigma,sharpening,alpha)

% s is the DWI signal, x is the directions vector, l is the spherical harmonics order and gamma is the
% coupling coefficient. Set gamma to 0.006.  sigma is the noise variance
% (Rician distribution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% add Rician noise %%%%%%

s=ricernd(s,sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phiX=acos(x(:,3));
thetaX=atan2(x(:,2),x(:,1));

phiY=acos(y(:,3));
thetaY=atan2(y(:,2),y(:,1));

%%%%%%%%%%%%% sharpening %%%%%%%%%%%%%%%%
idx=1;
for i=0:2:l,
    for j=1:2*i+1,
       SK(idx)=1+alpha*i*(i+1); % sharpening kernel
       idx=idx+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx=1;

B=CalcB(phiX,thetaX,l); %  calculate spherical harmonics, 0<phi<2pi and 0<theta<pi
By=CalcB(phiY,thetaY,l); 

for k=0:2:l,
    v=ones(1,2*k+1);
    Lv(idx:idx+size(v,2)-1)=k*v;
    idx=idx+size(v,2);
end

L=diag(Lv.^2.*(Lv+1).^2);
%P=diag(Pv);

T=(B'*B+gamma*L); %*inv(P);
S=B'*s';

C=T\S;

% extrapolate from 81 directions to 321 directions B->By

if sharpening
    Snew=By*(C.*SK');
else Snew=By*C;
end

figure,plot(s,'g');
hold on;
plot(Snew,'r');


