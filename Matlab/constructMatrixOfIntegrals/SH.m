function Yj=SH(k,theta,phi)

P=legendre(k,cos(theta))';

idx=1;

for m=-k:k;
    
    c=sqrt(((2*k+1).*factorial(k-abs(m)))./(4*pi.*factorial(k+abs(m))));
    
    if m==0,
        Yj(idx)=c*P(1);
        idx=idx+1;
    elseif m>=-k & m<0
         Yj(idx)=sqrt(2).*c.*P(abs(m)+1).*cos(m*phi);
        idx=idx+1;
    else Yj(idx)=sqrt(2).*c.*P(abs(m)+1).*sin(m*phi); 
        idx=idx+1;
    end
end
    

return

