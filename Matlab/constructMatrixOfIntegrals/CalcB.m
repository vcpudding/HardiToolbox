function B=CalcB(theta,phi,L)

for i=1:size(theta),
    n=1;
    for k=0:2:L,
        Yj=SH(k,theta(i),phi(i));
        B(i,n:n+size(Yj,2)-1)=Yj;
        n=n+size(Yj,2);
    end
end

    
return
