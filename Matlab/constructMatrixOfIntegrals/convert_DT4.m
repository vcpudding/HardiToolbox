function W=convert_DT4(D)
%converts a vectorized version D([1:15]) with unique coeeficients
%to a fully symmetric fourth order tensor W(:,:,:,:)
for i1=1:3
    for i2=1:3
        for i3=1:3
            for i4=1:3
                
                ix=0;
                iy=0;
                iz=0;
                
                if i1==1
                    ix=ix+1;
                elseif i1==2
                    iy=iy+1;
                elseif i1==3
                    iz=iz+1;
                end
                
                if i2==1
                    ix=ix+1;
                elseif i2==2
                    iy=iy+1;
                elseif i2==3
                    iz=iz+1;
                end
                
                if i3==1
                    ix=ix+1;
                elseif i3==2
                    iy=iy+1;
                elseif i3==3
                    iz=iz+1;
                end
                
                if i4==1
                    ix=ix+1;
                elseif i4==2
                    iy=iy+1;
                elseif i4==3
                    iz=iz+1;
                end
                
             
                if (ix==4)&(iy==0)&(iz==0)
                    W(i1,i2,i3,i4)=D(1)/computeFactor(4,0,0);
                elseif (ix==0)&(iy==4)&(iz==0)
                    W(i1,i2,i3,i4)=D(2)/computeFactor(0,4,0);
                elseif (ix==0)&(iy==0)&(iz==4)
                    W(i1,i2,i3,i4)=D(3)/computeFactor(0,0,4);
                elseif (ix==2)&(iy==2)&(iz==0)
                    W(i1,i2,i3,i4)=D(4)/computeFactor(2,2,0);
                elseif (ix==0)&(iy==2)&(iz==2)
                    W(i1,i2,i3,i4)=D(5)/computeFactor(0,2,2);
                elseif (ix==2)&(iy==0)&(iz==2)
                    W(i1,i2,i3,i4)=D(6)/computeFactor(2,0,2);
                elseif (ix==2)&(iy==1)&(iz==1)
                    W(i1,i2,i3,i4)=D(7)/computeFactor(2,1,1);
                elseif (ix==1)&(iy==2)&(iz==1)
                    W(i1,i2,i3,i4)=D(8)/computeFactor(1,2,1);
                elseif (ix==1)&(iy==1)&(iz==2)
                    W(i1,i2,i3,i4)=D(9)/computeFactor(1,1,2);
                elseif (ix==3)&(iy==1)&(iz==0)
                    W(i1,i2,i3,i4)=D(10)/computeFactor(3,1,0);
                elseif (ix==3)&(iy==0)&(iz==1)
                    W(i1,i2,i3,i4)=D(11)/computeFactor(3,0,1);
                elseif (ix==1)&(iy==3)&(iz==0)
                    W(i1,i2,i3,i4)=D(12)/computeFactor(1,3,0);
                elseif (ix==0)&(iy==3)&(iz==1)
                    W(i1,i2,i3,i4)=D(13)/computeFactor(0,3,1);
                elseif (ix==1)&(iy==0)&(iz==3)
                    W(i1,i2,i3,i4)=D(14)/computeFactor(1,0,3);
                elseif (ix==0)&(iy==1)&(iz==3)
                    W(i1,i2,i3,i4)=D(15)/computeFactor(0,1,3);
                end             
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function counter=computeFactor(x,y,z)
counter=0;
for i1=1:3
    for i2=1:3
        for i3=1:3
            for i4=1:3
                ix=0;
                iy=0;
                iz=0;
                if i1==1
                    ix=ix+1;
                elseif i1==2
                    iy=iy+1;
                elseif i1==3
                    iz=iz+1;
                end
                
                if i2==1
                    ix=ix+1;
                elseif i2==2
                    iy=iy+1;
                elseif i2==3
                    iz=iz+1;
                end
                
                if i3==1
                    ix=ix+1;
                elseif i3==2
                    iy=iy+1;
                elseif i3==3
                    iz=iz+1;
                end
                
                if i4==1
                    ix=ix+1;
                elseif i4==2
                    iy=iy+1;
                elseif i4==3
                    iz=iz+1;
                end
                if (ix==x)&(iy==y)&(iz==z)
                    counter=counter+1;
                end
            end
        end
    end
end
