function W=convertD2W_order6(D)
%converts a vectorized version D([1:15]) with unique coeeficients
%to a fully symmetric fourth order tensor W(:,:,:,:)
for i1=1:3
    for i2=1:3
        for i3=1:3
            for i4=1:3
                for i5=1:3
                    for i6=1:3
                
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
                
                if i5==1
                    ix=ix+1;
                elseif i5==2
                    iy=iy+1;
                elseif i5==3
                    iz=iz+1;
                end
                
                if i6==1
                    ix=ix+1;
                elseif i6==2
                    iy=iy+1;
                elseif i6==3
                    iz=iz+1;
                end
                
             
                if (ix==0)&(iy==0)&(iz==6)
                    W(i1,i2,i3,i4,i5,i6)=D(1)/computeFactor(0,0,6);
                elseif (ix==0)&(iy==1)&(iz==5)
                    W(i1,i2,i3,i4,i5,i6)=D(2)/computeFactor(0,1,5);
                elseif (ix==0)&(iy==2)&(iz==4)
                    W(i1,i2,i3,i4,i5,i6)=D(3)/computeFactor(0,2,4);
                elseif (ix==0)&(iy==3)&(iz==3)
                    W(i1,i2,i3,i4,i5,i6)=D(4)/computeFactor(0,3,3);
                elseif (ix==0)&(iy==4)&(iz==2)
                    W(i1,i2,i3,i4,i5,i6)=D(5)/computeFactor(0,4,2);
                elseif (ix==0)&(iy==5)&(iz==1)
                    W(i1,i2,i3,i4,i5,i6)=D(6)/computeFactor(0,5,1);
                elseif (ix==0)&(iy==6)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(7)/computeFactor(0,6,0);
                elseif (ix==1)&(iy==0)&(iz==5)
                    W(i1,i2,i3,i4,i5,i6)=D(8)/computeFactor(1,0,5);
                elseif (ix==1)&(iy==1)&(iz==4)
                    W(i1,i2,i3,i4,i5,i6)=D(9)/computeFactor(1,1,4);
                elseif (ix==1)&(iy==2)&(iz==3)
                    W(i1,i2,i3,i4,i5,i6)=D(10)/computeFactor(1,2,3);
                elseif (ix==1)&(iy==3)&(iz==2)
                    W(i1,i2,i3,i4,i5,i6)=D(11)/computeFactor(1,3,2);
                elseif (ix==1)&(iy==4)&(iz==1)
                    W(i1,i2,i3,i4,i5,i6)=D(12)/computeFactor(1,4,1);
                elseif (ix==1)&(iy==5)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(13)/computeFactor(1,5,0);
                elseif (ix==2)&(iy==0)&(iz==4)
                    W(i1,i2,i3,i4,i5,i6)=D(14)/computeFactor(2,0,4);
                elseif (ix==2)&(iy==1)&(iz==3)
                    W(i1,i2,i3,i4,i5,i6)=D(15)/computeFactor(2,1,3);
                elseif (ix==2)&(iy==2)&(iz==2)
                    W(i1,i2,i3,i4,i5,i6)=D(16)/computeFactor(2,2,2);
                elseif (ix==2)&(iy==3)&(iz==1)
                    W(i1,i2,i3,i4,i5,i6)=D(17)/computeFactor(2,3,1);
                elseif (ix==2)&(iy==4)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(18)/computeFactor(2,4,0);
                elseif (ix==3)&(iy==0)&(iz==3)
                    W(i1,i2,i3,i4,i5,i6)=D(19)/computeFactor(3,0,3);
                elseif (ix==3)&(iy==1)&(iz==2)
                    W(i1,i2,i3,i4,i5,i6)=D(20)/computeFactor(3,1,2);
                elseif (ix==3)&(iy==2)&(iz==1)
                    W(i1,i2,i3,i4,i5,i6)=D(21)/computeFactor(3,2,1);
                elseif (ix==3)&(iy==3)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(22)/computeFactor(3,3,0);
                elseif (ix==4)&(iy==0)&(iz==2)
                    W(i1,i2,i3,i4,i5,i6)=D(23)/computeFactor(4,0,2);
                elseif (ix==4)&(iy==1)&(iz==1)
                    W(i1,i2,i3,i4,i5,i6)=D(24)/computeFactor(4,1,1);
                elseif (ix==4)&(iy==2)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(25)/computeFactor(4,2,0);
                elseif (ix==5)&(iy==0)&(iz==1)
                    W(i1,i2,i3,i4,i5,i6)=D(26)/computeFactor(5,0,1);
                elseif (ix==5)&(iy==1)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(27)/computeFactor(5,1,0);
                elseif (ix==6)&(iy==0)&(iz==0)
                    W(i1,i2,i3,i4,i5,i6)=D(28)/computeFactor(6,0,0);
                end
            end
        end
    end
        end
    end
end


function counter=computeFactor(x,y,z)
counter=0;
for i1=1:3
    for i2=1:3
        for i3=1:3
            for i4=1:3
                for i5=1:3
                    for i6=1:3
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
                
                if i5==1
                    ix=ix+1;
                elseif i5==2
                    iy=iy+1;
                elseif i5==3
                    iz=iz+1;
                end
                
                if i6==1
                    ix=ix+1;
                elseif i6==2
                    iy=iy+1;
                elseif i6==3
                    iz=iz+1;
                end
                
                if (ix==x)&(iy==y)&(iz==z)
                    counter=counter+1;
                end
            end
        end
    end
   end
   end
end