function [OddMode,OddVecs,OddVals,EvenMode,EvenVecs,EvenVals] = EigVecsSort(EigVals,EigVecs,pointsx,dx, fullset10)
%EigVecsSort sorts eigenvectors and eigenvalues into odd and even modes,
%and numbers each mode. If we have a full set of eigenvalues (pointsx
%number of eigenvalues) calculated with eigs, then set fullset10 to 1,
%otherwise set it to 0



%Initialises the variables to hold the solutions
OddMode=[];
NewOddEigVals=[];
NewOddEigVecs=[];
EvenMode=[];
NewEvenEigVals=[];
NewEvenEigVecs=[];


        
        
          %First we sort the modes into odd and even
           for k=length(EigVals):-1:1
                    if abs(EigVecs(10*pointsx+10,k)-EigVecs(10*pointsx+pointsx-9,k))<  abs(EigVecs(10*pointsx+10,k)+EigVecs(10*pointsx+pointsx-9,k))
                           
                        NewEvenEigVals = [NewEvenEigVals; EigVals(k)];
                        NewEvenEigVecs = [NewEvenEigVecs, EigVecs(:,k)];
                        
                    else
                        NewOddEigVals = [NewOddEigVals; EigVals(k)];
                        NewOddEigVecs = [NewOddEigVecs, EigVecs(:,k)];
                        
                                
                    end
           end
            
            
            
           %If we have a full set of modes, the ordering of modes is
           %simplified
            if fullset10==1
                OddMode = 1:length(NewOddEigVals);
                 EvenMode = 1:length(NewEvenEigVals);
            end
            
            
            
    %If we don't have a full set of eigenvalues, we find out what number
    %eigenvalue we have by considering the number of turning points of the
    %surfactant concentration profile. The first block does odd
    %eigenvlaues, the second block does even eigenvlaues
 if fullset10==0           
            
    for k = 1:length(NewOddEigVals)
        
        SurfProf=zeros(pointsx+2,1);
           psizztopline=zeros(pointsx+2,1);

            for l=2:pointsx+1
                psizztopline(l) = -5*NewOddEigVecs(pointsx+l,k)+4*NewOddEigVecs(2*pointsx+l,k)-NewOddEigVecs(3*pointsx+l,k);
            end
        
        
        for l=2:pointsx+2
                SurfProf(l) = dx*trapz(psizztopline(1:l));
        end
        
        ZeroCross=0;
        
        if SurfProf(3)<SurfProf(2)
            SurfProf =  - SurfProf;
        end
        
        SurfProf = SurfProf -SurfProf(1)-0.2*max(SurfProf);
        
        for l=1:pointsx+1
            if sign(SurfProf(l+1))~= sign(SurfProf(l))
                
               ZeroCross = ZeroCross+1;
            end
        end
        
        
        if mod(ZeroCross,2)==0
                FoundMode=ZeroCross/2;
        else
            FoundMode = (ZeroCross-1)/2;
        end
        
        OddMode(k) = FoundMode;
        
    end






    
    for k = 1:length(NewEvenEigVals)
        
        SurfProf=zeros(pointsx+2,1);
           psizztopline=zeros(pointsx+2,1);

            for l=2:pointsx+1
                psizztopline(l) = -5*NewEvenEigVecs(pointsx+l,k)+4*NewEvenEigVecs(2*pointsx+l,k)-NewEvenEigVecs(3*pointsx+l,k);
            end
        
        
        for l=2:pointsx+2
                SurfProf(l) = dx*trapz(psizztopline(1:l));
        end
        
        ZeroCross=0;
        
         if SurfProf(3)<SurfProf(2)
            SurfProf =  - SurfProf;
        end
        
          SurfProf = SurfProf -SurfProf(1)-0.25*max(SurfProf);
        
        
        for l=1:pointsx+1
            if sign(SurfProf(l+1))~= sign(SurfProf(l))
                
               ZeroCross = ZeroCross+1;
            end
        end
        
        if mod(ZeroCross,2)==1
            FoundMode=(ZeroCross+1)/2;
        else
            FoundMode=(ZeroCross+2)/2;
        end
        
        EvenMode(k) = FoundMode;
        
       
    end
    
 end

 
 

OddVecs = NewOddEigVecs;
EvenVecs = NewEvenEigVecs;
OddVals = NewOddEigVals;
EvenVals = NewEvenEigVals;






end

