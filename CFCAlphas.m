function [ Alphas] = CFCAlphas(pointsx,pointsy,H, D)
%CFCAlphas calculates the first pointsx Alphas, and sorts them into even mode
%eigenvalues (column 1 of the output Alphas), and odd mode eigenvalues (column 2
%of the output Alphas). Note that only the 10 per cent or so smallest of the
%output Alphas are going to be particularly accurate.



lengthx = 2;



%Calculates the grid spacings dx and dy where we take into account the
%gridpoints used for the boundary conditions.
dx=lengthx/(pointsx+1);
dy=H/(pointsy+1);


%Calls a function which form the sparse Biharmonic operator matrix
[BiharmOp] = BiharmOpFD(pointsx,pointsy,lengthx,H);


  %Calculates the Q1 matrix as set out in the derivation
            Q1=sparse([],[],[],pointsx,pointsx,2);
            Q1(1,1)=1/dx^4;
            Q1(end,end)=1/dx^4;
            C=cell(1,pointsy);
            C(:)={Q1};
            SideBoundConds=blkdiag(C{:});

            
            
            %Calculates the Q2 matrix as set out in the derivation
            Q2=sparse([],[],[],pointsx,pointsx,3*pointsx);
            for i=1:pointsx
                Q2(i,i) = -2;
            end
            for i=1:pointsx-1
                Q2(i,i+1) = 1;
                Q2(i+1,i) = 1;
            end

            

            %Calculates the sparse matrix B as set out in the derivation
            B=sparse([],[],[],pointsx*pointsy,pointsx*pointsy,7*pointsx);
            B(1:pointsx,1:pointsx)=-2*Q2/(dx^2*dy^3);
            B(1:pointsx,pointsx+1:2*pointsx)=Q2/(2*dx^2*dy^3);
            
            
            
            
            %Alters the biharmonic operator so that the matrix BiharmOp is 
            %the matrix C as set out in the derivation after the following
            %lines of code
            BiharmOp(1:pointsx,1:pointsx)=BiharmOp(1:pointsx,1:pointsx)-eye(pointsx)/dy^4;
            BiharmOp(end-pointsx+1:end,end-pointsx+1:end)...
                =BiharmOp(end-pointsx+1:end,end-pointsx+1:end)+eye(pointsx)/dy^4;
            BiharmOp=BiharmOp+SideBoundConds;
            
            

            %Adds extra terms to the B matrix for non zero diffusion
            %coefficient D
            B(1:pointsx,1:pointsx) = B(1:pointsx,1:pointsx) -5*D*Q2/(dx^2*dy^4);
            B(1:pointsx,pointsx+1:2*pointsx) = B(1:pointsx,pointsx+1:2*pointsx)+ 4*D*Q2/(dx^2*dy^4);
            B(1:pointsx,2*pointsx+1:3*pointsx) =B(1:pointsx,2*pointsx+1:3*pointsx) -D*Q2/(dx^2*dy^4);




SparseIndexingDone=1

clearvars -except SecondOp BiharmOp pointsx pointsy dx dy  ...
                    H lengthx B D



%The following calculates a full set of Alphas.
 [EigVecs, EigVals]=eigs(B,BiharmOp,pointsx);
 

    
    %The following block of code sorts the eigenvalues into those
    %corresponding to even and odd modes of psi
 [~,~,OddVals,~,~,EvenVals] = EigVecsSort(diag(EigVals),EigVecs,pointsx,dx, 1);
 a = length(OddVals);
 b = length(EvenVals);
 Alphas = zeros(max(a,b),2);
 Alphas(1:b,1) = EvenVals;
 Alphas(1:a,2) = OddVals;
 
end
