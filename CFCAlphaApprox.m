function [ Alpha] = CFCAlphaApprox(pointsx,pointsy,H,OddEvenMode10,Mode,RefineAlpha, D)
%CFCAlphaApprox Solves the biharmonic equation for homogenous boundary
%conditions involving an eigenvalue alpha. The grid is spaced uniformly in
%both x (horizontal) and y (vertical) directions, but these spacings do not
%have to be equal to each other. The domain is lengthx units across and
%lengthy units tall, where the number of gridpoints for the unknowns is
%pointsx by pointsy. Two extra gridpoints exist in each row and column for
%the boundary conditions, and ghost rows and columns are used for the
%derivative boundary conditions. OddEvenMode10 should be selected as 1 if
%an odd mode is required, 0 if an even mode is required. Mode should be
%given the integer value of the mode you require from the selected odd or 
%even set. This is meant to be a shortened function which calculated
%approximations of the eigenvalues. If you have already an approximate
%value of Alpha enter this in as RefineAlpha, if not enter this as zero



lengthx = 2;


%Calculates the grid spacings dx and dy where we take into account the
%gridpoints used for the boundary conditions.
dx=lengthx/(pointsx+1);
dy=H/(pointsy+1);


%Calls the bihamonic operator
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




SparseIndexingForApproxDone=1

clearvars -except SecondOp BiharmOp pointsx pointsy dx dy RefineAlpha ...
                    Mode OddEvenMode10 H lengthx B D



                
                
%If we are not calculating a full set of eigenvalues, which initially calculate 
%one mode, and test whether or not it is the required mode, and then calculate more modes             
NumberOfEigVals=1;


%The following calculates either a full set of Alphas, or calculates from a
%given guess. The block then sorts into modes, and test whether or not the
%calculated alpha corresponds to the requested mode, and calculates more
%Alphas if not until the Alpha corresponding to the requested mode is
%obtained.

while exist('Alpha','var')==0

if RefineAlpha==0
     %Calculates pointsx Alphas
 [EigVecs, Eigvals]=eigs(B,BiharmOp,pointsx);
 [OddMode,~,OddVals,EvenMode,~,EvenVals] = EigVecsSort(diag(Eigvals),EigVecs,pointsx,dx,1);
else
    
    %Calculates only a limited number of eigenvalues
 [EigVecs, Eigvals]=eigs(B,BiharmOp,NumberOfEigVals,RefineAlpha);
 [OddMode,~ ,OddVals,EvenMode,~,EvenVals] = EigVecsSort(diag(Eigvals),EigVecs,pointsx,dx,0);
end
 


 

 if OddEvenMode10==1
     for k= 1:length(OddMode)
         if OddMode(k) == Mode
             Alpha = OddVals(k);
             break
         end
     end
 end
 
 
 if OddEvenMode10==0
     for k= 1:length(EvenMode)
         if EvenMode(k) == Mode
             Alpha = EvenVals(k);
             break
         end
     end
 end
 
         if exist('Alpha','var')==0
             %If we still havent found the requested mode we calculated
             %more modes for the next iteration of the loop
                NumberOfEigVals = NumberOfEigVals+2;
                    
                
               if OddEvenMode10==1
                   c=max(OddMode);
               else
                   c=max(EvenMode);
               end
               
               if isempty(c)
                   c=Mode;
               end
                %If we haven't found the requested mode which adjust the
                %guess for Alpha slightly to spped up the search process
                RefineAlpha = RefineAlpha +(Mode-c)*RefineAlpha/25;
        end
end


SortDone=1


clearvars -except Alpha
 
end
