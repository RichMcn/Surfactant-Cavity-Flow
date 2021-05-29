function [AlphaFromInt,Height,A_n2times4,SurfacePressure, SurfaceVelocity,...
                Vorticity,SurfactantProfile,StressProfile,Alpha,psi] = ...
                    CavityFlowCalculator(M,N,H,OddEvenMode10,Mode,D,AlphaApprox)
%CavityFlowCalculator Solves the biharmonic equation for the homogenous 
%boundary condition alpha psi_yy = - psi_xxy. The grid is spaced uniformly in
%both x (horizontal) and y (vertical) directions, but these spacings do not
%have to be equal to each other. The domain is 2 units across and
%H units tall, where the number of gridpoints for the solution is
%M (rows) by N (columns). Ghost rows and columns are used for the
%derivative boundary conditions. OddEvenMode10 should be selected as 1 if
%an odd mode is required, 0 if an even mode is required. Mode should be
%given the integer value of the mode you require from the selected odd or 
%even set. If you already know an approximate value for Alpha for the mode
%you require, enter it in as AlphaApprox. If not enter AlphasApprox as 0
%D sets the diffusion coefficient, which can be 0. For accurate computation
%of the first order film height deflection N should be an odd number.


%If you want to supress the plotting figures, turn this variable to 0
PlotFigures=1;


lengthx = 2;
pointsy = M-2;
pointsx = N-2;




%We run the code first for a coarse
%discretisation to approximate the eigenvalue, before running again
%specifically looking for the mode we want.
    if AlphaApprox ==0
        %We can use the lubrication theory limit if H is small
        if H<0.1
            if OddEvenMode10==1
                AlphaApprox = (Mode^2*pi^2*(H+4*D))/4;
            else
                AlphaApprox = ((2*Mode-1)^2*pi^2*(H+4*D))/16;
            end
        else
            %The following calculates an estimate for Alpha, which may take
            %a while for higher modes
                [AlphaApprox] = CFCAlphaApprox(150,150,H,OddEvenMode10,Mode,0,D);
                [AlphaApprox] = CFCAlphaApprox(250,250,H,OddEvenMode10,Mode,AlphaApprox,D);
                [AlphaApprox] = CFCAlphaApprox(400,400,H,OddEvenMode10,Mode,AlphaApprox,D);
                    if Mode>5
                          AlphaLoop = AlphaApprox;
                        [AlphaApprox] = CFCAlphaApprox(600,600,H,OddEvenMode10,Mode,AlphaApprox,D);
                            points=800;
                            if (Mode >10)&&(pointsx>600)
                                while abs(AlphaApprox-AlphaLoop)>0.1
                                    AlphaLoop = AlphaApprox;
                                    [AlphaApprox] = CFCAlphaApprox(points,points,H,OddEvenMode10,Mode,AlphaApprox,D);
                                    points=points+200;
                                        if points>=pointsx
                                            break;
                                        end
                                end
                            end
                    end
        end
    end
    AlphaApproxDone=AlphaApprox





%Calculates the grid spacings dx and dy where we take into account the
%gridpoints used for the boundary conditions.
dx=lengthx/(pointsx+1);
dy=H/(pointsy+1);







        %Calculates the Biharmonic finite difference operator using an external
        %function
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



%Clears unneeded variables to maximise available memory
clearvars  -except pointsx pointsy lengthx H ...
                    OddEvenMode10 Mode B ...
                        BiharmOp dx dy AlphaApprox D PlotFigures



 

%The following while loop is helpful for higher order modes and/or if the
%approximate value entered or calculated for Alpha is not accurate. It will
%loop until the requested mode is found.
Parameter=0;
OnLoop = 0;
while Parameter==0

    
    
        %Calculates the eigenvalue and eigenvector of
        %B\psi = Alpha C \psi for the selected mode using the approximate
        %value of Alpha calculated earlier.
        [usol, Alpha]=eigs(B,BiharmOp,1,AlphaApprox);
        
        CalculateEigenvalueAndVectorDone=1
        
      
        
%Calls a function which determines which mode has been calculated
[OddMode,OddVecs,OddVals,EvenMode,EvenVecs,EvenVals] = EigVecsSort(Alpha,usol,pointsx,dx, 0);




%The follwoing block of code determines whether or not the requested mode
%has been found, and if so breaks the loop.

if OddEvenMode10==1
        if OddMode==Mode
            Alpha = OddVals;
            usol= OddVecs;
            break;
        end
        ModeFound = OddMode;
else
        if EvenMode==Mode
             Alpha = EvenVals;
            usol= EvenVecs;
            break; 
        end
        ModeFound = EvenMode;
end

%If the requested mode has not been found, this block of code adjusts the
%value for Alpha and we recalculate again
if Parameter==0
    if (ModeFound<Mode)
        
          Alpha = Alpha+1;
    else 
            Alpha = Alpha-0.4;
    end
        
                if Alpha< 0 
                    Alpha = -Alpha;
                end
                
    CurrentGuessAlpha= Alpha
    OnLoop = OnLoop+1
    Loop = 'trying to find correct mode'
end
end


SortDone=1






%Creates a 2D grid for the solution
psi=zeros(pointsy+2,pointsx+2);






%The boundary conditions involve the two outermost columns and the
%bottom row and the top row as being all zeros. 
%This enters the rest of the solution from the vector usol to the grid
%usolgrid
for i=1:pointsy
    psi(pointsy+2-i,2:pointsx+1)=usol(1+(i-1)*pointsx:i*pointsx);
end



%Vectors for the two spatial dimensions
x=linspace(-lengthx/2,lengthx/2,pointsx+2);
y=linspace(-H,0,pointsy+2);





%The following block of code finds the normalisation factor NormScale

[psi_yy,~,~] = CalculateSecondDerivativesOfPsi(psi, pointsx,pointsy,dx,dy);
%This calculates the concentration profile
SurfactantProfile=zeros(pointsx+2,1);
for i=2:pointsx+2
SurfactantProfile(i)=trapz(dx,-psi_yy(end,1:i));
end
if SurfactantProfile (3)<SurfactantProfile(2)
    psi=-psi;
    SurfactantProfile = -SurfactantProfile;
end
%The following normalises the whole solution such that the amplitude of the
%surfactant concentration profile is 1.
NormScale=max(SurfactantProfile)-min(SurfactantProfile);
%This normalises the solution
psi=psi./NormScale;






%Calls functions which calculate the first, second and third derivatives of
%psi
[psi_y,psi_x] = CalculateFirstDerivativesOfPsi(psi, pointsx,pointsy,dx,dy);
[psi_yy,psi_xy,psi_xx] = CalculateSecondDerivativesOfPsi(psi, pointsx,pointsy,dx,dy);
[psi_yyy, ~, psi_xxy, ~] = CalculateThirdDerivativesOfPsi(psi, pointsx,pointsy,dx,dy);







%Calculate the concentration profile, via Gamma_x = - ps_yy, there is an
%alternate method which is Gamma = psi_xy/Alpha
SurfactantProfile=zeros(pointsx+2,1);
for i=2:pointsx+2
SurfactantProfile(i)=trapz(dx,-psi_yy(end,1:i));
end
%Repositions the surfactant concentration profile such that it has mean
%zero
SurfactantProfile = SurfactantProfile - sum(SurfactantProfile)/(pointsx+1);



%The calculates the surface stress profile
StressProfile=-psi_yy(end,:);



%This calculates the stress constant 4A_n2
A_n2times4=0;
for k=1:round(pointsx/2)
    PreviousSC=A_n2times4;
    A_n2times4 =(SurfactantProfile(k+2)-SurfactantProfile(k))/(2*dx);
                if ((SurfactantProfile(k+3)-SurfactantProfile(k+1))/(2*dx)- A_n2times4)>(A_n2times4- PreviousSC)
                    break;
                end
end





%Calculates the surface velocity
SurfaceVelocity=psi_y(end,:);



%Thid calculates the Vorticity
Vorticity = -(psi_yy+psi_xx);





%This calculates the surface pressure profile
SurfacePressure = zeros(pointsx+2,1);
for i=2:pointsx+2
    SurfacePressure(i) = dx*trapz(psi_yyy(end,1:i)+psi_xxy(end,1:i));
end
SurfacePressure = SurfacePressure - SurfacePressure(round(pointsx/2)+1);




%The following block of code calculates the height deflection 

%This calculates P - 2psi_xy
ModifiedPressure = SurfacePressure-2*psi_xy(end,:)';
%The following integrates this quantity across the surface once
H_r = zeros(pointsx+2,1);
for i=round(pointsx/2)+1:pointsx+2
    H_r(i) = -dx*trapz(ModifiedPressure(round(pointsx/2)+1:i));
end
for i= round(pointsx/2)+1:-1:1
    H_r(i) = -dx*trapz(ModifiedPressure(round(pointsx/2)+1:-1:i));
end
%This integrates this quantity again
He = zeros(pointsx+2,1);
for i=round(pointsx/2)+1:pointsx+2
    He(i) = dx*trapz(H_r(round(pointsx/2)+1:i));
end
for i= round(pointsx/2)+1:-1:1
    He(i) = dx*trapz(H_r(round(pointsx/2)+1:-1:i));
end
%The following calculates the constants needed for the hieght deflection
%and then calculates the height deflection
if OddEvenMode10==1
P1 = -dx*trapz(He-He(1))/(dx*trapz(x.^2-1));
Height = He +P1*(x'.^2-1) - He(1);
else
P1 = dx*trapz(-He-He(1)*x')/(dx*trapz(x.^2-1));
Height = He +P1*(x'.^2-1) +He(1)*x';
end
    




%The following block of code calculates the integrals used in the energy
%formula given in the notes. The output EnergyErrorNorm therefore gives a
%global measure of numerical error

GammaInt=Alpha*dx*trapz(SurfactantProfile.^2);
DiffusionInt = D*trapz(psi_yy(end,:).^2)*dy;
VorticityInt = dx*trapz(dy*trapz(Vorticity.^2,1));
EnergyErrorNorm = abs(GammaInt - DiffusionInt - VorticityInt)/(GammaInt-DiffusionInt)



%This calculates the value of Alpha from the integrals
AlphaFromInt = (VorticityInt+DiffusionInt)/(dx*trapz(SurfactantProfile.^2));

 

 %Plots figures
 if PlotFigures
CFCPlotFigures(x,y,psi,SurfactantProfile,psi_y,psi_x,SurfaceVelocity,StressProfile,Vorticity,SurfacePressure,Height)
 end


end