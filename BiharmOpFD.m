function [BiharmOp] = BiharmOpFD(pointsx,pointsy,lengthx,lengthy)
%BiharmOpFD produces a finite difference biharmonic operator


%Calculates the grid spacings dx and dy where we take into account the
%gridpoints used for the boundary conditions.
dx=lengthx/(pointsx+1);
dy=lengthy/(pointsy+1);

%Initialises matrix blocks
R=sparse([],[],[],pointsx,pointsx,3*pointsx);
T=sparse([],[],[],pointsx,pointsx,3*pointsx);

%Indexes matrix blocks
for i=1:pointsx
    R(i,i)=6/dx^4+8/(dx^2*dy^2)+6/dy^4;
end
for i=1:pointsx-1
    R(i,i+1)=-4/dx^4-4/(dx^2*dy^2);
    R(i+1,i)=-4/dx^4-4/(dx^2*dy^2);
end
for i=1:pointsx-2
    R(i,i+2)=1/dx^4;
    R(i+2,i)=1/dx^4;
end

for i=1:pointsx
    T(i,i)=-4/(dx^2*dy^2)-4/dy^4;
end
for i=1:pointsx-1
    T(i,i+1)=2/(dx^2*dy^2);
    T(i+1,i)=2/(dx^2*dy^2);
end



%Initialises the sparse matrix that will be the biharmonic operator as a 
%block matrix of R blocks
C=cell(1,pointsy);
C(:)={R};
BiharmOp=blkdiag(C{:});



%This creates a sparse matrix involving the T blocks which we add to the
%main operator
C=cell(1,pointsy-1);
C(:)={T};
Section = blkdiag(C{:});
Operator=[sparse([],[],[],pointsx,pointsx*pointsy,0) ; Section sparse([],[],[],pointsx*pointsy-pointsx,pointsx,0) ];

%Operator(pointsy+1:pointsy*pointsz,1:pointsy*pointsz-pointsy)=Section;




%This creates a sparse matrix involving the identity blocks which we add to the
%main operator
BiharmOp=BiharmOp+Operator+Operator';
G=sparse(1:pointsx,1:pointsx,ones(pointsx,1),pointsx,pointsx,pointsx);
C=cell(1,pointsy-2);
C(:)={G./dy^4};
Section = blkdiag(C{:});

Operator=[sparse([],[],[],2*pointsx,pointsx*pointsy,0) ; Section sparse([],[],[],pointsx*pointsy-2*pointsx,2*pointsx,0) ];

%Operator=sparse([],[],[],pointsy*pointsz,pointsy*pointsz,3*pointsy*pointsz);
% Operator(2*pointsy+1:pointsy*pointsz,1:pointsy*pointsz-2*pointsy)=Section./dz^4;
BiharmOp=BiharmOp+Operator+Operator';

clearvars -except BiharmOp


end

