function [psi_yy,psi_xy,psi_xx] = CalculateSecondDerivativesOfPsi(usolgrid, pointsx,pointsy,dx,dy)
%CalculateSecondDerivativesOfPsi Calculates the second derivatives of psi

psi_xy = zeros(pointsy+2,pointsx+2);

psi_yy =  ([usolgrid(2:end,:);zeros(1,pointsx+2)] - 2*usolgrid + [zeros(1,pointsx+2);usolgrid(1:end-1,:)])/(dy^2);
psi_yy(end,:) = (2*usolgrid(end,:)-5*usolgrid(end-1,:)+4*usolgrid(end-2,:)-usolgrid(end-3,:))/dy^2;
psi_yy(1,:) = (2*usolgrid(1,:)-5*usolgrid(2,:)+4*usolgrid(3,:)-usolgrid(4,:))/dy^2;

for i=2:pointsy+1
    for j=2:pointsx+1
        psi_xy(i,j) = (usolgrid(i+1,j+1) - usolgrid(i-1,j+1) - usolgrid(i+1,j-1) + usolgrid(i-1,j-1))/(4*dx*dy);
    end
end

for i=2:pointsy+1
    psi_xy(i,1) = (-3*(usolgrid(i+1,1)-usolgrid(i-1,1))/2+2*(usolgrid(i+1,2)-usolgrid(i-1,2))-0.5*(usolgrid(i+1,3)-usolgrid(i-1,3)))/(2*dy*dx);
    psi_xy(i,end) = (3*(usolgrid(i+1,end)-usolgrid(i-1,end))/2-2*(usolgrid(i+1,end-1)-usolgrid(i-1,end-1))+0.5*(usolgrid(i+1,end-2)-usolgrid(i-1,end-2)))/(2*dy*dx);
end

for i=2:pointsx+1
    psi_xy(1,i) = (-3*(usolgrid(1,i+1)-usolgrid(1,i-1))/2+2*(usolgrid(2,i+1)-usolgrid(2,i-1))-0.5*(usolgrid(3,i+1)-usolgrid(3,i-1)))/(2*dy*dx);
    psi_xy(end,i) = (3*(usolgrid(end,i+1)-usolgrid(end,i-1))/2-2*(usolgrid(end-1,i+1)-usolgrid(end-1,i-1))+0.5*(usolgrid(end-2,i+1)-usolgrid(end-2,i-1)))/(2*dy*dx);
end

psi_xy(end,1) =  - usolgrid(end-1,2)/(dx*dy);
psi_xy(end,end) =   usolgrid(end-1,end-1)/(dx*dy);
psi_xy(1,end) = -  usolgrid(2,end-1)/(dx*dy);
psi_xy(1,1) =   usolgrid(2,2)/(dx*dy);

psi_xx =  ([usolgrid(:,2:end),zeros(pointsy+2,1)] -2*usolgrid+ [zeros(pointsy+2,1),usolgrid(:,1:end-1)])/(dx^2);
psi_xx(:,end) = (2*usolgrid(:,end)-5*usolgrid(:,end-1)+4*usolgrid(:,end-2)-usolgrid(:,end-3))/(dx^2);
psi_xx(:,1) =(2*usolgrid(:,1)-5*usolgrid(:,2)+4*usolgrid(:,3)-usolgrid(:,4))/(dx^2);


end

