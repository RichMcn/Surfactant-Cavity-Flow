function [psi_y,psi_x] = CalculateFirstDerivativesOfPsi(usolgrid, pointsx,pointsy,dx,dy)
%CalculateFirstDerivativesOfPsi Calculates the first derivatives of psi

psi_y =  (  [usolgrid(2:end,:);zeros(1,pointsx+2)]-[zeros(1,pointsx+2);usolgrid(1:end-1,:)])/(2*dy);
psi_y(end,:) = (3*usolgrid(end,:)/2-2*usolgrid(end-1,:)+0.5*usolgrid(end-2,:))/dy;
psi_y(1,:) = (-3*usolgrid(1,:)/2+2*usolgrid(2,:)-0.5*usolgrid(3,:))/dy;


psi_x =  ([usolgrid(:,2:end),zeros(pointsy+2,1)] -   [zeros(pointsy+2,1),usolgrid(:,1:end-1)])/(2*dy);
psi_x(:,end) = (3*usolgrid(:,end)/2-2*usolgrid(:,end-1)+0.5*usolgrid(:,end-2))/dx;
psi_x(:,1) = (-3*usolgrid(:,1)/2+2*usolgrid(:,2)-0.5*usolgrid(:,3))/dx;


end

