function [psi_yyy, psi_xyy, psi_xxy, psi_xxx] = CalculateThirdDerivativesOfPsi(usolgrid, pointsx,pointsy,dx,dy)
%CalculateThirdDerivativesOfPsi calculates the third derivatives of psi

[psi_zz,~,psi_yy] = CalculateSecondDerivativesOfPsi(usolgrid, pointsx,pointsy,dx,dy);

[psi_yyy,psi_xyy] = CalculateFirstDerivativesOfPsi(psi_zz, pointsx,pointsy,dx,dy);

[psi_xxy,psi_xxx] = CalculateFirstDerivativesOfPsi(psi_yy, pointsx,pointsy,dx,dy);



end

