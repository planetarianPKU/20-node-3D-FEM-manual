function [shapeQ4,dhdrQ4,dhdsQ4,dhdtQ4]=shapefunctions_3d_20nodes(xi,eta,zeta)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric four-node Quadilateral shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapeQ4,dhdrQ4,dhdsQ4]=shapefunctions(xi,eta)  
%
%  Variable Description:
%     shapeQ4 - shape functions for four-node element
%     dhdrQ4 - derivatives of the shape functions w.r.t. r
%     dhdsQ4 - derivatives of the shape functions w.r.t. s
%     xi - r coordinate value of the selected point   
%     eta - s coordinate value of the selected point
%
%  Notes:
%     1st node at (-1,-1), 2nd node at (1,-1)
%     3rd node at (1,1), 4th node at (-1,1)
%------------------------------------------------------------------------

% shape functions

 shapeQ4(1)=0.125*(1+xi)*(1+eta)*(1+zeta)*(xi+eta+zeta-2);
 shapeQ4(2)=0.125*(1-xi)*(1+eta)*(1+zeta)*(-xi+eta+zeta-2);
 shapeQ4(3)=0.125*(1-xi)*(1-eta)*(1+zeta)*(-xi-eta+zeta-2);
 shapeQ4(4)=0.125*(1+xi)*(1-eta)*(1+zeta)*(xi-eta+zeta-2);
 shapeQ4(5)=0.125*(1+xi)*(1+eta)*(1-zeta)*(xi+eta-zeta-2);
 shapeQ4(6)=0.125*(1-xi)*(1+eta)*(1-zeta)*(-xi+eta-zeta-2);
 shapeQ4(7)=0.125*(1-xi)*(1-eta)*(1-zeta)*(-xi-eta-zeta-2);
 shapeQ4(8)=0.125*(1+xi)*(1-eta)*(1-zeta)*(xi-eta-zeta-2);
 shapeQ4(9)=0.25*(1-xi^2)*(1+eta)*(1+zeta);
 shapeQ4(10)=0.25*(1-eta^2)*(1-xi)*(1+zeta);
 shapeQ4(11)=0.25*(1-xi^2)*(1-eta)*(1+zeta);
 shapeQ4(12)=0.25*(1-eta^2)*(1+xi)*(1+zeta);
 shapeQ4(13)=0.25*(1-xi^2)*(1+eta)*(1-zeta);
 shapeQ4(14)=0.25*(1-eta^2)*(1-xi)*(1-zeta);
 shapeQ4(15)=0.25*(1-xi^2)*(1-eta)*(1-zeta);
 shapeQ4(16)=0.25*(1-eta^2)*(1+xi)*(1-zeta);
 shapeQ4(17)=0.25*(1-zeta^2)*(1+xi)*(1+eta);
 shapeQ4(18)=0.25*(1-zeta^2)*(1-xi)*(1+eta);
 shapeQ4(19)=0.25*(1-zeta^2)*(1-xi)*(1-eta);
 shapeQ4(20)=0.25*(1-zeta^2)*(1+xi)*(1-eta);
% derivatives
%xi
 dhdrQ4(1)=(xi/8 + 1/8)*(eta + 1)*(zeta + 1) + ((eta + 1)*(zeta + 1)*(eta + xi + zeta - 2))/8;
 dhdrQ4(2)=(xi/8 - 1/8)*(eta + 1)*(zeta + 1) - ((eta + 1)*(zeta + 1)*(eta - xi + zeta - 2))/8;
 dhdrQ4(3)=- ((eta - 1)*(zeta + 1)*(eta + xi - zeta + 2))/8 - (xi/8 - 1/8)*(eta - 1)*(zeta + 1);
 dhdrQ4(4)=((eta - 1)*(zeta + 1)*(eta - xi - zeta + 2))/8 - (xi/8 + 1/8)*(eta - 1)*(zeta + 1);
 dhdrQ4(5)=- ((eta + 1)*(zeta - 1)*(eta + xi - zeta - 2))/8 - (xi/8 + 1/8)*(eta + 1)*(zeta - 1);
 dhdrQ4(6)=- ((eta + 1)*(zeta - 1)*(xi - eta + zeta + 2))/8 - (xi/8 - 1/8)*(eta + 1)*(zeta - 1);
 dhdrQ4(7)=(xi/8 - 1/8)*(eta - 1)*(zeta - 1) + ((eta - 1)*(zeta - 1)*(eta + xi + zeta + 2))/8;
 dhdrQ4(8)=(xi/8 + 1/8)*(eta - 1)*(zeta - 1) - ((eta - 1)*(zeta - 1)*(eta - xi + zeta + 2))/8;
 dhdrQ4(9)=-(xi*(eta + 1)*(zeta + 1))/2;
 dhdrQ4(10)=(eta^2/4 - 1/4)*(zeta + 1);
 dhdrQ4(11)=(xi*(eta - 1)*(zeta + 1))/2;
 dhdrQ4(12)=-(eta^2/4 - 1/4)*(zeta + 1);
 dhdrQ4(13)=(xi*(eta + 1)*(zeta - 1))/2;
 dhdrQ4(14)=-(eta^2/4 - 1/4)*(zeta - 1);
 dhdrQ4(15)=-(xi*(eta - 1)*(zeta - 1))/2;
 dhdrQ4(16)=(eta^2/4 - 1/4)*(zeta - 1);
 dhdrQ4(17)=-(zeta^2/4 - 1/4)*(eta + 1);
 dhdrQ4(18)=(zeta^2/4 - 1/4)*(eta + 1);
 dhdrQ4(19)=-(zeta^2/4 - 1/4)*(eta - 1);
 dhdrQ4(20)=(zeta^2/4 - 1/4)*(eta - 1);
%eta

 dhdsQ4=[  (conj(zeta) + 1)*(conj(xi)/8 + 1/8)*(conj(eta) + conj(xi) + conj(zeta) - 2) + (conj(eta) + 1)*(conj(zeta) + 1)*(conj(xi)/8 + 1/8)
- (conj(eta) + 1)*(conj(zeta) + 1)*(conj(xi)/8 - 1/8) - (conj(zeta) + 1)*(conj(xi)/8 - 1/8)*(conj(eta) - conj(xi) + conj(zeta) - 2)
- (conj(eta) - 1)*(conj(zeta) + 1)*(conj(xi)/8 - 1/8) - (conj(zeta) + 1)*(conj(xi)/8 - 1/8)*(conj(eta) + conj(xi) - conj(zeta) + 2)
  (conj(eta) - 1)*(conj(zeta) + 1)*(conj(xi)/8 + 1/8) + (conj(zeta) + 1)*(conj(xi)/8 + 1/8)*(conj(eta) - conj(xi) - conj(zeta) + 2)
- (conj(eta) + 1)*(conj(zeta) - 1)*(conj(xi)/8 + 1/8) - (conj(zeta) - 1)*(conj(xi)/8 + 1/8)*(conj(eta) + conj(xi) - conj(zeta) - 2)
  (conj(eta) + 1)*(conj(zeta) - 1)*(conj(xi)/8 - 1/8) - (conj(zeta) - 1)*(conj(xi)/8 - 1/8)*(conj(xi) - conj(eta) + conj(zeta) + 2)
  (conj(zeta) - 1)*(conj(xi)/8 - 1/8)*(conj(eta) + conj(xi) + conj(zeta) + 2) + (conj(eta) - 1)*(conj(zeta) - 1)*(conj(xi)/8 - 1/8)
- (conj(eta) - 1)*(conj(zeta) - 1)*(conj(xi)/8 + 1/8) - (conj(zeta) - 1)*(conj(xi)/8 + 1/8)*(conj(eta) - conj(xi) + conj(zeta) + 2)
                                                                                             -(conj(zeta) + 1)*(conj(xi)^2/4 - 1/4)
                                                                                      (conj(eta)*(conj(xi) - 1)*(conj(zeta) + 1))/2
                                                                                              (conj(zeta) + 1)*(conj(xi)^2/4 - 1/4)
                                                                                     -(conj(eta)*(conj(xi) + 1)*(conj(zeta) + 1))/2
                                                                                              (conj(zeta) - 1)*(conj(xi)^2/4 - 1/4)
                                                                                     -(conj(eta)*(conj(xi) - 1)*(conj(zeta) - 1))/2
                                                                                             -(conj(zeta) - 1)*(conj(xi)^2/4 - 1/4)
                                                                                      (conj(eta)*(conj(xi) + 1)*(conj(zeta) - 1))/2
                                                                                             -(conj(xi) + 1)*(conj(zeta)^2/4 - 1/4)
                                                                                              (conj(xi) - 1)*(conj(zeta)^2/4 - 1/4)
                                                                                             -(conj(xi) - 1)*(conj(zeta)^2/4 - 1/4)
                                                                                              (conj(xi) + 1)*(conj(zeta)^2/4 - 1/4)];
dhdsQ4=dhdsQ4';                                                                                      
 %zeta
 
dhdtQ4=[  (conj(eta) + 1)*(conj(xi)/8 + 1/8)*(conj(eta) + conj(xi) + conj(zeta) - 2) + (conj(eta) + 1)*(conj(zeta) + 1)*(conj(xi)/8 + 1/8)
- (conj(eta) + 1)*(conj(zeta) + 1)*(conj(xi)/8 - 1/8) - (conj(eta) + 1)*(conj(xi)/8 - 1/8)*(conj(eta) - conj(xi) + conj(zeta) - 2)
  (conj(eta) - 1)*(conj(zeta) + 1)*(conj(xi)/8 - 1/8) - (conj(eta) - 1)*(conj(xi)/8 - 1/8)*(conj(eta) + conj(xi) - conj(zeta) + 2)
  (conj(eta) - 1)*(conj(xi)/8 + 1/8)*(conj(eta) - conj(xi) - conj(zeta) + 2) - (conj(eta) - 1)*(conj(zeta) + 1)*(conj(xi)/8 + 1/8)
  (conj(eta) + 1)*(conj(zeta) - 1)*(conj(xi)/8 + 1/8) - (conj(eta) + 1)*(conj(xi)/8 + 1/8)*(conj(eta) + conj(xi) - conj(zeta) - 2)
- (conj(eta) + 1)*(conj(zeta) - 1)*(conj(xi)/8 - 1/8) - (conj(eta) + 1)*(conj(xi)/8 - 1/8)*(conj(xi) - conj(eta) + conj(zeta) + 2)
  (conj(eta) - 1)*(conj(xi)/8 - 1/8)*(conj(eta) + conj(xi) + conj(zeta) + 2) + (conj(eta) - 1)*(conj(zeta) - 1)*(conj(xi)/8 - 1/8)
- (conj(eta) - 1)*(conj(zeta) - 1)*(conj(xi)/8 + 1/8) - (conj(eta) - 1)*(conj(xi)/8 + 1/8)*(conj(eta) - conj(xi) + conj(zeta) + 2)
                                                                                             -(conj(eta) + 1)*(conj(xi)^2/4 - 1/4)
                                                                                              (conj(xi) - 1)*(conj(eta)^2/4 - 1/4)
                                                                                              (conj(eta) - 1)*(conj(xi)^2/4 - 1/4)
                                                                                             -(conj(xi) + 1)*(conj(eta)^2/4 - 1/4)
                                                                                              (conj(eta) + 1)*(conj(xi)^2/4 - 1/4)
                                                                                             -(conj(xi) - 1)*(conj(eta)^2/4 - 1/4)
                                                                                             -(conj(eta) - 1)*(conj(xi)^2/4 - 1/4)
                                                                                              (conj(xi) + 1)*(conj(eta)^2/4 - 1/4)
                                                                                    -(conj(zeta)*(conj(eta) + 1)*(conj(xi) + 1))/2
                                                                                     (conj(zeta)*(conj(eta) + 1)*(conj(xi) - 1))/2
                                                                                    -(conj(zeta)*(conj(eta) - 1)*(conj(xi) - 1))/2
                                                                                     (conj(zeta)*(conj(eta) - 1)*(conj(xi) + 1))/2];
dhdtQ4=dhdtQ4';                                                                                      
 