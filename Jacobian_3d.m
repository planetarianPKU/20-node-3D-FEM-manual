function [jacobian]=Jacobian_3d(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord)

%------------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian for two-dimensional mapping
%
%  Synopsis:
%     [jacobian]=Jacobian(nnel,dhdr,dhds,xcoord,ycoord) 
%
%  Variable Description:
%     jacob2 - Jacobian for one-dimension
%     nnel - number of nodes per element   
%     dhdr - derivative of shape functions w.r.t. natural coordinate r
%     dhds - derivative of shape functions w.r.t. natural coordinate s
%     dhdt - derivative of shape functions w.r.t. natural coordinate t

%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
%     zcoord - z axis coordinate values of nodes
%------------------------------------------------------------------------

 jacobian=zeros(3,3);


 jacobian(1,1)=dhdr*xcoord';
 jacobian(1,2)=dhdr*ycoord';
 jacobian(1,3)=dhdr*zcoord';
 
 jacobian(2,1)=dhds*xcoord';
 jacobian(2,2)=dhds*ycoord';
 jacobian(2,3)=dhds*zcoord';

 jacobian(3,1)=dhdt*xcoord';
 jacobian(3,2)=dhdt*ycoord';
 jacobian(3,3)=dhdt*zcoord';