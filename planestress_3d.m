% Plane stress analysis of plates
% Plane stress analysis of a thin plate under tension at its extremes
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Warning : On running this the workspace memory will be deleted. Save any
% if any data present before running the code !!
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%--------------------------------------------------------------------------
% Code written by : Siva Srinivas Kolukula                                |
%                   Senior Research Fellow                                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |
%                   India                                                 |
% E-mail : allwayzitzme@gmail.com                                         |
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Variable descriptions                                                                                                 
%   k = element matrix for stiffness
%   f = element vector
%   stiffness = system matrix                                             
%   force = system vector                                                 
%   displacement = system nodal displacement vector
%   coordinates = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   gausspoint = matrix containing sampling points for bending term
%   gaussweight = matrix containing weighting coefficients for bending term
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%   B = matrix for kinematic equation for plane stress
%   D = matrix for material property for plane stress

%----------------------------------------------------------------------------            

%--------------------------------------------------------------------------
%  input data 
%--------------------------------------------------------------------------
clear 
clc
disp('Please wait Programme is under Run')
%--------------------------------------------------------------------------
% Input data for nodal coordinate values
%--------------------------------------------------------------------------
load  coordinates_3d.mat ;
%--------------------------------------------------------------------------
% Input data for nodal connectivity for each element
%--------------------------------------------------------------------------
load nodes_3d.mat ;
%
%coordinates=coordinates_3d;
nel = length(nodes) ;                  % number of elements
nnel=8;                                % number of nodes per element
ndof=3;                                % number of dofs per node (UX,UY,UZ)
nnode = length(coordinates) ;          % total number of nodes in system
sdof=nnode*ndof;                       % total system dofs  
edof=nnel*ndof;                        % degrees of freedom per element
% Units are in SI system
a = 1 ;                           % Length of the plate (along X-axes)
b = 1 ;                           % Length of the plate (along Y-axes)
c = 1 ;                           % Length of the plate (along Z-axes)
elementsalongX = 10 ;             % Number of elements along X-axes
elementsalongY = 10 ;             % Number of elements along Y-axes       
elementsalongZ = 10 ;             % Number of elements along Z-axes       

%
%PlotMesh(coordinates,nodes)

E = 2.1*10^11;                      % Youngs modulus
nu = 0.3;                           % Poisson's ratio
%t = 0.0254;                         % plate thickness
rho = 7840. ;                       % Density of the plate

nglx = 2; ngly = 2;   nglz = 2;       % 2x2x2 Gauss-Legendre quadrature 
nglxyz=nglx*ngly*nglz;            % number of sampling points per element

%--------------------------------------------------------------------------
% Input data for boundary conditions
%--------------------------------------------------------------------------
% (0,0) and (1,0) are fixed

  bcdof = [ 3*find(coordinates(:,3)==0)-2; 3*find(coordinates(:,3)==0)-1 ;3*find(coordinates(:,3)==0)] ; %需要修改
  bcval = zeros(1,length(bcdof)) ; %需要修改
%--------------------------------------------------------------------------
%  initialization of matrices and vectors
%--------------------------------------------------------------------------

force = zeros(sdof,1);                % system force vector
stiffness = sparse(sdof,sdof);         % system stiffness matrix
displacement = zeros(sdof,1);         % system displacement vector
eldepl = zeros(edof,1) ;              % element displacement vector
index = zeros(edof,1);                % index vector
B = zeros(6,edof);              % kinematic matrix for bending
D = zeros(6,6);                 % constitutive matrix for bending
% stressGP = zeros(nglxyz,3) ;             % Matrix containing stress components
% strainGP = zeros(nglx,3) ;              % Matrix containing strain components
% stress = zeros(nglx,3) ;
% stressG = zeros(nel,3) ;
% allstressGP = zeros(nnel*nel,3) ;
% allstressEXP = zeros(nnel*nel,3) ;
% allstressPR = zeros(nnel*nel,3) ;
%--------------------------------------------------------------------------
% force vector
%--------------------------------------------------------------------------
P = 1e5 ;       % Load
rightedge = find(coordinates(:,3)==a);
rightdof = 3*rightedge;
force(rightdof) = -P*b/(elementsalongY+1) ;
%leftedge = find(coordinates(:,2)==0);
%leftdof = 3*leftedge-1;
%UX索引是3*edge-2,UY是3*edge-1,UZ是3*edge
%force(leftdof) = P*b/(elementsalongY+1) ;

%--------------------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%--------------------------------------------------------------------------
[Gausspoint,Gaussweight]=GaussQuadrature(nglx);     % sampling points & weights
D = (E*(1-nu)/((1+nu)*(1-2*nu)*(1-nu)))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;  nu nu 1-nu 0 0 0; 0 0 0 0.5*(1-2*nu) 0 0;  0 0 0 0 0.5*(1-2*nu) 0; 0 0 0 0 0 0.5*(1-2*nu)] ;    % Constituent Matrix for Plane stress

for iel=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xx(i)=coordinates(nd(i),1);    % extract x value of the node
yy(i)=coordinates(nd(i),2);    % extract y value of the node
zz(i)=coordinates(nd(i),3);    % extract z value of the node

end

k = zeros(edof,edof);        % initialization of stiffness matrix

%--------------------------------------------------------------------------
%  numerical integration for stiffness matrix
%--------------------------------------------------------------------------

for intx=1:nglx
xi = Gausspoint(intx,1);                  % sampling point in x-axis
wtx = Gaussweight(intx,1);               % weight in x-axis
for inty=1:ngly
eta = Gausspoint(inty,1);                  % sampling point in y-axis
wty = Gaussweight(inty,1) ;              % weight in y-axis
 for intz=1:nglz
    zeta = Gausspoint(intz,1);                  % sampling point in z-axis
wtz = Gaussweight(intz,1) ;              % weight in z-axis
[shape,dhdr,dhds,dhdt] = shapefunctions_3d(xi,eta,zeta);     % compute shape functions and
                                    % derivatives at sampling point

jacobian = Jacobian_3d(nnel,dhdr,dhds,dhdt,xx,yy,zz);  % compute Jacobian

detjacob=det(jacobian);                 % determinant of Jacobian
invjacob=inv(jacobian);                 % inverse of Jacobian matrix

[dhdx,dhdy,dhdz]=shapefunctionderivatives_3d(nnel,dhdr,dhds,dhdt,invjacob); % derivatives w.r.t.
                                               % physical coordinate

B=fekineps_3d(nnel,dhdx,dhdy,dhdz);          % kinematic matrix for stiffness

%--------------------------------------------------------------------------
%  compute element stiffness matrix
%--------------------------------------------------------------------------

k = k+B'*D*B*wtx*wty*detjacob;
 
end
end                      % end of numerical integration loop for bending term
end

index = elementdof(nd,nnel,ndof); % extract system dofs associated with element

stiffness = assemble(stiffness,k,index);    % assemble element stiffness matrices 

end

%--------------------------------------------------------------------------
%   apply boundary conditions
%--------------------------------------------------------------------------

[stiffness,force] = constraints(stiffness,force,bcdof,bcval);

%--------------------------------------------------------------------------
% Solve the matrix equation 
%--------------------------------------------------------------------------
displacement = stiffness\force ;
UX = displacement(1:3:sdof) ;
UY = displacement(2:3:sdof) ;
UZ = displacement(3:3:sdof) ;
