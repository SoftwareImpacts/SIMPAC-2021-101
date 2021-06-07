function [data]  = defineEquation(data)
%DefineEquation converts a string, for the type of tpms to be generated, and outputs a symbolic function
%Inputs: implicit, a structure with fields
%   equation - string describing the TPMS desired, or the symbolic (gyroid is default)
%   type     - string describing the phase
%   v1       - isovalue 1 (geometric parameter)
%   v2       - isovalue 2 (geometric parameter)
%Outputs: 
%   u        - symbolic equation in terms of x,y,z

% Default Options 
if ~isfield(data,'equation')
    data.equation = "Gyroid";
end
if ~isfield(data,'type')
    data.type = "network";
end
if ~isfield(data,'v1')
    data.v2 = 0;
end
if ~isfield(data,'v2')
    data.v2 = 0;
end
if ~isfield(data,'bulksize')
    data.bulkSize = diag(data.tform.T(1:3,1:3));
end
if ~isfield(data,'res')
    data.res = 26;
end


% Create symbolic variables

%Definitions of Standard Implicit Surfaces
switch data.equation
    case "Diamond"
        u = @(x,y,z) (sin((2*pi)*x).*sin((2*pi)*y).*sin((2*pi)*z)+...
            sin((2*pi)*x).*cos((2*pi)*y).*cos((2*pi)*z)+...
            cos((2*pi)*x).*sin((2*pi)*y).*cos((2*pi)*z)+...
            cos((2*pi)*x).*cos((2*pi)*y).*sin((2*pi)*z));
        
    case "Primative"
        u = @(x,y,z) cos((2*pi).*x)+cos((2*pi).*y)+cos((2*pi).*z);
    case "Sinusoidal"
        u = @(x,y,z) 0.4.*sin((2*pi).*x)+...
            0.4.*sin((2*pi).*y)-5*(z-0.5);
    case "Sphere"
        u = @(x,y,z) (x-0.5).^2+(y-0.5).^2+(z-0.5).^2-1/4;
    case "P-norm 10 Cube"
        u = @(x,y,z) -0.5^10+(x-0.5).^10+(y-0.5).^10+(z-0.5).^10;
    case "Taurus"
        R = 0.4; r = 0.1;
        u = @(x,y,z) (sqrt((x-0.5).^2+(y-0.5).^2)-R).^2+(z-0.5).^2-r.^2;
    otherwise
        u = @(x,y,z) cos((2*pi).*x).*sin((2*pi).*y)+ ...
            cos((2*pi).*y).*sin((2*pi).*z)+cos((2*pi).*z).*sin((2*pi).*x);
end
data.implicit.u = u;

field.xq = linspace(0,data.bulkSize(1),data.res);
field.yq = linspace(0,data.bulkSize(2),data.res);
field.zq = linspace(0,data.bulkSize(2),data.res);
[field.X,field.Y,field.Z] = meshgrid(field.xq,field.yq,field.zq);


% passive (css) transform
[X,Y,Z] = transformPointsInverse(data.tform,field.X,field.Y,field.Z);

%Evaluate isosurface
switch data.type
    case "bounded"
        field.U = (u(X,Y,Z)-data.v1).*(u(X,Y,Z)-data.v2);
    case "surface"
        field.U = (u(X,Y,Z)-data.v1).^2-data.v2.^2;
    case "network"
        field.U = u(X,Y,Z)-data.v1;
    otherwise
        field.U = u(X,Y,Z);
end
data.field = field;
end
