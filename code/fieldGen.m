function F = fieldGen(F)
%FIELDGEN generates data from a field
%   Detailed explanation goes here
%Inputs:
%   - F a structure with the following fields
%       X,Y,Z,U - where U is the value of the field
F.solid = 1.0*(F.U<0);
[F.boundary,F.azimuth,elevation] = imgradient3(F.solid);
F.inclination = 90-elevation;
[F.shapeX,F.shapeY,F.shapeZ] = imgradientxyz(F.inclination);
[F.Nx, F.Ny, F.Nz] = imgradientxyz(F.solid);
F.mean = -1/2*divergence(F.X,F.Y,F.Z,F.Nx,F.Ny,F.Nz);
end
