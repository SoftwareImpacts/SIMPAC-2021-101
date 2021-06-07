function tform = transformObject(T)
%TRANSFORMOBJECT transformations for rotations and scale
%   Detailed explanation goes here
%Inputs:
%   - X object to be transformed
%   - T structure of the transform with fields,
%       scalings: a_x, a_y a_z
%       rotations: pitch, roll, yaw
%       offsets: o_x, o_y, o_z
%Outputs:
%   - Xout transformed object

%Handle missing fields and define default values
missing = ~isfield(T,{'a_x','a_y','a_z','pitch','roll','yaw','o_x','o_y','o_z'});
defaults = [1 1 1 0 0 0 1 1 1];
for i = 1:length(missing)
    if missing(i)
        T.(missing(i)) = defaults(i);
    end
end

rZ = [cosd(T.yaw) -sind(T.yaw) 0; sind(T.yaw) cosd(T.yaw) 0; 0 0 1];
rY = [cosd(T.roll) 0 -sind(T.roll); 0 1 0; sind(T.roll) 0 cosd(T.roll)];
rX = [1 0 0; 0 cosd(T.pitch) -sind(T.pitch); 0 sind(T.pitch) cosd(T.pitch)];
Axyz = [T.a_x 0 0; 0 T.a_y 0; 0 0 T.a_z];
trans = [T.o_x T.o_y T.o_z];
pad = [0; 0; 0];
rot = Axyz*rX*rY*rZ;
tform=affine3d([rot pad; trans 1]);
end