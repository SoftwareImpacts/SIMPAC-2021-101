
%% Main run function for testing functions and inputs
% Set directory
% Add to path the data, dependancies and results folder
function data = main(data)
tic
switch data.method
    case "implicit"
        %% If implicit generate the field first
        if ~isfield(data,'tform')
            data.tform = affine3d(); % if no transform is defined
        end
        [data] = defineEquation(data);
        data.field = fieldGen(data.field);
        
        %% Then generate the surface mesh (FV and properties)
        r = 0.7; % target mesh reduction
        data.FV = reducepatch(isosurface(data.field.X,data.field.Y,data.field.Z,data.field.U,0),r);
        data.FVcap = reducepatch(isocaps(data.field.X,data.field.Y,data.field.Z,data.field.U,0,'below'),r);
        data = triangulationProperties(data); %calculate the per-face properties
        
    case "field"
        data.field = fieldGen(data);
        
    case "triangulation"
        app.MeshData.R = max(app.FVmesh.vertices,[],1);
        xq = app.MeshData.R(1)*(0:c.n)/c.n; yq = app.MeshData.R(2)*(0:c.n)/c.n;
        zq = app.MeshData.R(3)*(0:c.n)/c.n; [X,Y,Z] = meshgrid(xq,yq,zq);
        app.MeshData.xq = xq(1:c.n+1); app.MeshData.yq = yq(1:c.n+1); app.MeshData.zq = zq(1:c.n+1);
        app.FVcap = []; app.FVcap.faces=[]; app.FVcap.vertices = [];
        app.FV = app.FVmesh;
        temp.vertices = 1+(c.n)*((app.FVmesh.vertices)./(app.MeshData.R+1));
        temp.faces = app.FVmesh.faces;
        U = imfill(1.0*polygon2voxel(temp,c.n+1,'wrap'));
        F = fieldGen(app,X,Y,Z,bwdist(U)-bwdist(1-U));
end
data.metrics.CPUtime = toc;
end

