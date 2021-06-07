function data = triangulationProperties(data)
%TRIANGULATIONPROPERTIES updates a triangulation(FV) object with more
%properties per-face
% Inputs:
%   - data structure with fields FV.faces and FV.vertices
% Outputs:
%   - augmented data structure
try
    FV = data.FV;
    Fin = data.FV.faces;
    Vin = data.FV.vertices;

% Compute properties for volume, normal etc...
d13= [(Vin(Fin(:,2),1)-Vin(Fin(:,3),1)), (Vin(Fin(:,2),2)-Vin(Fin(:,3),2)), (Vin(Fin(:,2),3)-Vin(Fin(:,3),3))];
d12= [(Vin(Fin(:,1),1)-Vin(Fin(:,2),1)), (Vin(Fin(:,1),2)-Vin(Fin(:,2),2)), (Vin(Fin(:,1),3)-Vin(Fin(:,2),3))];
cr = cross(d13,d12,2); % Edge 1-3x1-2 Cross-product
crNorm = sqrt(cr(:,1).^2+cr(:,2).^2+cr(:,3).^2);
FV.Area = 0.5*sqrt(cr(:,1).^2+cr(:,2).^2+cr(:,3).^2);% Area of each triangle
normal = -cr./crNorm;% Unit normal for each triangle
zMean = (Vin(Fin(:,1),3)+Vin(Fin(:,2),3)+Vin(Fin(:,3),3))/3;
FV.Nx = normal(:,1); FV.Ny = normal(:,2);  FV.Nz = normal(:,3);
FV.Volume = -FV.Area.*zMean.*FV.Nz; % Volume flux/contribution of each triangle


[FV.mean,FV.gaussian,FV.minK,FV.maxK] = patchcurvature(Fin,Vin,false); % Call patch curvature
FV.inclination = real(acosd(-FV.Nz));
FV.Azimuth = real(atan2d(FV.Ny,FV.Nx));

% L-PBF Surface Quality 0 is bad 1 is good.
FV.LPBFQuality = 0.5+0.5*tanh((FV.inclination-30)/20); %Low Elevation Error M2(Nz) [0 1]

% Manufacturability Paper Positional Error Function
c1 = 113.8; c2 = 260; c3 = 0.08119; c4 = 61.38; c5 = 0.001806;
FV.PositionalError = c1+c2./(1+exp(c3.*(FV.inclination+c4)+c5.*abs(FV.gaussian)));

% Surface Roughness Function.
c1 = 2.204; c2 = 63.76; c3 = 0.06736; c4 = 7.843;
FV.Ra = c1+c2./(1+exp(c3.*(FV.inclination+c4)));

%Calculate metrics based on the triangulation
data.metrics.errorFlag = 0;
data.metrics.surfaceArea = sum(FV.Area,'all','omitnan');
data.metrics.rmsMC = sum(sqrt(FV.mean.^2),'all','omitnan')/data.metrics.surfaceArea;
data.metrics.volume = abs(sum(FV.Volume,'all','omitnan'));
data.metrics.areaBelow30deg = sum(FV.Area(FV.inclination<30),'all','omitnan')/data.metrics.surfaceArea;
data.metrics.LpbfSimple = sum(FV.LPBFQuality.*FV.Area,'all','omitnan')/data.metrics.surfaceArea;
data.metrics.LpbfError_Mean = sum(FV.PositionalError.*FV.Area,'all','omitnan')/data.metrics.surfaceArea;
data.metrics.LpbfError_Max = prctile(FV.PositionalError,0.98,"all");
data.metrics.LpbfRa_Max = prctile(FV.Ra,0.98,"all");

data.FV = FV;

catch
    data.metrics.errorFlag = 1;
end
end

