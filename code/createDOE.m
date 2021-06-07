function [doeArray] = createDOE()
%CREATEDOE generates a DOE based so a single loop can be used for running
%   Detailed explanation goes here
%% Define variables to be investigated.


%Cell sizes - (in mm)
DOE.a_x = [10]; DOE.a_y = [10]; DOE.a_z = [10];

%Cell offsets - (in mm)
DOE.o_x = [0]; DOE.o_y = [0]; DOE.o_z = [0];

% Rotations - (in degrees)
DOE.pitch = [0]; DOE.roll = [0]; DOE.yaw = [0];

%Other inputs
DOE.equation = ["Gyroid"; "Diamond"; "Primative"];%  "Gyroid","Diamond"; "Primative"]; %equation type(s)
DOE.v1 = [0]; DOE.v2 = [0.5]; %isovalues
DOE.type = ["network"; "surface"]; % "network", "bounded", "surface"
DOE.res = 41;

% Convert into a table (design of experiments)
names = string(fieldnames(DOE));
levels = [];
for i = 1:length(names)
    levels = [levels, length(DOE.(names(i)))];
end

doeFF = fullfact(levels);
nel = size(doeFF,1);
doeArray = struct();

for i = 1:nel
    for j = 1:size(names,1)
        doeArray(i).(names(j)) = DOE.(names(j))(doeFF(i,j));
        doeArray(i).name = "DOE-Sample-"+int2str(1:nel);
    end
end
end
