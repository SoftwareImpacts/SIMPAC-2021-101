function [Vs,Fs,Cs] = poleFigure(N,Aw,n)
%% Function returns points for a heatmap mapped onto a unit sphere surface
% Inputs:   N  (n x 2or3)    - Orientation vector
%     if (n x 2) use spherical cs values in radians [azimuth,inclination]
%     if (n x 3) use cartesian values[Nx,Ny,Nz]
%           Aw  (n x 1)       - Data Weighting
%           n  (int)         - number of subdivisons for icosphere mesh
% Outputs:  Vs (nv x 3)      - verticies for patch
%           Fs (nf x 3)      - faces for patch
%           Cs (nv x 1)      - per vertex colour data
% Designed for outputs to be used in matlab patch graphic,
%   ex: patch('Faces',Fs,'Vertices',Vs,'CData',Cs,'FaceColor','interp')
%
%   Credit: Alistair Jones, 12/11/2020, RMIT University
if size(N,2)==2 %Convert from spherical cs (radians)
    N = [sin(N(2,:)).*cos(N(1,:)), sin(N(2,:)).*sin(N(1,:)), cos(N(2,:))];
end

totalA = sum(Aw,'all','omitnan'); %Total per face weighting
[Vs,Fs] = icosphere(n); %Generate an icosphere with n subdivisions
V1 = Vs(Fs(:,1),:); V2 = Vs(Fs(:,2),:); V3 = Vs(Fs(:,3),:); %Face vertices
Fc = (V1+V2+V3)/3; %Find face centroid
Fn = sqrt(Fc(:,1).^2+Fc(:,2).^2+Fc(:,3).^2); %Find length of face center
area = 2*pi*(1-Fn); % Approximate curved cap area of each triangle
ids = dsearchn(Fc./Fn,N); %Find which node of the icosphere is closest to the unit normal vector
Cs = 4*pi*accumarray(ids,Aw,[length(Fc) 1],@sum,0)./(area.*totalA); %Complete cumulative sum
end

function [vv,ff] = icosphere(n)
%ICOSPHERE Generate icosphere.
% Create a unit geodesic sphere created by subdividing a regular
% icosahedron with normalised vertices.
%
%   [V,F] = ICOSPHERE(N) generates to matrices containing vertex and face
%   data so that patch('Faces',F,'Vertices',V) produces a unit icosphere
%   with N subdivisions.
%   Based on C# code by Andres Kahler
%   http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
%   Adapted by Alistair Jones, 12/11/2020

% generate regular unit icosahedron (20 faced polyhedron)
[v,f] = icosahedron(); % size(v) = [12,3]; size(f) = [20,3];
% recursively subdivide triangle faces
for gen = 1:n
    f_ = zeros(size(f,1)*4,3);
    for i = 1:size(f,1) % for each triangle
        tri = f(i,:);
        % calculate mid points (add new points to v)
        [a,v] = getMidPoint(tri(1),tri(2),v);
        [b,v] = getMidPoint(tri(2),tri(3),v);
        [c,v] = getMidPoint(tri(3),tri(1),v);
        % generate new subdivision triangles
        nfc = [tri(1),a,c;
               tri(2),b,a;
               tri(3),c,b;
                    a,b,c];
        % replace triangle with subdivision
        idx = 4*(i-1)+1:4*i;
        f_(idx,:) = nfc;
    end
    f = f_; % update 
end
% remove duplicate vertices
[vv,~,ix] = unique(v,'rows');
% reassign faces to trimmed vertex list and remove any duplicate faces
ff = unique(ix(f),'rows');
end

function [i,v] = getMidPoint(t1,t2,v)
%GETMIDPOINT calculates point between two vertices
%   Calculate new vertex in sub-division and normalise to unit length
%   then find or add it to v and return index
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
% get vertice positions
p1 = v(t1,:); p2 = v(t2,:);
% calculate mid point (on unit sphere)
pm = (p1 + p2) ./ 2;
pm = pm./norm(pm);
% add to vertices list, return index
i = size(v,1) + 1;
v = [v;pm];
end

function [v,f] = icosahedron()
%ICOSAHEDRON creates unit regular icosahedron
%   Returns 12 vertex and 20 face values.
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
      1, t, 0; % v2
     -1,-t, 0; % v3
      1,-t, 0; % v4
      0,-1, t; % v5
      0, 1, t; % v6
      0,-1,-t; % v7
      0, 1,-t; % v8
      t, 0,-1; % v9
      t, 0, 1; % v10
     -t, 0,-1; % v11
     -t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));
% create faces
f = [ 1,12, 6; % f1
      1, 6, 2; % f2
      1, 2, 8; % f3
      1, 8,11; % f4
      1,11,12; % f5
      2, 6,10; % f6
      6,12, 5; % f7
     12,11, 3; % f8
     11, 8, 7; % f9
      8, 2, 9; % f10
      4,10, 5; % f11
      4, 5, 3; % f12
      4, 3, 7; % f13
      4, 7, 9; % f14
      4, 9,10; % f15
      5,10, 6; % f16
      3, 5,12; % f17
      7, 3,11; % f18
      9, 7, 8; % f19
     10, 9, 2];% f20
end