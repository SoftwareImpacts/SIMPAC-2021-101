function [Cmean,Cgaussian,k1,k2]=patchcurvature(F,V,usethird)
warning('off');
% This function calculates the principal curvature directions and values
% of a triangulated mesh. 
% [Cmean,Cgaussian,Dir1,Dir2,Lambda1,Lambda2,angle]=patchcurvature(F,V,usethird)
%
% inputs,
%   FV : A triangulated mesh (see Patch)
%   usethird : Use third order neighbour vertices for the curvature
%              fit, making it smoother but less local. true/ false (default)
% outputs,
%   Cmean : Mean Curvature
%   Cgaussian : Gaussian Curvature
%   Dir1 : XYZ Direction of first Principal component
%   Dir2 : XYZ Direction of second Principal component
%   Lambda1 : value of first Principal component
%   Lambda2 : value of second Principal component
% Function is written by D.Kroon University of Twente (August 2011)  
% Last Update, 15-1-2014 D.Kroon at Focal.

% Check inputs
if(nargin<3), usethird=false; end
% Number of vertices
nv=size(V,1);
% Calculate vertices normals
tri = triangulation(F,V);
N=tri.vertexNormal;

% Calculate Rotation matrices for the normals
M= zeros(3,3,nv);
Minv= zeros(3,3,nv);
for i=1:nv 
    [M(:,:,i),Minv(:,:,i)]=VectorRotationMatrix(N(i,:));
end

% Get neighbours of all vertices
Ne=vertex_neighbours(F,V);

% Loop through all vertices
k1=zeros(nv,1);
k2=zeros(nv,1);
Dir1=zeros(nv,3);
Dir2=zeros(nv,3);

for i=1:nv
   % Get first and second ring neighbours.
   if(~usethird)
       Nce=unique([Ne{Ne{i}}]);
   else
       % Get first, second and third ring neighbours
       Nce=unique([Ne{[Ne{Ne{i}}]}]);
   end
   
   Ve=V(Nce,:);

   % Rotate to make normal [-1 0 0]
   We=Ve*Minv(:,:,i);
   f=We(:,1); x=We(:,2); y=We(:,3); 
   
   % Fit patch
   % f(x,y) = ax^2 + by^2 + cxy + dx + ey + f
   FM=[x(:).^2 y(:).^2 x(:).*y(:) x(:) y(:) ones(numel(x),1)];
   abcdef=FM\f(:);
   a=abcdef(1); b=abcdef(2); c=abcdef(3);
   
   % Make Hessian matrix 
   % H =  [2*a c;c 2*b];
   Dxx = 2*a; Dxy=c; Dyy=2*b;
   
   [k1(i),k2(i),I1,I2]=eig2(Dxx,Dxy,Dyy);
   

   dir1=[0 I1(1) I1(2)]*M(:,:,i); 
   dir2=[0 I2(1) I2(2)]*M(:,:,i);
   Dir1(i,:)=dir1/sqrt(dir1(1)^2+dir1(2)^2+dir1(3)^2);
   Dir2(i,:)=dir2/sqrt(dir2(1)^2+dir2(2)^2+dir2(3)^2);
end

%Calculate face centre averages
k1 = min(max(k1,prctile(k1,5)),prctile(k1,95));
k2 = min(max(k2,prctile(k2,5)),prctile(k2,95));
k1 = (k1(F(:,1))+k1(F(:,2))+k1(F(:,3)))/3;
k2 = (k2(F(:,1))+k2(F(:,2))+k2(F(:,3)))/3;
Cmean=-(k1+k2)/2;
Cgaussian=k1.*k2;
end 


function [Lambda1,Lambda2,I1,I2]=eig2(Dxx,Dxy,Dyy)
% Compute the eigenvectors 
tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy; v2y = Dyy - Dxx + tmp;

% Normalize
mag = sqrt(v2x.^2 + v2y.^2); i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

% The eigenvectors are orthogonal
v1x = -v2y; v1y = v2x;

% Compute the eigenvalues
mu1 = (0.5*(Dxx + Dyy + tmp));
mu2 = (0.5*(Dxx + Dyy - tmp));

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
if(mu1<mu2)
    Lambda1=mu1;
    Lambda2=mu2;
    I2=[v1x v1y];
    I1=[v2x v2y];
else
    Lambda1=mu2;
    Lambda2=mu1;
    I2=[v2x v2y];
    I1=[v1x v1y];
end
end

function N=patchnormals(F,V)
% This function PATCHNORMALS calculates the normals of a triangulated
% mesh. PATCHNORMALS calls the patchnormal_double.c mex function which 
% first calculates the normals of all faces, and after that calculates 
% the vertice normals from the face normals weighted by the angles 
% of the faces.
[Nx,Ny,Nz]=patchnormals_double(double(F(:,1)),double(F(:,2)),double(F(:,3)),double(V(:,1)),double(V(:,2)),double(V(:,3)));
N=zeros(length(Nx),3);
N(:,1)=Nx; N(:,2)=Ny; N(:,3)=Nz;
end

function [M,Minv]=VectorRotationMatrix(v)
% [M,Minv]=VectorRotationMatrix(v,k)
v=(v(:)')/sqrt(sum(v.^2));
k=rand(1,3);
l = [k(2).*v(3)-k(3).*v(2), k(3).*v(1)-k(1).*v(3), k(1).*v(2)-k(2).*v(1)]; l=l/sqrt(sum(l.^2));
k = [l(2).*v(3)-l(3).*v(2), l(3).*v(1)-l(1).*v(3), l(1).*v(2)-l(2).*v(1)]; k=k/sqrt(sum(k.^2));
Minv=[v(:) l(:) k(:)];
M=inv(Minv);
end

function Ne=vertex_neighbours(F,V)
% Neighbourh cell array 
Ne=cell(1,size(V,1));

% Loop through all faces
for i=1:length(F)
    % Add the neighbors of each vertice of a face
    % to his neighbors list.
    Ne{F(i,1)}=[Ne{F(i,1)} [F(i,2) F(i,3)]];
    Ne{F(i,2)}=[Ne{F(i,2)} [F(i,3) F(i,1)]];
    Ne{F(i,3)}=[Ne{F(i,3)} [F(i,1) F(i,2)]];
end

% Loop through all neighbor arrays and sort them (Rotation same as faces)
for i=1:size(V,1)
    Pneighf=Ne{i};
    if(isempty(Pneighf))
        Pneig=[];
    else
        start=1;
        for index1=1:2:length(Pneighf)
            found=false;
            for index2=2:2:length(Pneighf)
                if(Pneighf(index1)==Pneighf(index2))
                    found=true; break
                end
            end
            if(~found)
                start=index1; break
            end
        end
        Pneig=[];
        Pneig(1)=Pneighf(start);
        Pneig(2)=Pneighf(start+1);
        
        % Add the neighbours with respect to original rotation
        for j=2+double(found):(length(Pneighf)/2)
            found = false;
            for index=1:2:length(Pneighf)
                if(Pneighf(index)==Pneig(end))
                    if(sum(Pneig==Pneighf(index+1))==0)
                        found =true;
                        Pneig=[Pneig Pneighf(index+1)];
                    end
                end
            end
            if(~found) % This only happens with weird edge vertices
                for index=1:2:length(Pneighf)
                    if(sum(Pneig==Pneighf(index))==0)
                        Pneig=[Pneig Pneighf(index)];
                        if(sum(Pneig==Pneighf(index+1))==0)
                            Pneig=[Pneig Pneighf(index+1)];
                        end
                    end
                end
            end
        end
        % Add forgotten neigbours
        if(length(Pneig)<length(Pneighf))
            for j=1:length(Pneighf)
                if(sum(Pneig==Pneighf(j))==0)
                    Pneig=[Pneig Pneighf(j)];
                end
            end
        end
    end
    Ne{i}=Pneig;
end
end
