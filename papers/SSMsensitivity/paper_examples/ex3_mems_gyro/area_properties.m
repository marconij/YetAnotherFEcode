function [xG, yG, A, Ixx, Iyy] = area_properties(X, Y)
% Calculates the (2D) area properties with respect to its centroid
%
% OUTPUTS
% xG,yG: center of mass
% Ixx, Iyy: moment of inertia (use to compute lumped mass inertia as:
%           Jx = (density*thickness)*Ixx;
%           Jy = (density*thickness)*Iyy;
%           Jz = Jx + Jy;
% A: volume and area
%
% INPUTS
% X,Y: contour points
%
% Original function authors: Damiano Milani, Daniele Giannini (Politecnico 
% di Milano, 2016)
%
% Adapted by Jacopo Marconi (Politecnico di Milano, 2024)

if length(X)~=length(Y) || length(X)<3 
    error('Error in polygon coordinates')
end

X(end+1)=X(1);
Y(end+1)=Y(1);

% area
A=0;
for i=1:length(X)-1
    A = A + 0.5*(X(i)*Y(i+1)-X(i+1)*Y(i));
end

% if X and Y coords are ccw then A>0, if cw A<0
if A<0 
    X = flipud(X);
    Y = flipud(Y);
    A = abs(A);
end

% centroid coords
xG=0;
yG=0;
for i=1:length(X)-1
    xG = xG + 1.0/6/A *(X(i)+X(i+1))*(X(i)*Y(i+1)-X(i+1)*Y(i));
    yG = yG + 1.0/6/A *(Y(i)+Y(i+1))*(X(i)*Y(i+1)-X(i+1)*Y(i));
end

% ref system in the centroid
x=X-xG;
y=Y-yG;

Ixx=0;
Iyy=0;

for i=1:length(x)-1
    Ixx = Ixx + 1.0/12 * (y(i)^2+y(i)*y(i+1)+y(i+1)^2)*(x(i)*y(i+1)-x(i+1)*y(i));
    Iyy = Iyy + 1.0/12 * (x(i)^2+x(i)*x(i+1)+x(i+1)^2)*(x(i)*y(i+1)-x(i+1)*y(i));
end

end


