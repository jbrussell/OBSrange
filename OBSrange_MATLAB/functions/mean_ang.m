function [ mu ] = mean_ang( theta)
% [ mu ] = mean_ang( theta )
% 
%  Calculates circular mean of angles in radians - deals with wraparound problem
%  
%  INPUTS
%    theta: vector or a matrix of angles, in radians
%           (if matrix, means will be calculated column-wise)
%  OUTPUTs
%    mu:    mean angles, in radians
% 
% Z. Eilon, 2011

if isrow(theta)
    theta = theta'; % flip if needed
end

n = size(theta,1);
m = size(theta,2);
mu   = zeros(1,m);

for im = 1:m
C=sum(cos(theta(:,im)));
S=sum(sin(theta(:,im)));
if C==0; error('angles sum to nothing'); end
if S==0; mu(im)=0; end
% fprintf('sample mean resultant length, mrlR, is %.3f \n',mrlR);
if S>0 && C>0
    mu(im)=atan(S/C); % in radians
elseif S<0 && C>0
    mu(im)=atan(S/C) + 2*pi; % in radians
elseif C<0
    mu(im)=atan(S/C)+pi; % in radians
end
end
end

