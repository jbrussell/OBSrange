function [xmatbs,ymatbs,zmatbs,amatbs,bmatbs,indxs] = bootstrap(x,y,z,a,b,niter)
% [xmatbs,ymatbs,zmatbs,amatbs,bmatbs,indxs]=bootstrap(x,y,z,a,b,niter)
% 
% Applies balanced random resampling, where each data point occurs niter times in
% the bootstraping scheme.
%
% The 1st column of each matrix is the full unscrambled dataset
% 
% J. Russell, 2018


[ xmatbs, indxs, ~ ] = balanced_resamp( x,niter );
xmatbs = [x,xmatbs];
indxs = [[1:length(x)]',indxs];
[ ~, ~, ymat ] = balanced_resamp( y,niter+1 );
ymatbs = ymat(indxs);
[ ~, ~, zmat ] = balanced_resamp( z,niter+1 );
zmatbs = zmat(indxs);
for i = 1:size(a,1)
    [ ~, ~, amat(i).amat ] = balanced_resamp( a(i,:),niter+1 );
    amatbs(i).amatbs = amat(i).amat(indxs);
end
[ ~, ~, bmat ] = balanced_resamp( b,niter+1 );
bmatbs = bmat(indxs);
