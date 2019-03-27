function [ dataout, indxs, datamat ] = balanced_resamp( data,niter )
%[ indxs ] = balanced_resampling( data,niter )
%   
%   Function to get indices for balanced resampling of data for
%   bootstrapping
%
% INPUTS:
%   data    = data vector that will be randomly resampled
%   niter   = number of iterations in bootstrap
%   
% OUTPUTS:
%  dataout  = Output data matrix with random resampling applied along the
%              columns. (Ndata x niter)
%   indxs   = Ndata x niter matrix, each column of which is a set of
%              indices from the original dataset. Each index appears a
%              total of precisely niter times in the whole matrix
%  datamat  = Data matrix before random resampling. Same size as indxs
%
% Written by Zach Eilon, 2014
% Adapted by Josh Russell 4/20/18

Ndata = length(data);

xx = [1:Ndata]';

XX = repmat(xx,niter,1);

ind = randperm(Ndata*niter);

YY = XX(ind);

indxs = reshape(YY,Ndata,niter);

datamat = reshape(repmat(data,1,niter),Ndata,niter);

dataout = datamat(indxs);



end

