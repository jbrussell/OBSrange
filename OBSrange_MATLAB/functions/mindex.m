function [ minX_ind ] = mindex( X,a )
% [ minX_ind ] = mindex( X,a )
% 
%   simple function to return the index of the minimum point in vector X
%   
%   if a second argument is given, the function outputs the index of the
%   point in X closest to the water level, a
%   basically just outputs the second output of the "min" function,
%   without giving you the magnitude of the minimum value.
% 
%   N.B. can be used as a zero finder if a==0
%
%   Intended for use when calling the value in one vector corresponding to
%   the minimum value in X - i.e more efficient than the clunkier:
%       Y(find(X==min(X)) or, more often, Y(find((X-a)==min(X-a)))
%   Instead, can now use
%       Y(mindex(X)) or Y(mindex(X,a))
% 
%  Z. Eilon, 2016

if nargin<2
[~,minX_ind] = min(X);
else
[~,minX_ind] = min(abs(X-a));
end
end

