function x = axpos(ax,val)
%  x = axpos(ax,val)
% 
%  Function to grab an element of the position vector for a particular set
%  of axes (or figure handle). 
%  I.e. 
%       pos = get(ax,'pos');
%       x = pos(val)
%  This is to allow one-line operations on the axes dimensions.
% 
% Z. Eilon 2017

if nargin < 2 || isempty(val)
    val = 1:4;
end


pos = get(ax,'pos');
x = pos(val);

end