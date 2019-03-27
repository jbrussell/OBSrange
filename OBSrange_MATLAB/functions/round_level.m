function [ xx ] = round_level( x, level )
%[xx] = round_level(x,level)
% 
% Rounds the number/vector x to the nearest "level". 
% E.g. round_level(5,4) = 4;
% E.g. round_level(5,3) = 6;
% 
% Z. Eilon, 2012

xx = round(x./level).*level;

end

