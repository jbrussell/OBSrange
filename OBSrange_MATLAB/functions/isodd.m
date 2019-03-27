function [ answ ] = isodd( x )
% [ answ ] = isodd( x )
% 
% determines if a number is odd or not
 
y=(x-1)/2;
yy=double(int8(y));
dy=y-yy;
dy=abs(dy-round(dy));
if dy==0
    answ=1;
elseif dy==0.5
    answ=0;
else 
    answ=0;
    fprintf('input is not an integer\n')
end

% mod(x,2)
end