function [ KMLEN, AZ ] = distance_km(LAT1,LON1,LAT2,LON2)
%[ KMLEN, AZ ] = distance_km(LAT1,LON1,LAT2,LON2)
%   
% This is just the distance function but with the Earth's ellipsoid input
%   by defult, such that the distance output is in km.
%   usage distance_km([LAT1,LON1],[LAT2,LON2]) is also ok
% 
% Z. Eilon, 2013

if nargin==2 % if the input is actually distance_km([LAT1,LON1],[LAT2,LON2])
    P1 = LAT1;
    P2 = LON1;
    LAT1 = P1(1);
    LON1 = P1(2);
    LAT2 = P2(1);
    LON2 = P2(2);  
end

[KMLEN,AZ] = distance(LAT1,LON1,LAT2,LON2,[6378.137, 0.08181919]);



end

