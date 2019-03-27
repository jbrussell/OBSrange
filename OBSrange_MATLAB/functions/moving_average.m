function [ B ] = moving_average( A,n )
%[ B ] = moving_average( A,n )
% 
% B is a tapered, n-point moving average of the vector A 
% n must be integer, and will be rounded up if even.
% 
%   Z. Eilon,  2017

if ~isodd(n)
    n = round_level(n,2)+1;
end

if isrow(A), A = A(:); end

dy = 0.5*(n-1);
L = size(A,1);

% mid full moving average section
C = toeplitz([1,zeros(1,L-n)],[ones(1,n),zeros(1,L-n)]);

% taper section % wish I could think of a way to avoid the loop!
D = zeros(dy,L);
f = 2*[1:dy]'-1;
for ii = 1:length(f)
    D(ii,1:f(ii)) = 1;
end

%complete moving average filter (not normalised)
G = [D;C;rot90(D,2)];

B = G*A./sum(G,2); % multiply and normalise
end

