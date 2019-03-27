function [ C ] = cov2cor( Cm )
%[ C ] = cov2cor( Cm )
%
% Convert model covariance matrix to correlation matrix (unit covariance
% matrix with ones along diagonal).

D = diag(sqrt(diag( Cm )));
C = D \ Cm / D;


end

