function [U, X0] = gen_GE_basis_T2B0(N, ETL, T0, dt, T2vals,B0vals)
% Generate basis
%
% Inputs:
%  N -- maximum number of T2 signals to simulate
%  ETL -- echo train length
%  T0 -- initial echoes time
%  TE (s) -- echo spacing
%  T2vals (s) -- array of T2 values to simulate
%  verbose -- print verbose output
%
% Outputs:
%  U -- temporal basis based on PCA
%  X -- [T, L] matrix of simulated signals

% randomly choose T2 values if more than N are given
if length(T2vals) > N
    idx = randperm(length(T2vals));
    T2vals = T2vals(idx(1:N));
end

TEs=(dt:dt:ETL*dt)+T0;

LT1 = length(T2vals);
LT2 = length(B0vals);

X0 = zeros(ETL, LT1, LT2);

for ii=1:LT1
    R2 = 1/T2vals(ii);
    for jj=1:LT2
    X0(:,ii,jj)=exp(-TEs(:)*R2).*exp(1i*2*pi*B0vals(jj)*TEs(:));
    end
end
X0=reshape(X0,ETL,[]);
[U, ~, ~] = svd(X0, 'econ');

end