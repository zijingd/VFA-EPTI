function [ im ] = temporal_forward(alpha, Phi, dims,flatten)
% temporal_forward orthognal back-projection according to:
%   im = Phi alpha
%
%    im: [I1, I2, ..., In, T] Fully sampled image at each echo time (oracle)
%    Phi: [T, K] orthogonal basis operator
%    alpha: [I1, I2, ..., In,  K] image coefficients for basis elements
%  I1, I2, ..., In are image/coil/other dimensions
%  T is echo train length (number of echo times)
%  K is the number of basis coefficients

if nargin < 4
    flatten = false;
end

alpha=reshape(alpha,dims);
% dims = size(alpha);
[T, K] = size(Phi);
L = prod(dims(1:end-1));

dims2 = dims;
dims2(end) = T;

% project onto basis
x = Phi * reshape(alpha, L, K).';

if flatten
    im = x;
else
    im = reshape(x.', dims2);
end
% im=real(im);
end
