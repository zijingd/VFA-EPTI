function [ y ] = Mask_adjoint( y, mask_sample )
%
y=reshape(y,size(mask_sample));
y = bsxfun(@times, y, mask_sample);

end

