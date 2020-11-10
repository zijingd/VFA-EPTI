function [ y ] = Mask_forward( y, mask_sample )
%
y = bsxfun(@times, y, mask_sample);
y=y(:);

end

