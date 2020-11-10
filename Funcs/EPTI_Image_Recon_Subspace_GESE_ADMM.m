function [im_recon,alpha,history]=EPTI_Image_Recon_Subspace_GESE_ADMM(kdata,mask_sample,smap,Phase_total,Phi,iter_ops,llr_ops,lsqr_ops)
% Tikhonov 
[nt,nx,ny,nc] = size(kdata);
kdata=permute(kdata,[2 3 4 1]);
mask_sample=permute(mask_sample,[2 3 4 1]);

K=size(Phi,2);

% Subspace basis projection
dims=[nx,ny,K];

T_for = @(a) temporal_forward(a, Phi, dims);
T_adj = @(x) temporal_adjoint(x, Phi);
% Image Phase
% Phase_total=Phase_T.*repmat(Phase0,[1,1,nt]);
P_for = @(x) bsxfun(@times, x, Phase_total);
P_adj = @(y) bsxfun(@times, y, conj(Phase_total));
% Sampling mask
% smap=repmat(smap,[1 1 1 nt]);
S_for = @(a) bsxfun(@times, smap, permute(a, [1, 2, 4, 3]));
S_adj = @(as) squeeze(sum(bsxfun(@times, conj(smap), as), 3));
% Sampling mask
F_for = @(x) fft2c_for(x);
F_adj = @(y) fft2c_adj(y);
% Sampling mask
M_for = @(y) Mask_forward(y, mask_sample);
M_adj = @(y) Mask_adjoint(y, mask_sample);

A_for = @(a) M_for(F_for(S_for(P_for(T_for(a)))));
A_adj = @(y) T_adj(P_adj(S_adj(F_adj(M_adj(y)))));

AHA = @(a) A_adj(A_for(a));
%% scaling
tmp = dimnorm(ifft2c(kdata), 3);
tmpnorm = dimnorm(tmp, 4);
tmpnorm2 = sort(tmpnorm(:), 'ascend');
% match convention used in BART
p100 = tmpnorm2(end);
p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
if (p100 - p90) < 2 * (p90 - p50)
    scaling = p90;
else
    scaling = p100;
end

% scaling=1;
fprintf('\nScaling: %f\n\n', scaling);

kdata = kdata ./ scaling;
ksp_adj = reshape(A_adj(kdata),dims);
%% ADMM
tic;

iter_ops.objfun = @(a, sv, lam) 0.5*norm_mat(kdata(:) - A_for(a))^2 + lam*sum(sv(:));

alpha_ref = RefValue;
alpha_ref.data = zeros(nx, ny, K);
history = iter_admm(alpha_ref, iter_ops, llr_ops, lsqr_ops, AHA, ksp_adj, @admm_callback);
toc;

disp(' ');
disp('Reconstruction Done');

%% Project and re-scale
alpha = alpha_ref.data;
% res_a=reshape(res,dims);
im_recon=temporal_forward(alpha(:), Phi, dims);
% im = T_for(alpha);

disp('Rescaling')
im_recon = im_recon * scaling;


