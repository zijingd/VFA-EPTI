function [im_recon,res_a]=EPTI_Image_Recon_Subspace_GESE_Tik(kdata,mask_sample,smap,Phase_total,Phi,a0,nIter,lambda)
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

total_size=nt*nx*ny*nc;
A_for = @(a) A_GESE_forward_Tik(a,lambda,T_for,P_for,S_for,F_for,M_for);
A_adj = @(y) A_GESE_adjoint_Tik(y,total_size,lambda,T_adj,P_adj,S_adj,F_adj,M_adj); 

a0 = a0(:);
y0 = [kdata(:); 0*a0(:)];
clear kdata;
[res,FLAG,RELRES,ITER,RESVEC] = lsqr(@aprod,y0,1e-6,nIter,[],[],a0(:),A_for,A_adj);
res_a=reshape(res,dims);
im_recon=temporal_forward(res, Phi, dims);

function [res,tflag] = aprod(a,A_for,A_adj,tflag)	
	if strcmp(tflag,'transp')
        res = A_adj(a);
    else
        res = A_for(a);
    end

