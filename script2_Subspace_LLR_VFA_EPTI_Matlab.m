% Example script for VFA-EPTI reconstruction using low-rank subspace
% related to the following works
% 1) "Variable flip angle Echo Planar Time-Resolved Imaging (vFA-EPTI) for fast high-resolution gradient echo myelin water imaging", Zijing Dong et al.
% 2) "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al.,MRM,2019
% 3) "EPTI with Subspace Reconstruction and Optimized Spatiotemporal Encoding", Zijing Dong et al.2020
% Zijing Dong <zijingd@mit.edu> Nov/2020 

clear;
close all;
addpath(genpath('Funcs'));
%% load data and basis
directory_data = 'Data/';
load([directory_data,'kdata.mat'])
load([directory_data,'Basis_GRE_withPhase.mat'])
U_single = U;
load([directory_data,'Basis_GE_3cx_T1T2s_Phase_B1_EPGX.mat'])
load([directory_data,'B0_Sens.mat'])
%% Pre calcuation of phase evolution
TEs_GRE=param.TEs_GRE(:);

mask_sample=kdata~=0;
disp(['Undersampling Rate = :',num2str(1./mean(mask_sample(:)))]);
[nt,ny,nz,nc,nFA]=size(kdata);
MSK_extended = P_dB_tmp~=0;
% --------------Generate Phase--------------    
Phase_T=zeros(ny,nz,nt);
for t=1:nt
    Phase_T(:,:,t)=exp(1i*2*pi*P_dB_tmp*TEs_GRE(t)).*Phase0_tmp;  
end
%% -------------- Pre-reconstruction for each FA to update B0 field --------------
P_dB_rec=zeros(ny,nz,nFA);
Phase0_rec=zeros(ny,nz,nFA);
K = 5;
Phi = U_single(:,1:K);
a0=zeros([ny,nz,K]);
N_iteration=50;
lambda = 0;   % tikhonov regularization
nt_to_fit = 10:30; % echo used to estimate B0

for FA = 1:nFA
    mask_tmp=squeeze(mask_sample(:,:,:,:,FA));
    kdata_tmp=squeeze(kdata(:,:,:,:,FA));
    disp('Reconstruction Start');
    tic;
    im_recon=EPTI_Image_Recon_Subspace_GESE_Tik(kdata_tmp,mask_tmp,sens_map_tmp,Phase_T,Phi,a0,N_iteration,lambda);
    toc;
    disp('Reconstruction Done');
    im_recon=(im_recon).*Phase_T;
    close all;
    figure; imshow3(permute(cat(3,abs(im_recon(end:-1:1,:,4:8:42))),[1 2 3]),[0 4],[1 5]);
    PHS=angle(im_recon);
    P_dB_rec(:,:,FA) = dB_fitting(PHS(:,:,nt_to_fit),TEs_GRE(nt_to_fit),logical(MSK_extended),1);
    figure; imshow(P_dB_rec(end:-1:1,:,FA),[-100 100]);
    Phase=[];
    for t = 1:length(nt_to_fit)
        Phase(:,:,t)=(im_recon(:,:,nt_to_fit(t)).*exp(-1i*2*pi*P_dB_rec(:,:,FA)*TEs_GRE(nt_to_fit(t))));
    end
    Phase0_rec(:,:,FA)=angle(mean(Phase,3));
end

Phase0_rec=exp(1i*Phase0_rec);

Phase_T=zeros(ny,nz,nt,nFA);
for fa=1:nFA
for t=1:nt
    Phase_T(:,:,t,fa)=exp(1i*2*pi*P_dB_rec(:,:,fa)*TEs_GRE(t)).*Phase0_rec(:,:,fa);  
end
end

Phase_T = reshape(Phase_T,[ny nz nFA*nt]);
kdata = permute(kdata,[1 5 2 3 4]);
kdata = reshape(kdata,[nt*nFA,ny,nz,nc]);
mask_sample = kdata~=0;
%% ----------Joint Reconstruction--------------
K=18;
Phi = U(:,1:K);
im_rec=zeros(ny,nz,nt*nFA);
iter_ops.max_iter = 20;
llr_ops.lambda = .005;
iter_ops.rho = 0.1;
llr_ops.block_dim = [8, 4];
lsqr_ops.max_iter = 10;
lsqr_ops.tol = 1e-4;
disp('Reconstruction Start');
tic;
[im_recon,a]=EPTI_Image_Recon_Subspace_GESE_ADMM(kdata,mask_sample,sens_map_tmp,Phase_T,Phi,iter_ops,llr_ops,lsqr_ops);
toc;
disp('Reconstruction Done');
im_recon=(im_recon).*Phase_T;
close all;
echo_show = 5:nt:nt*nFA;
figure; imshow3(permute(abs(im_recon(end:-1:1,:,echo_show)),[1 2 3]),[0 6],[2 length(echo_show)/2]);
im_recon=single(im_recon);
save([directory_data,'Recon_Subspace_LLR_K',num2str(K),'.mat'],'im_recon','-v7.3');
