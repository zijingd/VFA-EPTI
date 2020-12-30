% Example script to generate VFA-EPTI sampling pattern
% related to the following works
% 1) "Variable flip angle Echo Planar Time-Resolved Imaging (vFA-EPTI) for fast high-resolution gradient echo myelin water imaging", Zijing Dong et al.
% 2) "Echo planar time-resolved imaging (EPTI)", Fuyixue Wang et al.,MRM,2019
% 3) "EPTI with Subspace Reconstruction and Optimized Spatiotemporal Encoding", Zijing Dong et al.2020
% Zijing Dong <zijingd@mit.edu> Nov/2020 

clear;
clc;
data_size = [50 88 72 8];  % nt ny nz nc
Block_size_y=8;
Block_size_z=4;

Nblock_central = [5,8]; % number of blocks for oversampled central k-space

shift_y_odd = [0,2,4,6,0,2,4,6]; % shift along ky from the initial pattern of the odd echo sections, 8 values for 8 FAs
shift_z_odd = [0,0,0,0,2,2,2,2]; % shift along kz from the initial pattern of the odd echo sections
shift_y_even = [0,2,4,6,0,2,4,6]; % shift along ky from the initial pattern of the even echo sections
shift_z_even = [2,2,2,2,0,0,0,0]; % shift along kz from the initial pattern of the even echo sections

% shift parameters for central k-space
shift_y_odd_central = [4,6,0,2,4,6,0,2]; % shift along ky from the initial pattern of the odd echo sections, 8 values for 8 FAs
shift_z_odd_central = [0,0,0,0,2,2,2,2]; % shift along kz from the initial pattern of the odd echo sections
shift_y_even_central = [4,6,0,2,4,6,0,2]; % shift along ky from the initial pattern of the even echo sections
shift_z_even_central = [2,2,2,2,0,0,0,0]; % shift along kz from the initial pattern of the even echo sections

mask_sample_FA=[]; % sampling mask of different FAs [nt,ny,nz,nc]
for nFA = 1:8
    mask_tmp = EPTI_sampling_mask_Shift(data_size,Block_size_y,Block_size_z,shift_y_odd(nFA),shift_z_odd(nFA),shift_y_even(nFA),shift_z_even(nFA));
    mask_tmp2 = EPTI_sampling_mask_Shift(data_size,Block_size_y,Block_size_z,shift_y_odd_central(nFA),shift_z_odd_central(nFA),shift_y_even_central(nFA),shift_z_even_central(nFA));
    mask_Center = Gen_EPTI_sampling_Square_mask(data_size,Block_size_y,Block_size_z,Nblock_central(1),Nblock_central(2));
    mask_sample_FA{nFA} = mask_tmp+mask_tmp2.*mask_Center;
end

disp(['Undersampling Rate = :',num2str(1./mean(mask_sample_FA{1}(:)))]);

figure; imshow(squeeze(sum(mask_sample_FA{1}(1:8,:,:,1),1)),[])
figure; imshow(squeeze(sum(mask_sample_FA{2}(1:8,:,:,1),1)),[])
