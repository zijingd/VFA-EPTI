function mask_sample = Gen_EPTI_sampling_Square_mask(data_size,Block_size_y,Block_size_z,Sample_ky,Sample_kz)

nt = data_size(1);
ny = data_size(2);
nz = data_size(3);
nc = data_size(4);

N_block_y=ceil(ny/Block_size_y);
N_block_z=ceil(nz/Block_size_z);

mask=ones(Sample_ky,Sample_kz);
mask=zpad(mask,N_block_y,N_block_z);
% figure; imshow(mask,[]);

mask_sample = zeros(ny,nz,nt);
for y=1:N_block_y
    for z=1:N_block_z
        idy = Block_size_y*(y-1)+1:Block_size_y*y;
        idz = Block_size_z*(z-1)+1:Block_size_z*z;
        if mask(y,z)==1
            mask_sample(idy,idz,:)=1;
        end
    end
end
% figure; imshow(mask_sample(:,:,1),[]);
mask_sample=repmat(mask_sample,[1,1,1,nc]);
mask_sample=permute(mask_sample,[3 1 2 4]);
end

