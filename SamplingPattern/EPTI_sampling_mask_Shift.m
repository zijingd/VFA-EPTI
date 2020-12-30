function [mask_sample] = EPTI_sampling_mask_Shift(data_size,Block_size_y,Block_size_z,shift_y_odd,shift_z_odd,shift_y_even,shift_z_even)
if nargin<2
    Block_size_y=8;
    Block_size_z=4;
end
nt = data_size(1);
ny = data_size(2);
nz = data_size(3);
nc = data_size(4);

%% Generate sampling mask
Block_size_t=4;
Jump_y = 4;
Jump_z = 1;
N_dt_Group=ceil(nt/Block_size_t);
%% MASK_1
mask_sample1 = zeros(Block_size_y,Block_size_z,Block_size_t);

index_y=1;
index_z=1;
for t=1:Block_size_t/2
    mask_sample1(index_y,index_z,t)=1;
    index_y=mod(index_y-1+Jump_y,Block_size_y)+1;
    index_z=mod(index_z-1+Jump_z,Block_size_z)+1;
end
index_y=Jump_y/2+1;
for t=Block_size_t/2+1:Block_size_t
    mask_sample1(index_y,index_z,t)=1;
    index_y=mod(index_y-1+Jump_y,Block_size_y)+1;
    index_z=mod(index_z-1+Jump_z,Block_size_z)+1;
end
%% MASK_2
mask_sample2 = mask_sample1;

mask_sample11=circshift(mask_sample1,[shift_y_odd shift_z_odd 0]);
mask_sample22=circshift(mask_sample2,[shift_y_odd shift_z_odd 0]);

mask_sample=cat(1,mask_sample11,mask_sample22);
mask_sample=repmat(mask_sample,ceil(ny/Block_size_y/2),ceil(nz/Block_size_z),1,nc);

mask_sample_p1=mask_sample;
mask_sample1=circshift(mask_sample1,[shift_y_even shift_z_even 0]);
mask_sample2=circshift(mask_sample2,[shift_y_even shift_z_even 0]);
mask_sample=cat(1,mask_sample1,mask_sample2);
mask_sample=repmat(mask_sample,ceil(ny/Block_size_y/2),ceil(nz/Block_size_z),1,nc);
mask_sample_p2=mask_sample;

mask_sample_temp=[];
for i=1:N_dt_Group
    if mod(i,2)==1
        mask_sample_temp=cat(3,mask_sample_temp,mask_sample_p1);
    else
        mask_sample_temp=cat(3,mask_sample_temp,mask_sample_p2);
    end
end

mask_sample_temp=permute(mask_sample_temp,[3 1 2 4]);
mask_sample=mask_sample_temp(1:nt,1:ny,1:nz,1:nc);

end
