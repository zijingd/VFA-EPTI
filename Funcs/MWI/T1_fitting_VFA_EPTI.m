function [T1,PD] = T1_fitting_VFA_EPTI(imgParam,Echo_use)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
FA = imgParam.fa;
TR = imgParam.tr;
mask = imgParam.mask;
img = abs(imgParam.img(:,:,:,Echo_use,:));
[nx,ny,nz,nt,nFA] = size(img);
T1 = zeros(nx,ny,nz);
PD = zeros(nx,ny,nz);
t1 = zeros([nt,1]);
m0 = zeros([nt,1]);

for x = 1:nx
%   disp(x);
  for y = 1:ny
    for z = 1:nz
    if mask(x,y,z)>0
        for t = 1:nt
            S = squeeze(img(x,y,z,t,:));
            [t1(t), m0(t)] = DESPOT1(S,FA,TR);
        end
        T1(x,y,z) = mean(t1);
        PD(x,y,z) = mean(m0);
    end
    end
  end
end

end

