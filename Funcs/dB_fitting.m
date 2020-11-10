function [ dB,fit_error,phase ] = dB_fitting( phase,TEs,mask,unwrap_flag)
% Fitting of dB
if(unwrap_flag==1)
    phase=unwrap(phase,[],3);
end
phase=phase./(2*pi);
[nx,ny,nt]=size(phase);
A = [ones(nt,1),TEs(1:nt)];   
dB=zeros(nx,ny);
fit_error=zeros(nx,ny);
for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
           % fit voxel
           signal=double(squeeze(phase(ii,jj,:)));
           param = A\signal;
           phase_fitting=A*param;
           error = abs(phase_fitting-signal);
           dB(ii,jj) = param(2);
           fit_error(ii,jj)=sqrt(sum(error.^2)/size(error,1));   
        end
    end % end jj 
end % end ii
dB(isnan(dB)) = 0;
% dB=dB-mean(dB(mask(:)));
dB = dB.*mask;
phase=phase.*(2*pi);
end

