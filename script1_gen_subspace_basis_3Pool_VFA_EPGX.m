clear;
close all;
addpath(genpath('Funcs'));
%% 
directory_data = 'Data/';
load([directory_data,'Prior_Estimates.mat'])
load([directory_data,'Acq_params.mat'])
%% Simulate and generate basis
Global_B0 = 0; % 
B1 = 0.8:0.1:1.2;
Num=length(S0_EW);     % different B0 offsets to correct inaccurate calibration
N_B0 = length(Global_B0);
N_B1 = length(B1);
npulse = 50;
phiCycle = RF_phase_cycle(npulse,50);
N_pick = 1000;
X0 = [];
for i = 1:N_B0
    disp(['Simulating ',num2str(i),'/',num2str(N_B0)]);
    for j = 1:N_B1
    disp(['Simulating ',num2str(j),'/',num2str(N_B1)]);
    
    idx = randperm(Num);                   % random pick values
    Freq_EW_tmp = Freq_EW(idx(1:N_pick));
    Freq_IW_tmp = Freq_IW(idx(1:N_pick));  
    Freq_MW_tmp = Freq_MW(idx(1:N_pick));  

    Freq_EW_tmp = Freq_EW_tmp+Global_B0(i); 
    Freq_IW_tmp = Freq_IW_tmp+Global_B0(i);  
    Freq_MW_tmp = Freq_MW_tmp+Global_B0(i);

    S0_EW_tmp = S0_EW(idx(1:N_pick));
    S0_IW_tmp = S0_IW(idx(1:N_pick));  
    S0_MW_tmp = S0_MW(idx(1:N_pick));  
    
    T1_mw_tmp = T1_mw(idx(1:N_pick));  
    T1_ow_tmp = T1_ow(idx(1:N_pick));
    X0_tmp = zeros(nt_GE*nFA,N_pick);

    T2vals_EW_tmp = T2vals_EW(idx(1:N_pick));
    T2vals_IW_tmp = T2vals_IW(idx(1:N_pick));
    T2vals_MW_tmp = T2vals_MW(idx(1:N_pick));
    X0_tmpX = zeros(nt_GE*nFA,N_pick);
    b1 = B1(j);
    for ii=1:N_pick
        for kfa=1:length(FA)
           T3D_all{kfa} = PrecomputeT(phiCycle,d2r(FA(kfa)*b1));
        end
        sHat = mwi_model_2T13T2scc_epgx_EPTI(FA,TEs,TR,S0_MW_tmp(ii),S0_IW_tmp(ii),S0_EW_tmp(ii),...
        T2vals_MW_tmp(ii),T2vals_IW_tmp(ii),T2vals_EW_tmp(ii),T1_mw_tmp(ii),T1_ow_tmp(ii),...
        Freq_MW_tmp(ii),Freq_IW_tmp(ii),Freq_EW_tmp(ii),Global_B0(i),0,b1,2,npulse,T3D_all);
        sHat=permute(sHat,[2 1]);
        X0_tmpX(:,ii) = sHat(:);
    end
    X0 = cat(2,X0,X0_tmpX);
    end
end
X0=single(X0);
disp('start SVD');
[U, ~, ~] = svd(X0, 'econ');
disp('SVD done');
X = X0;
clear X0;
%%
K = 18; % subspace size
Phi = U(:,1:K);
Z = Phi*Phi'*X;
err = norm(X(:) - Z(:)) / norm(X(:));
fprintf('Relative norm of error: %.6f\n', err);
figure;
plot(real(X(:,1:1000:end)), 'linewidth', 2);
plot(abs(X(:,1:1000:end)), 'linewidth', 2);
xlabel('Virtual echo time (ms)');
ylabel('Signal value');
ftitle('Signal evolutions for distribution of T2 values', 24)
faxis;
figure;
plot(real(Phi), 'linewidth', 3);
ftitle('Subspace curves', 24)
faxis;
%%
param.dt=TEs(2)-TEs(1);
param.t0=t0;
param.TEs_GRE=TEs;
save([directory_data,'Basis_GE_3cx_T1T2s_Phase_B1_EPGX','.mat'],'Phi','U','param','-v7.3');