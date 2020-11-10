clear;
close all;
addpath(genpath('Funcs'));
%% 
directory_data = 'Data/';
load([directory_data,'Prior_Estimates.mat'])
load([directory_data,'Acq_params.mat'])
%%
T2vals=[5:1:100,102:3:200,210:10:300];
B0vals=[-50:1:50]; % unit Hz
K = 8; % subspace size
N = 256; % maximum number of unique T2s values for training
%% pulse parameters
param.TEs_GRE = TEs(:);   % echo time for readout unit: s
dt =TEs(2)-TEs(1);
T2vals = T2vals/1000;
%%
[U, X] = gen_GE_basis_T2B0(N, nt_GE, t0, dt, T2vals, B0vals);
Phi = U(:,1:K);
Z = Phi*Phi'*X;
err = norm(X(:) - Z(:)) / norm(X(:));
fprintf('Relative norm of error: %.6f\n', err);

show_fig=1;
if show_fig==1
    figure;
    plot(TEs*1000, real(X(:,1:1:end)), 'linewidth', 2); hold on;
    xlabel('Virtual echo time (ms)');
    ylabel('Signal value');
    ftitle('Signal evolutions for distribution of T2 T2* and T1 values', 24)
    faxis;

    figure;
    subplot(1,2,1); plot(real(Phi), 'linewidth', 3);
    ftitle('Subspace curves Real Part', 24)
    faxis;
    subplot(1,2,2); plot(imag(Phi), 'linewidth', 3);
    ftitle('Subspace curves Imag Part', 24)
    faxis;
    %% Project the signal evolutions onto the subspace
    figure;
    plot(TEs*1000, real(Z(:,1:1:end)), 'linewidth', 2); hold on;
    xlabel('Virtual echo time (ms)');
    ylabel('Signal value');
    ftitle('Projected signal evolutions', 24)
    faxis;
end

param.dt=dt;
param.t0=t0;
save([directory_data,'Basis_GRE_withPhase','.mat'],'Phi','X','U','param','-v7.3');