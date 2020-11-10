%% Test 1: Steady state gradient echo vs EPG
load_module_epg_epgx;
% addpath(genpath('lib'));
% addpath(genpath('EPGX-src'));

%% SPGR (RF spoiling)
% Use simple parameters

%%% Sequences
% TR = 5;
TR = 30;    % ms
alpha = 30; % degree
phi0 = 50;  % RF spoiling initial phase

%%% Relaxation parameters: single pool, copy MT model
T1=1000;    %ms
T2=100;     %ms

%%% Relaxation parameters: MT
T1_MT = [1000 1000];    %ms
f_MT = 0.117;           %
k_MT = 4.3e-3;          % ms^-1
T2_MT = 100;            % ms

%%% Relaxation parameters: exchange
T1x = [1000 500]; % s
T2x = [100 20];  % s
ka = 2e-3;                 % s^-1
fb = 0.2;               % fraction
freqoffset = 10e-3; % hz

%%% RF saturation factor for MT
G = 15.1;% us 
b1 = 13; % uT
gam = 267.5221 *1e-3; %< rad /ms /uT
trf = d2r(alpha)/(gam*b1);% ms
b1sqrdtau = b1^2*trf;

% EPG simulations
npulse=200;
phi = RF_phase_cycle(npulse,phi0);

%% Test 1: Effect of T1 on steady-state excitation
phi0 = 50;   % RF phase
alpha = 30;  % FA
TR = 30;

T1a = 3000;
T2a = 100;
npulse = floor(5*T1a/TR);    % number of excitation pulses
AA = d2r(alpha)*ones(npulse,1);
phia = RF_phase_cycle(npulse,phi0);
tmpa = abs(EPG_GRE(AA,phia,TR,T1a,T2a));

T1b = 500;
T2b = 100;
phib = RF_phase_cycle(npulse,phi0);
tmpb = abs(EPG_GRE(AA,phib,TR,T1b,T2b));

figure(11)
subplot(211)
plot(tmpa);hold on
plot(tmpb); legend('Long T1','Short T1');hold off;title('Signal of each excitation');
subplot(212)
plot(abs(diff(tmpa)));hold on
plot(abs(diff(tmpb)));hold off
legend('Long T1','Short T1');title('Abs difference between successive excitation');
inda = find(abs(diff(tmpa))<1e-7);
indb = find(abs(diff(tmpb))<1e-7);
inda(1)
indb(1)

%% Test 2: Effect of T2 on T1w curve
phi0 = 50;   % RF phase
alpha_arr = 1:90 ;  % flip angle array
TR = 30;

T1 = 1000;
T2a = 100;
T2b = 20;
npulse = floor(5*T1/TR);    % number of excitation pulses

Sig = zeros(length(alpha_arr),3);
for kfa=1:length(alpha_arr)
    alpha = alpha_arr(kfa);
    AA = d2r(alpha)*ones(npulse,1);
    phi = RF_phase_cycle(npulse,phi0);
    [F0,Fn,Zn,F] = EPG_GRE(AA,phi,TR,T1,T2a);
%     tmpb = abs(EPG_GRE(AA,phi,TR,T1,T2b));
%     tmpc = abs(EPG_GRE(AA,phi,TR,T1,0));
    
%     Sig(kfa,1) = tmpa(end);
%     Sig(kfa,2) = tmpb(end);
%     Sig(kfa,3) = tmpc(end);
end

% figure(12)
% plot(Sig);hold on
% legend('Long T2','Short T2','Zero T2');hold off;title('Signal of different T2');

%% Test 3: Effect of frequency offset between 2 pools on T1w curve
phi0 = 50;   % RF phase
alpha_arr = 1:90 ;  % flip angle array
TR = 30;

T1x = [1000 500]; % ms
T2x = [100 20];  % ms
ka = 2e-3;                 % ms^-1
fx = 0.2;               % fraction
fba = 0;
fbb = 10/1000;
fbc = 50/1000;
fbd = -50/1000;
npulse = floor(5*T1x(1)/TR);    % number of excitation pulses

Sig = zeros(length(alpha_arr),4);
for kfa=1:length(alpha_arr)
    alpha = alpha_arr(kfa);
    AA = d2r(alpha)*ones(npulse,1);
    phi = RF_phase_cycle(npulse,phi0);
    tmpa = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,ka,'delta',fba);
    tmpb = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,ka,'delta',fbb);
    tmpc = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,ka,'delta',fbc);
    tmpd = EPGX_GRE_BM(AA,phi,TR,T1x,T2x,fx,ka,'delta',fbd);
    
    Sig(kfa,1) = tmpa(end);
    Sig(kfa,2) = tmpb(end);
    Sig(kfa,3) = tmpc(end);
    Sig(kfa,4) = tmpd(end);
end

figure(13)
plot(abs(Sig));hold on
legend('No offset','10Hz','50Hz','-50Hz');hold off;title('Signal of different frequency offset');



%% Compute whole T1w curve for different epg models
phi = 50;
alpha_arr = 1:90 ;  % flip angle array
nfa=length(alpha_arr);
npulse = floor(5*T1/TR);    % number of excitation pulses
phi = RF_phase_cycle(npulse,phi);

% if 0 % load stored result if 0
if 1
Sig = zeros(nfa,3);
% figure(4)
figure
clf
for ii=1:nfa
    alpha = alpha_arr(ii);
    AA = d2r(alpha)*ones(npulse,1);
   
    % Compute RF spoling phase cycles

    % single pool
    % standard Ernst equation
    tmp1 = EPG_GRE(AA,phi,TR,T1,T2*0);
    % account for RF spoiling
    tmp2 = EPG_GRE(AA,phi,TR,T1,T2);


    % MT
    tmp3 = EPGX_GRE_MT(AA,phi,b1sqrdtau*ones(npulse,1),...
    TR,T1_MT,T2_MT,f_MT,k_MT,G);

    % BM, 2-pool exchange model, no frequency offset between 2 pools
    tmp4 = EPGX_GRE_BMsplit(AA,phi,TR,T1x,T2x,fx,kx);
    % BM, 2-pool exchange model, with frequency offset between 2 pools 
    tmp5 = EPGX_GRE_BMsplit(AA,phi,TR,T1x,T2x,fx,kx,'delta',freqoffset);
    
    % Assuming stead-state is reached in the last excitation
    Sig(ii,1) = abs(tmp1(end));     % STD
    Sig(ii,2) = abs(tmp2(end));     % STD RF spoiling
    Sig(ii,3) = abs(tmp3(end));     % MT
    Sig(ii,4) = abs(tmp4{1}(end));  % BM, Xfreq, long T1
    Sig(ii,5) = abs(tmp4{2}(end));  % BM, Xfreq, short T1
    Sig(ii,6) = abs(tmp5{1}(end));  % BM, freq, long T1
    Sig(ii,7) = abs(tmp5{2}(end));  % BM, freq, short T1
    
    disp([ii nfa])
    %save Sig Sig
    plot(alpha_arr,Sig,'-')
    drawnow
    pause(0.0001)
    
end
    save bin/SigT1models Sig
else
    load bin/SigT1models
end
legend('standard Ernst','standard with rfspoiling','MT','exchange large pool','exchange small pool','exchange fs Fwater','exchange fs Mwater')
%    subplot(122)
%         plot(alpha_arr,bsxfun(@times,Sig,1./max(Sig,[],1)),'-')