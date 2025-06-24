
% adaptive FSST toolbox v1.0
% this is an example for the two-component LFM signal

% Lin Li, March 2018
% Copyright (c) 2018 Xidian University
% 
% 
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following papers:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn
% or jiangq@umsl.edu. 




clear;
close all;
clc;

load Data\LinearChirp.mat signal
N = 256;

Fs = 256;
Ts = 1/Fs;
t = [0:N-1]/Fs; % s
fi = t.*Fs;    % Hz  

% c1 = 12;
% c2 =64;
% 
% r1 = 50;
% r2 =-34; 
% 
% x1 = cos(2*pi*(c1*t+r1/2*t.^2));
% x2 = cos(2*pi*(c2*t+r2/2*t.^2));
% f1 = c1+r1*t;
% f2 = c2+r2*t;
% s = x1+x2;

[s, f, r, c] = signal_gen('n');
r1 = r(1);
r2 = r(2);
c1 = c(1);
c2 = c(2);
f1 = f(1, :);
f2 = f(2, :);


figure;
plot(t,s,'k-','linewidth',1);
title("Signal in Time domain")

figure;
plot(t,f1,'r-','linewidth',2);
hold on;
plot(t,f2,'b-','linewidth',2);
title("Instantenous frequencies")

%%
gamma = 0.001;

lmd_1 = 0.2:0.01:0.8;    %%duration
sgm_est2 = zeros(len(lmd_1));
for i = 1:len(lmd_1)
    lmd = lmd(i);
    arfa = 1/(2*pi)*sqrt(2*log(1/lmd)); 
    ci_1 = [ci_1 arfa./((f2-f1)/N)];        % Eq.(47)
    tic
    Ak = 2*pi*arfa*(abs(r1/(N^2)) + abs (r2/(N^2)));  % normalized
    if length(Ak)==1
        Ak = ones(1,N)*Ak;
    end
    Bk = (f2-f1)/N;
    Ck = sqrt(Bk.^2-8*arfa*Ak);
    for bb = 1:N
        if Ak(bb)~=0 && real(Ck(bb)) == Ck(bb)
            ci_2(bb) = max(1,(Bk(bb)-Ck(bb))/(2*Ak(bb)));
        else
            ci_2(bb) = ci_1(bb);
        end
    end
       
    ci_1 = ci_1/N;
    ci_2 = ci_2/N;
    
    %% the Renyi Entropy method
    sgm_1 = 0.005;
    d_sgm = 0.001;
    sgm = sgm_1:d_sgm:0.1;
    al = 2.5;
    t1 = 4;
    [sgm_u,sgm_R1,sgm_R2] = Renyi_STFT_SST(s,sgm,al,t1);         % Entropy
    ci_tv2=smooth(sgm_R1,20,'rlowess');                          % \sigma_{Re}
    ci_tv3=smooth(sgm_R2,20,'rlowess');                          % \sigma_{Re2}
    
    gamma1 = 0.3;
    ci_est = tv_para_fast(s,sgm_1,sgm_u,d_sgm,lmd,gamma1);  
    ci_tv1=smooth(ci_est,20,'rlowess'); 
    sgm_est2(i) = ci_tv1;
end

figure;
plot(t,ci_1,'b--','linewidth',2);
hold on;
plot(t,ci_2,'k--','linewidth',2);
% plot(t,sgm_u,'m-','linewidth',2);
% plot(t,ci_tv1,'r-','linewidth',2);
% plot(t,ci_tv2,'g-','linewidth',2);
% plot(t,ci_tv3,'c-','linewidth',2);
% legend('\sigma_1(t)', '\sigma_2(t)','\sigma_u(t)','\sigma_{est}(t)','\sigma_{Re}(t)','\sigma_{Re2}(t)','Location','north')
% axis([0 1 0 0.11]);
% set(gca,'FontSize',20)
% xlabel('Time (s)','FontSize',20);
% ylabel('\sigma','FontSize',20);
% title('Various time-varying parameters of \sigma(t)');

%%
sigma_tv = ci_2;
[Vs_tv Ts1_tv Ts2_tv] = adap_stft_sst(s,gamma,sigma_tv*Fs);
figure;
imageSQ(t,fi,abs(Vs_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive STFT with \sigma_2(t)');

figure;
imageSQ(t,fi,abs(Ts1_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive FSST with \sigma_2(t)');

figure;
imageSQ(t,fi,abs(Ts2_tv));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Adaptive 2nd-order FSST with \sigma_2(t)');

%% estimated sigma for fsst
sigma_tv = ci_tv1;
[Vs_tv Ts1_tv Ts2_tv] = adap_stft_sst(s,gamma,sigma_tv*Fs);
figure;
subplot(131);
imagesc(t,fi,abs(Vs_tv));
% axis([0 1 0 40]);
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Adaptive STFT with the estimated \sigma(t)');

subplot(132);
imagesc(t,fi,abs(Ts1_tv));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Adaptive FSST with the estimated \sigma(t)');

subplot(133);
imagesc(t,fi,abs(Ts2_tv));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Adaptive 2nd-order FSST with the estimated \sigma(t)');

%%
% regular-PT adaptive SST with Re1
FSST_rpt = regular_pt_adap_sst(s,gamma,ci_tv2); 
figure;
subplot(121);
imageSQ(t,fi,abs(FSST_rpt));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('Regular-PT adaptive SST with \sigma_{Re1}(t)');


% regular-PT adaptive 2nd-order SST with Re2
FSST2_rpt = regular_pt_adap_sst2(s,gamma,ci_tv3);
subplot(122);
imageSQ(t,fi,abs(FSST2_rpt));
axis xy;
xlabel('Time (s)','FontSize',20);
ylabel('Frequency (Hz)','FontSize',20);
set(gca,'FontSize',20);
title('2nd-order regular-PT adaptive SST with \sigma_{Re2}(t)');

%% conventional stft and fsst
[STFT, SST, VSST] = stft_sst2(s,gamma,0.02);

figure;
imagesc(t, fi, abs(STFT))
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Conventional STFT when \sigma=0.02');


figure;
imagesc(t,fi,abs(SST));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Conventional SST when \sigma=0.02');

figure;
imagesc(t,fi,abs(VSST));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Conventional 2nd-order SST when \sigma=0.02');

toc
