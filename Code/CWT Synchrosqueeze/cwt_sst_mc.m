function [Ws, Rs, Ts, Ts2, aj] = cwt_sst_mc(s,gamma,mu,ci,nv)
% computes the CWT Ws, the phase transform of CWT Rs, and the
% synchrosqueezing transform SST Ts of the signal s.
% Input:
%   s: signal
%   gamma: threshold
%   mu,ci: window parameter
% Output:
%   Ws, Rs, Ts
% By Dr. Lin Li, 11/15/2016
% Copyright (c) 2017 Xidian University
% 
% This program can not be used for commercialization without the authorization of its author(s). 

% If you use this toolbox for research, please must cite the following paper:
% [1] Lin Li, Haiyan Cai and Qingtang Jiang, "Adaptive synchrosqueezing transform with a time-varying parameter for non-stationary signal separation," preprint, 2018, arXiv:1812.11364.
% [2] Lin Li, Haiyan Cai, Hongxia Han, Qingtang Jiang, and Hongbing Ji, "Adaptive short-time Fourier transform and synchrosqueezing transform for non-stationary signal separation," preprint, 2018, arXiv:1812.11292.
% [3] Haiyan Cai, Qingtang Jiang, Lin Li, and Bruce W. Suter, "Analysis of adaptive short-time Fourier transform-based synchrosqueezing", preprint, 2018, arXiv:1812.11033.

% For any comment or bug report, please send e-mail to lilin@xidian.edu.cn or jiangq@umsl.edu. 


ol = length(s);  % original length
s = s(:);      % Turn into column vector
s = s - mean(s);      %zero mean
n = length(s);     % checking length of signal
p2 = ceil(log2(n));
s = [s;zeros(2^p2-n,1)];     % pading the signal by zeros
n = length(s);
    
L = log2(n);
na = L*nv;    % number of scales
j = [1:na];
m = [0:n-1];   % time index
aj = 2.^(j/nv);    % scales

% Padding, symmetric
sleft = flipud(s(2:n/2+1));
sright = flipud(s(end-n/2:end-1));
x = [sleft; s; sright];
n1 = length(sleft);
clear xleft xright;

N = length(x);
Ws = zeros(na,N);
dWs = Ws;

x = x(:).';
xh = fft(x);

%%% computes CWT with "Morlet"
cs = pi^(-1/4)*ci^(-1/2)*(1+exp(-mu^2*ci^2)-2*exp(-3/4*mu^2*ci^2))^(-1/2);
psihfn = @(w)cs*(2*pi*ci^2)^(1/2)*(exp(-1/2*ci^2*(w-mu).^2)-exp(-1/2*ci^2*(w.^2+mu^2)));

% psihfn = @(w) exp(-2*pi^2*ci^2*(w-mu).^2);  % lin li 1/15/2018

k = 0:(N-1);
xi = zeros(1,N);
xi(1:N/2+1) = 2*pi/N*[0:N/2];
xi(N/2+2:end) = 2*pi/N*[-N/2+1:-1];

for ai = 1:na
    a = aj(ai);
    psih = psihfn(a*xi);
    dpsih = (i*xi) .* psih;
    
    xcpsi = (ifft(psih .* xh));
    Ws(ai, :) = xcpsi;
    dxcpsi = (ifft(dpsih .* xh));
    dWs(ai, :) = dxcpsi;
end

Ws = Ws(:, n1+1:n1+n);
dWs = dWs(:, n1+1:n1+n);

hWs = (abs(Ws)>gamma).*(1./Ws);
Rs = 1/(2*pi)*dWs.*hWs;    % phase transform

wl = [0:2:n-1]/n/2;          % linear frequency bin
Ts = zeros(n/2,n);
Ts2 = zeros(n/2, n);

dw = wl(2)-wl(1);    % frequency step
for bi = 1:n
    for wi = 2:n/2
        w = wl(wi);
        tmp1 = find(imag(Rs(:,bi))>w-dw/2 & imag(Rs(:,bi))<=w+dw/2 );
        if length(tmp1)>=1
            Ts(wi,bi) = Ts(wi,bi) + log(2)/nv*ones(1,length(tmp1))*Ws(tmp1,bi);  
            Ts2(wi,bi) = Ts2(wi,bi) + log(2)/nv * aj(tmp1).^(-1/2)*Ws(tmp1,bi);    
        end
    end
end

Ws = Ws(:,1:ol);
Rs = Rs(:,1:ol);
Ts = Ts(:,1:ol);
Ts2 = Ts2(:,1:ol);

end



