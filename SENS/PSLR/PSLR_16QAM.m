% 计算全重叠和半重叠PSLR
clc; clear; close all;

%%规模
flag = 3;
switch flag
    case 1
        n_fft = 512;            % 512点FFT
        n_sc = 200;            % 有效子载波数
        n_sym = 40;            % OFDM符号数
    case 2
        n_fft = 1024;           % 2048点FFT
        n_sc = 600;             % 600个有效子载波
        n_sym = 100;            % 100个OFDM符号
    case 3
        n_fft = 2048;           % 2048点FFT
        n_sc = 600;             % 600个有效子载波
        n_sym = 100;            % 100个OFDM符号
end

%% 信号参数
fc = 10e9;              % 载波频率 10GHz
delta_f = 120e3;        % 子载波间隔 120kHz
T_O = 1 / delta_f;      % 有效符号时长 = 1/120kHz = 8.33μs
T_cp = 0.58e-6;         % CP时长 0.58μs 
T_sym = T_O + T_cp;     % 总符号时长 = 8.91μs



fs = n_fft * delta_f;   % 采样率

% 计算CP长度（按FFT点数）
n_cp = round(T_cp * fs);  % CP采样点数


rho = 0.5;
% 设置延迟-多普勒范围
c = 3e8;  % 光速
lambda = c / fc;  % 波长 ≈ 0.03 m

% 延迟范围 100m
max_range = 100;  
tau_max =  max_range / c; 
N_tau = 201;
tau_vec = linspace(-tau_max, tau_max, N_tau);

% 多普勒范围 ±100m/s
max_velocity = 100;  
fd_max = 2 * max_velocity / lambda;
N_fd = 201;
fd_vec = linspace(-fd_max, fd_max, N_fd);

%% 随机16qam信号
qam_symbols = qammod(randi([0 15], n_sc, n_sym), 16, 'UnitAveragePower', true);
S_fd1 = zeros(n_fft, n_sym);

% 使用n_sc个子载波
sc_start = 2;
sc_end =  n_sc + 1;
S_fd1(sc_start:sc_end, :) = qam_symbols;


% 计算模糊函数和PSLR
draw_flag = 1;
tic
[PSLR_weighted1, PSLR_delay1, PSLR_doppler1] = calc_PSLR_weighted(S_fd1, tau_vec, fd_vec, fs, draw_flag);
toc
%% 
tic
[PSLR_weighted2, PSLR_delay2, PSLR_doppler2] = calc_PSLR_fast(S_fd1, tau_vec, fd_vec, fs, rho);
toc