clc; clear; close all;

%% ==============================
%% 参数初始化
%% ==============================
SystemSetup_sens;
currentFolder = fileparts(mfilename('fullpath'));
load(fullfile(currentFolder, 'parameter_sens.mat'));

SNR_dB        = 0;              % 仿真信噪比 (dB)
range_view    = [20 50];        % 距离显示范围 (m)
velocity_view = [-200 200];     % 速度显示范围 (m/s)

%% ==============================
%% 发射信号生成
%% ==============================
rng(42);
tx_data        = randi([0 1], 1, N_b_data);     % 随机比特流
[txSig, S_eff] = OFDM_sens_tx(tx_data);
% txSig : 时域列向量，长度 = (N_fft+N_cp) × N_sym + 前导零 + PSS
% S_eff : N_sc × N_sym 频域发射矩阵（导频+数据，与接收端对齐）

%% ==============================
%% 多目标回波仿真
%% ==============================
tic;
c      = 3e8;
lambda = c / fc;
N_sig  = length(txSig);
t      = (0 : N_sig-1).' / fs;     % 时间向量 (s)

% --- 目标定义 ---
targets(1).range    = 32;   targets(1).velocity = 10;  targets(1).rcs = 0.08;
targets(2).range    = 35;   targets(2).velocity = 10;  targets(2).rcs = 0.07;
n_targets = length(targets);

% 计算每个目标的时延和多普勒频移
for i = 1:n_targets
    targets(i).tau        = 2 * targets(i).range / c;          % 双程时延 (s)
    targets(i).delay_samp = round(targets(i).tau * fs);         % 时延样本数
    targets(i).fd         = 2 * targets(i).velocity * fc / c;  % 多普勒频移 (Hz)
end

fprintf('========== 目标参数 ==========\n');
for i = 1:n_targets
    fprintf('目标%d: 距离=%5.1f m, 速度=%6.1f m/s, 多普勒=%7.1f Hz, 时延=%.2f μs\n', ...
        i, targets(i).range, targets(i).velocity, ...
        targets(i).fd, targets(i).tau * 1e6);
end
fprintf('==============================\n\n');

% 叠加多目标回波
rxSig = zeros(N_sig, 1);
for i = 1:n_targets
    d   = targets(i).delay_samp;
    fd  = targets(i).fd;
    att = targets(i).rcs;
    if d < N_sig
        rxSig(1+d : end) = rxSig(1+d : end) + ...
            att * txSig(1 : end-d) .* exp(1j*2*pi*fd * t(1 : end-d));
    end
end

% 添加 AWGN 噪声
sig_pwr   = mean(abs(rxSig).^2);
noise_pwr = sig_pwr / 10^(SNR_dB/10);
rxSig     = rxSig + sqrt(noise_pwr/2) * (randn(N_sig,1) + 1j*randn(N_sig,1));
fprintf('信噪比: %.1f dB\n\n', SNR_dB);

%% ==============================
%% 雷达接收：直接提取频域矩阵
%% ==============================
% 绕开通信链路（CFO补偿/均衡会破坏雷达相位），
% 用PSS定位后直接去CP+FFT，保留完整的信道相位信息。

% 1. 数据OFDM符号起始样本（跳过前导零和PSS）
pss_len    = 2^PSS_stage - 1;
data_start = N_zeros + pss_len + 1;

% 2. 逐符号去CP + FFT
n_sym_total = data_col * spread_ratio;          % 总符号数
rx_fft_mat  = zeros(N_fft, n_sym_total);
for i = 1:n_sym_total
    idx_start          = data_start + (i-1)*N_symbol + N_cp;
    rx_fft_mat(:, i)   = fft(rxSig(idx_start : idx_start + N_fft - 1), N_fft);
end

% 3. 提取有效子载波（与 S_eff 子载波索引完全对齐）
Y_eff = rx_fft_mat(pos_row : pos_row + N_sc - 1, :);   % N_sc × N_sym

%% ==============================
%% 信道估计
%% ==============================
H_est = Y_eff ./ S_eff;     % N_sc × N_sym，每个时频格点的信道响应

%% ==============================
%% Range-Doppler Map
%% ==============================
% 零填充至数据点数的8倍以上，保证曲线光滑（不改变分辨率）
Nfft_range = 2^nextpow2(N_sc   * 8);   % 距离维FFT点数
Nfft_dopp  = 2^nextpow2(n_sym  * 8);   % 多普勒维FFT点数

% 距离维：沿子载波方向做IFFT（时延域）
range_profile = ifft(H_est, Nfft_range, 1);

% 多普勒维：沿符号方向做FFT（多普勒域），fftshift使零频居中
RDM = fftshift(fft(range_profile, Nfft_dopp, 2), 2);

RDM = abs(RDM);
RDM = RDM / max(RDM(:));           % 归一化

%% ==============================
%% 坐标轴计算
%% ==============================
% --- 距离轴 ---
% IFFT bin间距对应时延 Δτ = 1/(N_sc×Δf) = 1/B
% 距离 R = c×τ/2，故每bin对应 c/(2B) × (Nfft_range/N_sc) 的距离间隔
% 化简后：range_axis(n) = n × c / (2 × Nfft_range × Δf)
range_axis = (0 : Nfft_range-1) * c / (2 * Nfft_range * delta_f);   % (m)

% --- 多普勒轴 ---
% 每个OFDM符号（含CP）对应一个"脉冲"，重复周期为 T_sym = N_symbol/fs
% 符号重复频率（等效PRF）= 1/T_sym
% 多普勒bin间距 = PRF/Nfft_dopp，fftshift后范围为 [-PRF/2, +PRF/2)
PRF          = 1 / T_sym;                                             % 符号重复频率 (Hz)
doppler_axis = (-Nfft_dopp/2 : Nfft_dopp/2-1) * (PRF / Nfft_dopp);  % 多普勒频率 (Hz)
velocity_axis = doppler_axis * lambda / 2;                            % 径向速度 (m/s)

% --- 有效带宽（用于分辨率显示）---
B = N_sc * delta_f;

%% ==============================
%% 分辨率信息
%% ==============================
range_res    = c / (2 * B);                         % 距离分辨率 (m)
velocity_res = lambda / (2 * n_sym * T_sym);        % 速度分辨率 (m/s)
fprintf('距离分辨率:  %.2f m\n',  range_res);
fprintf('速度分辨率:  %.2f m/s\n', velocity_res);
fprintf('最大不模糊距离: %.1f m\n',  c / (2 * delta_f));
fprintf('最大不模糊速度: ±%.1f m/s\n\n', lambda * PRF / 4);

%% ==============================
%% 峰值检测（局部最大值 + 阈值）
%% ==============================
threshold = 0.2;    % 相对阈值（相对于全局最大值）
BW_tau    = 1;      % 距离维搜索半径（bin）
BW_fd     = 1;      % 多普勒维搜索半径（bin）

[N_R, N_V] = size(RDM);
peak_r_idx = [];
peak_v_idx = [];

for i = 1+BW_tau : N_R-BW_tau
    for j = 1+BW_fd : N_V-BW_fd
        patch    = RDM(i-BW_tau:i+BW_tau, j-BW_fd:j+BW_fd);
        center   = RDM(i, j);
        neighbor = patch(:);
        neighbor(neighbor == center) = [];
        if all(center > neighbor) && center > threshold
            peak_r_idx(end+1) = i;
            peak_v_idx(end+1) = j;
        end
    end
end

% 按峰值大小降序排列
peak_vals          = RDM(sub2ind(size(RDM), peak_r_idx, peak_v_idx));
[~, sort_idx]      = sort(peak_vals, 'descend');
peak_r_idx         = peak_r_idx(sort_idx);
peak_v_idx         = peak_v_idx(sort_idx);

%% ==============================
%% 目标参数估计
%% ==============================
n_det = min(length(peak_r_idx), 5);
detected(n_det) = struct('range', 0, 'velocity', 0, 'peak', 0, ...
                          'r_idx', 0, 'v_idx', 0);

fprintf('检测到目标数量: %d\n', n_det);
for i = 1:n_det
    detected(i).range    = range_axis(peak_r_idx(i));
    detected(i).velocity = velocity_axis(peak_v_idx(i));
    detected(i).peak     = RDM(peak_r_idx(i), peak_v_idx(i));
    % 在坐标格上找最近索引（用于切片）
    [~, detected(i).r_idx] = min(abs(range_axis    - detected(i).range));
    [~, detected(i).v_idx] = min(abs(velocity_axis - detected(i).velocity));
    fprintf('目标%d: 距离 %5.1f m,  速度 %6.1f m/s\n', ...
        i, detected(i).range, detected(i).velocity);
end
fprintf('\n');

%% ==============================
%% 绘图
%% ==============================

% --- 三维 RDM ---
figure;
mesh(range_axis, velocity_axis, RDM');
xlabel('距离 (m)');  ylabel('速度 (m/s)');  zlabel('归一化幅度');
title('Range-Doppler Map（3D）');
view(45, 30);
axis([range_view velocity_view 0 1]);
hold on;
for i = 1:n_det
    plot3(detected(i).range, detected(i).velocity, detected(i).peak, ...
        'r*', 'MarkerSize', 15, 'LineWidth', 2);
    text(detected(i).range, detected(i).velocity, detected(i).peak * 1.1, ...
        sprintf('目标%d\n(%.1fm, %.1fm/s)', i, detected(i).range, detected(i).velocity), ...
        'Color', 'r', 'FontSize', 10, 'HorizontalAlignment', 'center');
end

% --- 二维 RDM ---
figure;
imagesc(range_axis, velocity_axis, RDM');
set(gca, 'YDir', 'normal');
xlabel('距离 (m)');  ylabel('速度 (m/s)');
title('Range-Doppler Map（2D）');
axis([range_view velocity_view]);
colorbar;
hold on;
for i = 1:n_det
    plot(detected(i).range, detected(i).velocity, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
    text(detected(i).range, detected(i).velocity, ...
        sprintf('  检测%d', i), 'Color', 'r', 'FontSize', 10);
end

% --- 距离切片（取第一目标的速度bin）---
figure;
plot(range_axis, RDM(:, detected(1).v_idx));
xlabel('距离 (m)');  ylabel('归一化幅度');
title('Range Slice');
xlim(range_view);  grid on;

% --- 速度切片（取第一目标的距离bin）---
figure;
plot(velocity_axis, RDM(detected(1).r_idx, :));
xlabel('速度 (m/s)');  ylabel('归一化幅度');
title('Velocity Slice');
xlim(velocity_view);  grid on;

toc;