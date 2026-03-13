function [PSLR_weighted, PSLR_delay, PSLR_doppler] = calc_PSLR_weighted(S_fd, tau_vec, fd_vec, fs, draw_flag, rho)
% 计算加权PSLR (论文公式41)
% PSLR = rho * PSLR_delay + (1-rho) * PSLR_doppler
%
% 其中:
%   PSLR_delay   = 10*log10(P_p(χ(τ,0)) / P_s(χ(τ,0)))  论文公式(41)上半部分
%   PSLR_doppler = 10*log10(P_p(χ(0,f_d)) / P_s(χ(0,f_d)))  论文公式(41)下半部分
%   P_p(·) = 主瓣峰值
%   P_s(·) = 最高旁瓣峰值
%
% 输入:

%   S_fd: 频域符号矩阵 [n_fft x n_sym]
%   fs: 采样率
%   n_fft, n_cp: OFDM参数
%   rho: 权重系数 (0~1)
%        rho=1   -> 只考虑距离性能
%        rho=0   -> 只考虑速度性能
%        rho=0.5 -> 均衡
%
% 输出:
%   PSLR_weighted: 加权PSLR (dB)
%   PSLR_delay, PSLR_doppler: 各自的PSLR (dB)
%   chi, chi_0delay, chi_0doppler: 模糊函数及其切片(用于可视化)

if nargin < 6
    rho = 0.5;  % 默认均衡
end
fprintf('rho = %.2f (%.0f%%距离 + %.0f%%速度)\n', rho, rho*100, (1-rho)*100);
[n_fft, ~] = size(S_fd);
%s_t = ifft(S_fd, n_fft);
s_t = reshape(ifft(S_fd, n_fft), [], 1);  % 注意reshape
%% 计算模糊函数
fprintf('\n计算模糊函数...\n');


if draw_flag
    % 画图时才算完整模糊函数（慢）
    [chi, chi_0delay, chi_0doppler] = calc_ambig(s_t, tau_vec, fd_vec, fs);
else
    % 优化时只算两个切片（快）
    [chi_0delay, chi_0doppler] = calc_ambig_slice(s_t, tau_vec, fd_vec, fs);
    chi = [];
end
%% 计算 PSLR for Zero-Doppler切片 χ(τ, 0) - 论文公式(41)上半部分
fprintf('\n--- 计算 Delay PSLR ---\n');

% 找主瓣峰值
[P_peak_delay, peak_idx_delay] = max(chi_0doppler);
fprintf('主瓣峰值位置: τ = %.3f us\n', tau_vec(peak_idx_delay)*10^6);
mainlobe_indices_delay = find_mainlobe_derivative(chi_0doppler, peak_idx_delay);

% 旁瓣区域: 除主瓣外的所有点
chi_sidelobe_delay = chi_0doppler;
chi_sidelobe_delay(mainlobe_indices_delay) = 0;

% 找最高旁瓣
[P_sidelobe_delay, sidelobe_idx_delay] = max(chi_sidelobe_delay);
fprintf('最高旁瓣位置: τ = %.3f us\n', tau_vec(sidelobe_idx_delay)*10^6);

% 计算PSLR (公式41)
if P_sidelobe_delay > 0
    PSLR_delay = 10 * log10(P_peak_delay / P_sidelobe_delay);
else
    PSLR_delay = Inf;
    warning('旁瓣为0, PSLR = Inf');
end

fprintf('Delay PSLR = %.2f dB\n', PSLR_delay);

%% 计算 PSLR for Zero-Delay切片 χ(0, f_d) - 论文公式(41)下半部分
fprintf('\n--- 计算 Doppler PSLR ---\n');

% 找主瓣峰值
[P_peak_doppler, peak_idx_doppler] = max(chi_0delay);
fprintf('主瓣峰值位置: f_d = %.2f Hz\n', fd_vec(peak_idx_doppler))
mainlobe_indices_doppler = find_mainlobe_derivative(chi_0delay, peak_idx_doppler);

% 旁瓣区域
chi_sidelobe_doppler = chi_0delay;
chi_sidelobe_doppler(mainlobe_indices_doppler) = 0;

% 找最高旁瓣
[P_sidelobe_doppler, sidelobe_idx_doppler] = max(chi_sidelobe_doppler);
fprintf('最高旁瓣位置: f_d = %.2f Hz\n', fd_vec(sidelobe_idx_doppler));

% 计算PSLR
if P_sidelobe_doppler > 0
    PSLR_doppler = 10 * log10(P_peak_doppler / P_sidelobe_doppler);
else
    PSLR_doppler = Inf;
    warning('旁瓣为0, PSLR = Inf');
end

fprintf('Doppler PSLR = %.2f dB\n', PSLR_doppler);

%% 计算加权PSLR - 论文公式(41)
PSLR_weighted = rho * PSLR_delay + (1 - rho) * PSLR_doppler;

fprintf('\n========== PSLR结果 ==========\n');
fprintf('Delay PSLR:    %.2f dB (权重 %.0f%%)\n', PSLR_delay, rho*100);
fprintf('Doppler PSLR:  %.2f dB (权重 %.0f%%)\n', PSLR_doppler, (1-rho)*100);
fprintf('Weighted PSLR: %.2f dB\n', PSLR_weighted);
fprintf('==============================\n\n');

%% 可视化(可选)
% 如果draw_flag=1,自动绘图
if draw_flag
    visualize_PSLR(tau_vec, fd_vec, chi, chi_0delay, chi_0doppler, ...
                   mainlobe_indices_delay, mainlobe_indices_doppler, ...
                   peak_idx_delay, peak_idx_doppler, ...
                   sidelobe_idx_delay, sidelobe_idx_doppler);
end

end
