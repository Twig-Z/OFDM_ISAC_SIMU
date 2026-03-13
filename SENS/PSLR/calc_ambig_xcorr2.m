function [chi, chi_zero_delay, chi_zero_doppler] = calc_ambig_xcorr2(s_t, tau_vec, fd_vec, fs)

%% ========= 信号整理 =========
s_t = reshape(s_t, [], 1);
N_total = length(s_t);

t = (0:N_total-1).' / fs;

fprintf('生成时域信号: N_total=%d, T_total=%.3f ms\n', N_total, N_total/fs*1e3);

%% ========= 参数 =========
N_tau = length(tau_vec);
N_fd  = length(fd_vec);

chi = zeros(N_tau, N_fd);

%% ========= 找零点索引 =========
[~, zero_tau_idx] = min(abs(tau_vec));
[~, zero_fd_idx]  = min(abs(fd_vec));

%% ========= 多普勒调制 =========
fprintf('生成多普勒调制信号...\n');

doppler_signal = zeros(N_total, N_fd);

for j = 1:N_fd
    doppler_signal(:,j) = s_t .* exp(-1j*2*pi*fd_vec(j)*t);
end

%% ========= 计算相关 =========
fprintf('计算模糊函数...\n');

for j = 1:N_fd

    % 与原信号做互相关
    r = xcorr(doppler_signal(:,j), s_t);

    % xcorr对应的延迟
    lags = (-N_total+1:N_total-1).' / fs;

    % 插值到 tau_vec
    chi(:,j) = abs(interp1(lags, r, tau_vec, 'linear', 0));

end

%% ========= 提取切片 =========
chi_zero_delay   = chi(zero_tau_idx, :);
chi_zero_doppler = chi(:, zero_fd_idx);

%% ========= 归一化 =========
max_chi = max(chi(:));

if max_chi > 0
    chi = chi / max_chi;
    chi_zero_delay = chi_zero_delay / max_chi;
    chi_zero_doppler = chi_zero_doppler / max_chi;
else
    warning('模糊函数全为零！');
end

%% ========= 峰值位置 =========
[max_val, max_idx] = max(chi(:));
[row, col] = ind2sub(size(chi), max_idx);

peak_tau = tau_vec(row);
peak_fd  = fd_vec(col);

fprintf('模糊函数计算完成!\n');
fprintf('峰值位置: τ=%.3e s, fd=%.3f Hz\n', peak_tau, peak_fd);
fprintf('峰值大小: %.3e\n', max_val);

end