function [chi, chi_zero_delay, chi_zero_doppler] = calc_ambig(s_t, tau_vec, fd_vec, fs)
s_t = reshape(s_t, [], 1);
N_total = length(s_t);
t = (0:N_total-1).' / fs;

fprintf('生成时域信号: N_total=%d, T_total=%.3f ms\n', N_total, N_total/fs*1e3);

%% ========== 参数准备 ==========
N_tau = length(tau_vec);
N_fd = length(fd_vec);
chi = zeros(N_tau, N_fd);

% 找到零时延和多普勒的索引（用于切片提取）
[~, zero_tau_idx] = min(abs(tau_vec));
[~, zero_fd_idx] = min(abs(fd_vec));

% 预计算exp(-j2πf_d*t)项
exp_fd = zeros(N_total, N_fd);
for j = 1:N_fd
    exp_fd(:, j) = exp(-1j * 2 * pi * fd_vec(j) * t);
end

fprintf('计算模糊函数: %d x %d 网格...\n', N_tau, N_fd);

%% ========== 计算完整模糊函数 ==========
for i = 1:N_tau
    tau = tau_vec(i);
    delay_samples = tau * fs;     
    n = (0:N_total-1).';
    s_delayed = interp1(n, s_t, n - delay_samples, 'spline', 0);

    % 计算该时延下所有多普勒点的模糊函数
    for j = 1:N_fd
        % 公式(38): χ(τ, f_d) = |∫ s(t) · s*(t-τ) · exp(-j2πf_d·t) dt|
        integrand = s_t .* conj(s_delayed) .* exp_fd(:, j);
        
        % 数值积分
        integral_value = trapz(t, integrand);
        chi(i, j) = abs(integral_value);
    end
    
    % 进度显示
    if mod(i, max(1, floor(N_tau/10))) == 0
        %fprintf('  进度: %d/%d (%.1f%%)\n', i, N_tau, 100*i/N_tau);
    end
end

%% ========== 提取切片（基于计算好的chi矩阵） ==========
% 直接从已计算的chi矩阵中提取，确保一致性
chi_zero_delay = chi(zero_tau_idx, :);      % τ=0 行
chi_zero_doppler = chi(:, zero_fd_idx);     % f_d=0 列

%% ========== 统一归一化 ==========
% 以完整模糊函数的峰值为基准进行归一化
max_chi = max(chi(:));
if max_chi > 0
    chi = chi / max_chi;
    chi_zero_delay = chi_zero_delay / max_chi;
    chi_zero_doppler = chi_zero_doppler / max_chi;
else
    warning('模糊函数全为零！');
end

%% ========== 峰值位置分析 ==========
[max_val, max_idx] = max(chi(:));
[max_row, max_col] = ind2sub(size(chi), max_idx);
peak_tau = tau_vec(max_row);
peak_fd = fd_vec(max_col);

fprintf('模糊函数计算完成!\n');
fprintf('  峰值位置: [τ, f_d] = [%.3e s, %.3f Hz]\n', peak_tau, peak_fd);
fprintf('  峰值大小: %.3e\n', max_val);

% 检查主峰是否在原点（对于自模糊函数应该如此）
if abs(peak_tau) > 1e-9 || abs(peak_fd) > 1e-9
    fprintf('  注意：主峰不在原点！τ=%.2e s, f_d=%.2f Hz\n', peak_tau, peak_fd);
    fprintf('  这可能表明：\n');
    fprintf('    1. 信号s(t)本身包含固有延迟/频偏\n');
    fprintf('    2. 计算的是互模糊函数\n');
    fprintf('    3. 数值误差\n');
end

end