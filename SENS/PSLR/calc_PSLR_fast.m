function [PSLR_weighted, PSLR_delay, PSLR_doppler] = calc_PSLR_fast(S_fd, tau_vec, fd_vec, fs, rho)
% 计算加权PSLR 
% PSLR = rho * PSLR_delay + (1-rho) * PSLR_doppler

[n_fft, ~] = size(S_fd);
s_t = reshape(ifft(S_fd, n_fft), [], 1);  % 注意reshape
[chi_0delay, chi_0doppler] = calc_ambig_slice(s_t, tau_vec, fd_vec, fs);

%% 计算 PSLR for Zero-Doppler切片 χ(τ, 0)
fprintf('\n--- 计算 Delay PSLR ---\n');

% 找主瓣峰值
[P_peak_delay, peak_idx_delay] = max(chi_0doppler);
mainlobe_indices_delay = find_mainlobe_derivative(chi_0doppler, peak_idx_delay);

% 旁瓣区域: 除主瓣外的所有点
chi_sidelobe_delay = chi_0doppler;
chi_sidelobe_delay(mainlobe_indices_delay) = 0;

% 找最高旁瓣
[P_sidelobe_delay, sidelobe_idx_delay] = max(chi_sidelobe_delay);
fprintf('最高旁瓣位置: τ = %.3f us\n', tau_vec(sidelobe_idx_delay)*10^6);

% 计算PSLR
if P_sidelobe_delay > 0
    PSLR_delay = 10 * log10(P_peak_delay / P_sidelobe_delay);
else
    PSLR_delay = Inf;
    warning('旁瓣为0, PSLR = Inf');
end

%% 计算 PSLR for Zero-Delay切片 χ(0, f_d)
fprintf('\n--- 计算 Doppler PSLR ---\n');

% 找主瓣峰值
[P_peak_doppler, peak_idx_doppler] = max(chi_0delay);
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

%% 计算加权PSLR - 论文公式(41)
PSLR_weighted = rho * PSLR_delay + (1 - rho) * PSLR_doppler;
fprintf('\n========== PSLR结果 ==========\n');
fprintf('Delay PSLR:    %.2f dB (权重 %.0f%%)\n', PSLR_delay, rho*100);
fprintf('Doppler PSLR:  %.2f dB (权重 %.0f%%)\n', PSLR_doppler, (1-rho)*100);
fprintf('Weighted PSLR: %.2f dB\n', PSLR_weighted);
fprintf('==============================\n\n');
end
