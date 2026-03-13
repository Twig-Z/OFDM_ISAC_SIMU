function [chi_zero_delay, chi_zero_doppler] = calc_ambig_slice(s_t, tau_vec, fd_vec, fs)

s_t     = reshape(s_t, [], 1);
N_total = length(s_t);
t       = (0:N_total-1).' / fs;   % 列向量 [N x 1]

%% ================= Zero-Doppler 切片 χ(τ,0) =================
% χ(τ,0) = | ∫ s(t) · s*(t-τ) dt |
%
% 实现方法：
% 使用向量化分数时延（fractional delay）实现 s(t-τ) 的计算。
% 对所有 τ 一次性构造插值查询矩阵，通过 interp1 完成分数延迟，
% 然后进行向量化积分计算。

n = (0:N_total-1).';                 % 时间采样索引 [N x 1]

delay_samps = tau_vec(:).' * fs;    % τ对应的延迟采样数 [1 x N_tau]

% 查询索引矩阵：每列对应一个 τ
% query_idx(n,i) = n - delay_samps(i)
query_idx = n - delay_samps;         % [N x N_tau]

% 对 s(t) 进行分数延迟插值
s_delayed_mat = interp1(n, s_t, query_idx, 'spline', 0);   % [N x N_tau]

% 向量化积分
% sum( s(t) · conj(s(t-τ)) ) · dt
chi_zero_doppler = abs(sum(s_t .* conj(s_delayed_mat), 1)).' / fs;   % [N_tau x 1]


%% ================= Zero-Delay 切片 χ(0,f_d) =================
% χ(0,f_d) = | ∫ s(t)·s*(t) · e^{-j2πf_d t} dt |
%          = | ∫ |s(t)|^2 · e^{-j2πf_d t} dt |
% 通过构造指数矩阵实现向量化DFT计算。
s_power = abs(s_t).^2;                               % 信号功率 |s(t)|² [N x 1]

exp_matrix = exp(-1j * 2*pi * t * fd_vec(:).');      % 指数矩阵 [N x N_fd]

chi_zero_delay = abs(s_power.' * exp_matrix) / fs;   % [1 x N_fd]
%% ================= 统一归一化 =================
% 以 τ=0 处的自相关峰值为基准（xcorr在lag=0处最大），
% 与 Method2 以完整2D矩阵峰值归一化等价（自模糊函数峰值在原点）
peak_val = abs(sum(s_t .* conj(s_t))) / fs;   % = xcorr(τ=0)/fs
if peak_val > 0
    chi_zero_doppler = chi_zero_doppler / peak_val;
    chi_zero_delay   = chi_zero_delay   / peak_val;
else
    warning('calc_ambig_fast: 切片全为零！');
end

end