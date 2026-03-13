function Y_eff = OFDM_sens_rx(rxSig)
% OFDM_radar_rx.m
% 雷达感知专用接收函数
% 与 OFDM_proc_rx.m 共用同一套同步/解调思路，但在FFT后直接输出原始
% 频域矩阵 Y_eff，不做任何CFO补偿或信道均衡，以保留雷达所需的完整相位信息。
%
% 输入：
%   rxSig   - 接收时域信号（列向量）
% 输出：
%   Y_eff   - N_sc × (data_col*spread_ratio) 原始频域矩阵
%             与发射端 S_eff 子载波索引完全对齐，可直接做 H_est = Y_eff ./ S_eff

currentFolder = fileparts(mfilename('fullpath'));
load(fullfile(currentFolder, 'parameter_sens.mat'));

len_rx = length(rxSig);

%% ================= PSS序列 =================
% 与 OFDM_proc_rx 完全一致
rx_pss_code = spread_mseq(PSS_stage, PSS_ptap1, PSS_regi1, 1);
rx_pss_code = rx_pss_code * 2 - 1;
rx_pss_row  = rx_pss_code.';

%% ================= PSS同步 =================
% 与 OFDM_proc_rx 完全一致
one_frame_len = N_zeros + length(rx_pss_code) + N_symbol * (data_col * spread_ratio);
pss_index     = 1:length(rx_pss_code);
search_range  = min(2 * one_frame_len, len_rx - length(pss_index));

rx_corr = zeros(search_range, 1);
for delay = 0:search_range-1
    rx_window       = delay + pss_index;
    rx_corr(delay+1) = rx_pss_row.' * rxSig(rx_window);
end

[~, p1] = max(abs(rx_corr));
rx_temp      = rx_corr;
rx_temp(p1)  = 0;
[~, p2]      = max(abs(rx_temp));
p            = min(p1, p2);

%% ================= OFDM符号定位 =================
% 与 OFDM_proc_rx 完全一致
rx_pos_cp_start = p + length(rx_pss_code) - 1;
rx_pos_cp = rx_pos_cp_start : N_symbol : ...
    (rx_pos_cp_start + N_symbol * (data_col * spread_ratio - 1));

max_needed = rx_pos_cp(end) + N_cp + N_fft - 1;
if max_needed > len_rx
    rxSig  = [rxSig; zeros(max_needed - len_rx, 1)];
end

%% ================= 去CP =================
% 与 OFDM_proc_rx 完全一致
rx_dcp = zeros(N_fft, data_col * spread_ratio);
for i = 1:data_col * spread_ratio
    rx_dcp(:, i) = rxSig( ...
        rx_pos_cp(i) + N_cp : ...
        rx_pos_cp(i) + N_cp + N_fft - 1);
end

%% ================= FFT =================
% 与 OFDM_proc_rx 完全一致
rx_fft  = fft(rx_dcp, N_fft);

Y_eff = rx_fft(pos_row : pos_row + N_sc - 1, :);

end