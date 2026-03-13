%% ============================SNR_TIME_DOMAIN==================================
%  OFDM_proc_rx.m
%  功能：
%  接收带扩频与信道编码的 OFDM 信号，并完成解调与 BER 计算。
%
%  主要流程：
%  3. PSS 前导相关同步
%  4. OFDM 符号定位
%  5. 去循环前缀
%  6. FFT
%  7. 导频信道估计
%  8. 频域均衡
%  9. 解扩
% 10. QPSK 解调
% 11. 维特比译码
% ==============================================================

function [rx_data,data_stream,Y] = OFDM_proc_rx(rxSig)


currentFolder = fileparts(mfilename('fullpath'));
load(fullfile(currentFolder, 'parameter_ttt.mat'));

len_rx = length(rxSig);

%%
rng(1);

pilot_seq = 2 * (randi([0 1], pilot_num, data_col * spread_ratio)) - 1;

%% ================= 扩频码 =================

code = spread_mseq(stage, ptap1, regi1, data_row);

code = code * 2 - 1;


%% ================= 卷积编码 =================

trellis = poly2trellis(L, [133 171]);


%% ================= PSS序列 =================

PSS_code = spread_mseq(PSS_stage, PSS_ptap1, PSS_regi1, 1);

PSS_code2 = PSS_code * 2 - 1;

PSS_transpose = PSS_code2.';


%% ================= PSS同步 =================


one_frame_len = N_zeros + length(PSS_code2) + N_symbol * (data_col * spread_ratio);
pss_vector = 1:length(PSS_code2);
search_range = min(2 * one_frame_len, len_rx - length(pss_vector));
auto_correlation_vector = zeros(search_range, 1);

for delay = 0:len_rx-length(pss_vector)

    window = delay + pss_vector;

    auto_correlation_vector(delay+1) = ...
        PSS_transpose.' * (rxSig(window));

end

[~,p1] = max(abs(auto_correlation_vector));

temp = auto_correlation_vector;

temp(p1) = 0;

[~,p2] = max(abs(temp));

p = min(p1,p2);

%% ================= OFDM符号定位 =================

pos_cp1 = p + length(PSS_code2) - 1;

pos_cp = pos_cp1 : N_symbol : ...
    (pos_cp1 + N_symbol * (data_col * spread_ratio - 1));
max_needed = pos_cp(end) + N_cp + N_fft - 1;
if max_needed > len_rx
    rxSig = [rxSig; zeros(max_needed - len_rx, 1)];
    len_rx1 = length(rxSig);
end

%% ================= CFO估计（利用CP相关）=================
% 利用CP与OFDM符号尾部的重复结构估计频偏
cfo_est = zeros(data_col * spread_ratio, 1);
for i = 1:data_col * spread_ratio
    % 取CP段和对应符号尾部
    cp_part  = rxSig(pos_cp(i) : pos_cp(i) + N_cp - 1);
    tail_part = rxSig(pos_cp(i) + N_fft : pos_cp(i) + N_fft + N_cp - 1);
    % 相位差估计频偏
    cfo_est(i) = angle(cp_part' * tail_part) / (2 * pi * N_fft);
end
cfo_mean = mean(cfo_est);

%% ================= CFO补偿（时域） =================
n = (0 : len_rx1 - 1).';
rxSig = rxSig .* exp(-1j * 2 * pi * cfo_mean * n);
%% ================= 去CP =================

dcp = zeros(N_fft, data_col * spread_ratio);

for i = 1:data_col * spread_ratio

    dcp_temp = rxSig( ...
        pos_cp(i) + N_cp : ...
        pos_cp(i) + N_cp + N_fft - 1);

    dcp(:,i) = dcp_temp;

end

%% ================= FFT =================

fft_data = fft(dcp, N_fft);

data3 = fft_data(pos_row : pos_row + data_row + pilot_num - 1, :);
Y = fft_data;
%% ================= 信道估计 =================

Rx_pilot = data3(P_f_station,:);

h = Rx_pilot ./ pilot_seq;

H = interp1(P_f_station', h, data_station', 'linear','extrap');

%% ================= 均衡 =================

data_aftereq = data3(data_station,:) ./ H;

data_aftereq0 = data_aftereq;

%% ================= 并串 =================

data_stream = reshape(data_aftereq,[],1);

%% ================= 解扩 =================

demspread_data = despread(data_aftereq0, code);

demspread_data_row = reshape(demspread_data, [], 1);

%% ================= QPSK解调 =================

demodulation_data = pskdemod(demspread_data_row, M, pi/M);

De_data1 = reshape(demodulation_data, [], 1);

De_data2 = de2bi(De_data1);

De_Bit = reshape(De_data2', 1, []);

%% ================= 解交织 =================

y_deintlv = matdeintrlv(De_Bit, log2(M), length(De_Bit) / log2(M));

%% ================= 维特比译码 =================

rx_c_de = vitdec(y_deintlv, trellis, tblen, 'trunc', 'hard');
rx_data = rx_c_de;

end