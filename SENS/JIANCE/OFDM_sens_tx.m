function [txSig, X_eff] = OFDM_sens_tx(tx_bits)

currentFolder = fileparts(mfilename('fullpath'));
load(fullfile(currentFolder, 'parameter_sens.mat'));

%% ================= 卷积编码 =================

trellis = poly2trellis(L, [133 171]);

tx_code = convenc(tx_bits, trellis);

% 交织编码，深度log2(M)
tx_intlv = matintrlv(tx_code, log2(M), length(tx_code) / log2(M));

%% ================= QPSK调制 =================

tx_bit_mat = reshape(tx_intlv, log2(M), [])';

tx_sym_idx = bi2de(tx_bit_mat);

tx_symbol = pskmod(tx_sym_idx, M, pi/M);

%% ================= 扩频 =================

tx_spread_code = spread_mseq(stage, ptap1, regi1, data_row);

tx_spread_code = tx_spread_code * 2 - 1;

tx_symbol = reshape(tx_symbol, data_row, []);

tx_spread_signal = spread(tx_symbol, tx_spread_code);

tx_spread_signal = reshape(tx_spread_signal, [], 1);

tx_symbol = tx_spread_signal;

%% ================= 插入导频 =================

rng(1);

tx_pilot = 2 * (randi([0 1], pilot_num, data_col * spread_ratio)) - 1;

tx_grid = zeros(N_sc, data_col * spread_ratio);

tx_grid(P_f_station, :) = tx_pilot;

if data_row * data_col > length(tx_symbol)

    tx_symbol_pad = [tx_symbol; zeros(data_row * data_col - length(tx_symbol), 1)];

else

    tx_symbol_pad = tx_symbol;

end

tx_data_mat = reshape(tx_symbol_pad, data_row, data_col * spread_ratio);

tx_grid(data_station, :) = tx_data_mat;

%% ================= 子载波映射 =================

tx_grid = [zeros(pos_row - 1, data_col * spread_ratio); tx_grid];

tx_ifft_input = [tx_grid; zeros(N_fft - N_sc - pos_row + 1, data_col * spread_ratio)];

X_eff = tx_ifft_input(pos_row : pos_row + N_sc - 1 , :);

%% ================= OFDM调制 =================

tx_ifft = ifft(tx_ifft_input);

tx_ofdm_cp = [tx_ifft(N_fft - N_cp + 1:end, :); tx_ifft];

%% ================= PSS同步 =================

tx_pss_code = spread_mseq(PSS_stage, PSS_ptap1, PSS_regi1, 1);

tx_pss_code = tx_pss_code * 2 - 1;

tx_pss = sqrt(1/N_fft) * tx_pss_code(:);

%% ================= 串并转换 =================

tx_serial = reshape(tx_ofdm_cp, [], 1);

tx_zeros = zeros(N_zeros, 1);

txSig = [tx_zeros; tx_pss; tx_serial];

txSig = txSig / max(abs(txSig));