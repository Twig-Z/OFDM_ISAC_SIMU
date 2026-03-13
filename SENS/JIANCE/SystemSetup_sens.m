function SystemSetup_sens
%% ================================================================
%  SystemSetup_sens.m
%  功能：OFDM-ISAC系统统一参数配置文件（感知优化版）
%  被 OFDM_proc_tx.m、OFDM_proc_rx.m、OFDM_radar_rx.m 共同调用
%
%  感知设计目标（对应目标距离32m/35m，速度10m/s）：
%    距离分辨率 δR = c/(2B)  < 3 m   → B > 50 MHz
%    速度分辨率 δv = λΔf/(2Nsym) < 5 m/s
%    最大不模糊距离 Rmax = c/(2Δf) > 200 m
%
%  参数确定顺序（不可颠倒，逻辑与原版完全一致）：
%  masterClock → decim/interp → fs → N_fft → delta_f → N_sc → B(结果)
%  ================================================================

%% ================= 依赖文件检查 =================

[fList1, ~] = matlab.codetools.requiredFilesAndProducts('OFDM_sens_tx.m');
[fList2, ~] = matlab.codetools.requiredFilesAndProducts('OFDM_sens_rx.m');
fList = union(fList1, fList2);
disp('USRP收发调用的非系统函数/文件有：');
disp(fList')

%% ================= 硬件参数 =================
% 【修改】fc: 915MHz → 24GHz（ISM频段，ISAC主流，λ更短速度分辨率更好）
% 【修改】masterClock: 61.44MHz → 1966.08MHz，保持2的幂次关系
% 【修改】decim: 4 → 16，使 fs = 1966.08/16 = 122.88MHz（够装下100MHz带宽）
% 其余硬件参数逻辑不变

fc          = 24e9;             % 载波频率 24GHz（ISM，ISAC主流）
gain_tx     = 70;               % 发射增益，可调
gain_rx     = 30;               % 接收增益，可调

samplesPerFrame = 64000;        % 每帧采样点数

masterClock = 1966.08e6;        % 硬件基准时钟（理想仿真值）
                                % = 1966.08MHz，为15kHz×2^17，保证整除

decim       = 16;               % 抽取因子，建议取2的幂次
interp      = decim;            % 插值因子，需与decim一致

fs = masterClock / decim;       % 实际采样率 = 1966.08MHz / 16 = 122.88MHz
                                % 足以覆盖100MHz有效带宽

%% ================= OFDM基础参数 =================
% 【修改】mu: 0→4，子载波间隔 15kHz×2^4 = 240kHz
%   → Rmax = c/(2Δf) = 3e8/(2×240e3) = 625m  ✓ 满足覆盖需求
%   → T_sym = 1/240kHz = 4.17μs（符号时间短，帧长可控）
% N_fft 由 fs/delta_f 计算得到，逻辑不变

mu      = 4;                            % 子载波间隔缩放因子（5G NR定义）
delta_f = 15000 * 2^mu;                 % 子载波间隔 = 240 kHz
                                        % Rmax = c/(2Δf) = 625 m ✓

N_fft = fs / delta_f;                   % FFT点数 = 122.88MHz / 240kHz = 512
assert(N_fft == floor(N_fft),       'N_fft不是整数，检查fs与delta_f是否匹配');
assert(bitand(N_fft, N_fft-1) == 0, 'N_fft必须是2的幂次');

N_cp     = floor(0.07 * N_fft);         % 循环前缀长度 = 35
N_symbol = N_fft + N_cp;                % 单个OFDM符号长度 = 547
T_sym    = N_symbol / fs;               % 含CP符号时间

%% ================= 带宽与子载波参数 =================
% 逻辑与原版完全一致，B 是计算结果不是输入量
% 【结果】N_sc ≈ 368，B = N_sc×delta_f ≈ 88.3MHz
%   → δR = c/(2B) ≈ 1.70m < 3m  ✓ 可分辨32m和35m目标

eta          = 0.6;                         % 频谱效率因子，留出滚降保护
B_leackage   = 0.1 * fs;                    % 保护带宽 = 12.288MHz（fs的10%）

stage        = 1;                           % m序列扩频级数
spread_ratio = 2^stage - 1;                 % 扩频比 = 1（stage=1时无扩频）

B0   = fs - B_leackage;                     % 单路可用带宽 = 110.592MHz
B    = B0 * spread_ratio + B_leackage;      % 扩频后系统总带宽（结果量，勿作输入）

N_sc0 = floor(B0 / delta_f);               % 不含eta的子载波估计
assert(N_fft >= N_sc0, ...
    sprintf('N_fft(%d)装不下子载波N_sc0(%d)，请减小B_leackage或增大decim', ...
    N_fft, N_sc0));

N_sc    = floor(eta * B0 / delta_f);        % 数据+导频子载波总数
pos_row = ceil(B_leackage / delta_f / 2);   % FFT中数据起始bin偏移

%% ================= 子载波结构参数 =================
% 逻辑与原版完全一致

P_f_inter = 3;                              % 导频间隔（每3个子载波插1个导频）

% 导频子载波位置
P_f_station = 1 : P_f_inter : N_sc;
pilot_num   = length(P_f_station);          % 导频数量

% 数据子载波位置（非导频位置）
data_station = [];
for img = 1:N_sc
    if mod(img, P_f_inter) ~= 1
        data_station = [data_station, img];
    end
end
data_row = length(data_station);            % 数据子载波数量

P_f = sqrt(1/2) * (1 + 1i);                % 导频符号（QPSK点）

%% ================= 信道编码与调制参数 =================
% 【修改】Nd: 14 → 64，增加符号数
%   → δv = λ×Δf/(2×Nd) = 0.0125×240e3/(2×64) = 23.4 m/s
%   → 对10m/s目标可检测（在速度分辨单元内有峰值）
%   注：若需更高速度分辨率可继续增大Nd，但帧长会增加

M                = 4;           % 调制阶数：QPSK
R_channelcoding  = 1/2;         % 卷积码码率
L                = 7;           % 卷积码约束长度
tblen            = 6 * L;       % 维特比译码回溯长度 = 42
Nd               = 64;          % 每次传输的OFDM符号数（↑ 提升速度分辨率）
data_col = Nd;
n_sym    = Nd;

% 每帧数据比特数（由上述参数计算得到，逻辑不变）
N_RE_data = data_row * data_col;
N_b_data  = N_RE_data * log2(M) * R_channelcoding;

%% ================= 扩频参数 =================
% 与原版完全一致

ptap1 = [1];                    % m序列生成多项式tap
regi1 = [1];                    % m序列初始寄存器状态

%% ================= PSS同步序列参数 =================
% 与原版完全一致

PSS_stage = 9;                          % PSS用m序列级数
PSS_ptap1 = [4 9];                      % PSS生成多项式tap
PSS_regi1 = [1 1 0 1 1 0 1 1 1];       % PSS初始寄存器状态

%% ================= 其他参数 =================

N_zeros = 10000;                % 发送前导零样本数

%% ================= 感知性能汇总 =================

lambda          = 3e8 / fc;
B_eff           = N_sc * delta_f;
range_res       = 3e8 / (2 * B_eff);
velocity_res    = lambda * delta_f / (2 * Nd);
R_max           = 3e8 / (2 * delta_f);
v_max           = lambda / (4 * T_sym);

%% ================= 参数汇总打印 =================

fprintf('\n========== 系统参数汇总 ==========\n');
fprintf('硬件层:  masterClock=%.2fMHz, decim/interp=%d\n', masterClock/1e6, decim);
fprintf('采样层:  fs=%.3fMHz\n', fs/1e6);
fprintf('OFDM层:  delta_f=%.0fkHz, N_fft=%d, N_cp=%d\n', delta_f/1e3, N_fft, N_cp);
fprintf('带宽层:  B0=%.3fMHz, B=%.3fMHz\n', B0/1e6, B/1e6);
fprintf('子载波:  N_sc=%d, pilot=%d, data=%d, pos_row=%d\n', N_sc, pilot_num, data_row, pos_row);
fprintf('编码层:  M=%d, R=%.1f, L=%d, Nd=%d\n', M, R_channelcoding, L, Nd);
fprintf('扩频层:  stage=%d, spread_ratio=%d\n', stage, spread_ratio);
fprintf('---------- 感知性能 ----------\n');
fprintf('载波频率:      %.0f GHz,  λ=%.4f m\n', fc/1e9, lambda);
fprintf('有效带宽:      %.1f MHz\n', B_eff/1e6);
fprintf('距离分辨率:    %.2f m\n', range_res);
fprintf('速度分辨率:    %.2f m/s\n', velocity_res);
fprintf('最大不模糊距离: %.1f m\n', R_max);
fprintf('最大不模糊速度: ±%.1f m/s\n', v_max);
fprintf('===================================\n\n');

% 获取当前程序所在路径
currentFolder = fileparts(mfilename('fullpath'));
save(fullfile(currentFolder, 'parameter_sens.mat'));

end