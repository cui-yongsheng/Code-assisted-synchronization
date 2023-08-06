%% 主执行脚本：LDPC编码(CCSDS推荐码字，不考虑打孔)-QPSK(格雷码调制)
clc
clear
close all
rng default
addpath(genpath('../../functions'));
addpath(genpath('../../dependences'));
addpath(genpath('../../datas'));
addpath(genpath('../../mex'));

%% 测试控制参数
if_FreqPhaseErr = 1;                        % 是否引入相偏及频偏
if_timingErr = 1;                           % 是否引入定时偏差
FreqErr= 2e-4;                              % 测试频率误差
PhaseErr= 160;                               % 测试相位误差
timingErr = 2;                              % 附加时延（样本点数）
EbN0_scale = 2;                   % EbN0测试范围  (理想同步时，码率1/2译码门限2dB)
nFrame = 100;                                % 测试帧数

%% 信号调制及编码参数
M=4;                            % 调制阶数 (QPSK)
sps = 4;                        % 单位符号采样数
rolloff = 0.35;                 % 根升余弦滤波器滚降系数
span = 6;                       % 根升余弦滤波器阶数
code_r=1/6;                     % LDPC编码效率 (1/2,1/4)
code_K=1024;                    % LDPC编码信息位数量
code_M = 512;                   % CCSDS标准规定的校验矩阵参数（非必要不修改）
H_src=load('LDPC_非规则码.mat');
H=H_src.H;                      % LDPC校验矩阵
code_N=size(H,2);               % LDPC码字长度
code_r_real=code_K/code_N;      % 实际LDPC编码效率（保留打孔位导致）
L = code_N/log2(M);             % 调制符号数量

%% 通信工具箱对象实例化
ldpcEnc = comm.LDPCEncoder(H);                                  % 发端-LDPC编码器
qpskMod = comm.QPSKModulator('BitInput',true);                  % 发端-QPSK调制器
txfilter = comm.RaisedCosineTransmitFilter...
    ('RolloffFactor',rolloff,'FilterSpanInSymbols',....
    span,'OutputSamplesPerSymbol',sps);                         % 发端-根升余弦滚降滤波器
fixedDelay = dsp.Delay(timingErr);                              % 发端-附加时延模块
freqOffset = comm.PhaseFrequencyOffset...
    ('FrequencyOffset',FreqErr,....
    'PhaseOffset',PhaseErr,'SampleRate',4);                    % 发端-附加频偏、相位模块
rxfilter = comm.RaisedCosineReceiveFilter...
    ('RolloffFactor',rolloff,'FilterSpanInSymbols',span,...
    'DecimationFactor',1,'InputSamplesPerSymbol',sps);          % 收端-根升余弦滚降滤波器
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % 收端-QPSK解调器(软解调，并根据输入的噪声功率，
    'DecisionMethod','Approximate log-likelihood ratio', ...    % 在假定复信号功率为1的条件下进行LLR归一化)
    'VarianceSource','Input port');
ldpcDec = comm.LDPCDecoder('ParityCheckMatrix',H,'DecisionMethod','Soft decision','MaximumIterationCount',20,'OutputValue','Whole codeword');

%% 不同EbN0下循环测试
BER_test=zeros(1,length(EbN0_scale));              % 译码误码率
for i=1:length(EbN0_scale)
    disp(['####### 当前 EbN0 ' num2str(EbN0_scale(i)) ' dB  #######']);
    %% 发端生成信号
    QPSK_frame_length=code_N/log2(M);                       % QPSK符号帧长度
    QPSK_symbol_frames=zeros(QPSK_frame_length,nFrame);     % 复数形式的QPSK符号帧集合
    src_bit_frames=logical(randi([0 1],code_K,nFrame));     % 随机生成信源帧集合
    code_bit_frames=zeros(code_N,nFrame);                   % LDPC编码后的信源帧集合
    for frame_index=1:nFrame
        code_bit_frames(:,frame_index)=ldpcEnc(src_bit_frames(:,frame_index));      % LDPC编码后的信源帧 (ToDo: 目前暂不考虑帧同步的实现)
        QPSK_symbol_frames(:,frame_index)=qpskMod(code_bit_frames(:,frame_index));  % 编码信源帧映射为复数形式的QPSK符号帧
    end
    QPSK_symbol_stream=reshape(QPSK_symbol_frames,[],1);    % 复数形式的QPSK符号流
    QPSK_symbol_stream=[QPSK_symbol_stream;zeros(span,1)];  % 符号流中加入全零尾片段（克服升余弦滤波产生的尾符号丢失）
    txSig=txfilter(QPSK_symbol_stream);                     % 脉冲成形，产生发端基带信号
    disp('发送端基带信号已生成...')
    
    %% 引入定时偏差及多普勒
    if(if_FreqPhaseErr==1)
        txSig = freqOffset(txSig);
    end
    if(if_timingErr==1)
        txSig = TEmodel(txSig,timingErr);
    end
    
    %% 引入信道噪声
    EbN0=EbN0_scale(i);                           % 归一化编码前比特能量Eb/N0
    EbN0_code=EbN0+10*log10(code_r_real);         % 归一化编码后比特能量Eb_code/N0
    EsN0=10*log10(2)+EbN0_code;                   % 归一化调制符号能量Es/N0 ( QPSK调制 )
    SNR=EsN0-10*log10(sps);                       % 位同步前的采样序列信噪比
    rxSig=awgn(txSig,SNR,'measured');
    symbol_noise_var=10^(-1*EsN0/10);
    disp('信道噪声已引入...')
    
    %% 收端测试(matlab工具箱载波同步和定时同步模块信噪比工作范围要求QPSK的EsN0高于3dB，此外需要导频克服相位模糊)
    rxSig=rxfilter(rxSig);                                  % 匹配滤波                   
    rxSig=rxSig(span*sps+1:end);                                % 去除匹配滤波与定时同步拖尾                     
    rxQPSK_symbol_stream=reshape(rxSig,[],nFrame);

    %% 收端测试DNN+EM
    BER=zeros(1,nFrame);
    timePre=ones(1,nFrame);
    phasePre=ones(1,nFrame);
    freqPre=ones(1,nFrame);
    Est_phase_d = ones(2,nFrame);
    Est_freq_d = ones(2,nFrame);
    Est_time_d = ones(2,nFrame);
    Est_phase_v = ones(21,nFrame);
    Est_freq_v = ones(21,nFrame);
    Est_time_v = ones(21,nFrame);
    tic
    for frame_id=1:nFrame
        if(mod(frame_id,1000)==0)
            disp(['####### 第 ' num2str(frame_id) ' 帧  #######']);
        end
        rxFrame_checked_sample=rxQPSK_symbol_stream(:,frame_id);
%         %添加相频偏移
%         Error_init_Iterator=comm.PhaseFrequencyOffset('PhaseOffset',PhaseErr,'FrequencyOffset',FreqErr,'SampleRate',sps);
%         rxFrame_checked=Error_init_Iterator(rxFrame_checked_sample);
        % 相频时延预搜索
        % Todo:增加候选列表机制
        [phasePre(frame_id),freqPre(frame_id),timePre(frame_id)]=phaseFreqTimeSearch(rxFrame_checked_sample,symbol_noise_var,H,sps);
        %使用预搜索结果修正
        Error_pre_Iterator=comm.PhaseFrequencyOffset('PhaseOffset',phasePre(frame_id),'FrequencyOffset',freqPre(frame_id),'SampleRate',sps);
        rxFrame_checked_sample=Error_pre_Iterator(rxFrame_checked_sample);
        rxFrame_checked_sample=TEmodel(rxFrame_checked_sample,timePre(frame_id));
        %   
        [Est_phase_d(:,frame_id), Est_freq_d(:,frame_id) ,Est_time_d(:,frame_id) ]=single_DNN_estimate_joint(rxFrame_checked_sample,symbol_noise_var,H,sps);
        %使用预测结果修正
        Error_amend_Iterator=comm.PhaseFrequencyOffset('PhaseOffset',-Est_phase_d(end,frame_id),'FrequencyOffset',-Est_freq_d(end,frame_id),'SampleRate',sps);
        rxFrame_checked_sample=Error_amend_Iterator(rxFrame_checked_sample);
        rxFrame_checked_sample=TEmodel(rxFrame_checked_sample,-Est_time_d(end,frame_id));
        %EM预测
        [Est_phase_v(:,frame_id),Est_freq_v(:,frame_id),Est_time_v(:,frame_id),LLR_vector]=sequential_EM_estimate_joint(rxFrame_checked_sample,symbol_noise_var,H,sps);
        %使用预测结果进行译码
        decode_bit_frame = (LLR_vector<0);
        [~,BER(frame_id)]=biterr(decode_bit_frame(1:code_K),src_bit_frames(:,frame_id));
    end
    toc
    BER_test(i)=mean(BER);

    %% 重置通信工具箱对象
    release(ldpcEnc)
    release(ldpcDec)
    release(qpskMod)
    release(qpskDemod)
    release(txfilter)
    release(rxfilter)

end

%% 绘图输出
figure(1)
semilogy(EbN0_scale,BER_test,'k*-')
xlabel('Eb/N0(dB)')
ylabel('BER')
grid on;

% load chirp;
% sound(y,Fs)

