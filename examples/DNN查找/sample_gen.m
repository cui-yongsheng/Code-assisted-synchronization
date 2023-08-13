%% 生成样本数据集
clc
clear
close all
rng default
addpath(genpath('../../functions'));
addpath(genpath('../../dependences'));
addpath(genpath('../../datas'));

%% 测试控制参数
EbN0_scale=2;                                       % EbN0测试范围  (译码门限 1.4)
BP_times=2;                                         % BP译码迭代次数
run_times=1e5;                                      % 生成帧数

%% 信号调制及编码参数
sps = 4;                        % 单位符号采样数
rolloff = 0.35;                 % 根升余弦滤波器滚降系数
span = 6;                       % 根升余弦滤波器阶数
code_r=1/6;                     % LDPC编码效率 (1/2,1/4)
code_K=1024;                    % LDPC编码信息位数量
H_src=load('LDPC_非规则码.mat');
H=H_src.H;                      % LDPC校验矩阵
code_N=size(H,2);               % LDPC码字长度
M=4;                            % 调制阶数 (QPSK)
code_r_real=code_K/code_N;      % 实际LDPC编码效率

%% 通信工具箱对象实例化
ldpcEnc = comm.LDPCEncoder(H);     % LDPC编码器
ldpcDec_for_LLR = comm.LDPCDecoder('ParityCheckMatrix',H,'DecisionMethod','Soft decision','MaximumIterationCount',BP_times,'OutputValue','Whole codeword','FinalParityChecksOutputPort',true);
qpskMod = comm.QPSKModulator('BitInput',true);                  % QPSK调制器
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK解调器(软解调，并根据输入的噪声功率，
    'DecisionMethod','Approximate log-likelihood ratio', ...    % 在假定复信号功率为1的条件下进行LLR归一化
    'VarianceSource','Input port');
txfilter = comm.RaisedCosineTransmitFilter...
    ('RolloffFactor',rolloff,'FilterSpanInSymbols',....
    span,'OutputSamplesPerSymbol',sps);                         % 发端-根升余弦滚降滤波器
rxfilter = comm.RaisedCosineReceiveFilter...
    ('RolloffFactor',rolloff,'FilterSpanInSymbols',span,...
    'DecimationFactor',1,'InputSamplesPerSymbol',sps);          % 收端-根升余弦滚降滤波器(不进行降采样 Dec=1)

%% 预置统计量
label_val_p=zeros(run_times,1);
label_val_f=zeros(run_times,1);
label_val_t=zeros(run_times,1);
Fun_d0=zeros(run_times,1);
Fun_d0_p=zeros(run_times,1);
Fun_d0_f=zeros(run_times,1);
Fun_d0_t=zeros(run_times,1);
Fun_d1=zeros(run_times,1);
Fun_d1_p=zeros(run_times,1);
Fun_d1_f=zeros(run_times,1);
Fun_d1_t=zeros(run_times,1);
Fun_d2=zeros(run_times,1);
Fun_d2_p=zeros(run_times,1);
Fun_d2_f=zeros(run_times,1);
Fun_d2_t=zeros(run_times,1);
Fun_d3=zeros(run_times,1);
Fun_d3_p=zeros(run_times,1);
Fun_d3_f=zeros(run_times,1);
Fun_d3_t=zeros(run_times,1);


%% 帧数
parfor i=1:run_times
    if(mod(i,1000)==0)
        disp(['####### 第 ' num2str(i) ' 帧  #######']);
    end
    %% 发端生成信号
    QPSK_frame_length=code_N/log2(M);                       % QPSK符号帧长度
    src_bit_frame=logical(randi([0 1],code_K,1));           % 随机生成信源帧集合
    code_bit_frame=ldpcEnc(src_bit_frame);                  % LDPC编码后的信源帧 (ToDo: 目前暂不考虑帧同步的实现)
    QPSK_symbol_frame=qpskMod(code_bit_frame);              % 编码信源帧映射为复数形式的QPSK符号帧
    QPSK_symbol_stream=[QPSK_symbol_frame;zeros(span,1)];   % 符号流中加入全零尾片段（克服升余弦滤波产生的尾符号丢失）
    txSig=txfilter(QPSK_symbol_stream);                     % 脉冲成形，产生发端基带信号

    %% 当前位置
    label_val_t(i)=rand()*2-1;
    label_val_p(i)=rand()*60-30;
    label_val_f(i)=rand()*1e-4-5e-5;
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',label_val_f(i),....
    'PhaseOffset',label_val_p(i),'SampleRate',sps);
    txSig=ErrIterator(txSig);
    txSig=TEmodel(txSig,label_val_t(i));

    %% 引入信道噪声
    EbN0=EbN0_scale;                              % 归一化编码前比特能量Eb/N0
    EbN0_code=EbN0+10*log10(code_r_real);         % 归一化编码后比特能量Eb_code/N0
    EsN0=10*log10(2)+EbN0_code;                   % 归一化调制符号能量Es/N0 ( QPSK调制 )
    SNR=EsN0-10*log10(sps);                       % 位同步前的采样序列信噪比
    rxSig=awgn(txSig,SNR,'measured');
    symbol_noise_var=10^(-1*EsN0/10);
    
    %% 收端处理
    rxSig=rxfilter(rxSig);                                  % 匹配滤波
    rxSig=rxSig((span)*sps+1:end);                          % 去除匹配滤波拖尾

    %% 收端处理0 （当前位置）
    rxSig_checked_sample=rxSig;
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);
    
    %% 计算目标函数
    LLR_vector_self=ldpcDec_for_LLR(LLR);
    Fun_d0(i)=0.5*mean(log(1+cosh(LLR_vector_self)));
    BPSK_rec=tanh(LLR_vector_self/2);
    QPSK_rec=(1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end))/sqrt(2);
    Fun_d0_p(i)=imag(sum(rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    Fun_d0_f(i)=imag(sum((2*pi.*(0:(length(rxSig_checked)-1)).').*rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    sample_time = 1:length(rxSig_checked_sample);
    rxFrame_adjust_sample=(interp1(sample_time,rxSig_checked_sample',sample_time+0.01))';
    rxFrame_gradient_sample=(rxFrame_adjust_sample-rxSig_checked_sample)./0.01;
    rxFrame_gradient_trail=rxFrame_gradient_sample(1:sps:end);
    Fun_d0_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q函数一阶导

    %% 收端处理1 （相邻位置）
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',1e-5,....
    'PhaseOffset',5,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(rxSig);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,0);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);

    %% 计算目标函数1
    LLR_vector_self=ldpcDec_for_LLR(LLR);
    Fun_d1(i)=0.5*mean(log(1+cosh(LLR_vector_self)));
    BPSK_rec=tanh(LLR_vector_self/2);
    QPSK_rec=(1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end))/sqrt(2);
    Fun_d1_p(i)=imag(sum(rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    Fun_d1_f(i)=imag(sum((2*pi.*(0:(length(rxSig_checked)-1)).').*rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    sample_time = 1:length(rxSig_checked_sample);
    rxFrame_adjust_sample=(interp1(sample_time,rxSig_checked_sample',sample_time+0.01))';
    rxFrame_gradient_sample=(rxFrame_adjust_sample-rxSig_checked_sample)./0.01;
    rxFrame_gradient_trail=rxFrame_gradient_sample(1:sps:end);
    Fun_d1_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q函数一阶导

    %% 收端处理2 （相邻位置）
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',0,....
    'PhaseOffset',0,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(rxSig);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,1);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);

    %% 计算目标函数2
    LLR_vector_self=ldpcDec_for_LLR(LLR);
    Fun_d2(i)=0.5*mean(log(1+cosh(LLR_vector_self)));
    BPSK_rec=tanh(LLR_vector_self/2);
    QPSK_rec=(1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end))/sqrt(2);
    Fun_d2_p(i)=imag(sum(rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    Fun_d2_f(i)=imag(sum((2*pi.*(0:(length(rxSig_checked)-1)).').*rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    sample_time = 1:length(rxSig_checked_sample);
    rxFrame_adjust_sample=(interp1(sample_time,rxSig_checked_sample',sample_time+0.01))';
    rxFrame_gradient_sample=(rxFrame_adjust_sample-rxSig_checked_sample)./0.01;
    rxFrame_gradient_trail=rxFrame_gradient_sample(1:sps:end);
    Fun_d2_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q函数一阶导

    %% 收端处理3 （相邻位置）
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',1e-5,....
    'PhaseOffset',5,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(rxSig);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,1);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);

    %% 计算目标函数3
    LLR_vector_self=ldpcDec_for_LLR(LLR);
    Fun_d3(i)=0.5*mean(log(1+cosh(LLR_vector_self)));
    BPSK_rec=tanh(LLR_vector_self/2);
    QPSK_rec=(1i*BPSK_rec(1:2:end)+BPSK_rec(2:2:end))/sqrt(2);
    Fun_d3_p(i)=imag(sum(rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    Fun_d3_f(i)=imag(sum((2*pi.*(0:(length(rxSig_checked)-1)).').*rxSig_checked.*conj(QPSK_rec)))/symbol_noise_var;
    sample_time = 1:length(rxSig_checked_sample);
    rxFrame_adjust_sample=(interp1(sample_time,rxSig_checked_sample',sample_time+0.01))';
    rxFrame_gradient_sample=(rxFrame_adjust_sample-rxSig_checked_sample)./0.01;
    rxFrame_gradient_trail=rxFrame_gradient_sample(1:sps:end);
    Fun_d3_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q函数一阶导
end

%% 特征工程构造查找函数
look_fun_table_t = (Fun_d3_t-Fun_d0_t)/2;
look_fun_table_p = (Fun_d3_p-Fun_d0_p)/2;
look_fun_table_f = (Fun_d3_f-Fun_d0_f)/2;
look_fun_table_t1 = Fun_d0_t./Fun_d3_t;
look_fun_table_p1 = Fun_d0_p./Fun_d3_p;
look_fun_table_f1 = Fun_d0_f./Fun_d3_f;
look_fun_table = [look_fun_table_t look_fun_table_p look_fun_table_f];
look_fun_table_1 = [look_fun_table_t1 look_fun_table_p1 look_fun_table_f1];

%% 数据打包
point0 = [Fun_d0 Fun_d0_p Fun_d0_f Fun_d0_t];
point1 = [Fun_d1 Fun_d1_p Fun_d1_f Fun_d1_t];
point2 = [Fun_d2 Fun_d2_p Fun_d2_f Fun_d2_t];
point3 = [Fun_d3 Fun_d3_p Fun_d3_f Fun_d3_t];
label_val_f = label_val_f.*10e5;
df=[point0 point1 point2 point3 look_fun_table look_fun_table_1 label_val_p label_val_f label_val_t];

