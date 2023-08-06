%% 主执行脚本：目标函数曲线仿真。
clc
clear
close all
rng default
addpath(genpath('..\..\functions'));
addpath(genpath('..\..\dependences'));
addpath(genpath('..\..\datas'));

%% 测试控制参数
phaseErrVector=-180:2:180;
freqErrVector=(-8:0.4:8)*1e-4;
timeErrVector=-8:1:8;                   % 附加时延（样本点数）
EbN0_scale = 2;                         % EbN0测试范围  (理想同步时，码率1/2译码门限2dB)
nFrame = 1;                             % 测试帧数
BP_times=1;                             % LDPC译码迭代次数

%% 信号调制及编码参数
M=4;                            % 调制阶数 (QPSK)
sps = 16;                       % 单位符号采样数
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
%% 通信工具箱对象实例化
ldpcEnc = comm.LDPCEncoder(H);     % LDPC编码器
ldpcDec_for_LLR = comm.LDPCDecoder('ParityCheckMatrix',H,'DecisionMethod','Soft decision','MaximumIterationCount',BP_times,'OutputValue','Whole codeword','FinalParityChecksOutputPort',true);
qpskMod = comm.QPSKModulator('BitInput',true);                  % QPSK调制器
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK解调器(软解调，并根据输入的噪声功率，
    'DecisionMethod','Approximate log-likelihood ratio', ...    % 在假定复信号功率为1的条件下进行LLR归一化
    'VarianceSource','Input port');
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff,'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',sps);   % 发端根升余弦滚降滤波器
rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff,'FilterSpanInSymbols',span,'DecimationFactor',1,'InputSamplesPerSymbol',sps);        % 收端根升余弦滚降滤波器


%% 在不同EbN0下循环测试
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

    %% 引入信道噪声
    EbN0=EbN0_scale;                              % 归一化编码前比特能量Eb/N0
    EbN0_code=EbN0+10*log10(code_r_real);         % 归一化编码后比特能量Eb_code/N0
    EsN0=10*log10(2)+EbN0_code;                   % 归一化调制符号能量Es/N0 ( QPSK调制 )
    SNR=EsN0-10*log10(sps);                       % 位同步前的采样序列信噪比
    rxSig=awgn(txSig,SNR,'measured');
    symbol_noise_var=10^(-1*EsN0/10);

    %% 收端测试
    Object_Fun_1=zeros(length(phaseErrVector),length(freqErrVector),length(timeErrVector));       % LLR目标函数
    Object_Fun_2=zeros(length(phaseErrVector),length(freqErrVector),length(timeErrVector));       % 校验子范数
    for timeErr_index=1:length(timeErrVector)
        disp([num2str(timeErr_index) '/' num2str(length(timeErrVector))])
        for freqErr_index=1:length(freqErrVector)
            for phaseErr_index=1:length(phaseErrVector)
                freqErrIterator=comm.PhaseFrequencyOffset('SampleRate',sps,'PhaseOffset',-1*phaseErrVector(phaseErr_index),'FrequencyOffset',-1*freqErrVector(freqErr_index));
                rxSig_checked=TEmodel(rxSig,timeErrVector(timeErr_index));
                rxSig_checked=freqErrIterator(rxSig_checked);
                rxFiltSig=rxfilter(rxSig_checked);
                rxFiltSig=rxFiltSig(1:sps:end);
                rxQPSK_symbol_stream=rxFiltSig(span+1:end);
                rxLLR_bit_stream=qpskDemod(rxQPSK_symbol_stream,symbol_noise_var);
                rxLLR_frames=reshape(rxLLR_bit_stream,code_N,nFrame);
                %% 以LLR为目标函数
                LLR_vector_self=ldpcDec_for_LLR(rxLLR_frames(:,1));
                Object_Fun_1(phaseErr_index,freqErr_index,timeErr_index)=mean(abs(LLR_vector_self));
                %% 校验子范数
                Hard_vector=(LLR_vector_self<=0);   
                Object_Fun_2(phaseErr_index,freqErr_index,timeErr_index)=sum(1-double(rem(sum(H(:,Hard_vector),2),2)));
            end
        end
    end
end
%% 目标函数归一化
Object_Fun_1 = (Object_Fun_1-min(Object_Fun_1,[],'all'))./(max(Object_Fun_1,[],'all')-min(Object_Fun_1,[],'all'));
Object_Fun_2 = (Object_Fun_2-min(Object_Fun_2,[],'all'))./(max(Object_Fun_2,[],'all')-min(Object_Fun_2,[],'all'));

%% 保存数据
save('../../result/Objtest.mat','Object_Fun_1','Object_Fun_2','phaseErrVector','freqErrVector','timeErrVector')
%% 绘图输出
[X,Y]=meshgrid(phaseErrVector,freqErrVector);
figure(1)
mesh(X,Y,Object_Fun_1(:,:,10).')
xlim([-180 180])
xticks(-180:90:180)
xlabel('载波相位 (deg)')
ylabel('载波频率 (NFO)')
zlabel('归一化函数值')
view(63,15)

figure(2)
mesh(X,Y,Object_Fun_2(:,:,10).')
xlim([-180 180])
xticks(-180:90:180)
xlabel('载波相位 (deg)')
ylabel('载波频率 (NFO)')
zlabel('归一化函数值')
% set(gcf, 'renderer', 'opengl'); % 输出时建议使用导出设置，选择300dpi的矢量输出格式
view(63,15)
% [X,Y,Z]=meshgrid(freqErrVector,phaseErrVector,timeErrVector);
% x = reshape(X,[],1);
% y = reshape(Y,[],1);
% z = reshape(Z,[],1);
% object_Fun_1=reshape(Object_Fun_1,[],1);
% object_Fun_2=reshape(Object_Fun_2,[],1);
% figure(1)
% scatter3(x,y,z,40,object_Fun_1,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% xlabel('载波频率 (NFO)')
% ylabel('载波相位 (deg)')
% zlabel('时延 (point)')
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = '归一化函数值';
% 
% figure(2)
% scatter3(x,y,z,40,object_Fun_2,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% xlabel('载波频率 (NFO)')
% ylabel('载波相位 (deg)')
% zlabel('时延 (point)')
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = '归一化函数值';
