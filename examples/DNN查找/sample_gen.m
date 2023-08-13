%% �����������ݼ�
clc
clear
close all
rng default
addpath(genpath('../../functions'));
addpath(genpath('../../dependences'));
addpath(genpath('../../datas'));

%% ���Կ��Ʋ���
EbN0_scale=2;                                       % EbN0���Է�Χ  (�������� 1.4)
BP_times=2;                                         % BP�����������
run_times=1e5;                                      % ����֡��

%% �źŵ��Ƽ��������
sps = 4;                        % ��λ���Ų�����
rolloff = 0.35;                 % ���������˲�������ϵ��
span = 6;                       % ���������˲�������
code_r=1/6;                     % LDPC����Ч�� (1/2,1/4)
code_K=1024;                    % LDPC������Ϣλ����
H_src=load('LDPC_�ǹ�����.mat');
H=H_src.H;                      % LDPCУ�����
code_N=size(H,2);               % LDPC���ֳ���
M=4;                            % ���ƽ��� (QPSK)
code_r_real=code_K/code_N;      % ʵ��LDPC����Ч��

%% ͨ�Ź��������ʵ����
ldpcEnc = comm.LDPCEncoder(H);     % LDPC������
ldpcDec_for_LLR = comm.LDPCDecoder('ParityCheckMatrix',H,'DecisionMethod','Soft decision','MaximumIterationCount',BP_times,'OutputValue','Whole codeword','FinalParityChecksOutputPort',true);
qpskMod = comm.QPSKModulator('BitInput',true);                  % QPSK������
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...           % QPSK�����(������������������������ʣ�
    'DecisionMethod','Approximate log-likelihood ratio', ...    % �ڼٶ����źŹ���Ϊ1�������½���LLR��һ��
    'VarianceSource','Input port');
txfilter = comm.RaisedCosineTransmitFilter...
    ('RolloffFactor',rolloff,'FilterSpanInSymbols',....
    span,'OutputSamplesPerSymbol',sps);                         % ����-�������ҹ����˲���
rxfilter = comm.RaisedCosineReceiveFilter...
    ('RolloffFactor',rolloff,'FilterSpanInSymbols',span,...
    'DecimationFactor',1,'InputSamplesPerSymbol',sps);          % �ն�-�������ҹ����˲���(�����н����� Dec=1)

%% Ԥ��ͳ����
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


%% ֡��
parfor i=1:run_times
    if(mod(i,1000)==0)
        disp(['####### �� ' num2str(i) ' ֡  #######']);
    end
    %% ���������ź�
    QPSK_frame_length=code_N/log2(M);                       % QPSK����֡����
    src_bit_frame=logical(randi([0 1],code_K,1));           % ���������Դ֡����
    code_bit_frame=ldpcEnc(src_bit_frame);                  % LDPC��������Դ֡ (ToDo: Ŀǰ�ݲ�����֡ͬ����ʵ��)
    QPSK_symbol_frame=qpskMod(code_bit_frame);              % ������Դ֡ӳ��Ϊ������ʽ��QPSK����֡
    QPSK_symbol_stream=[QPSK_symbol_frame;zeros(span,1)];   % �������м���ȫ��βƬ�Σ��˷��������˲�������β���Ŷ�ʧ��
    txSig=txfilter(QPSK_symbol_stream);                     % ������Σ��������˻����ź�

    %% ��ǰλ��
    label_val_t(i)=rand()*2-1;
    label_val_p(i)=rand()*60-30;
    label_val_f(i)=rand()*1e-4-5e-5;
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',label_val_f(i),....
    'PhaseOffset',label_val_p(i),'SampleRate',sps);
    txSig=ErrIterator(txSig);
    txSig=TEmodel(txSig,label_val_t(i));

    %% �����ŵ�����
    EbN0=EbN0_scale;                              % ��һ������ǰ��������Eb/N0
    EbN0_code=EbN0+10*log10(code_r_real);         % ��һ��������������Eb_code/N0
    EsN0=10*log10(2)+EbN0_code;                   % ��һ�����Ʒ�������Es/N0 ( QPSK���� )
    SNR=EsN0-10*log10(sps);                       % λͬ��ǰ�Ĳ������������
    rxSig=awgn(txSig,SNR,'measured');
    symbol_noise_var=10^(-1*EsN0/10);
    
    %% �ն˴���
    rxSig=rxfilter(rxSig);                                  % ƥ���˲�
    rxSig=rxSig((span)*sps+1:end);                          % ȥ��ƥ���˲���β

    %% �ն˴���0 ����ǰλ�ã�
    rxSig_checked_sample=rxSig;
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);
    
    %% ����Ŀ�꺯��
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
    Fun_d0_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q����һ�׵�

    %% �ն˴���1 ������λ�ã�
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',1e-5,....
    'PhaseOffset',5,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(rxSig);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,0);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);

    %% ����Ŀ�꺯��1
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
    Fun_d1_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q����һ�׵�

    %% �ն˴���2 ������λ�ã�
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',0,....
    'PhaseOffset',0,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(rxSig);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,1);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);

    %% ����Ŀ�꺯��2
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
    Fun_d2_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q����һ�׵�

    %% �ն˴���3 ������λ�ã�
    ErrIterator = comm.PhaseFrequencyOffset('FrequencyOffset',1e-5,....
    'PhaseOffset',5,'SampleRate',sps);
    rxSig_checked_sample=ErrIterator(rxSig);
    rxSig_checked_sample=TEmodel(rxSig_checked_sample,1);
    rxSig_checked=rxSig_checked_sample(1:sps:end);
    LLR=qpskDemod(rxSig_checked,symbol_noise_var);

    %% ����Ŀ�꺯��3
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
    Fun_d3_t(i)=real(sum(rxFrame_gradient_trail.*conj(QPSK_rec)))/symbol_noise_var;   %Q����һ�׵�
end

%% �������̹�����Һ���
look_fun_table_t = (Fun_d3_t-Fun_d0_t)/2;
look_fun_table_p = (Fun_d3_p-Fun_d0_p)/2;
look_fun_table_f = (Fun_d3_f-Fun_d0_f)/2;
look_fun_table_t1 = Fun_d0_t./Fun_d3_t;
look_fun_table_p1 = Fun_d0_p./Fun_d3_p;
look_fun_table_f1 = Fun_d0_f./Fun_d3_f;
look_fun_table = [look_fun_table_t look_fun_table_p look_fun_table_f];
look_fun_table_1 = [look_fun_table_t1 look_fun_table_p1 look_fun_table_f1];

%% ���ݴ��
point0 = [Fun_d0 Fun_d0_p Fun_d0_f Fun_d0_t];
point1 = [Fun_d1 Fun_d1_p Fun_d1_f Fun_d1_t];
point2 = [Fun_d2 Fun_d2_p Fun_d2_f Fun_d2_t];
point3 = [Fun_d3 Fun_d3_p Fun_d3_f Fun_d3_t];
label_val_f = label_val_f.*10e5;
df=[point0 point1 point2 point3 look_fun_table look_fun_table_1 label_val_p label_val_f label_val_t];

