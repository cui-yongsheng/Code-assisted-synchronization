function [F_max,Angle_max] = give_CZT_graph(Fs,Signal_sequence,f1,f2)
%GIVE_CZT_GRAPH 输入时域信号，绘制CZT频谱 (非系统文件! 仅用于调试)
%   CZT点数同输入点数相同,但应尽量使点数为偶数,否则频率轴可能出现频率分辨率尺度的偏移
%   F_max 3*1数组 幅值最高的前K个频率值（Hz）
%   f_accurate 当前的频率分辨率
%% 函数内置的超参数
K=1; %标记幅值前K的频谱点
H=0; %峰值文字浮动高度
%% 获取输入信号的相关信息
N_sample=length(Signal_sequence);
f_accurate=Fs/N_sample;
m=N_sample; %CZT点数
w=exp(-1j*2*pi*(f2-f1)/(N_sample*Fs)); %细化频段步长
a=exp(1j*2*pi*f1/Fs);                  %细化频段起始点
%% CZT变换获取频域采样
Signal_CZT=czt(Signal_sequence,m,w,a);
Signal_CZT_abs=abs(Signal_CZT);
Signal_CZT_angle = angle(Signal_CZT);
f_axis=((0:length(Signal_CZT)-1)*(f2-f1)/length(Signal_CZT)) + f1;
%% 搜索幅值前K的频率值，按从大到小排列
[~,F_Index]=maxk(Signal_CZT_abs,K);
F_max=f_axis(F_Index);
Angle_max=Signal_CZT_angle(F_Index);

end

