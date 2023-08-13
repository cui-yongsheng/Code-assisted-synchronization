barker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];	%13位巴克码序列

[y1, x1] = xcorr(barker);	%计算13位巴克码序列的自相关值
        
figure(1);
plot(x1, abs(y1), '-o');
% axis tight;
axis([-13 13 0 14]);
xlabel('\itj', 'FontName', 'Times New Roman');
ylabel('\it\gamma(j)', 'FontName', 'Times New Roman');

mseq=idinput(13,'prbs');                           % pn序列
[y1, x1] = xcorr(mseq);	%计算13位巴克码序列的自相关值
        
figure(2);
plot(x1, abs(y1), '-o');
% axis tight;
axis([-13 13 0 14]);
xlabel('\itj', 'FontName', 'Times New Roman');
ylabel('\it\gamma(j)', 'FontName', 'Times New Roman');