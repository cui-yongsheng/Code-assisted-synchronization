load('Objtest.mat')
%% 绘图输出
[X,Y]=meshgrid(phaseErrVector,freqErrVector);
figure(1)
mesh(X,Y,Object_Fun_1(:,:,9).')
xlim([-180 180])
xticks(-180:90:180)
xlabel('载波相位 (deg)')
ylabel('载波频率 (NFO)')
zlabel('归一化函数值')
view(63,15)

figure(2)
mesh(X,Y,Object_Fun_2(:,:,9).')
xlim([-180 180])
xticks(-180:90:180)
xlabel('载波相位 (deg)')
ylabel('载波频率 (NFO)')
zlabel('归一化函数值')
set(gcf, 'renderer', 'opengl'); % 输出时建议使用导出设置，选择300dpi的矢量输出格式
view(63,15)