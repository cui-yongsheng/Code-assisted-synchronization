function H_info = gen_H_info(H)
%GEN_H_INFO 校验矩阵预处理函数，仅需执行一次，生成的校验矩阵信息可供译码器反复调用
%   输入:     1.  校验矩阵                        H           sparse(M*N)
%
%   输出:     1.  校验矩阵预处理结果              H_info      cell(1*3)
%             

%% 校验矩阵信息提取
H_weight = nnz(H);              % 校验矩阵重量
[r_index,~] = find(H~=0);       % 非零项行位置向量
ColWeight = full(sum(H,1));     % 校验矩阵列重
RowWeight = full(sum(H,2));     % 校验矩阵行重
RowIndex = ones(size(H,1),max(RowWeight))*(H_weight+1);  % 每行非零项索引
for row = 1:size(H,1)
    RowIndex(row,1:RowWeight(row)) = find(r_index==row).';
end

%% 结果输出
H_info.H_weight = H_weight;
H_info.RowIndex = RowIndex;
H_info.ColWeight = ColWeight;

end
