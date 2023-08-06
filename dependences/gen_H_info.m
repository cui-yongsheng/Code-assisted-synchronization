function H_info = gen_H_info(H)
%GEN_H_INFO У�����Ԥ������������ִ��һ�Σ����ɵ�У�������Ϣ�ɹ���������������
%   ����:     1.  У�����                        H           sparse(M*N)
%
%   ���:     1.  У�����Ԥ������              H_info      cell(1*3)
%             

%% У�������Ϣ��ȡ
H_weight = nnz(H);              % У���������
[r_index,~] = find(H~=0);       % ��������λ������
ColWeight = full(sum(H,1));     % У���������
RowWeight = full(sum(H,2));     % У���������
RowIndex = ones(size(H,1),max(RowWeight))*(H_weight+1);  % ÿ�з���������
for row = 1:size(H,1)
    RowIndex(row,1:RowWeight(row)) = find(r_index==row).';
end

%% ������
H_info.H_weight = H_weight;
H_info.RowIndex = RowIndex;
H_info.ColWeight = ColWeight;

end
