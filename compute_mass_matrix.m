function M = compute_mass_matrix(w, J)
    % COMPUTE_MASS_MATRIX 构造每个单元的局部质量矩阵对角项
    %
    % 输入：
    %   w  — Q×1 列向量，GLL 权重
    %   J  — Ncells×Q×Q 数组，第 c 层是单元 c 上每个 GLL 点的雅可比行列式
    %
    % 输出：
    %   M  — Ncells×Q×Q 数组，第 c 层的 (i,j) 分量等于 J(c,i,j)*w(i)*w(j)
    
    [Ncells, Q, ~] = size(J);
    
    % 构造 w_i * w_j 矩阵
    WiWj = w(:) * w(:).';    % Np×Np
    
    % 直接按元素相乘
    % M(c,i,j) = J(c,i,j) * WiWj(i,j)
    M = J .* reshape(WiWj, [1, Q, Q]);
end