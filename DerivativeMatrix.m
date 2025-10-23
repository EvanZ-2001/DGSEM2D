%% Gauss-Lobatto-Legendre型数值微分的微分矩阵

% Q是积分点的个数,nodes是积分点的数组

function DDX = DerivativeMatrix(Q,nodes)
    % 计算非对角线元素
    nline = 1:1:Q;
    % legendre_nline = LegendrePolynomial(nodes(nline),Q-1);% 计算出积分点向量的Q-1阶Legendre多项式
    legendre_nline = legendreP(Q-1,nodes(nline));% 计算出积分点向量的Q-1阶Legendre多项式
    legendre_nline_i = repmat(legendre_nline,1,Q); % 列向量横向扩展成矩阵
    legendre_nline_j = repmat(legendre_nline',Q,1); % 行向量纵向扩展成矩阵
    matrix_A = legendre_nline_i./legendre_nline_j;
    
    xi_nline_i = repmat(nodes,1,Q);
    xi_nline_j = repmat(nodes',Q,1);
    matrix_B = 1./(xi_nline_i-xi_nline_j);

    DDX = matrix_A.*matrix_B; % 数值微分矩阵拆分为2个矩阵逐项相乘

    %对角线置零
    for i=1:Q 
        DDX(i,i)=0;
    end
    
    %对角线两个端点单独赋值
    DDX(1,1) = -Q*(Q-1)/4;
    DDX(Q,Q) = Q*(Q-1)/4;
end
