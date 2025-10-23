function numF = numerical_flux_all(cell_id, U, nodes, Nx, Ny, w, g, x_xi, x_eta, y_xi, y_eta)
% NUMERICAL_FLUX_ALL 计算第 cell_id 单元所有 GLL 点的界面数值通量贡献
%
% 输入：
%   cell_id – 当前单元编号
%   U       – Q×Q×3×Ncells 数组，守恒量 [h, hu, hv]
%   nodes   – Q×Q×2×Ncells 数组，nodes(i,j,1,c)=x, nodes(i,j,2,c)=y（仅用于边界处理）
%   w       – Q×1 向量，GLL 权重
%   g       – 重力加速度
%   x_xi, x_eta, y_xi, y_eta  - 映射导数
%
% 输出：
%   numF    – Q×Q×3 张量，累计的数值通量项（已含插值权重）
%
% 注意：
%   单元边界面编号对应关系：1-下边界，2-右边界，3-上边界，4-左边界

    % 假设全局定义 Nx, Ny
    % global Nx Ny

    [Q, ~, ~, ~] = size(U);
    numF = zeros(Q, Q, 3);

    % 当前单元在网格中的 (ix,iy) 位置
    ix = mod(cell_id-1, Nx) + 1;
    iy = floor((cell_id-1)/Nx) + 1;

    % 当前单元所有积分点的坐标
    nodes_ele = nodes(:,:,:,cell_id);

    % 四个顶点坐标
    xA = nodes(1,1,1,cell_id); yA = nodes(1,1,2,cell_id);
    xB = nodes(Q,1,1,cell_id); yB = nodes(Q,1,2,cell_id);
    xC = nodes(Q,Q,1,cell_id); yC = nodes(Q,Q,2,cell_id);
    xD = nodes(1,Q,1,cell_id); yD = nodes(1,Q,2,cell_id);

    % 左边界 (xi = -1), 面编号 faceN = 4 
    faceN = 4;
    n = FaceNormal(nodes_ele,faceN); % 单位外法向
    L4 = sqrt( (xA-xD)^2 + (yA-yD)^2 ) ; % 该边界长度
    for j = 1:Q
        U_in = reshape(U(1, j, :, cell_id), [], 1); % Q×1
        if ix > 1 % 不是物理边界
            neighbor = cell_id - 1;
            U_out  = reshape(U(Q, j, :, neighbor), [], 1); % Q×1
        else % 是物理边界
            U_out  = BoundaryState(U_in, faceN);
        end
        Fn = numerical_flux_LF_Roe(U_in, U_out, n, g);  % 3×1 向量
        numF(1, j, :) = reshape(numF(1, j, :), [], 1) + (L4/2) * w(j) * Fn;
    end

    % 右边界 (xi = +1), 面编号 faceN = 2 
    faceN = 2;
    n = FaceNormal(nodes_ele,faceN);% 单位外法向
    L2 = sqrt( (xB-xC)^2 + (yB-yC)^2 ) ; % 该边界长度
    for j = 1:Q
        U_in = reshape(U(Q, j, :, cell_id), [], 1); % Q×1
        if ix < Nx % 不是物理边界
            neighbor = cell_id + 1;
            U_out = reshape(U(1, j, :, neighbor), [], 1); % Q×1
        else % 是物理边界
            U_out  = BoundaryState(U_in, faceN);
        end
        Fn = numerical_flux_LF_Roe(U_in, U_out, n, g);
        numF(Q, j, :) = reshape(numF(Q, j, :), [], 1) + (L2/2) * w(j) * Fn;
    end

    % 下边界 (eta = -1), 面编号 faceN = 1
    faceN = 1;
    n = FaceNormal(nodes_ele,faceN);% 单位外法向
    L1 = sqrt( (xA-xB)^2 + (yA-yB)^2 ) ; % 该边界长度
    for i = 1:Q
        U_in = reshape(U(i, 1, :, cell_id), [], 1); % Q×1
        if iy > 1
            neighbor = cell_id - Nx;
            U_out = reshape(U(i, Q, :, neighbor), [], 1); % Q×1
        else
            U_out  = BoundaryState(U_in, faceN);
        end
        Fn = numerical_flux_LF_Roe(U_in, U_out, n, g);
        numF(i, 1, :) = reshape(numF(i, 1, :), [], 1) + (L1/2) * w(i) * Fn;
    end

    % 上边界 (eta = +1), 面编号 faceN = 3
    faceN = 3;
    n = FaceNormal(nodes_ele,faceN);% 单位外法向
    L3 = sqrt( (xC-xD)^2 + (yC-yD)^2 ) ; % 该边界长度
    for i = 1:Q
        U_in = reshape(U(i, Q, :, cell_id), [], 1); % Q×1
        if iy < Ny
            neighbor = cell_id + Nx;
            U_out = reshape(U(i, 1, :, neighbor), [], 1); % Q×1
        else
            U_out  = BoundaryState(U_in, faceN);
        end
        Fn = numerical_flux_LF_Roe(U_in, U_out, n, g);
        numF(i, Q, :) = reshape(numF(i, Q, :), [], 1) + (L3/2) * w(i) * Fn;
    end
end
