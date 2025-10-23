function  n = FaceNormal(nodes_ele,faceN)
% 输入：
%     node_ele - 当前单元所有积分点的坐标，维度Q*Q*2
%     faceN - 面编号，1-下边界，2-右边界，3-上边界，4-左边界
% 输出：
%     n - 当前单元指定面编号的边界的单位外法向

    % 节点数
    Q = size(nodes_ele,1);
    
    % 四个顶点坐标
    xA = nodes_ele(1,1,1); yA = nodes_ele(1,1,2);
    xB = nodes_ele(Q,1,1); yB = nodes_ele(Q,1,2);
    xC = nodes_ele(Q,Q,1); yC = nodes_ele(Q,Q,2);
    xD = nodes_ele(1,Q,1); yD = nodes_ele(1,Q,2);
    % switch face
    %     case 1, n = [0, -1];  % 下
    %     case 2, n = [1, 0];   % 右
    %     case 3, n = [0, 1];  % 上
    %     case 4, n = [-1, 0];   % 左
    % end
      % 根据 faceN 选取对应边的切向量（端点差分）
    switch faceN
      case 1  % 下边界 AB
        dx = xB - xA;  
        dy = yB - yA;
      case 2  % 右边界 BC
        dx = xC - xB;
        dy = yC - yB;
      case 3  % 上边界 CD
        dx = xD - xC;
        dy = yD - yC;
      case 4  % 左边界 DA
        dx = xA - xD;
        dy = yA - yD;
      otherwise
        error('faceN 必须是 1–4 之间的整数');
    end

    % 切向量 [dx,dy] 顺时针旋转 90° 得到外法向 [dy, -dx]，再归一化
    n = [ dy, -dx ];
    n = n / norm(n);
end