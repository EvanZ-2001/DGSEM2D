function [J, xi_x, xi_y, eta_x, eta_y] = mapping(nodes, xi_ref)
% MAPPING 计算每个单元的几何映射雅可比及映射导数（基于双线性顶点映射）
%
% 输入：
%   nodes   – Q×Q×2×Ncells 数组，nodes(i,j,1,c)=x, nodes(i,j,2,c)=y
%   xi_ref  – (Q×1) 列向量，参考域 GLL 节点 xi 坐标（Q = P+1）
%
% 输出：
%   J       – Q×Q×Ncells数组，雅可比行列式 J(c,i,j)
%   xi_x    – Q×Q×Ncells 数组，∂xi/∂x
%   xi_y    – Q×Q×Ncells 数组，∂xi/∂y
%   eta_x   – Q×Q×Ncells 数组，∂eta/∂x
%   eta_y   – Q×Q×Ncells 数组，∂eta/∂y

[Q, ~, ~, Ncells] = size(nodes);

% 构造 (i,j) 对应的参考坐标 xi, eta
% [Xi, Eta] = meshgrid(xi_ref, xi_ref);
[Eta, Xi] = meshgrid(xi_ref, xi_ref); % ij对调改动
% Xi(i,j)=xi_ref(i)、Eta(i,j)=xi_ref(j)

% 预分配，全部赋值为0
J     = zeros(Q, Q, Ncells);
xi_x  = zeros(Q, Q, Ncells); 
xi_y  = zeros(Q, Q, Ncells);
eta_x = zeros(Q, Q, Ncells);
eta_y = zeros(Q, Q, Ncells);

% 对每个单元逐点计算
for c = 1:Ncells
    % 四个顶点坐标
    % xA = nodes(c,1,1,1); yA = nodes(c,1,1,2);
    % xB = nodes(c,Q,1,1); yB = nodes(c,Q,1,2);
    % xC = nodes(c,Q,Q,1); yC = nodes(c,Q,Q,2);
    % xD = nodes(c,1,Q,1); yD = nodes(c,1,Q,2);
    xA = nodes(1,1,1,c); yA = nodes(1,1,2,c);
    xB = nodes(Q,1,1,c); yB = nodes(Q,1,2,c);
    xC = nodes(Q,Q,1,c); yC = nodes(Q,Q,2,c);
    xD = nodes(1,Q,1,c); yD = nodes(1,Q,2,c);

    for i = 1:Q
        for j = 1:Q
            xi = Xi(i,j);
            eta = Eta(i,j);
            % 双线性形函数及偏导
            N1   = (1-xi)*(1-eta)/4;   dN1_dxi  = -(1-eta)/4;  dN1_deta  = -(1-xi)/4;
            N2   = (1+xi)*(1-eta)/4;   dN2_dxi  =  (1-eta)/4;  dN2_deta  = -(1+xi)/4;
            N3   = (1+xi)*(1+eta)/4;   dN3_dxi  =  (1+eta)/4;  dN3_deta  =  (1+xi)/4;
            N4   = (1-xi)*(1+eta)/4;   dN4_dxi  = -(1+eta)/4;  dN4_deta  =  (1-xi)/4;
            % 映射导数
            x_xi  = dN1_dxi*xA + dN2_dxi*xB + dN3_dxi*xC + dN4_dxi*xD;
            x_eta = dN1_deta*xA + dN2_deta*xB + dN3_deta*xC + dN4_deta*xD;
            y_xi  = dN1_dxi*yA + dN2_dxi*yB + dN3_dxi*yC + dN4_dxi*yD;
            y_eta = dN1_deta*yA + dN2_deta*yB + dN3_deta*yC + dN4_deta*yD;
            % 雅可比行列式
            J(i,j,c) = x_xi * y_eta - x_eta * y_xi;
            % 逆映射导数
            xi_x(i,j,c)  =  y_eta / J(i,j,c);
            xi_y(i,j,c)  = -x_eta / J(i,j,c);
            eta_x(i,j,c) = -y_xi / J(i,j,c);
            eta_y(i,j,c) =  x_xi / J(i,j,c);
        end
    end
end
end
