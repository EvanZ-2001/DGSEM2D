% mesh_generate.m
% 将任意四边形域均匀划分为 Nx×Ny 个单元，
% 并在每个单元内部按 GLL 点生成局部网格坐标

function nodes = mesh_generate(domain, Nx, Ny, Q, xi_ref)
    % 输入：
    %   domain – 4×2 矩阵，四边形域四个顶点 [xA,yA; xB,yB; xC,yC; xD,yD]
    %   Nx, Ny – 将域在 ξ 和 η 方向上均匀划分的单元数
    %   Q     – DGSEM 每方向 GLL 点数
    % 输出：
    %   nodes  – Q×Q×2×Ncells 数组，
    %            nodes(i,j,1:2,c) = [x,y] 是第 c 个单元序号坐标为 (i,j) 的 GLL 点的物理坐标
    %            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

    % 1. 参考域 GLL 节点

    [Eta0, Xi0] = meshgrid(xi_ref, xi_ref); % ij对调改动
    % Xi0(i,j)=xi_ref(i)、Eta0(i,j)=xi_ref(j)

    Xi0  = Xi0(:);
    Eta0 = Eta0(:);

    Ncells = Nx * Ny;

    % 2. 预分配
    nodes = zeros(Q, Q, 2, Ncells);

    % 3. 遍历每个单元
    cell_id = 0;
    for j = 1:Ny
        for i = 1:Nx
            cell_id = cell_id + 1;

           % N1_c, N2_c, N3_c, N4_c 是对应于P_A, P_B, P_C, P_D四个顶点的双线性基函数
            % u 和 v 的取值范围均为 [0, 1]
            
            N1_c = @(u, v) (1 - u) .* (1 - v);
            N2_c = @(u, v) u .* (1 - v);
            N3_c = @(u, v) u .* v;
            N4_c = @(u, v) (1 - u) .* v;

            P_A = domain(1,:);
            P_B = domain(2,:);
            P_C = domain(3,:);
            P_D = domain(4,:);

            u = (i-1)/Nx;
            v = (j-1)/Ny;
            domain_cell(1,:) = N1_c(u,v)*P_A + N2_c(u,v)*P_B + N3_c(u,v)*P_C + N4_c(u,v)*P_D;

            u = i/Nx;
            v = (j-1)/Ny;
            domain_cell(2,:) = N1_c(u,v)*P_A + N2_c(u,v)*P_B + N3_c(u,v)*P_C + N4_c(u,v)*P_D;

            u = i/Nx;
            v = j/Ny;
            domain_cell(3,:) = N1_c(u,v)*P_A + N2_c(u,v)*P_B + N3_c(u,v)*P_C + N4_c(u,v)*P_D;

            u = (i-1)/Nx;
            v = j/Ny;
            domain_cell(4,:) = N1_c(u,v)*P_A + N2_c(u,v)*P_B + N3_c(u,v)*P_C + N4_c(u,v)*P_D;

            % 双线性映射到物理域
            phys_pts = bilinear_map_global(Xi0, Eta0, domain_cell);  % Nloc×2

            % 重塑并存入 nodes
            phys_pts = reshape(phys_pts, [Q, Q, 2]);
            nodes(:, :, :, cell_id) = phys_pts;
        end
    end
end

function phys = bilinear_map_global(xi, eta, domain_cell)
    % xi: Q^2×1，每个点对应的xi坐标
    % eta : Q^2×1，每个点对应的eta坐标
    % domain_cell : 4×2，当前单元顶点 [x1,y1; x2,y2; x3,y3; x4,y4]，按逆时针分布

    x1  = domain_cell(1,1); y1 = domain_cell(1,2);
    x2  = domain_cell(2,1); y2 = domain_cell(2,2);
    x3  = domain_cell(3,1); y3 = domain_cell(3,2);
    x4  = domain_cell(4,1); y4 = domain_cell(4,2);

    % 双线性形函数
    N1 = (1 - xi).*(1 - eta)/4;
    N2 = (1 + xi).*(1 - eta)/4;
    N3 = (1 + xi).*(1 + eta)/4;
    N4 = (1 - xi).*(1 + eta)/4;

    phys_x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
    phys_y = N1*y1 + N2*y2 + N3*y3 + N4*y4;
    phys   = [phys_x, phys_y];
end
