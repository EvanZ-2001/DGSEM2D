% 实现思路：物理通量和边界数值通量的计算是h、hu、hv三项一起计算的，但物理通量右端项、除以质量矩阵需要h、hu、hv分别计算

function RHS = compute_RHS(U, nodes, Nx, Ny, w, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g, zb)
    % compute_RHS 按矩阵对应项相乘的离散公式构造 RHS
    %
    % 输入：
    %   U       – Q×Q×3×Ncells 数组，守恒量 [h, hu, hv]
    %   nodes   – Q×Q×2×Ncells 数组，nodes(i,j,1,c)=x, nodes(i,j,2,c)=y（仅用于数值通量）
    %   w       – (Q×1) GLL 权重
    %   J       – Q×Q×Ncells 雅可比行列式
    %   dL_ref      – Q×Q 参考域一维基函数导数矩阵
    %   xi_x, xi_y, eta_x, eta_y – Q×Q×Ncells 映射导数
    %   g       – 重力加速度
    %   zb    – Q×Q×Ncells 底面高程数组
    %
    % 输出：
    %   RHS     – Q×Q×3×Ncells，每个分量的 RHS
    
    [Q, ~, ~, Ncells] = size(U);
    RHS = zeros(size(U));
    
    % 质量对角：w_i * w_j, 同样用于所有单元
    WiWj = w(:) * w(:).';

    for c = 1:Ncells
        % 当前单元解
        Uc = U(:,:,:,c);               % Q×Q×3

        % 本单元的Jacobi行列式
        J_c = J(:,:,c);

        % 本单元的映射导数
        xi_x_c = xi_x(:,:,c);
        eta_x_c = eta_x(:,:,c);
        xi_y_c = xi_y(:,:,c);
        eta_y_c = eta_y(:,:,c);

        % 物理通量
        [F, G] = physical_flux(Uc, g);  % Q×Q×3

        % 界面数值通量
        numF = numerical_flux_all(c, U, nodes, Nx, Ny, w, g); % Q×Q×3
        %function numF = numerical_flux_all(cell_id, U, nodes, Nx, Ny, w, g)

        % 计算源项
        % 本单元的地形高程
        zb_c = zb(:,:,c); % QxQ 矩阵
        
        T1_zb = dL_ref * zb_c;
        T2_zb = zb_c * dL_ref.';
        T3_zb = dL_ref * zb_c;
        T4_zb = zb_c * dL_ref.';
        
        % 这就是离散的 g*h*∇zb
        grad_zb_x_discrete = g * J_c .* WiWj .* Uc(:,:,1) .* (xi_x_c .* T1_zb + eta_x_c .* T2_zb);
        grad_zb_y_discrete = g * J_c .* WiWj .* Uc(:,:,1) .* (xi_y_c .* T3_zb + eta_y_c .* T4_zb);
        
        % 构造新的源项 S_c (注意，它现在是一个QxQ矩阵，不再是QxQx3)
        Source = zeros(Q, Q, 3);
        Source(:,:,2) = grad_zb_x_discrete; % -g*h*dz/dx 变为 +∫(g*h*zb*φ_x)
        Source(:,:,3) = grad_zb_y_discrete; % -g*h*dz/dy 变为 +∫(g*h*zb*φ_y)
        % --- 结束新的源项计算 ---
        
        % 对每个守恒量分量构造离散公式
        for comp = 1:3

            % 体积分项
            F_c = F(:,:,comp);
            G_c = G(:,:,comp);

            W = diag(w);  % Q×Q 权重对角阵
            dLw = W * dL_ref;  % (W)×dL_ref ，用于 η 项的权重 w_q

            % —— 4. 四项求和 via 矩阵乘 —— 
            WF = W * F_c;        % 每行乘 w_p，用于 ξ 项
            WG = W * G_c;
            
            % T1(i,j) = w_j ∑_p w_p f(U_pj) h_i'(ξ_p)
            T1 = dL_ref.' * WF * W;    
            
            % T2(i,j) = w_i ∑_q w_q f(U_iq) h_j'(η_q)
            T2 = W * (F_c * dLw);       
            
            % T3(i,j) = w_j ∑_p w_p g(U_pj) h_i'(ξ_p)
            T3 = dL_ref.' * WG * W;    
            
            % T4(i,j) = w_i ∑_q w_q g(U_iq) h_j'(η_q)
            T4 = W * (G_c * dLw);       

            %  对应公式中的 |J| ∂ξ/∂x·T1 + |J| ∂η/∂x·T2 + |J| ∂ξ/∂y·T3 + |J| ∂η/∂y·T4
            divV = J_c .* ( xi_x_c .* T1 ...
                        + eta_x_c .* T2 ...
                        + xi_y_c .* T3 ...
                        + eta_y_c .* T4 );

            % 面数值通量（已包含积分权重）
            numF_c = numF(:,:,comp);

            % 源项
            S_c = Source(:,:,comp);
            
            % 合成 RHS 并除以质量对角
            % 合成 RHS
            % 核心变化在这里：用 Vtot = divV - S_c - numF_c
            Vtot = divV - S_c - numF_c; % 注意 S_c 前面的符号是负号
            Mdiag = J(:,:,c) .* WiWj;     % Q×Q
            RHS(:,:,comp,c) = Vtot ./ Mdiag;
        end
    end
    
end

