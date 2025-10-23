% main.m
% DGSEM 求解二维浅水方程主程序

clear; clc; close all;

%% 物理与数值参数
g        = 9.81;             % 重力加速度
domain   = [0, 0;            % 物理域四角坐标 [x,y]
            25, 0;
            25, 1;
            0, 1];

Nx       =  20;               % x 方向单元数
Ny       =  4;               % y 方向单元数
P        =  3;               % 多项式阶数 (GLL 点数 = P+1)
Q = P+1;                    % GLL 点数 = P+1
CFL      =  0.2;             % CFL 数
dt = 0.001; % 时间步长
% dt = 0.0001; % 时间步长
% dt = 0.00005; % 时间步长
t_final  =  10;             % 终止时间
% output_dt=  0.005;             % 可视化输出时间间隔
% output_dt=  0.005;             % 可视化输出时间间隔
output_dt=  dt;             % 可视化输出时间间隔

%% 网格生成及基函数预处理
% [nodes, elems] = mesh_generate(domain, Nx, Ny, P);
[xi_ref, w_ref] = GLLNodesAndWeights(Q); 

nodes = mesh_generate(domain, Nx, Ny, Q ,xi_ref); % nodes是全局积分点的坐标
% nodes  – Ncells×Q×Q×2 数组，
%            nodes(c,i,j,1:2) = [x,y] 是第 c 个单元序号坐标为 (i,j) 的 GLL 点的物理坐标
%            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

% [L_ref, dL_ref] = lagrange_basis(xi_ref);
dL_ref = DerivativeMatrix(Q,xi_ref);
%    dL_ref — N×N 矩阵，dL_ref(i,j) = ℓ_j'(xi(i))，基函数导数

%% 映射及局部算子计算
[J, xi_x, xi_y, eta_x, eta_y] = mapping(nodes,xi_ref);

%% 生成底面地形
zb = zb_generate(nodes);
% 三者均为Ncells×Q×Q 数组
% zb(c,i,j) 第 c 个单元序号坐标为 (i,j) 的 GLL 点的底面高度
%            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

%% 初始化守恒量 U(c,i,j,comp)
Ncells = size(nodes,1);
% Np     = P + 1;
U      = zeros(Ncells, Q, Q, 3);  % 最后维度 1→h, 2→hu, 3→hv
% U(c,i,j,1:3) = 是第 c 个单元序号坐标为 (i,j) 的 GLL 点的三个守恒变量
%            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

h0 = 2.0;  % 初始静水水面高度
for c = 1:Ncells
    U(c,:,:,1) = h0 - squeeze(zb(c,:,:));
    % hu, hv 初值设为 0
end

%% 时间积分主循环 (RK3)
t            = 0;
next_output  = 0;

%% 调试选项
step = 0;

%%
while t < t_final
    
    % output_data = false;%调试参数

    %  输出 & 可视化
    if t >= next_output - 1e-8
        fprintf('t = %.2f / %.2f\n', t, t_final);
        postprocess(U, nodes, P, t,'output.gif',Nx, Ny, zb)
        next_output = next_output + output_dt;

        % output_data = true;%调试参数

    end

    % % 单步测试
    % RHS1 = compute_RHS(U, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb_x,zb_y);
    % U   = U + dt * RHS1;

    %  三阶 TVD RK3
    % step = 1; %调试参数
    % RHS1 = compute_RHS_debug(U, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb_x,zb_y,step,t,'diag_out',output_data); % TODO：调试
    RHS1 = compute_RHS(U, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g, zb);
    U1   = U + dt * RHS1;
    U1   = limiter_weno(U1, Nx, Ny, w_ref, dL_ref);

    % step = 2; %调试参数
    % RHS2 = compute_RHS_debug(U1, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb_x,zb_y,step,t,'diag_out',output_data); % TODO：调试
    RHS2 = compute_RHS(U1, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g, zb);
    U2   = 0.75*U + 0.25*(U1 + dt * RHS2);
    U2   = limiter_weno(U2, Nx, Ny, w_ref, dL_ref);

    % step = 3; %调试参数
    % RHS3 = compute_RHS_debug(U2, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb_x,zb_y,step,t,'diag_out',output_data); % TODO：调试
    RHS3 = compute_RHS(U2, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g, zb);
    U    = (1/3)*U + (2/3)*(U2 + dt * RHS3);
    U    = limiter_weno(U, Nx, Ny, w_ref, dL_ref);

    t = t + dt;



    % 调试显示
    % h_ele1 = squeeze(U(1,:,:,1));
    % disp(h_ele1);
end

fprintf('计算完成，t = %.2f\n', t_final);
