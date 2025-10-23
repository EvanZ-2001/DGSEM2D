% main.m
% DGSEM 求解二维浅水方程主程序

clear; clc; close all;

%% 物理与数值参数
g        = 9.81;             % 重力加速度
domain   = [0, 0;            % 物理域四角坐标 [x,y]
            25, 0;
            25, 1;
            0, 1];

Nx       =  50;               % x 方向单元数
Ny       =  3;               % y 方向单元数
P        =  3;               % 多项式阶数 (GLL 点数 = P+1)
Q = P+1;                    % GLL 点数 = P+1
CFL      =  0.2;             % CFL 数

% dt = 0.001; % 时间步长
dt = 0.005; % 时间步长
% dt = 0.01; % 时间步长
% dt = 0.0001; % 时间步长
% dt = 0.00005; % 时间步长

% t_final  =  0.02;           % 终止时间
% t_final  =  10;           % 终止时间
t_final  =  100;           % 终止时间

% output_dt=  0.01;             % 可视化输出时间间隔
% output_dt=  0.02;             % 可视化输出时间间隔
output_dt=  0.05;             % 可视化输出时间间隔
% output_dt=  dt;             % 可视化输出时间间隔

Qin = 4.42;
% Qin = 0; 
%%TODO:把Qin的设置放在主程序，修改子程序，使得子程序通过调用参数获取Qin
% Warning:目前主程序main和子程序BoundaryState中都需要设置Qin，注意保持一致

%% 网格生成及基函数预处理
[xi_ref, w_ref] = GLLNodesAndWeights(Q); 

nodes = mesh_generate(domain, Nx, Ny, Q ,xi_ref); % nodes是全局积分点的坐标
%   nodes  – Q×Q×2×Ncells 数组，
%            nodes(i,j,1:2,c) = [x,y] 是第 c 个单元序号坐标为 (i,j) 的 GLL 点的物理坐标
%            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

dL_ref = DerivativeMatrix(Q,xi_ref);
%    dL_ref — N×N 矩阵，dL(i,j) = ℓ_j'(xi(i))，基函数导数

%% 映射及局部算子计算
[J, xi_x, xi_y, eta_x, eta_y] = mapping(nodes,xi_ref);
% 输出变量均为Q×Q×Ncells数组

%% 生成底面地形
zb = zb_generate(nodes);
% zb    – Q×Q×Ncells 底面高程数组
% zb(i,j,c) 第 c 个单元序号坐标为 (i,j) 的 GLL 点的底面高度
%            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

%% 初始化守恒量 U(i,j,comp,c)
Ncells = Nx*Ny;
% Np     = P + 1;
U      = zeros(Q, Q, 3,Ncells);  % 第三个维度 1→h, 2→hu, 3→hv
% U       – Q×Q×3×Ncells 数组，守恒量 [h, hu, hv]
% U(i,j,1:3,c) = 是第 c 个单元序号坐标为 (i,j) 的 GLL 点的三个守恒变量
%            注意：序号坐标为 (i,j)表示x坐标第i个、y坐标第j个

h0 = 2.0;  % 初始静水水面高度
for c = 1:Ncells
    % h初值设置为h0-zb，水平水面
    U(:,:,1,c) = h0 - zb(:,:,c);
    % hu, hv 初值设为 0

    % hu在上游边界设置为Q0，具体做法为寻找上游单元的第一列点
    ix = mod(c-1, Nx) + 1;
    iy = floor((c-1)/Nx) + 1;
    if ix == 1
        U(1,:,2,c) = Qin*ones(Q,1);
    end
end


%% 时间积分主循环 (RK3)
t            = 0;
next_output  = 0;

%% 调试选项
% step = 0;

%%
while t < t_final
    
    % output_data = false;%调试参数

     %% 输出 & 可视化
    if t >= next_output - 1e-8
        fprintf('t = %.2f / %.2f\n', t, t_final);
        postprocess(U, nodes, P, t,'output.gif',Nx, Ny,zb);
        next_output = next_output + output_dt;
    end

    output_data = true;%调试参数

    %  三阶 TVD RK3
    % step = 1; %调试参数
    % RHS1 = compute_RHS_debug(U, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb,step,t,'diag_out',output_data); % TODO：调试
    RHS1 = compute_RHS(U, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb);
    U1   = U + dt * RHS1;
    U1   = limiter_weno(U1, Nx, Ny, w_ref, dL_ref);

    % step = 2; %调试参数
    % RHS2 = compute_RHS_debug(U1, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb,step,t,'diag_out',output_data); % TODO：调试
    RHS2 = compute_RHS(U1, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb);
    U2   = 0.75*U + 0.25*(U1 + dt * RHS2);
    U2   = limiter_weno(U2, Nx, Ny, w_ref, dL_ref);

    % step = 3; %调试参数
    % RHS3 = compute_RHS_debug(U2, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb,step,t,'diag_out',output_data); % TODO：调试
    RHS3 = compute_RHS(U2, nodes, Nx, Ny, w_ref, J, dL_ref, xi_x, xi_y, eta_x, eta_y, g,zb);
    U    = (1/3)*U + (2/3)*(U2 + dt * RHS3);
    U    = limiter_weno(U, Nx, Ny, w_ref, dL_ref);

    t = t + dt;



    % 调试显示
    % h_ele1 = squeeze(U(1,:,:,1));
    % disp(h_ele1);
end

fprintf('计算完成，t = %.2f\n', t_final);
