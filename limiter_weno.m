function U_lim = limiter_weno(U, Nx, Ny, w, dL_ref)
% LIMITER_WENO  基于 Zhong & Shu (2013) 的 WENO 限制器，先用 minmod 识别问题单元
%
%   U_lim = limiter_weno(U, nodes, w, L, elems)
%
% 输入：
%   U     – Q×Q×3×Ncells 数组，守恒量 [h, hu, hv]
%   Nx, Ny – x方向和y方向的单元数
%   w     – Q×1 列向量，GLL 权重
%   dL_ref    – Q×Q 参考域一维基函数导数矩阵
%
% 输出：
%   U_lim – 限制后的解，同尺寸 U
%
% 方法概述：
%   1. 对每个单元，计算 cell-average，与左右、上下邻居的平均值做 minmod 判断，
%      若在 ξ 或 η 方向上出现局部极值（minmod ≤ 0），则该单元为问题单元。
%   2. 仅对问题单元做五模板 WENO 重构（Zhong & Shu 简单 WENO 限制器），保持平均不变。

[Q, ~, nVar, Ncells] = size(U);
U_lim = U;

troubled_count = 0;        % 统计本步被标记为troubled的单元数

% 一维 GLL 权重矩阵
WiWj = w * w.';
Wsum = sum(WiWj(:));

% 一维 GLL 权重
W1   = w(:).';                  % 1xQ
Wsum1= sum(W1);

% 线性权重（中间大、邻居小）
if Q==3
    gamma = [0.996;0.001;0.001;0.001;0.001];
else
    gamma = [0.9996;0.0001;0.0001;0.0001;0.0001];
end
eps = 1e-6;

for c = 1:Ncells
    % 单元在网格中的位置
    ix = mod(c-1, Nx) + 1;
    iy = floor((c-1)/Nx) + 1;
    % 邻居 ID（边界时重复自己）
    idL = (ix>1)  *(c-1) + (ix==1)  *c;
    idR = (ix<Nx) *(c+1) + (ix==Nx) *c;
    idB = (iy>1)  *(c-Nx) + (iy==1)  *c;
    idT = (iy<Ny) *(c+Nx) + (iy==Ny) *c;
    
    % 计算本单元及邻居平均
    avg0 = zeros(nVar,1);
    avgL = avg0; avgR = avg0; avgB = avg0; avgT = avg0;
    for v = 1:nVar
      Uc   = reshape(U(:,:,v,c),[Q,Q]);
      avg0(v) = sum(sum(WiWj .* Uc))/ Wsum;
      avgL(v) = sum(sum(WiWj .* reshape(U(:,:,v,idL),[Q,Q])))/ Wsum;
      avgR(v) = sum(sum(WiWj .* reshape(U(:,:,v,idR),[Q,Q])))/ Wsum;
      avgB(v) = sum(sum(WiWj .* reshape(U(:,:,v,idB),[Q,Q])))/ Wsum;
      avgT(v) = sum(sum(WiWj .* reshape(U(:,:,v,idT),[Q,Q])))/ Wsum;
    end
    
    % % minmod 判别：若在任一方向出现局部极值，则标记问题单元
    % troubled = false;
    % for v = 1:nVar
    %   dL = avg0(v) - avgL(v);
    %   dR = avgR(v) - avg0(v);
    %   dB = avg0(v) - avgB(v);
    %   dT = avgT(v) - avg0(v);
    %   if dL*dR <= 0 || dB*dT <= 0
    %     troubled = true;
    %     break;
    %   end
    % end

    %% 修改后的TVB minmod判别
    % ---------- TVB troubled-cell indicator ----------
    % 面平均（用 GLL 权重做 1D 积分），以水深 h 为主检测量；需要可对每个分量或特征空间分别检测
    
    % 方向尺度（均匀网格下用 1/Nx, 1/Ny；非均匀网格建议传入每单元的 hchar(c)）
    hx = 25.0 / Nx; 
    hy = 1.0 / Ny;
    hchar_x = hx; 
    hchar_y = hy;

    
    % TVB 参数（可调）
    % M_TVB = 50; 
    % M_TVB = 20; 
    M_TVB = 0.5;
    % M_TVB = 0.1; 
    % M_TVB = 0.01;
    
    % 面平均：左/右边（沿 y 做 1D 加权平均）；下/上边（沿 x 做 1D 加权平均）
    % 注意：sum_i w_i * U(1,i) / sum_i w_i 是左边面平均
    uL = zeros(1,nVar); uR = uL; uB = uL; uT = uL;
    for v = 1:nVar
        % Uc = squeeze(U(c,:,:,v));        % QxQ
        Uc = reshape(U(:,:,v,c),[Q,Q]);  % QxQ
        % uL(v) = (W1 * squeeze(Uc(1,:)).') / Wsum1;    % i=1 列，沿 j 加权
        % uR(v) = (W1 * squeeze(Uc(Q,:)).') / Wsum1;    % i=Q 列
        % uB(v) = (W1 * squeeze(Uc(:,1)))   / Wsum1;    % j=1 行，沿 i 加权
        % uT(v) = (W1 * squeeze(Uc(:,Q)))   / Wsum1;    % j=Q 行
        uL(v) = (W1 * reshape(Uc(1,:),[],1)) / Wsum1;    % i=1 列，沿 j 加权
        uR(v) = (W1 * reshape(Uc(Q,:),[],1)) / Wsum1;    % i=Q 列
        uB(v) = (W1 * reshape(Uc(:,1),[],1)) / Wsum1;    % j=1 行，沿 i 加权
        uT(v) = (W1 * reshape(Uc(:,Q),[],1)) / Wsum1;    % j=Q 行
    end
    
    troubled = false;
    for v = 1:nVar
        % 与相邻单元平均值的“跳跃”作参照
        % x 方向：用右/左边界处(单元内)的面平均与单元平均之差，和相邻单元平均差做 TVB 测试
        aR = (uR(v) - avg0(v));                   % 右边界“斜率”
        aL = (avg0(v) - uL(v));                   % 左边界“斜率”
        bR = (avgR(v) - avg0(v));                 % 右邻与本单元平均的差
        bL = (avg0(v) - avgL(v));                 % 本单元与左邻平均的差
    
        thR = mmTVB(aR, [bR, bL], hchar_x, M_TVB);
        thL = mmTVB(aL, [bR, bL], hchar_x, M_TVB);
    
        % y 方向同理
        aT = (uT(v) - avg0(v));
        aB = (avg0(v) - uB(v));
        bT = (avgT(v) - avg0(v));
        bB = (avg0(v) - avgB(v));

        thT = mmTVB(aT, [bT, bB], hchar_y, M_TVB);
        thB = mmTVB(aB, [bT, bB], hchar_y, M_TVB);
    
        % 只要有任意一边“会被 TVB-minmod 修改”，就标记为 troubled
        if (thR ~= aR) || (thL ~= aL) || (thT ~= aT) || (thB ~= aB)
            troubled = true; 
            break;
        end
    end
    % ---------- 结束 TVB 指标 ----------

    
    if ~troubled
      continue;
    else
        troubled_count = troubled_count + 1;   % 记一次
    end
    
    % 对问题单元应用五模板 WENO 重构
    % 五模板：中心、左、右、下、上
    nbrs = [c, idL, idR, idB, idT];
    US   = cell(5,1);
    for l = 1:5
      US{l} = zeros(Q,Q,nVar);
      Uc = reshape(U(:,:,:,nbrs(l)),[Q,Q,nVar]);         % Q×Q×nVar
      % 在本单元参考点插值(模板多项式)
      for v = 1:nVar
        US{l}(:,:,v) = Uc(:,:,v);
      end
      % 保持平均
      for v = 1:nVar
        avg_l = sum(sum(WiWj .* US{l}(:,:,v)))/ Wsum;
        US{l}(:,:,v) = US{l}(:,:,v) + (avg0(v)-avg_l);
      end
    end
    
    % 计算平滑性指标 β_l
    beta = zeros(5,1);
    for l = 1:5
      for v = 1:nVar
        Uc = US{l}(:,:,v);
        % 求两个方向的导数
        dU_dxi  = dL_ref * Uc;          % Q×Q
        dU_deta = Uc * dL_ref.';        % Q×Q
        beta(l) = beta(l) + sum(sum(WiWj .* (dU_dxi.^2 + dU_deta.^2)));
      end
    end
    
    % 非线性权重
    alpha = gamma ./ ( (eps + beta).^2 );
    w_nl  = alpha / sum(alpha);
    
    % 重构
    for v = 1:nVar
      Unew = zeros(Q,Q);
      for l = 1:5
        Unew = Unew + w_nl(l) * US{l}(:,:,v);
      end
      U_lim(:,:,v,c) = Unew;
    end
end

troubled_ratio = troubled_count / Ncells;
fprintf('[WENO] troubled cells: %.2f%%  (%d / %d)\n',100*troubled_ratio, troubled_count, Ncells);

end

%% TVB判别函数
%% 无量纲M TVB判别函数
function theta = mmTVB(a, vec, h, M)
    if abs(a) <= M*h*h
        theta = a;
    else
        theta = minmod([a, vec(:).']);
    end
end

%% 固定值M TVB判别函数
% function theta = mmTVB(a, vec, h, M)
%     if abs(a) <= M
%         theta = a;
%     else
%         theta = minmod([a, vec(:).']);
%     end
% end

%% minmod函数
function m = minmod(vals)
    if all(vals > 0)
        m = min(abs(vals));          % 等价于 min(vals) 对于正数
    elseif all(vals < 0)
        m = -min(abs(vals));
    else
        m = 0;
    end
end
