% TODO:检查是否正确
function [F, G] = physical_flux(U, g)
    % PHYSICAL_FLUX 计算二维浅水方程在 x 和 y 方向的物理通量
    %
    %   [F, G] = physical_flux(U, g)
    %
    % 输入：
    %   U(:,:,1) = h  — 水深
    %   U(:,:,2) = hu — x 方向动量
    %   U(:,:,3) = hv — y 方向动量
    %   g         — 重力加速度
    %
    % 输出：
    %   F(:,:,1) = hu
    %   F(:,:,2) = hu.^2./h + 0.5*g*h.^2
    %   F(:,:,3) = hu.*hv./h
    %
    %   G(:,:,1) = hv
    %   G(:,:,2) = hu.*hv./h
    %   G(:,:,3) = hv.^2./h + 0.5*g*h.^2
    
    % 提取守恒量
    h  = U(:,:,1);
    hu = U(:,:,2);
    hv = U(:,:,3);
    
    % 计算速度分量
    u = hu ./ h;
    v = hv ./ h;
    
    % 预分配通量张量
    [np_i, np_j, ~] = size(U);
    F = zeros(np_i, np_j, 3);
    G = zeros(np_i, np_j, 3);
    
    % x 方向通量 F
    F(:,:,1) = hu;
    F(:,:,2) = hu .* u + 0.5 * g * h.^2;
    F(:,:,3) = hu .* v;
    
    % y 方向通量 G
    G(:,:,1) = hv;
    G(:,:,2) = hu .* v;
    G(:,:,3) = hv .* v + 0.5 * g * h.^2;
end
