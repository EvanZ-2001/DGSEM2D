function zb = zb_generate(nodes)
% ZB_GENERATE  根据给定的底面表达式计算底面高程
%
%   zb = zb_generate(nodes)
%
% 输入：
%   nodes – Q×Q×2×Ncells 数组，nodes(i,j,1,c)=x，nodes(i,j,2,c)=y
%
% 输出：
%   zb    – Q×Q×Ncells 底面高程数组
%

    [Q, ~, ~, Ncells] = size(nodes);

    % 提取所有单元内的物理坐标（保持原维度）
    X = reshape(nodes(:,:,1,:), Q, Q, Ncells);
    Y = reshape(nodes(:,:,2,:), Q, Q, Ncells);

    % 选择需要的底面表达式
    zb_fun = @Zb_Expression5;  % 例如选择第5种底面类型
    % zb_fun = @Zb_Expression1;
    % zb_fun = @Zb_Expression2;
    % zb_fun = @Zb_Expression3;

    % 直接向量化调用底面函数
    zb = zb_fun(X, Y);

    % 尺寸检查（可选）
    if ~isequal(size(zb), [Q, Q, Ncells])
        error('zb_fun 输出尺寸必须与输入 (Q×Q×Ncells) 一致');
    end
end

%% ------------------------------------------------------------------------
%% 第1种底面情况：平底
function zb = Zb_Expression1(x, y)
    zb = zeros(size(x));
end

%% 第2种底面情况：平面倾斜底
function zb = Zb_Expression2(x, y)
    zb = 0.01*x + 0.02*y;
end

%% 第3种底面情况：三角函数底
function zb = Zb_Expression3(x, y)
    zb = 0.25 + 0.25 * cos(pi * (x - 12.5) / 12.5);
end

%% 第5种底面情况：无穷阶光滑二次凸起
function zb = Zb_Expression5(x, y) %#ok<INUSD>
    x1 = 7.5;
    x2 = 12.5;
    a  = 9.37;

    zb = zeros(size(x));
    mask = (x > x1) & (x < x2);

    % 仅在中间区间计算指数项
    if any(mask(:))
        c0 = 0.2 * exp((4*a)/(x2 - x1)^2);
        xm = x(mask);
        zb(mask) = c0 .* exp(-a ./ ((xm - x1) .* (x2 - xm)));
    end
end
