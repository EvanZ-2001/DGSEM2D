function postprocess(U, nodes, P, t, outputGif, Nx, Ny, zb)
% POSTPROCESS  绘制若干条 y=const 的水平剖面曲线（每条对应一个 y 方向 GLL 点）
% 在每个时间步调用以生成动图帧。曲线显示 η = h + z_b 随 x 的分布，并随时间更新。
%
% 输入：
%   U       – Q×Q×3×Ncells 数组，守恒量 [h, hu, hv]
%   nodes   – Q×Q×2×Ncells 数组，nodes(i,j,1,c)=x, nodes(i,j,2,c)=y
%   P     – 多项式阶数（Q = P+1）
%   t     – 当前时间
%   outputGif – 输出 GIF 文件名
%   Nx, Ny – 单元数（x/y 方向），要求 Ncells = Nx*Ny
%   zb    – Q×Q×Ncells 底面高程数组

    % —— 动图初始化：首次（或极早）时间步删旧文件 ——
    if t <= 1e-5
        if exist(outputGif, 'file')
            delete(outputGif);
        end
    end

    % 尺寸
    [Q, ~, ~, Ncells] = size(U);

    % 选取 y 方向的哪一排单元（Ny=1：唯一一排；Ny>1：取中间那一排）
    iy_sel = min(max(1, ceil(Ny/2)), Ny);

    % 统一纵轴范围（可按需调整）
    hmin_assum = 0.0;
    hmax_assum = 3.0;

    % 每条曲线的全局采样点数（沿 x 串接各单元）
    Ncols_line = Nx*(Q-1) + 1;

    % ——— 绘图 ———
    fig = figure(100);
    % set(fig, 'Units','pixels', 'Position',[100,100,1600,600]);
    set(fig, 'Units','pixels', 'Position',[100,100,1200,900]);
    clf; hold on; grid on; box on;

    % x 轴范围：改为全域 nodes 的最小/最大，避免使用单条探针可能的截断
    x_all = nodes(:,:,1,:);
    xlim([min(x_all(:)), max(x_all(:))]);

    % 逐个 y 向 GLL 点 j0 生成曲线
    curveHandles = gobjects(Q,1);
    % zb_legend_handle = [];   % 新增：用于图例的底面句柄（只保留一个）
    for j0 = 1:Q
        x_line   = zeros(Ncols_line, 1);
        eta_line = zeros(Ncols_line, 1);   % η = h + z_b
        zb_line  = zeros(Ncols_line, 1);   % 新增：底面 z_b
        k = 1;

        for ix = 1:Nx
            c = (iy_sel-1)*Nx + ix;  % 当前单元编号（该排里第 ix 个单元）

            % 取该单元在第 j0 行上的 x、η、z_b
            x_cell   = reshape(nodes(:,j0,1,c),[],1);
            h_cell   = reshape(U(:,j0,1,c),[],1);
            zb_cell  = reshape(zb(:,j0,c),[],1);
            eta_cell = h_cell + zb_cell;

            % 若该单元 x 序列反向，做同步反转，保证从左到右拼接
            if x_cell(1) > x_cell(end)
                x_cell   = flipud(x_cell);
                eta_cell = flipud(eta_cell);
                zb_cell  = flipud(zb_cell);
            end

            if ix > 1
                % 从第二个单元起，跳过 i=1（与前一单元右边界重复点）
                x_line(k:(k+Q-2))    = x_cell(2:Q);
                eta_line(k:(k+Q-2))  = eta_cell(2:Q);
                zb_line(k:(k+Q-2))   = zb_cell(2:Q);   % 新增
                k = k + (Q-1);
            else
                % 第一个单元保留全部 Q 个点
                x_line(k:(k+Q-1))    = x_cell(1:Q);
                eta_line(k:(k+Q-1))  = eta_cell(1:Q);
                zb_line(k:(k+Q-1))   = zb_cell(1:Q);   % 新增
                k = k + Q;
            end
        end

        % —— 先画底面，再画水面（让水面盖在上面、便于辨识）——
        % 仅把第一条底面线加入图例，其余底面线隐藏以避免图例过长
        % if j0 == 1
        %     zb_legend_handle = plot(x_line, zb_line, 'k-', 'LineWidth', 1.0, ...
        %                             'DisplayName', 'z_b (底面)');
        % else
        %     plot(x_line, zb_line, 'k-', 'LineWidth', 1.0, ...
        %          'HandleVisibility', 'off');
        % end

        % 底面曲线不给图例
        plot(x_line, zb_line, 'k-', 'LineWidth', 1.0, ...
                 'HandleVisibility', 'off');

        % 用该条曲线的 y 值做标签（取该行 y 的中位数以避免轻微数值抖动）
        y_line_all = zeros(Nx*Q - (Nx-1), 1);
        kk = 1;
        for ix = 1:Nx
            c = (iy_sel-1)*Nx + ix;
            y_cell = reshape(nodes(:,j0,2,c),[],1);
            if y_cell(1) > y_cell(end)
                y_cell = flipud(y_cell);
            end
            if ix > 1
                y_line_all(kk:(kk+Q-2)) = y_cell(2:Q);
                kk = kk + (Q-1);
            else
                y_line_all(kk:(kk+Q-1)) = y_cell(1:Q);
                kk = kk + Q;
            end
        end
        y_tag = median(y_line_all,'omitnan');

        % 水面曲线（带图例 y 标签）
        curveHandles(j0) = plot(x_line, eta_line, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('y = %.6g', y_tag));
        % plot(x_line, eta_line, 'LineWidth', 1.5, ...
        % 'DisplayName', sprintf('y = %.6g', y_tag));
    end


    % 统一坐标与标题
    xlabel('x');
    ylabel('自由液面高程 \eta = h + z_b');
    ylim([hmin_assum, hmax_assum]);

    % if ~isempty(zb_legend_handle)
    %     legend([curveHandles; zb_legend_handle], 'Location', 'northeastoutside');
    % else
    %     legend(curveHandles, 'Location', 'northeastoutside');
    % end
    legend(curveHandles, 'Location', 'northeastoutside');

    title(sprintf('固定 y 的水平剖面：\\eta 随 x 变化  (t = %.4f s, 选用 iy = %d)', t, iy_sel));
    drawnow;

    % —— 捕捉当前帧并写入 GIF ——
    frame = getframe(fig);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);

    if ~exist(outputGif, 'file')
        imwrite(A, map, outputGif, 'gif', 'LoopCount', Inf, 'DelayTime', 0.05);
    else
        imwrite(A, map, outputGif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end
