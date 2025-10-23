function UR = BoundaryState(UL, faceN)
% 构造边界外的虚拟状态 UR

    UR = UL;  % 默认镜像状态
    h = UL(1); hu = UL(2); hv = UL(3);
    u = hu / h; v = hv / h;

    BCNum = 1; %选择实施的边界条件，1或2

    if BCNum == 1 %上游给定流量，下游给定高度的边界条件
        switch faceN
            case 4  % x=0：上游流量
                % 1）正常情形，Q=4.42
                Qin = 4.42; 
                % 2）测试情形，Q=0
                % Qin = 0;
                
                % h_in = 2;
                % u = Qin / h_in; v = 0;
                % UR(:) = [h_in, h_in*u, h_in*v];
                UR(:) = [h, Qin, 0];  %% DEBUG:上游边界条件的高度和y方向流量如何设置
                
            case 2  % x=25：下游水深
                h_out = 2;
                % UR(:) = [h_out, h_out*u, h_out*v]; %% DEBUG:下游出流条件是设置流量与内点相同还是设置速度与内点相同
                UR(:) = [h_out, hu, hv];
            case {1,3}  % y边界：滑移
                UR(:) = [h, hu, -hv];
        end
    end

    if BCNum == 2 %各个避免都是反射边界条件，用于测试静止算例
        % 所有边界都是滑移墙
        if faceN == 1 || faceN == 3 % 上下边界
            UR = [h; hu; -hv];
        else % 左右边界
            UR = [h; -hu; hv];
        end
    end
end
