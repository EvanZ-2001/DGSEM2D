function F = numerical_flux_LF_Roe(UL, UR, n, g)
% NUMERICAL_FLUX_LF_ROE  Roe‐平均速度下的 Lax–Friedrichs 数值通量
%   输入：
%     UL, UR — 长度 3 列向量 [h; hu; hv]
%     n      — 2×1 单位法向量
%     g      — 重力加速度
%   输出：
%     F      — 3×1 列向量，数值通量

  % 提取左/右态
  hL  = UL(1); huL = UL(2); hvL = UL(3);
  hR  = UR(1); huR = UR(2); hvR = UR(3);
  uL  = huL / hL; vL = hvL / hL;
  uR  = huR / hR; vR = hvR / hR;

  % 法向速度分量
  u_barL = uL * n(1) + vL * n(2);
  u_barR = uR * n(1) + vR * n(2);

  % 切向速度分量
  v_barL = -uL*n(2) + vL*n(1);
  v_barR = -uR*n(2) + vR*n(1);

  % >>> 修正 ：将左右状态向量 UL, UR 旋转到 (n, t) 坐标系 TODO：调试
  UL_rot = [hL; hL * u_barL; hL * v_barL];
  UR_rot = [hR; hR * u_barR; hR * v_barR];

  % 物理通量在法向方向上的投影
  FL = [ hL * u_barL;
         hL * u_barL^2 + 0.5*g*hL^2;
         hL * u_barL * v_barL ];
  FR = [ hR * u_barR;
         hR * u_barR^2 + 0.5*g*hR^2;
         hR * u_barR * v_barR ];

  % Roe 平均态
  sqrt_hL = sqrt(hL);
  sqrt_hR = sqrt(hR);
  hRoe    = sqrt_hL * sqrt_hR;
  u_barRoe    = (sqrt_hL*u_barL + sqrt_hR*u_barR) / (sqrt_hL + sqrt_hR);
  v_barRoe    = (sqrt_hL*v_barL + sqrt_hR*v_barR) / (sqrt_hL + sqrt_hR);
  % unRoe   = u_Roe*n(1) + vRoe*n(2);
  cL = sqrt(g * hL);
  cR = sqrt(g * hR);
  cRoe    = sqrt(0.5 * (cL^2 + cR^2));

  % Lax–Friedrichs 粘滞系数
  alpha = abs(u_barRoe) + cRoe;

  % LF 数值通量
  % Fn = 0.5*(FL + FR) - 0.5*alpha*(UR - UL);
  Fn = 0.5*(FL + FR) - 0.5*alpha*(UR_rot - UL_rot);

  % 将法向/切向通量还原到全局坐标 (x,y)
  % t    = [-n(2); n(1)];  % 切向单位矢量
  F    = zeros(3,1);
  F(1) = Fn(1);
  F(2) = Fn(2)*n(1) - Fn(3)*n(2);
  F(3) = Fn(2)*n(2) + Fn(3)*n(1);
end
