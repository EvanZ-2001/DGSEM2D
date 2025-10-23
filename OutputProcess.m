
Ncells = meta.Ncells;
Nx = meta.Nx;
Ny = meta.Ny;
Q = meta.Q;

for c = 1:Ncells
    % 计算单元索引
    ix = mod(c-1, Nx) + 1;
    iy = floor((c-1)/Nx) + 1;
    for i = 1:Q
      for j = 1:Q
        row = (iy-1)*Q + j;
        col = (ix-1)*Q + i;
        xi_x_g(row,col) = xi_x(c,i,j);
        xi_y_g(row,col) = xi_y(c,i,j);
        eta_x_g(row,col) = eta_x(c,i,j);
        eta_y_g(row,col) = eta_y(c,i,j);
        % nodes_g_x(row,col) = nodes(c,i,j,1);
        % nodes_g_y(row,col) = nodes(c,i,j,2);
        for k=1:3
            T1_g(row,col,k) = T1(c,i,j,k);
            T2_g(row,col,k) = T2(c,i,j,k);
            T3_g(row,col,k) = T3(c,i,j,k);
            T4_g(row,col,k) = T4(c,i,j,k);
            RHS_g(row,col,k) = RHS(c,i,j,k);
            Vtot_g(row,col,k) = Vtot(c,i,j,k);
            divV_g(row,col,k) = divV(c,i,j,k);
            numF_g(row,col,k) = numF(c,i,j,k);
            S_g(row,col,k) = S(c,i,j,k);
        end
        numF_g_1(row,col) = numF(c,i,j,1);
        numF_g_2(row,col) = numF(c,i,j,2);
        numF_g_3(row,col) = numF(c,i,j,3);
        divV_g_1(row,col) = divV(c,i,j,1);
        divV_g_2(row,col) = divV(c,i,j,2);
        divV_g_3(row,col) = divV(c,i,j,3);
        Vtot_g_1(row,col) = Vtot(c,i,j,1);
        Vtot_g_2(row,col) = Vtot(c,i,j,2);
        Vtot_g_3(row,col) = Vtot(c,i,j,3);
        RHS_g_1(row,col) = RHS(c,i,j,1);
        RHS_g_2(row,col) = RHS(c,i,j,2);
        RHS_g_3(row,col) = RHS(c,i,j,3);
      end
    end
end