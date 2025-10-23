%% alpha=beta=1的P阶Jacobi多项式的导数值

function ddx = DxJacobiPolynomial(x,P)
    if P ==0
        ddx=0;
        return
    end
    b1 = 2*(P+1)*(1-x^2);
    b2 = -2*P*(P+1)*x;
    b3 = 2*(P+1)^2;
    ddx = (b2*JacobiPolynomial(x,P) + b3*JacobiPolynomial(x,P-1) )/b1;
end
