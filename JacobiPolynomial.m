%% alpha=beta=1的P阶Jacobi多项式的值

function y = JacobiPolynomial(x,P)
    if P == 0
        y=1;
        return;
    end
    if P == 1
        y=2*x;
        return;
    end
    yVec = zeros(P+1,1);
    yVec(1) = 1;
    yVec(2) = 2*x;
    for i = 2:P
        n = i-1;
        a1 = 4*(n+1)^2*(n+3);
        a3 = 4*(n+1)*(2*n+3)*(n+2);
        a4 = 4*(n+1)^2*(n+2);
        yVec(i+1) = (a3*x*yVec(i) - a4*yVec(i-1) )/a1;
    end
    y = yVec(P+1);
end