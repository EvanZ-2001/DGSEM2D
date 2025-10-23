%% 确定谱元法的Gauss-Lobatto-Legendre积分点坐标和积分权系数
% 积分点为alpha=beta=1的Q-2阶Jacobi多项式零点再加上-1和1
% 积分点同时也是数值微分的插值点
% Q为积分点的个数

function [xi,w] = GLLNodesAndWeights(Q)
    xi = zeros(Q,1);
    x = zeros(Q-2,1);  
    err = 1*10^(-8);%容许的误差
    xi(1) = -1;
    xi(Q) = 1;
    P = Q-2;%Jacobi多项式阶数
    for k = 1:P
        r = -cos(((2*k-1)*pi)/(2*P)); %Chebychev 多项式零点
        if k>1
            r = (r+x(k-1))/2;
        end
        j = 0;
        while(1)
            s = 0 ;
            if k > 1
                for i=1:k-1
                    s = s + 1/(r-x(i));
                end
            end
            delta = - JacobiPolynomial(r,P)/(DxJacobiPolynomial(r,P)-JacobiPolynomial(r,P)*s);
            r = r + delta;
            if delta < err %判断迭代变化值是否小于容许误差
                break;
            end
            j = j+1;
            if(j>50) 
                disp("无法计算出零点");
                disp(k);
                break;
            end
        end
        if abs(r)<err %校正x=0的积分点
            r=0;
        end
        x(k)=r;
    end
    xi(2:Q-1) = x;

    w = 2./(Q*(Q-1).*(legendreP(Q-1,xi)).^2);%积分权系数
end