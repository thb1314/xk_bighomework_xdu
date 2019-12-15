clear all;
A = [1 0;0 2];
B = [1;1];
C = [1 2];

u_k = 1;
AA = [A zeros(2,1);C 0];
BB = [B;0];

[G0_z,G0_s] = ss2tf(A,B,C,0);

P_MAT = [
    -1,-1,-5;
    -1+1j,-1-1j,-5;
    -5,-5,-25;
    -10+10j,-10-10j,-50;
];

s_value = inf;
s_max = s_value;
s_min = -s_value;
len_p_mat = length(P_MAT);
sub_plot_rows = ceil(len_p_mat / 2);

AA = [A zeros(2,1);C 0];
BB = [B;0];
syms s;

for p_idx=1:1:len_p_mat
    
    P = P_MAT(p_idx,:);
    K = acker(AA,BB,P);
    % 系统反馈矩阵
    K1 = K(1:2);
    K2 = K(3);

    Q = 3 * P(1:2);
    % 状态观测器反馈矩阵
    L = (acker(A',C',Q))';

    G0 = ([s^2 s^1 s^0] * G0_z') / ([s^2 s^1 s^0] * G0_s');
    G1 = simplify( K1 * inv(s*eye(2) - (A - L * C)) * B );
    [n,d] = numden(G1);
    G1_z = sym2poly(n);
    G1_s = sym2poly(d);
    G2 = simplify( K1 * inv(s*eye(2) - (A - L * C)) * L );
    [n,d] = numden(G2);
    G2_z = sym2poly(n);
    G2_s = sym2poly(d);

    G3 = simplify( inv(1+G1) );
    [n,d] = numden(G3);
    G3_z = sym2poly(n);
    G3_s = sym2poly(d);
    [t,x,y] = sim('xiankong05_nolinear.mdl', [0:0.01:50]);

    % 求取系统开环传递函数
    G = (1/s) * K2 * (G3 * G0) / (1 + G3 * G0 * G2);
    G = simplify(G);
    % pretty(G);
    %% 绘制图像
    if imag(P) == [0 0 0]
        sub_title = sprintf('\\lambda_1=%.02f,\\lambda_2=%.02f，\\lambda_3=%.02f',P);
    else
        sub_title = sprintf('\\lambda_1=%.02f+%.02fj,',real(P(1)),imag(P(1)));
        sub_title = strcat(sub_title, sprintf('\\lambda_2=%.02f+%.02fj',real(P(2)),imag(P(2))));
        sub_title = strcat(sub_title, sprintf('\\lambda_3=%.02f+%.02fj',real(P(3)),imag(P(3))));
    end

    subplot(sub_plot_rows,2,p_idx);
    plot(t,simout0);
    hold on;
    plot(t,simout3);
    % plot(t,simout2);
    title(sub_title);
    xlabel('t');
    ylabel('y');
    legend('输入信号','误差信号','Location','SouthEast');
end
% s_title = sprintf('状态反馈模型,输入为阶跃信号,饱和非线性系数为%.02f',s_value);
s_title = sprintf('串并联等价模型,输入为阶跃信号,饱和非线性系数为%.02f',s_value);
f=figure(1);
set(f,'name',s_title);
suptitle(s_title);

