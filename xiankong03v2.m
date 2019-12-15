clear all;
A = [1 0;0 2];
B = [1;1];
C = [1 2];

% 先判断系统是否能控
M = ctrb(A,B);
rank_M = rank(M);
new_M = [A B;C 0];
rank_new_M = rank(new_M);

if rank_M == 2 && rank_new_M == 3
    disp('增广系统可控，可以通过状态反馈任意配置极点')
else
    disp('增广系统不可控，不能通过状态反馈任意配置极点，程序终止')
    return;
end


P_MAT = [
    -1,-1,-5;
    -1+1j,-1-1j,-5;
    -5,-5,-25;
    -10+10j,-10-10j,-50;
];
s_value = 10;
s_max = s_value;
s_min = -s_value;
len_p_mat = length(P_MAT);
sub_plot_rows = ceil(len_p_mat / 2);

AA = [A zeros(2,1);C 0];
BB = [B;0];

for p_idx=1:1:len_p_mat
    P = P_MAT(p_idx,:);
    K = acker(AA,BB,P);
    K1 = K(1:2);
    K2 = K(3);
    % disp(K1);
    % disp(K2);
    [t,x,y] = sim('xiankong03_nolinearv3.mdl', [0:0.1:50]);

    
    %% 绘制图像
    if imag(P) == [0 0 0]
        sub_title = sprintf('\\lambda_1=%.02f,\\lambda_2=%.02f，\\lambda_3=%.02f',P);
    else
        sub_title = sprintf('\\lambda_1=%.02f+%.02fj,',real(P(1)),imag(P(1)));
        sub_title = strcat(sub_title, sprintf('\\lambda_2=%.02f+%.02fj',real(P(2)),imag(P(2))));
        sub_title = strcat(sub_title, sprintf('\\lambda_3=%.02f+%.02fj',real(P(3)),imag(P(3))));
    end
    
    subplot(sub_plot_rows,2,p_idx);
    plot(t,simout3);
    hold on;
    % plot(t,simout2);
    plot(t,simout1,'--');
    % plot(t,simout0);
    title(sub_title);
    xlabel('t');
    ylabel('y');
    legend('误差信号','干扰信号','Location','SouthEast');
end
s_title = sprintf('干扰信号为脉冲信号加白噪声，输入为阶跃信号');
f=figure(1);
set(f,'name',s_title);
suptitle(s_title);

