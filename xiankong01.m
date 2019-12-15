clear all
A = [1 0;0 2];
B = [1;1];
C = [1 2];
% 先判断是否可以任意配置极点
M = ctrb(A,B);
rank_M = rank(M);
if rank_M == 2
    disp('系统可控，可以通过状态反馈任意配置极点')
else
    disp('系统不可控，不能通过状态反馈任意配置极点，程序终止')
    return;
end
C = -C;

P_MAT = [
    -1,-1;
    -1+1j,-1-1j;
    -5,-5;
    -10+10j,-10-10j;
    % -30,-30;
    % -30+30j,-30-30j;
    ];
s_value = 10;
s_max = s_value;
s_min = -s_value;
len_p_mat = length(P_MAT);
sub_plot_rows = ceil(len_p_mat / 2);
% 调整单位阶跃信号的K
u_k = 10;
for p_idx=1:1:len_p_mat
    P = P_MAT(p_idx,:);
    K = acker(A,B,P);
    % step(A-B*K, B, C, 0)
    % [a,b] = ss2tf(A-B*K, B, C, 0)
    % step(sys)
    [t,x,y] = sim('xiankong01_nolinear.mdl', [0:0.01:10]);

    %% 求取ts tp tr
    y_steady_value = dcgain(A-B*K, B, C, 0);
    % 求取tp
    [ymax,tp] = max(simout2);
    % disp(ymax);
    max_iter_step = length(simout2);
    if ymax > y_steady_value+1e-6
        tp = t(tp);
    else
        tp = inf;
    end

    % 求取 tr 
    % tr = 0.0;
    % 对应系统没有超调的情况
    if ymax < y_steady_value+1e-6
        left_value = 0.1 * y_steady_value;
        idx = 1;
        t1 = 0.0;
        while idx <= max_iter_step && simout2(idx) < left_value
            idx = idx + 1;
        end
        if idx < max_iter_step
            t1 = t(idx);
        end
        right_value = 0.9 * y_steady_value;
        t2 = 0.0;
        while idx <= max_iter_step && simout2(idx) < right_value
            idx = idx + 1;
        end
        if idx < max_iter_step
            t2 = t(idx);
        end
        tr = t2 - t1;
        % disp(tr);
    else
    % 对应系统有超调的情况
        idx = 1;
        while simout2(idx) < y_steady_value
            idx = idx + 1;
        end 
        tr = t(idx);
    end  

    % 求解 ts
    idx = max_iter_step;
    while simout2(idx) > 0.98 * y_steady_value && simout2(idx) < 1.02 * y_steady_value
        idx = idx - 1;
    end
    % disp(idx+1);
    ts = inf;
    if idx + 1 <= max_iter_step
        ts = t(idx+1);
    end
    
    % ts,tp,tr

    %% 绘制图像
    if imag(P) == [0 0]
        sub_title = sprintf('\\lambda_1=%.02f,\\lambda_2=%.02f',P);
    else
        sub_title = sprintf('\\lambda_1=%.02f+%.02fj,',real(P(1)),imag(P(1)));
        sub_title = strcat(sub_title, sprintf('\\lambda_2=%.02f+%.02fj',real(P(2)),imag(P(2))));
    end
    sub_title = strcat(sub_title, sprintf(',steady value=%.02f',y_steady_value));
    
    subplot(sub_plot_rows,2,p_idx);
    plot(t,simout2);
    title(sub_title);
    xlabel('t');
    ylabel('y');
    legend(sprintf('t_s=%.02f\nt_r=%.02f\nt_p=%.02f',ts,tr,tp),'Location','SouthEast');
end
s_title = sprintf('单位阶跃K=%d,非线性u_{max}=%.02f', u_k, s_value);
f=figure(1);
set(f,'name',s_title);
suptitle(s_title);
