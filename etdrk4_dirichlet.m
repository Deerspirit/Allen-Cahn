% ETDRK4方法解具有Dirichlet边值条件的Allen-Cahn方程
% 添加非线性项g(u)=λu^5

% 参数设置
gamma = 0.1;      % 扩散系数
epsilon = 0.1;     % 界面宽度参数
lambda = 0.01;    % 非线性项系数
N = 280;          % 空间离散点数
dt = 1/10000;     % 时间步长
tmax = 0.01;      % 最大计算时间
nt = tmax/dt;     % 时间步数

% 空间离散
a = -1;           % 计算区域左端点
b = 1;            % 计算区域右端点

% 生成Chebyshev点和微分矩阵
[x, D] = cheb(N);
D2 = D^2;         % 二阶微分矩阵

% 初始条件
u = 0.53*x + 0.47*sin(-1.5*pi*x);

% 边界条件处理
% 设置边界条件：u(-1,t) = -1, u(1,t) = 1
bc_left = -1;
bc_right = 1;

% 修改微分矩阵以应用边界条件
D2(1,:) = 0; D2(1,1) = 1;     % 左边界条件
D2(N+1,:) = 0; D2(N+1,N+1) = 1; % 右边界条件

% 存储结果用于绘图
plot_times = [0 0.001 0.01]; % 需要绘图的时间点
plot_indices = round(plot_times/dt) + 1;
plot_indices = plot_indices(plot_indices <= nt+1);
results = zeros(N+1, length(plot_indices));
result_index = 1;

if plot_indices(1) == 1
    results(:, result_index) = u;
    result_index = result_index + 1;
end

% 绘制初始条件
figure(1);
plot(x, u, 'LineWidth', 1.5);
title(['t = 0']);
xlabel('x');
ylabel('u');
drawnow;

% 线性算子
L = gamma * D2;

% 应用边界条件到线性算子
L(1,:) = 0; L(1,1) = 1;     % 左边界条件
L(N+1,:) = 0; L(N+1,N+1) = 1; % 右边界条件

% 计算ETDRK4所需的系数
E = exp(dt*L);
E2 = exp(dt*L/2);

% 使用contour积分计算ETDRK4系数
M = 32; % 积分点数
r = exp(1i*pi*((1:M)-0.5)/M); % 积分路径上的点

% 计算ETDRK4系数
Lr = dt*diag(L) + r; % 只使用L的对角线元素，而不是整个矩阵
Q = dt*mean((exp(Lr/2)-1)./Lr, 2);
f1 = dt*mean((-4-Lr+exp(Lr).*(4-3*Lr+Lr.^2))./Lr.^3, 2);
f2 = dt*mean((2+Lr+exp(Lr).*(-2+Lr))./Lr.^3, 2);
f3 = dt*mean((-4-3*Lr-Lr.^2+exp(Lr).*(4-Lr))./Lr.^3, 2);

% 确保所有系数都是实数
Q = real(Q);
f1 = real(f1);
f2 = real(f2);
f3 = real(f3);

% 改为直接计算对角线元素
E = real(exp(dt*diag(L)));
E2 = real(exp(dt*diag(L)/2));

% 时间推进
for n = 1:nt
    % 计算非线性项 N(u) = -γ/ε² * (u³-u+λu⁵)
    Nu = -gamma/epsilon^2 * ((u.^3-u) + lambda*u.^5);
    
    % 应用边界条件到非线性项
    Nu(1) = 0;      % 左边界条件
    Nu(N+1) = 0;    % 右边界条件
    
    % ETDRK4方法的四个阶段 - 使用元素级乘法
    a = E2.*u + Q.*Nu;
    
    % 应用边界条件
    a(1) = bc_left;
    a(N+1) = bc_right;
    
    % 计算第二阶段的非线性项
    Na = -gamma/epsilon^2 * ((a.^3-a) + lambda*a.^5);
    Na(1) = 0;      % 左边界条件
    Na(N+1) = 0;    % 右边界条件
    
    % 计算第三阶段
    b = E2.*u + Q.*Na;
    b(1) = bc_left;
    b(N+1) = bc_right;
    
    % 计算第三阶段的非线性项
    Nb = -gamma/epsilon^2 * ((b.^3-b) + lambda*b.^5);
    Nb(1) = 0;      % 左边界条件
    Nb(N+1) = 0;    % 右边界条件
    
    % 计算第四阶段
    c = E2.*a + Q.*(2*Nb-Nu);
    c(1) = bc_left;
    c(N+1) = bc_right;
    
    % 计算第四阶段的非线性项
    Nc = -gamma/epsilon^2 * ((c.^3-c) + lambda*c.^5);
    Nc(1) = 0;      % 左边界条件
    Nc(N+1) = 0;    % 右边界条件
    
    % 最终更新
    u = E.*u + f1.*Nu + f2.*(Na+Nb) + f3.*Nc;
    
    % 应用边界条件
    u(1) = bc_left;
    u(N+1) = bc_right;
    
    % 存储特定时间点的结果
    if ismember(n+1, plot_indices)
        results(:, result_index) = u;
        result_index = result_index + 1;
    end
    
    % 绘制特定时间点的结果
    if ismember(n, round(plot_times/dt))
        figure(1);
        plot(x, u, 'LineWidth', 1.5);
        title(['t = ', num2str(n*dt)]);
        xlabel('x');
        ylabel('u');
        drawnow;
    end
end

% 绘制最终结果
figure(2);
plot_count = length(plot_indices);
for i = 1:plot_count
    subplot(1, plot_count, i);
    plot(x, real(results(:,i)), 'LineWidth', 1.5);
    title(['t = ', num2str(plot_times(i))]);
    xlabel('x');
    ylabel('u');
    axis([-1 1 -1.2 1.2]);
end

% 绘制时空图
figure(3);
[X, T] = meshgrid(x, plot_times);
imagesc(x, plot_times, real(results'));
colorbar;
title('Allen-Cahn方程的时空演化');
xlabel('x');
ylabel('t');

% Chebyshev微分矩阵函数
function [x, D] = cheb(N)
    % 计算Chebyshev点
    x = cos(pi*(0:N)/N)';
    
    % 初始化微分矩阵
    D = zeros(N+1, N+1);
    
    % 计算对角线外的元素
    c = [2; ones(N-1, 1); 2] .* (-1).^(0:N)';
    for i = 0:N
        for j = 0:N
            if i ~= j
                D(i+1, j+1) = c(i+1) / (c(j+1) * (x(i+1) - x(j+1)));
            end
        end
    end
    
    % 计算对角线元素
    D(1, 1) = (2*N^2 + 1) / 6;
    D(N+1, N+1) = -D(1, 1);
    for i = 2:N
        D(i, i) = -x(i) / (2 * (1 - x(i)^2));
    end
end