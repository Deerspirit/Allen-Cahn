% Crank-Nicolson方法解具有Dirichlet边值条件的Allen-Cahn方程
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

% 构建Crank-Nicolson方法的系数矩阵
W = (gamma*dt/2) * D2;
A = eye(N+1) - W;
B = eye(N+1) + W + (gamma*dt/epsilon^2)*eye(N+1);

% 时间推进
for n = 1:nt
    % 计算非线性项 f(u) = (u^2-1)u + λu^5
    f = (u.^2 - 1).*u + lambda*u.^5;
    
    % 构建右侧向量
    b = B*u - (gamma*dt/epsilon^2)*f;
    
    % 应用边界条件
    b(1) = bc_left;
    b(N+1) = bc_right;
    
    % 求解线性系统
    u = A\b;
    
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
    plot(x, results(:,i), 'LineWidth', 1.5);
    title(['t = ', num2str(plot_times(i))]);
    xlabel('x');
    ylabel('u');
end

% 绘制时空图
figure(3);
[X, T] = meshgrid(x, plot_times);
Z = results';
surf(X, T, Z);
shading interp;
colormap(jet);
xlabel('x');
ylabel('t');
zlabel('u');
title('Allen-Cahn方程的时空演化');

% Chebyshev微分矩阵生成函数
function [x, D] = cheb(N)
    % 生成Chebyshev点和微分矩阵
    if N == 0
        x = 1;
        D = 0;
        return
    end
    
    % Chebyshev点
    x = cos(pi*(0:N)/N)';
    
    % 初始化微分矩阵
    D = zeros(N+1, N+1);
    
    % 计算非对角元素
    c = [2; ones(N-1, 1); 2] .* (-1).^(0:N)';
    for i = 0:N
        for j = 0:N
            if i ~= j
                D(i+1, j+1) = c(i+1) / (c(j+1) * (x(i+1) - x(j+1)));
            end
        end
    end
    
    % 计算对角元素
    D = D - diag(sum(D, 2));
end