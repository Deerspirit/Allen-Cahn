% ETDRK4方法解具有周期边值条件的Allen-Cahn方程
% 添加非线性项g(u)=λu^5

% 参数设置
gamma = 0.01;     % 扩散系数
epsilon = 0.1;     % 界面宽度参数
lambda = 0.01;    % 非线性项系数
N = 128;          % 空间离散点数
dt = 0.1;         % 时间步长
tmax = 50;        % 最大计算时间
nt = tmax/dt;     % 时间步数

% 空间离散
L = 2*pi;         % 计算区域长度
dx = L/N;         % 空间步长
x = (0:N-1)'*dx; % 空间网格点

% 波数
k = [0:N/2-1 0 -N/2+1:-1]'; % 波数向量

% 初始条件（随机分块）
u = zeros(N,1);
nb = 8;  % 分块数
for j = 0:N/nb-1
    u(j*nb+1:(j+1)*nb) = (rand(1)-0.5)*ones(nb,1);
end

% 存储结果用于绘图
plot_times = [0 1 5 10 20 50]; % 需要绘图的时间点
plot_indices = round(plot_times/dt) + 1;
plot_indices = plot_indices(plot_indices <= nt+1);
results = zeros(N, length(plot_indices));
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
L = -gamma * k.^2;

% 计算ETDRK4所需的系数
E = exp(dt*L);
E2 = exp(dt*L/2);

% 使用contour积分计算ETDRK4系数
M = 32; % 积分点数
r = exp(1i*pi*((1:M)-0.5)/M); % 积分路径上的点

% 计算ETDRK4系数
Lr = dt*L(:) + r; % 扩展L到积分路径上的点
Q = dt*mean((exp(Lr/2)-1)./Lr, 2);
f1 = dt*mean((-4-Lr+exp(Lr).*(4-3*Lr+Lr.^2))./Lr.^3, 2);
f2 = dt*mean((2+Lr+exp(Lr).*(-2+Lr))./Lr.^3, 2);
f3 = dt*mean((-4-3*Lr-Lr.^2+exp(Lr).*(4-Lr))./Lr.^3, 2);

% 时间推进
for n = 1:nt
    % 计算非线性项 N(u) = -γ/ε² * (u³-u+λu⁵)
    Nuk = -gamma/epsilon^2 * fft((u.^3-u) + lambda*u.^5);
    
    % ETDRK4方法的四个阶段
    a = E2.*fft(u) + Q.*Nuk;
    
    % 计算第二阶段的非线性项
    Na = -gamma/epsilon^2 * fft((real(ifft(a)).^3-real(ifft(a))) + lambda*real(ifft(a)).^5);
    
    % 计算第三阶段
    b = E2.*fft(u) + Q.*Na;
    
    % 计算第三阶段的非线性项
    Nb = -gamma/epsilon^2 * fft((real(ifft(b)).^3-real(ifft(b))) + lambda*real(ifft(b)).^5);
    
    % 计算第四阶段
    c = E2.*a + Q.*(2*Nb-Nuk);
    
    % 计算第四阶段的非线性项
    Nc = -gamma/epsilon^2 * fft((real(ifft(c)).^3-real(ifft(c))) + lambda*real(ifft(c)).^5);
    
    % 最终更新
    fft_u = E.*fft(u) + f1.*Nuk + f2.*(Na+Nb) + f3.*Nc;
    u = real(ifft(fft_u));
    
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
    subplot(2, ceil(plot_count/2), i);
    plot(x, results(:,i), 'LineWidth', 1.5);
    title(['t = ', num2str(plot_times(i))]);
    xlabel('x');
    ylabel('u');
end

% 计算误差和收敛阶（如果需要）
function [error, order] = compute_error(dt_values)
    % 参数设置
    gamma = 0.01;
    epsilon = 0.1;
    lambda = 0.01;
    N = 128;
    tmax = 50;
    
    % 计算参考解（使用非常小的时间步长）
    dt_ref = 1/640;
    nt_ref = tmax/dt_ref;
    
    % 初始条件
    u_init = zeros(N,1);
    nb = 8;
    for j = 0:N/nb-1
        u_init(j*nb+1:(j+1)*nb) = (rand(1)-0.5)*ones(nb,1);
    end
    
    % 计算参考解
    u_ref = solve_allen_cahn_etdrk4(u_init, gamma, epsilon, lambda, N, dt_ref, nt_ref);
    
    % 计算不同时间步长的误差
    errors = zeros(size(dt_values));
    for i = 1:length(dt_values)
        dt = dt_values(i);
        nt = tmax/dt;
        u_approx = solve_allen_cahn_etdrk4(u_init, gamma, epsilon, lambda, N, dt, nt);
        errors(i) = max(abs(u_approx - u_ref));
    end
    
    % 计算收敛阶
    order = zeros(length(dt_values)-1, 1);
    for i = 1:length(dt_values)-1
        order(i) = log2(errors(i)/errors(i+1));
    end
    
    error = errors;
end

% 求解Allen-Cahn方程的函数（使用ETDRK4方法）
function u = solve_allen_cahn_etdrk4(u_init, gamma, epsilon, lambda, N, dt, nt)
    % 波数
    k = [0:N/2-1 0 -N/2+1:-1]';
    
    % 线性算子
    L = -gamma * k.^2;
    
    % 计算ETDRK4所需的系数
    E = exp(dt*L);
    E2 = exp(dt*L/2);
    
    % 使用contour积分计算ETDRK4系数
    M = 32; % 积分点数
    r = exp(1i*pi*((1:M)-0.5)/M); % 积分路径上的点
    
    % 计算ETDRK4系数
    Lr = dt*L(:) + r; % 扩展L到积分路径上的点
    Q = dt*mean((exp(Lr/2)-1)./Lr, 2);
    f1 = dt*mean((-4-Lr+exp(Lr).*(4-3*Lr+Lr.^2))./Lr.^3, 2);
    f2 = dt*mean((2+Lr+exp(Lr).*(-2+Lr))./Lr.^3, 2);
    f3 = dt*mean((-4-3*Lr-Lr.^2+exp(Lr).*(4-Lr))./Lr.^3, 2);
    
    % 初始条件
    u = u_init;
    
    % 时间推进
    for n = 1:nt
        % 计算非线性项 N(u) = -γ/ε² * (u³-u+λu⁵)
        Nuk = -gamma/epsilon^2 * fft((u.^3-u) + lambda*u.^5);
        
        % ETDRK4方法的四个阶段
        a = E2.*fft(u) + Q.*Nuk;
        
        % 计算第二阶段的非线性项
        Na = -gamma/epsilon^2 * fft((real(ifft(a)).^3-real(ifft(a))) + lambda*real(ifft(a)).^5);
        
        % 计算第三阶段
        b = E2.*fft(u) + Q.*Na;
        
        % 计算第三阶段的非线性项
        Nb = -gamma/epsilon^2 * fft((real(ifft(b)).^3-real(ifft(b))) + lambda*real(ifft(b)).^5);
        
        % 计算第四阶段
        c = E2.*a + Q.*(2*Nb-Nuk);
        
        % 计算第四阶段的非线性项
        Nc = -gamma/epsilon^2 * fft((real(ifft(c)).^3-real(ifft(c))) + lambda*real(ifft(c)).^5);
        
        % 最终更新
        fft_u = E.*fft(u) + f1.*Nuk + f2.*(Na+Nb) + f3.*Nc;
        u = real(ifft(fft_u));
    end
end