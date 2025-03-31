% ETDRK4方法解二维的Allen-Cahn方程
% 添加非线性项g(u)=λu^5

% 参数设置
gamma = 0.01;     % 扩散系数
epsilon = 0.1;     % 界面宽度参数
lambda = 0.01;    % 非线性项系数
N = 128;          % 空间离散点数
dt = 0.1;         % 时间步长
tmax = 20;        % 最大计算时间
nt = tmax/dt;     % 时间步数

% 空间离散
L = 2*pi;         % 计算区域长度
x = (0:N-1)*L/N; % x方向网格点
y = (0:N-1)*L/N; % y方向网格点
[X, Y] = meshgrid(x, y);

% 波数
k = [0:N/2-1 0 -N/2+1:-1]; % 波数向量
[KX, KY] = meshgrid(k, k);  % 二维波数网格

% 初始条件（随机分块）
u = zeros(N, N);
nb = 8;  % 分块数
for j = 0:N/nb-1
    for j2 = 0:N/nb-1
        u(j*nb+1:(j+1)*nb, j2*nb+1:(j2+1)*nb) = 0.5*(rand(1)-0.5)*ones(nb, nb);
    end
end

% 绘制初始条件
figure(1);
subplot(2, 2, 1);
imagesc(x, y, u);
axis square;
colorbar;
title('初始条件 t = 0');
xlabel('x');
ylabel('y');
drawnow;

% 线性算子
L = -gamma * (KX.^2 + KY.^2);

% 计算ETDRK4所需的系数
E = exp(dt*L);
E2 = exp(dt*L/2);

% 使用contour积分计算ETDRK4系数
M = 16; % 积分点数
r = exp(1i*pi*((1:M)-0.5)/M); % 积分路径上的点

% 计算ETDRK4系数
Lr = dt*L(:) + r; % 扩展L到积分路径上的点
Q = dt*mean((exp(Lr/2)-1)./Lr, 2);
f1 = dt*mean((-4-Lr+exp(Lr).*(4-3*Lr+Lr.^2))./Lr.^3, 2);
f2 = dt*mean((2+Lr+exp(Lr).*(-2+Lr))./Lr.^3, 2);
f3 = dt*mean((-4-3*Lr-Lr.^2+exp(Lr).*(4-Lr))./Lr.^3, 2);

% 重塑系数为二维
Q = reshape(Q, size(L));
f1 = reshape(f1, size(L));
f2 = reshape(f2, size(L));
f3 = reshape(f3, size(L));

% 时间推进
plot_times = [0, 5, 10, 20]; % 需要绘图的时间点
plot_indices = round(plot_times/dt) + 1;
plot_indices = plot_indices(plot_indices <= nt+1);

for n = 1:nt
    % 计算非线性项 N(u) = -γ/ε² * (u³-u+λu⁵)
    Nu = -gamma/epsilon^2 * ((u.^3-u) + lambda*u.^5);
    Nuk = fft2(Nu);
    
    % ETDRK4方法的四个阶段
    a = E2.*fft2(u) + Q.*Nuk;
    
    % 计算第二阶段的非线性项
    ua = real(ifft2(a));
    Na = -gamma/epsilon^2 * ((ua.^3-ua) + lambda*ua.^5);
    
    % 计算第三阶段
    b = E2.*fft2(u) + Q.*fft2(Na);
    
    % 计算第三阶段的非线性项
    ub = real(ifft2(b));
    Nb = -gamma/epsilon^2 * ((ub.^3-ub) + lambda*ub.^5);
    
    % 计算第四阶段
    c = E2.*a + Q.*(2*fft2(Nb)-Nuk);
    
    % 计算第四阶段的非线性项
    uc = real(ifft2(c));
    Nc = -gamma/epsilon^2 * ((uc.^3-uc) + lambda*uc.^5);
    
    % 最终更新
    u_hat = E.*fft2(u) + f1.*Nuk + f2.*(fft2(Na)+fft2(Nb)) + f3.*fft2(Nc);
    u = real(ifft2(u_hat));
    
    % 绘制特定时间点的结果
    if ismember(n, round(plot_times(2:end)/dt))
        figure(1);
        subplot(2, 2, find(round(plot_times/dt) == n));
        imagesc(x, y, u);
        axis square;
        colorbar;
        title(['t = ', num2str(n*dt)]);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
end

% 绘制最终结果的等高线图
figure(2);
contourf(X, Y, u, 20);
axis square;
colorbar;
title(['Allen-Cahn方程在t = ', num2str(tmax), '时的等高线图']);
xlabel('x');
ylabel('y');

% 绘制最终结果的3D表面图
figure(3);
surf(X, Y, u);
shading interp;
colormap(jet);
title(['Allen-Cahn方程在t = ', num2str(tmax), '时的3D表面图']);
xlabel('x');
ylabel('y');
zlabel('u');
axis tight;

% 分析不同非线性项系数的影响
function analyze_nonlinear_effect()
    % 参数设置
    gamma = 0.01;
    epsilon = 0.1;
    N = 128;
    dt = 0.1;
    tmax = 20;
    nt = tmax/dt;
    
    % 测试不同的非线性项系数
    lambda_values = [0, 0.01, 0.05, 0.1];
    
    % 初始条件
    u_init = zeros(N, N);
    nb = 8;
    for j = 0:N/nb-1
        for j2 = 0:N/nb-1
            u_init(j*nb+1:(j+1)*nb, j2*nb+1:(j2+1)*nb) = 0.5*(rand(1)-0.5)*ones(nb, nb);
        end
    end
    
    % 波数
    k = [0:N/2-1 0 -N/2+1:-1];
    [KX, KY] = meshgrid(k, k);
    
    % 线性算子
    L = -gamma * (KX.^2 + KY.^2);
    
    % 计算ETDRK4所需的系数
    E = exp(dt*L);
    E2 = exp(dt*L/2);
    
    % 使用contour积分计算ETDRK4系数
    M = 16;
    r = exp(1i*pi*((1:M)-0.5)/M);
    Lr = dt*L(:) + r;
    Q = dt*mean((exp(Lr/2)-1)./Lr, 2);
    f1 = dt*mean((-4-Lr+exp(Lr).*(4-3*Lr+Lr.^2))./Lr.^3, 2);
    f2 = dt*mean((2+Lr+exp(Lr).*(-2+Lr))./Lr.^3, 2);
    f3 = dt*mean((-4-3*Lr-Lr.^2+exp(Lr).*(4-Lr))./Lr.^3, 2);
    
    % 重塑系数为二维
    Q = reshape(Q, size(L));
    f1 = reshape(f1, size(L));
    f2 = reshape(f2, size(L));
    f3 = reshape(f3, size(L));
    
    % 创建图形窗口
    figure(4);
    
    % 对每个非线性项系数进行计算
    for lambda_idx = 1:length(lambda_values)
        lambda = lambda_values(lambda_idx);
        
        % 复制初始条件
        u = u_init;
        
        % 时间推进
        for n = 1:nt
            % 计算非线性项 N(u) = -γ/ε² * (u³-u+λu⁵)
            Nu = -gamma/epsilon^2 * ((u.^3-u) + lambda*u.^5);
            Nuk = fft2(Nu);
            
            % ETDRK4方法的四个阶段
            a = E2.*fft2(u) + Q.*Nuk;
            
            % 计算第二阶段的非线性项
            ua = real(ifft2(a));
            Na = -gamma/epsilon^2 * ((ua.^3-ua) + lambda*ua.^5);
            
            % 计算第三阶段
            b = E2.*fft2(u) + Q.*fft2(Na);
            
            % 计算第三阶段的非线性项
            ub = real(ifft2(b));
            Nb = -gamma/epsilon^2 * ((ub.^3-ub) + lambda*ub.^5);
            
            % 计算第四阶段
            c = E2.*a + Q.*(2*fft2(Nb)-Nuk);
            
            % 计算第四阶段的非线性项
            uc = real(ifft2(c));
            Nc = -gamma/epsilon^2 * ((uc.^3-uc) + lambda*uc.^5);
            
            % 最终更新
            u_hat = E.*fft2(u) + f1.*Nuk + f2.*(fft2(Na)+fft2(Nb)) + f3.*fft2(Nc);
            u = real(ifft2(u_hat));
        end
        
        % 绘制结果
        subplot(2, 2, lambda_idx);
        imagesc(x, y, u);
        axis square;
        colorbar;
        title(['\lambda = ', num2str(lambda)]);
        xlabel('x');
        ylabel('y');
    end
end