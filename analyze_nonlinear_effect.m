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
    
    % 空间离散
    L = 2*pi;         % 计算区域长度
    x = (0:N-1)*L/N; % x方向网格点
    y = (0:N-1)*L/N; % y方向网格点
    [X, Y] = meshgrid(x, y);
    
    % 波数
    k = [0:N/2-1 0 -N/2+1:-1]; % 波数向量
    [KX, KY] = meshgrid(k, k);  % 二维波数网格
    
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
    
    % 存储最终结果
    final_results = zeros(N, N, length(lambda_values));
    
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
        
        % 存储最终结果
        final_results(:,:,lambda_idx) = u;
        
        % 绘制结果
        subplot(2, 2, lambda_idx);
        imagesc(x, y, u);
        axis square;
        colorbar;
        title(['\lambda = ', num2str(lambda)]);
        xlabel('x');
        ylabel('y');
    end
    
    % 计算和绘制不同lambda值下的能量
    figure(5);
    energies = zeros(length(lambda_values), 1);
    % 重新定义波数，确保维度正确
    k = [0:N/2-1 0 -N/2+1:-1]; % 波数向量
    [KX, KY] = meshgrid(k, k);  % 二维波数网格
    
    % 定义空间离散步长
    L = 2*pi;         % 计算区域长度
    dx = L/N;
    
    for i = 1:length(lambda_values)
        u = final_results(:,:,i);
        lambda = lambda_values(i);
        % 计算能量泛函 E(u) = ∫(γ/2|∇u|² + 1/(4ε²)(u²-1)² + λ/(6ε²)u⁶)dx
        % 这里使用简化的离散计算
        u_x = fft2(u) .* (1i*KX);
        u_y = fft2(u) .* (1i*KY);
        grad_u_squared = abs(ifft2(u_x)).^2 + abs(ifft2(u_y)).^2;
        potential = (1/(4*epsilon^2))*(u.^2-1).^2 + (lambda/(6*epsilon^2))*u.^6;
        energy = sum(sum(gamma/2*grad_u_squared + potential)) * (dx)^2;
        energies(i) = energy;
    end
    
    % 绘制能量随lambda变化的曲线
    plot(lambda_values, energies, 'o-', 'LineWidth', 2);
    title('能量随非线性项系数\lambda的变化');
    xlabel('\lambda');
    ylabel('能量');
    grid on;
    
    % 绘制不同lambda值下的相图
    figure(6);
    for i = 1:length(lambda_values)
        subplot(2, 2, i);
        u = final_results(:,:,i);
        contourf(X, Y, u, 20);
        axis square;
        colorbar;
        title(['\lambda = ', num2str(lambda_values(i)), ' 的相图']);
        xlabel('x');
        ylabel('y');
    end
end