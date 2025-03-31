% 测试所有六个Allen-Cahn方程求解文件是否能正常运行

disp('开始测试所有Allen-Cahn方程求解文件...');

% 测试etdrk4_2d.m
try
    disp('测试etdrk4_2d.m...');
    run('etdrk4_2d.m');
    disp('etdrk4_2d.m 运行成功!');
    
    % 测试analyze_nonlinear_effect函数
    disp('测试analyze_nonlinear_effect.m...');
    analyze_nonlinear_effect();
    disp('analyze_nonlinear_effect.m 运行成功!');
catch e
    disp(['etdrk4_2d.m 运行失败: ' e.message]);
end

% 测试etdrk4_periodic.m
try
    disp('测试etdrk4_periodic.m...');
    run('etdrk4_periodic.m');
    disp('etdrk4_periodic.m 运行成功!');
catch e
    disp(['etdrk4_periodic.m 运行失败: ' e.message]);
end

% 测试etdrk4_dirichlet.m
try
    disp('测试etdrk4_dirichlet.m...');
    run('etdrk4_dirichlet.m');
    disp('etdrk4_dirichlet.m 运行成功!');
catch e
    disp(['etdrk4_dirichlet.m 运行失败: ' e.message]);
end

% 测试semi_implicit_periodic.m
try
    disp('测试semi_implicit_periodic.m...');
    run('semi_implicit_periodic.m');
    disp('semi_implicit_periodic.m 运行成功!');
catch e
    disp(['semi_implicit_periodic.m 运行失败: ' e.message]);
end

% 测试semi_implicit_neumann.m
try
    disp('测试semi_implicit_neumann.m...');
    run('semi_implicit_neumann.m');
    disp('semi_implicit_neumann.m 运行成功!');
catch e
    disp(['semi_implicit_neumann.m 运行失败: ' e.message]);
end

% 测试crank_nicolson_dirichlet.m
try
    disp('测试crank_nicolson_dirichlet.m...');
    run('crank_nicolson_dirichlet.m');
    disp('crank_nicolson_dirichlet.m 运行成功!');
catch e
    disp(['crank_nicolson_dirichlet.m 运行失败: ' e.message]);
end

disp('所有测试完成!');