# Allen-Cahn方程数值求解

本项目基于论文《谱方法求解Allen-Cahn方程与Cahn-Hilliard方程》，实现了六种不同的数值方法来求解Allen-Cahn方程，并添加了非线性项。

## 项目结构

1. `semi_implicit_periodic.m` - 半隐式方法解具有周期边值条件的Allen-Cahn方程
2. `etdrk4_periodic.m` - ETDRK4方法解具有周期边值条件的Allen-Cahn方程
3. `crank_nicolson_dirichlet.m` - Crank-Nicolson方法解具有Dirichlet边值条件的Allen-Cahn方程
4. `etdrk4_dirichlet.m` - ETDRK4方法解具有Dirichlet边值条件的Allen-Cahn方程
5. `stabilized_semi_implicit.m` - Allen-Cahn方程稳定的一阶半隐式格式
6. `etdrk4_2d.m` - ETDRK4方法解二维的Allen-Cahn方程

## 方程说明

Allen-Cahn方程的标准形式为：

```
∂_t u = γ[Δu - (1/ε²)f(u)]
```

其中，f(u) = (u² - 1)u 是原始的非线性项。在本实现中，我们添加了额外的非线性项 g(u) = λu⁵。

## 运行环境

- MATLAB R2018b或更高版本
- 需要Signal Processing Toolbox用于FFT操作
