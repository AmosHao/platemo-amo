# 订单配送问题数据与可视化

## 文件说明

### 1. `order_data.m` - 数据文件
包含所有问题数据：
- **dotss**: 坐标矩阵（配送中心、商家1-10、客户点11-20）
- **forbidden_zones**: 矩形禁飞区 [x_min, y_min, x_max, y_max]
- **crowded_zones**: 圆形人流密集区 [cx, cy, radius]
- **noise_zones**: 矩形噪音敏感区 [x_min, y_min, x_max, y_max]
- **penalty_forbidden**: 禁飞区惩罚系数
- **penalty_crowded**: 人流密集区惩罚系数
- **penalty_noise**: 噪音敏感区惩罚系数

### 2. `plot_order_map.m` - 地图绘制脚本
绘制包含所有元素的地图：
- 配送中心（紫色大圆点）
- 商家（蓝色圆点，编号1-10）
- 客户点（绿色圆点，编号11-20）
- 禁飞区（红色矩形，半透明）
- 人流密集区（淡绿色圆形，半透明）
- 噪音敏感区（橙色矩形，半透明）

## 使用方法

### 修改数据
直接编辑 `order_data.m` 文件，修改相应的数据即可。

### 绘制地图
在MATLAB中运行：
```matlab
cd forOrderNew26
plot_order_map
```

### 在问题代码中使用
`planner_order_v3.m` 会自动从 `order_data.m` 引入数据。如果找不到数据文件，会使用默认数据并显示警告。

## 数据格式说明

### 坐标矩阵 (dotss)
- 第1行：配送中心
- 第2-11行：商家1-10
- 第12-21行：客户点11-20
- 格式：[x, y, z]，单位：米

### 禁飞区 (forbidden_zones)
- 格式：[x_min, y_min, x_max, y_max]
- 矩形区域，红色显示

### 人流密集区 (crowded_zones)
- 格式：[cx, cy, radius]
- 圆形区域，淡绿色显示

### 噪音敏感区 (noise_zones)
- 格式：[x_min, y_min, x_max, y_max]
- 矩形区域，橙色显示

