# NSGA-II-Order 算法

## 说明

这是专门为`planner_order_v3`问题设计的NSGA-II变体算法。

## 问题

标准的NSGA-II使用`OperatorGA`，其中`GApermutation`函数假设所有元素都是唯一的（标准permutation编码）。但`planner_order_v3`问题中：
- 序列包含重复的0（分隔符）
- 点的编号是1到20，但序列长度是22（包含2个0）

这导致`GApermutation`中的`setdiff`操作失败。

## 解决方案

1. **修改问题编码**：将`planner_order_v3.m`中的`encoding`从5（permutation）改为2（integer）
2. **创建自定义操作符**：使用A-amos算法中的交叉和变异函数：
   - `SegmentBasedCrossover`：基于段的交叉
   - `DiscreteMutation`：离散序列变异
   - `OrderConstraintRepair`：顺序约束修复

## 使用方法

1. 确保`planner_order_v3.m`中的`encoding`已改为2
2. 在PlatEMO界面中选择`NSGA-II-Order`算法
3. 运行优化

## 参数

- 交叉概率：0.8
- 变异概率：0.1

这些参数可以在`NSGAIIOrder.m`中的`OperatorOrder`函数中调整。
