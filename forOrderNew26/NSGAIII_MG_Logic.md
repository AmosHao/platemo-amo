# NSGAIII\_MG 算法逻辑说明（用于 M / G 消融）

本文档解释 `PlatEMO/Algorithms/Multi-objective optimization/NSGA-III-amo/NSGAIII_MG.m` 中 **MG（同时注入 MatingSelectionAMO + OperatorAMO）** 的具体流程，并重点说明：

- **M（`MatingSelectionAMO`）**：如何在目标空间聚类的基础上，为每个父代1匹配一个“尽量跨簇”的父代2（严格索引对应）。
- **G（`OperatorAMO`）**：如何在“可选的严格配对父代2”存在时进行交叉；以及当不提供父代2时，如何在父代集合内随机选择一个不同的父代作为交叉对象。

> 对照消融：Baseline 为 `PlatEMO/Algorithms/Multi-objective optimization/NSGA-III/NSGAIII.m`。  
> 只加 M：`NSGAIII_M.m`；只加 G：`NSGAIII_G.m`；加 M+G：`NSGAIII_MG.m`。

---

## 1. `NSGAIII_MG.m`：每一代在做什么

核心循环（省略环境选择细节）逻辑如下：

1. **父代1选择（p1）**  
   仍然沿用 NSGA-III 的约束锦标赛（constraint tournament），得到长度为 \(N\) 的索引向量 `p1`：
   - `p1 = TournamentSelection(2, Problem.N, CV)`  
   - 其中 \(CV = \sum \max(0, cons)\)（每个个体的总约束违反度）
   - 输出 `p1` **长度为 \(N\)**，允许重复索引（同一父代可多次进入交配池）

2. **父代2匹配（p2）：注入 M（MatingSelectionAMO）**  
   调用 `MatingSelectionAMO`，在目标空间聚类结果上，为每个 `p1(k)` 匹配一个 `p2(k)`：
   - `[p1,p2,state] = MatingSelectionAMO(Population, N, 'Kmax', Kmax, 'p1', p1, 'state', state);`
   - 这里 **`p1` 和 `p2` 是严格一一对应的**：第 \(k\) 对父代就是 `(p1(k), p2(k))`
   - `state` 跨代复用，避免每代从零开始聚类初始化

3. **将索引转为父代个体数组**  
   - `Parents1 = Population(p1(:));`
   - `Parents2 = Population(p2(:));`

4. **子代生成：注入 G（OperatorAMO）并严格配对**  
   - `Offspring = OperatorAMO(Problem, Parents1, 'P2', Parents2, 'pc', pc, 'er', er);`
   - 这一步会生成 **\(N\) 个子代**（每个 `Parents1(i)` 产生 1 个子代）
   - **当发生交叉时**，第二父代固定使用 `Parents2(i)`（严格索引匹配，不随机）

5. **环境选择（NSGA-III 的 EnvironmentalSelection）**  
   将父代种群 `Population` 与子代 `Offspring` 合并后，使用 NSGA-III 环境选择保留 \(N\) 个进入下一代：
   - `Population = EnvironmentalSelection([Population,Offspring], N, Z, Zmin);`

---

## 2. M：`MatingSelectionAMO` 的逻辑（聚类 + 跨簇配偶）

文件：`PlatEMO/Algorithms/Multi-objective optimization/A-amos/MatingSelectionAMO.m`

### 2.1 输入输出

- **输入**
  - `Population`：个体数组，要求每个个体有 `objs`
  - `N`：需要配对的数量（通常等于 `Problem.N`）
  - 可选参数：
    - `'Kmax'`：聚类最多簇数
    - `'state'`：跨代复用的聚类状态（`clustTag/clustName/centroid`）
    - `'p1'`：外部给定的父代1索引（长度 N）；不传则内部随机生成

- **输出**
  - `p1`：父代1索引，1×N
  - `p2`：父代2索引，1×N
  - `state`：更新后的聚类状态（供下一代复用）

### 2.2 聚类是“对整个种群一次性完成”

每一代会对 `Population.objs` 更新一次聚类：

- 若 `state` 无效或长度不匹配，会初始化：
  - `state.clustTag  = (1:popSize)'`
  - `state.centroid  = Population.objs`
- 然后调用：
  - `[clustTag,clustName,centroid] = ObjectiveSpaceClustering(Population.objs, ...)`

其中 `clustTag(i)` 表示第 `i` 个个体所属簇。

### 2.3 为每个 `p1(k)` 选择 `p2(k)`：优先跨簇

对 `k = 1..N`：

- `idx1 = p1(k)`，`tag1 = clustTag(idx1)`
- 候选集合优先取 **不同簇**：
  - `otherIdx = find(clustTag ~= tag1)`，再排除自身 `idx1`
- 如果不同簇为空（全在同簇），退化为“全局随机（除自身）”
- `p2(k)` 从候选集合中随机选一个

因此：

- **`p2(k)` 尽量不与 `p1(k)` 同簇**（跨簇交配增强多样性）
- `p1` 与 `p2` **严格按 k 对应**（配对关系明确）

---

## 3. G：`OperatorAMO` 的逻辑（交叉/变异 + 严格配对父代2）

文件：`PlatEMO/Algorithms/Multi-objective optimization/A-amos/OperatorAMO.m`

### 3.1 输入形式与两种模式

`OperatorAMO(Problem, Parents, ...)` 里：

- `Parents`（本文对应 `Parents1`）：长度为 \(N\) 的父代数组（内部取 `Parents.decs`）
- 可选参数：
  - `'pc'`：交叉概率（默认 0.5）
  - `'er'`：变异探索参数（默认 0.25）
  - `'P2'`：额外父代（本文对应 `Parents2`），用于交叉时严格配对

两种模式：

- **模式 A（严格配对）**：传入 `'P2', Parents2`  
  - 交叉时第二父代固定为 `P2(i,:)`，与 `Parents(i,:)` 严格索引匹配
- **模式 B（无 P2）**：不传 `'P2'`  
  - 交叉时第二父代从 `Parents` 中随机选一个 **且保证不是自己**

### 3.2 严格索引匹配的校验

当传入 `P2` 时，算子会检查：

- `size(P2Dec,1) == size(ParentDec,1)`（行数相同，确保一一配对）
- `size(P2Dec,2) == size(ParentDec,2)`（决策维度相同）

不满足则直接报错，避免“看似运行但父代配错”的隐蔽问题。

### 3.3 每个父代生成 1 个子代（N 个子代）

主体循环 `for i = 1:N`：

1. 设 `p1 = ParentDec(i,:)`
2. 决定是否交叉：
   - `useCrossover = (rand < pc) && (N >= 2)`
3. **如果交叉**
   - 若提供 `P2`：`p2 = P2Dec(i,:)`（严格配对）
   - 否则：从 `1..N` 随机选 `idx2 != i`，取 `p2 = ParentDec(idx2,:)`
   - 生成子代：`child = SegmentBasedCrossover(p1,p2,...)`
4. **否则变异**
   - `child = OrderPreservingMutation(p1, ..., er)`
5. 统一修复段结构：
   - `child = OrderSegmentRepair(child, ...)`
6. 写入 `OffDec(i,:)`

最后：

- `Offspring = Problem.Evaluation(OffDec);`
- **返回子代数量 = N**

---

## 4. 为什么该实现适合做 M / G 消融

- **只加 M**：只改变“父代如何配对”（交配选择），算子仍可用 GA（并可用 `OperatorGAhalf` 统一每代评估为 N）
- **只加 G**：只改变“给定父代集合时如何产生子代”（算子），交配选择仍是原锦标赛
- **加 M+G**：父代2来自聚类跨簇匹配，并且在算子中按索引严格使用该父代2，从而验证 M 与 G 的协同效应

关键是：这些版本都可以做到**每代评估子代数量一致（N）**，保证消融对比不因“评估预算不同”而偏置。

