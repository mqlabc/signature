# 复杂系统signature算法的开源实现
## sig2.2
sig2.2包含了常用复杂系统signature算法的实现，其目录结构如下：
```
sig2.2
├── DESCRIPTION
├── NAMESPACE
├── R
│   ├── auto_judge.R
│   ├── fra.R
│   ├── get_sig_boland.R
│   ├── get_sig_da.R
│   ├── get_sig_def.R
│   ├── load_packages.R
│   ├── minimalVertexPairCutSets.R
│   └── simulate_sig.R
├── man
│   └── hello.Rd
└── sig2.2.Rproj
```
其中，get_sig_def.R使用定义计算signature，get_sig_boland.R实现boland算法，get_sig_da.R实现了达高峰老师提出的算法，simulate_sig.R则使用随机方法模拟构造失效情况计算signature。
以上文件可能使用了minimalVertexPairCutSets.R计算网络结构的最小割集，使用fra.R将小数形式的signature分量转换为分数形式。
auto_judge.R通过人为设置的规则自动判断输入系统结构适用于哪种算法（判断标准是用时长短）。
## shinyapp.R
shinyapp.R使用R提供的shiny框架，封装以上算法，以直观的方式绘出系统的网络结构，并给出系统的signature，并且方便用户对各种方法进行比较。

