# MultiSTE-Filter-for-Sensor-Grid

This repository provides a reference implementation of a **hybrid Bayesian particle filter** for **multi-source term estimation (MultiSTE)** using measurements from a **static sensor grid**. The filter estimates the locations and emission rates of multiple airborne release sources with an **unknown and bounded number of sources**, under a superposition-based measurement model.

The implementation follows the inference framework proposed in the paper:

**Mr.MSTE: Multi-robot Multi-Source Term Estimation with Wind-Aware Coverage Control**  
arXiv: https://arxiv.org/abs/2512.17001

Although the full framework in the paper addresses **mobile multi-robot sensing**, this repository focuses on the **sensor-grid benchmark configuration**, which replaces mobile agents with a fixed spatial sensor network. The same hybrid particle filter and physics-informed state transition model are used, enabling direct comparison between static and mobile sensing strategies.

### Key features
- Hybrid particle filter for joint estimation of **source parameters and source count**
- Physics-informed **birth, death, and merge** transition model
- Permutation-invariant representation of multi-source beliefs
- Superposition-based likelihood model with detection thresholds
- Benchmark setup for evaluating static sensor networks

### Purpose
This code is released to:
- facilitate reproducible research,
- support benchmarking against mobile sensing strategies,
- encourage extensions to alternative sensing layouts and inference models.

If you use this code in your research, please cite the associated paper.
