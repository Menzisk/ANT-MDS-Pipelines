<div align="center">

# 🧬 ANT-MDS-Pipelines

### Trajectory processing and analysis pipelines for the ANT(3″)-Ia computational biology study

![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?style=flat-square&logo=python&logoColor=white)
![Field](https://img.shields.io/badge/Field-Computational%20Biology-orange?style=flat-square)
![Domain](https://img.shields.io/badge/Domain-Molecular%20Dynamics-6A5ACD?style=flat-square)
![Target](https://img.shields.io/badge/Target-ANT(3%E2%80%B3)--Ia-2E8B57?style=flat-square)

<br>

**dPCA** · **FEL** · **DCCM** · **Trajectory format conversion**

</div>

---

## Overview

This repository holds the core pipelines used to process and analyse molecular dynamics (MD) trajectories for the **ANT(3″)-Ia** aminoglycoside-modifying enzyme study. It brings together the analyses needed to move from raw trajectories to interpretable descriptions of conformational dynamics.

## What's Inside

| Pipeline | Purpose |
|---|---|
| **dPCA** | Dihedral principal component analysis of backbone conformational motion |
| **FEL** | Free energy landscapes projected onto principal components |
| **DCCM** | Dynamic cross-correlation matrices showing correlated residue motions |
| **Format conversion** | Converting trajectories between formats for downstream tools |

## Repository Structure

This project is part of a **three-repository ecosystem**:

```
ANT-MDS
├── 🔧 ANT-MDS-Pipelines   ← you are here (core pipelines)
├── 📊 ANT-MDS-Analysis    (downstream analysis scripts and results)
└── 🗄️ ANT-MDS-Database    (trajectory data)
```

| Repository | Description | Visibility |
|---|---|---|
| 🔧 **ANT-MDS-Pipelines** | Core processing and analysis pipelines (this repo) | Public |
| 📊 [**ANT-MDS-Analysis**](https://github.com/Menzisk/ANT-MDS-Analysis) | Downstream analysis scripts and results | Public |
| 🗄️ [**ANT-MDS-Database**](https://github.com/Menzisk/ANT-MDS-Database) | Trajectory data | Public |

## 🔄 Workflow

```mermaid
flowchart LR
    A[🗄️ Raw trajectories<br/>ANT-MDS-Database] --> B[🔁 Format conversion]
    B --> C[📐 dPCA]
    B --> D[🔗 DCCM]
    C --> E[🏔️ FEL]
    C --> F[📊 ANT-MDS-Analysis]
    D --> F
    E --> F
```

## Getting Started

```bash
git clone https://github.com/Menzisk/ANT-MDS-Pipelines.git

```

> Requires **Python 3.10+**.

## Citation

Manuscript in preparation. If you use this work, please cite:

> Sikakane M. *et al.*, 2026.

## Author

**Menzi Sikakane**
[![GitHub](https://img.shields.io/badge/GitHub-Menzisk-181717?style=flat-square&logo=github&logoColor=white)](https://github.com/Menzisk)

---

<div align="center">
<sub>Built for the ANT(3″)-Ia computational biology study 🧬</sub>
</div>
