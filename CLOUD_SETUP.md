# Cloud Execution Guide

This guide describes how to pull the repository to your cloud instance, install dependencies, and run the parameter recovery scripts.

## 1. Pull Repository
Access your cloud terminal and clone/pull the repository:

```bash
# If not cloned yet
git clone https://github.com/BayesianBoi/decisionmaking-igt.git
cd decisionmaking-igt

# If already cloned
cd decisionmaking-igt
git pull origin main
```

## 2. Install System Dependencies (JAGS)

### Option A: Conda (Recommended for Cloud)
If you are in a Conda environment (likely):
```bash
conda install -c conda-forge jags
```

### Option B: Ubuntu/Debian (apt)
If `apt-get` cannot find `jags`:
```bash
sudo add-apt-repository universe
sudo apt-get update
sudo apt-get install -y jags
```

### Option C: Fedora/CentOS/RHEL
```bash
sudo dnf install jags
```

## 3. Install R Packages
Run this command in R to install all required packages (`R2jags`, `parallel`, `ggpubr`, `extraDistr`, `truncnorm`, `pacman`):

```bash
Rscript -e 'if(!require("pacman")) install.packages("pacman"); pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm, rjags, coda)'
```

## 4. Run Parameter Recovery
Execute the recovery scripts. Each script will:
1.  Simulate 48 subjects (matching Ahn 2014 HC group).
2.  Fit the JAGS model.
3.  Save plots to `analysis/outputs/recovery/`.

It is recommended to run them in the background (using `nohup` or `screen`/`tmux`) as they may take time.

### Option A: Sequential Run
```bash
Rscript analysis/1_analysis/2_Recovery/recovery_pvl_delta.R
Rscript analysis/1_analysis/2_Recovery/recovery_orl.R
Rscript analysis/1_analysis/2_Recovery/recovery_eef.R
```

### Option B: Parallel / Background
```bash
nohup Rscript analysis/1_analysis/2_Recovery/recovery_pvl_delta.R > pvl.log 2>&1 &
nohup Rscript analysis/1_analysis/2_Recovery/recovery_orl.R > orl.log 2>&1 &
nohup Rscript analysis/1_analysis/2_Recovery/recovery_eef.R > eef.log 2>&1 &
```

## 5. Check Results
Results and plots will be saved in `analysis/outputs/recovery/` (or similar output directory defined in scripts).
Check the logs (`pvl.log`, etc.) for progress.
