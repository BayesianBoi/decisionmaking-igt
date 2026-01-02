# GitHub Deployment Guide

## Pre-Deployment Verification ✅

All checks passed:
- ✅ 3 v2 model files exist (pvl_delta_v2, vse_v2, orl_v2)
- ✅ fit_models.R configured to use v2 models
- ✅ 5 data files found in data/raw
- ✅ All utility functions present
- ✅ Documentation complete
- ✅ .gitignore configured

## Step-by-Step GitHub Setup

### 1. Create GitHub Repository

**Option A: Via GitHub Website (Recommended)**

1. Go to https://github.com/new
2. Repository settings:
   - Name: `igt-decision-models` (or your choice)
   - Description: `Hierarchical Bayesian models for Iowa Gambling Task data`
   - Visibility: **Private** (recommended, contains research data)
   - **DO NOT** initialize with README, .gitignore, or license
3. Click "Create repository"
4. Keep the page open - you'll need the URLs

**Option B: Via GitHub CLI**

```bash
gh repo create igt-decision-models --private --source=. --remote=origin
```

### 2. Initialize Local Git Repository

```bash
cd /Users/nielsvaerbak/Desktop/decision_making

# Initialize git
git init

# Add all files
git add .

# Check what will be committed
git status

# Create initial commit
git commit -m "Initial commit: Complete IGT analysis pipeline

- PVL-Delta, VSE, and ORL models in JAGS (v2 versions)
- Data loading and harmonization utilities
- Full diagnostic and validation pipeline
- Cloud deployment documentation
- All models tested and working"
```

### 3. Connect to GitHub Remote

```bash
# Add GitHub remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/igt-decision-models.git

# Or if using SSH:
git remote add origin git@github.com:YOUR_USERNAME/igt-decision-models.git

# Verify remote
git remote -v
```

### 4. Push to GitHub

```bash
# Push to GitHub
git branch -M main
git push -u origin main
```

You should see output like:
```
Enumerating objects: XX, done.
Counting objects: 100% (XX/XX), done.
...
To github.com:YOUR_USERNAME/igt-decision-models.git
 * [new branch]      main -> main
```

### 5. Verify on GitHub

1. Go to https://github.com/YOUR_USERNAME/igt-decision-models
2. You should see:
   - `analysis/` folder with all scripts and models
   - `data/` folder with raw data
   - `R/` folder (original reference code)
   - Documentation files (README.md, etc.)

## Deploying to Cloud from GitHub

### Method 1: Direct Clone (Recommended)

On your cloud machine:

```bash
# Install git
sudo apt-get update
sudo apt-get install -y git

# Clone repository
git clone https://github.com/YOUR_USERNAME/igt-decision-models.git
cd igt-decision-models

# If private repo, you'll need authentication:
# Option A: Use personal access token
# Option B: Use SSH key (setup first)
```

**Setting up authentication for private repos:**

```bash
# Option 1: HTTPS with Personal Access Token
# 1. Create token at: https://github.com/settings/tokens
# 2. Clone with token:
git clone https://YOUR_TOKEN@github.com/YOUR_USERNAME/igt-decision-models.git

# Option 2: SSH (More secure)
# 1. Generate SSH key on cloud machine:
ssh-keygen -t ed25519 -C "your_email@example.com"

# 2. Add to GitHub:
cat ~/.ssh/id_ed25519.pub
# Copy output and paste at: https://github.com/settings/keys

# 3. Clone with SSH:
git clone git@github.com:YOUR_USERNAME/igt-decision-models.git
```

### Method 2: Deploy Key (For automation)

If you want the cloud machine to auto-pull updates:

```bash
# On cloud machine, generate deploy key
ssh-keygen -t ed25519 -C "cloud-deploy" -f ~/.ssh/cloud_deploy
cat ~/.ssh/cloud_deploy.pub

# Add this key to your repo's deploy keys:
# GitHub repo → Settings → Deploy keys → Add deploy key
# Paste the public key, optionally allow write access

# Configure git to use this key
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/cloud_deploy

# Clone
git clone git@github.com:YOUR_USERNAME/igt-decision-models.git
```

## Complete Cloud Setup After Clone

```bash
cd igt-decision-models

# Install dependencies
sudo apt-get install -y r-base r-base-dev jags
sudo R -e "install.packages(c('rjags', 'coda', 'dplyr'), repos='https://cloud.r-project.org/')"

# Verify setup
Rscript analysis/quick_test.R

# Run pipeline
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &
```

## Updating Code

If you make changes locally and want to update the cloud:

```bash
# On local machine
git add .
git commit -m "Description of changes"
git push

# On cloud machine
git pull origin main

# Re-run if needed
Rscript analysis/fit_models.R
```

## GitHub Repository Structure

Your repo will look like:

```
igt-decision-models/
├── .gitignore                      # Ignores outputs, system files
├── GITHUB_DEPLOYMENT.md            # This file
├── README.md                       # Project overview
├── analysis/                       # Main analysis pipeline
│   ├── README.md
│   ├── CLOUD_SETUP.md
│   ├── DEPLOYMENT_CHECKLIST.md
│   ├── FINAL_SUMMARY.md
│   ├── IMPLEMENTATION_NOTES.md
│   ├── TROUBLESHOOTING.md
│   ├── fit_models.R
│   ├── quick_test.R
│   ├── test_data_loading.R
│   ├── models/
│   │   ├── pvl_delta_v2.jags
│   │   ├── vse_v2.jags
│   │   ├── orl_v2.jags
│   │   └── ...
│   ├── utils/
│   │   ├── load_data.R
│   │   ├── prepare_jags_data.R
│   │   ├── diagnostics.R
│   │   ├── ppc.R
│   │   └── ...
│   └── outputs/                    # Not tracked (in .gitignore)
├── data/                           # Research data
│   └── raw/
│       ├── Ahn_2014/
│       ├── Fridberg_2010/
│       └── Steingroever_2014/
└── R/                              # Reference code
    ├── 0_preprocess/
    ├── 1_analysis/
    └── ...
```

## Best Practices

### Commits

Make meaningful commits:
```bash
# Good
git commit -m "Fix numerical stability in VSE model"

# Less good
git commit -m "update"
```

### Branches (Optional)

For experimental changes:
```bash
# Create branch for testing
git checkout -b test-new-priors

# Make changes, test
...

# If good, merge back
git checkout main
git merge test-new-priors

# Push
git push origin main
```

### .gitignore

Already configured to ignore:
- Analysis outputs (*.rds, *.pdf)
- System files (.DS_Store)
- R history files
- Log files
- Temporary files

Output files are NOT tracked (you'll download them separately after cloud run).

## Troubleshooting

### Issue: "Permission denied (publickey)"

**Solution:** Set up SSH keys or use HTTPS with token

```bash
# Use HTTPS instead
git remote set-url origin https://github.com/YOUR_USERNAME/igt-decision-models.git
```

### Issue: "Large files warning"

If data files are very large (>100MB):

**Solution 1:** Use Git LFS
```bash
git lfs install
git lfs track "*.rdata"
git add .gitattributes
git commit -m "Track large files with LFS"
```

**Solution 2:** Store data separately
- Upload data to cloud storage (AWS S3, Google Drive, etc.)
- Add download script to setup

### Issue: "Authentication failed"

**Solution:** Create personal access token
1. Go to https://github.com/settings/tokens
2. Generate new token (classic)
3. Select scopes: `repo` (full control)
4. Copy token
5. Use token as password when pushing

## Security Considerations

### For Private Research Data

1. **Use private repository** ✓
2. **Review .gitignore** to ensure no outputs are committed ✓
3. **Consider these files:**
   - Keep: Analysis code, models, documentation
   - Review: Raw data (is it anonymized?)
   - Exclude: Any personally identifiable information

### For Public Repository

If you later want to make it public:

```bash
# Review what would be public
git ls-tree -r main --name-only

# Remove sensitive data if needed
git filter-branch --tree-filter 'rm -rf data/sensitive' HEAD

# Make public via GitHub settings
```

## Quick Reference Commands

```bash
# Clone on cloud
git clone https://github.com/YOUR_USERNAME/igt-decision-models.git

# Update from cloud
git pull origin main

# Commit and push changes
git add .
git commit -m "Description"
git push origin main

# Check status
git status

# View history
git log --oneline
```

## Cost Estimate

**GitHub:**
- Public repos: Free unlimited
- Private repos: Free (up to 500MB, 1GB bandwidth/month)
- This repo: ~10-50 MB (without outputs)
- **Cost: $0** ✓

**Cloud Compute:**
- As documented in CLOUD_SETUP.md
- ~$0.50-$1.70 per full run

---

**Next Steps:**
1. Create GitHub repository
2. Push code following steps above
3. Clone on cloud machine
4. Follow DEPLOYMENT_CHECKLIST.md

Good luck! 🚀
