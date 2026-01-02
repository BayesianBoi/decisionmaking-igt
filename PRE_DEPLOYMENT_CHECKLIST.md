# Pre-Deployment Checklist

Complete this checklist before deploying to GitHub and cloud.

## ✅ Code Verification

- [x] All v2 models exist (pvl_delta_v2, vse_v2, orl_v2)
- [x] fit_models.R uses v2 models
- [x] Data files present (5 IGT datasets)
- [x] Utilities load without errors
- [x] quick_test.R runs successfully
- [x] .gitignore configured
- [x] Documentation complete

## ✅ File Cleanup

Run this to verify no unnecessary files:

```bash
cd /Users/nielsvaerbak/Desktop/decision_making

# Check for backup/temp files
find . -name "*.bak" -o -name "*~" -o -name ".DS_Store"

# Should return nothing or minimal output
```

If files found:
```bash
# Clean them up
find . -name "*.bak" -delete
find . -name "*~" -delete
find . -name ".DS_Store" -delete
```

## ✅ Final Verification Test

```bash
cd /Users/nielsvaerbak/Desktop/decision_making

# This should complete successfully in ~2 minutes
Rscript analysis/quick_test.R
```

Expected output:
- ✓ Data loads (48 subjects)
- ✓ Model compiles without errors
- ✓ R-hat values ~1.0-1.05
- ✓ Parameter estimates reasonable:
  - mu_A: ~0.01-0.10
  - mu_alpha: ~0.1-0.5
  - mu_cons: ~1-3
  - mu_lambda: ~0.1-1

## 📦 GitHub Setup Steps

### Step 1: Create GitHub Repository

Go to: https://github.com/new

Settings:
- **Name**: `igt-decision-models` (or your choice)
- **Description**: `Hierarchical Bayesian models for Iowa Gambling Task`
- **Visibility**: Private ✓ (recommended for research data)
- **Initialize**: Leave all unchecked ✓

Click "Create repository"

### Step 2: Initialize Local Repository

```bash
cd /Users/nielsvaerbak/Desktop/decision_making

# Initialize git
git init

# Add all files
git add .

# Verify what will be committed
git status
```

Review the output - should see:
- ✓ analysis/ directory
- ✓ data/ directory
- ✓ docs/ directory
- ✓ R/ directory
- ✓ README.md, .gitignore, etc.
- ✗ No *.log files
- ✗ No *.bak files
- ✗ No .DS_Store files

### Step 3: Create Initial Commit

```bash
git commit -m "Initial commit: IGT analysis pipeline

- Complete JAGS models (PVL-Delta, VSE, ORL)
- Data harmonization utilities
- Full diagnostic pipeline
- Cloud deployment documentation
- Tested and validated"
```

### Step 4: Connect to GitHub

Replace `YOUR_USERNAME` with your GitHub username:

```bash
# Using HTTPS
git remote add origin https://github.com/YOUR_USERNAME/igt-decision-models.git

# Or using SSH (if you have SSH keys set up)
git remote add origin git@github.com:YOUR_USERNAME/igt-decision-models.git

# Verify
git remote -v
```

### Step 5: Push to GitHub

```bash
git branch -M main
git push -u origin main
```

You should see:
```
Enumerating objects: ...
Counting objects: 100% ...
Writing objects: 100% ...
To github.com:YOUR_USERNAME/igt-decision-models.git
 * [new branch]      main -> main
```

### Step 6: Verify on GitHub

1. Go to: https://github.com/YOUR_USERNAME/igt-decision-models
2. Check that you see:
   - ✓ README.md displays nicely
   - ✓ analysis/ folder
   - ✓ data/ folder
   - ✓ All documentation files

## 🚀 Cloud Deployment

### Step 1: Provision Cloud Instance

Recommended specs:
- **Provider**: AWS, GCP, or DigitalOcean
- **Instance type**: 8 vCPUs, 16GB RAM
- **OS**: Ubuntu 20.04 or 22.04
- **Storage**: 20GB minimum

Example AWS: `c6i.2xlarge`
Example GCP: `c2-standard-8`

### Step 2: Clone Repository on Cloud

```bash
# SSH into your cloud instance
ssh user@your-cloud-ip

# Clone repository
git clone https://github.com/YOUR_USERNAME/igt-decision-models.git
cd igt-decision-models

# If private repo, you'll need to authenticate
# See GITHUB_DEPLOYMENT.md for authentication options
```

### Step 3: Setup Environment

```bash
# Install R and JAGS
sudo apt-get update
sudo apt-get install -y r-base r-base-dev jags

# Install R packages
sudo R -e "install.packages(c('rjags', 'coda', 'dplyr'), repos='https://cloud.r-project.org/')"
```

### Step 4: Verify Setup

```bash
# Quick test (should complete in ~2 minutes)
Rscript analysis/quick_test.R
```

If this passes ✓, you're ready to run the full pipeline!

### Step 5: Run Full Pipeline

```bash
# Start in background
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &

# Monitor progress
tail -f fitting.log
```

Expected runtime: **3-5 hours**

## 📋 Troubleshooting Before Deployment

### Issue: quick_test.R fails

**Check:**
```bash
# Are v2 models being used?
grep "v2.jags" analysis/quick_test.R

# Do all v2 models exist?
ls analysis/models/*_v2.jags
```

Should see all three: pvl_delta_v2, vse_v2, orl_v2

### Issue: Data not loading

**Check:**
```bash
# Do data files exist?
ls data/raw/*/IGTdata*.txt

# Test data loading
Rscript -e "source('analysis/utils/load_data.R'); load_ahn_2014('HC')"
```

### Issue: Git push fails

**Solutions:**

For authentication error:
```bash
# Create personal access token at:
# https://github.com/settings/tokens

# Use token as password when pushing
```

For large file warning:
```bash
# Check file sizes
du -sh data/raw/*

# If > 100MB, consider Git LFS or separate storage
```

## ✅ Final Checklist

Before proceeding:

- [ ] quick_test.R passes locally
- [ ] All files cleaned up (no .bak, .DS_Store)
- [ ] GitHub repository created
- [ ] Code pushed to GitHub successfully
- [ ] Repository visible on GitHub
- [ ] Cloud instance provisioned (or ready to provision)
- [ ] GITHUB_DEPLOYMENT.md reviewed
- [ ] CLOUD_SETUP.md reviewed
- [ ] DEPLOYMENT_CHECKLIST.md ready to use

## 🎯 Success Criteria

After deployment, you should have:

1. **Local**: Clean git repository
2. **GitHub**: Code safely backed up and version controlled
3. **Cloud**: Ready to clone and run
4. **Documentation**: Clear guides for all steps

## 📊 Expected Results

After cloud run completes:

- 3 .rds files (one per model)
- Convergence diagnostics (R-hat < 1.1)
- Model comparison results
- Total cost: ~$0.50-$1.70

## 🆘 Need Help?

1. Check `TROUBLESHOOTING.md` for common issues
2. Verify with `quick_test.R` on smaller data
3. Review error messages in logs
4. Ensure using v2 models (not original versions)

---

**Ready to deploy?** ✅

Follow the steps above, then proceed with:
1. `GITHUB_DEPLOYMENT.md` for GitHub setup
2. `analysis/DEPLOYMENT_CHECKLIST.md` for cloud deployment

Good luck! 🚀
