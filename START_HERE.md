# 🚀 START HERE - Deployment Guide

This is your complete guide to deploying the IGT analysis pipeline.

## ✅ What You Have

A complete, tested, and documented pipeline for fitting hierarchical Bayesian models to Iowa Gambling Task data.

**Status**: Ready for deployment ✅

## 📚 Documentation Overview

1. **`PRE_DEPLOYMENT_CHECKLIST.md`** ← Start here for deployment prep
2. **`GITHUB_DEPLOYMENT.md`** ← GitHub setup guide
3. **`analysis/CLOUD_SETUP.md`** ← Cloud instance setup
4. **`analysis/DEPLOYMENT_CHECKLIST.md`** ← Step-by-step cloud deployment
5. **`README.md`** ← Project overview

## 🎯 Quick Deployment Path

### Option 1: Deploy via GitHub (Recommended)

```bash
# 1. Verify locally (2 minutes)
cd /Users/nielsvaerbak/Desktop/decision_making
Rscript analysis/quick_test.R

# 2. Push to GitHub (follow PRE_DEPLOYMENT_CHECKLIST.md)
git init
git add .
git commit -m "Initial commit: IGT analysis pipeline"
git remote add origin https://github.com/YOUR_USERNAME/igt-decision-models.git
git push -u origin main

# 3. Clone on cloud (follow GITHUB_DEPLOYMENT.md)
ssh user@cloud-ip
git clone https://github.com/YOUR_USERNAME/igt-decision-models.git

# 4. Run on cloud (follow analysis/DEPLOYMENT_CHECKLIST.md)
cd igt-decision-models
sudo apt-get install r-base jags
Rscript analysis/quick_test.R  # Verify
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &
```

### Option 2: Direct Transfer (Alternative)

```bash
# 1. Package code
tar -czf decision_making.tar.gz analysis/ data/ R/

# 2. Transfer to cloud
scp decision_making.tar.gz user@cloud-ip:~/

# 3. Follow analysis/CLOUD_SETUP.md for setup
```

## 📋 What Gets Deployed

✅ **Analysis Pipeline**
- 3 JAGS models (pvl_delta_v2, vse_v2, orl_v2)
- Data loading and harmonization
- Fitting scripts
- Diagnostic utilities
- Model comparison tools

✅ **Data**
- 173 subjects from 3 published studies
- 81,831 total trials

✅ **Documentation**
- Complete setup guides
- Troubleshooting documentation
- Usage examples

## ⚡ What to Expect

**Local Testing**: 2 minutes
- Tests 48 subjects
- Verifies models compile
- Checks convergence

**Cloud Full Run**: 3-5 hours
- Fits all 173 subjects
- 3 models × (2000 burn-in + 5000 samples)
- Cost: ~$0.50-$1.70

**Results**:
- 3 model fit files (.rds)
- Convergence diagnostics
- Parameter estimates
- Model comparison

## 🔍 Pre-Flight Check

Run this now:

```bash
cd /Users/nielsvaerbak/Desktop/decision_making

# Verify everything works
Rscript analysis/quick_test.R

# Expected: Completes in ~2 minutes with R-hat ~1.0
```

If this passes ✅, you're ready to deploy!

## 📖 Detailed Guides

### For GitHub Setup
→ Read: `GITHUB_DEPLOYMENT.md`
- Create repository
- Push code
- Clone on cloud
- Authentication options

### For Cloud Deployment
→ Read: `analysis/CLOUD_SETUP.md`
- Instance requirements
- Installation steps
- Running the pipeline
- Monitoring progress
- Downloading results

### For Step-by-Step Process
→ Read: `analysis/DEPLOYMENT_CHECKLIST.md`
- Complete checklist format
- Every step detailed
- Troubleshooting tips
- Success criteria

## ⚠️ Important Notes

1. **Use v2 models**: All scripts automatically use the stable v2 versions
2. **Private repo recommended**: Data contains research subjects
3. **Cost awareness**: Cloud run ~$0.50-$1.70, remember to stop instance after
4. **Runtime**: Budget 3-5 hours for full pipeline

## 🆘 If Something Goes Wrong

1. **Model won't initialize**
   - Check using v2 models (not original versions)
   - See `analysis/TROUBLESHOOTING.md`

2. **Out of memory**
   - Use larger instance (16GB+ RAM)
   - Or reduce to single study: `fit_all_studies = FALSE`

3. **Authentication errors**
   - See `GITHUB_DEPLOYMENT.md` for token/SSH setup

## ✨ You're Ready!

Everything is tested and validated. Just follow the guides:

1. **Now**: Read `PRE_DEPLOYMENT_CHECKLIST.md`
2. **Next**: Follow `GITHUB_DEPLOYMENT.md`
3. **Then**: Use `analysis/DEPLOYMENT_CHECKLIST.md` on cloud

Good luck with your analysis! 🎉

---

**Questions?**
- Technical details → `analysis/IMPLEMENTATION_NOTES.md`
- Model fixes → `analysis/TROUBLESHOOTING.md`
- Project overview → `README.md`

## 📊 For Manuscript Preparation

After cloud run completes and you've downloaded the fitted models:

```bash
# Generate all publication outputs
Rscript analysis/generate_paper_outputs.R
```

This creates:
- **Trace plots** (convergence checking)
- **Posterior density plots** (parameter distributions)
- **Forest plots** (model comparison)
- **CSV tables** (ready for manuscript)
- **APA-formatted results** (copy-paste ready)

See `analysis/PAPER_WORKFLOW.md` for complete manuscript preparation guide.

