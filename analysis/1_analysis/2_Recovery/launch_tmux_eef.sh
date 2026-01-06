#!/bin/bash

# ==============================================================================
# EEF Parameter Recovery - Tmux Launcher
# ==============================================================================
# This script launches multiple R processes in a tmux session to run 
# parameter recovery in parallel.
#
# Usage: ./launch_tmux_eef.sh
# ==============================================================================

SESSION_NAME="EEF_Recovery"
NUM_JOBS=10             # Number of tmux panes (workers)
ITERS_PER_JOB=10        # Iterations per worker (Total = NUM_JOBS * ITERS_PER_JOB = 100)
OUTPUT_DIR="analysis/outputs/recovery/parts_eef"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Clean up any existing partial files? (Optional, maybe ask user)
# rm -f "$OUTPUT_DIR"/*.rds

# Check if session exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching..."
    tmux attach-session -t $SESSION_NAME
    exit 0
fi

# Create new session
echo "Creating new session: $SESSION_NAME"
tmux new-session -d -s $SESSION_NAME

# Loop to create panes and run jobs
# We start from 1 to NUM_JOBS
for i in $(seq 1 $NUM_JOBS); do
    SEED=$((1000 + i))
    OUT_FILE="${OUTPUT_DIR}/recovery_batch_${SEED}.rds"
    
    CMD="Rscript analysis/1_analysis/2_Recovery/recovery_eef_batch.R --seed $SEED --iter $ITERS_PER_JOB --output $OUT_FILE; echo 'Job $i Done. Press Enter to close.'; read"

    if [ "$i" -eq 1 ]; then
        # First pane (already created by new-session)
        tmux send-keys -t $SESSION_NAME "$CMD" C-m
    else
        # Split window and run command
        tmux split-window -t $SESSION_NAME
        tmux select-layout -t $SESSION_NAME tiled
        tmux send-keys -t $SESSION_NAME "$CMD" C-m
    fi
    
    echo "Started Job $i (Seed $SEED)"
done

# Attach to session
tmux attach-session -t $SESSION_NAME
