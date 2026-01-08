#!/bin/bash

# EEF Parameter Recovery - Tmux Launcher
#
# This script launches multiple R processes in a tmux session to run
# parameter recovery in parallel. Each pane runs a batch with different
# random seeds so we can aggregate results afterwards.
#
# Usage from project root
#   ./scripts/recovery/launch_tmux_eef.sh

SESSION_NAME="EEF_Recovery"
NUM_JOBS=10
ITERS_PER_JOB=10
OUTPUT_DIR="outputs/recovery/parts_eef"

mkdir -p "$OUTPUT_DIR"

# Check if session exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching..."
    tmux attach-session -t $SESSION_NAME
    exit 0
fi

echo "Creating new session $SESSION_NAME"
tmux new-session -d -s $SESSION_NAME

for i in $(seq 1 $NUM_JOBS); do
    SEED=$((1000 + i))
    OUT_FILE="${OUTPUT_DIR}/recovery_eef_batch_${SEED}.rds"
    
    CMD="Rscript scripts/recovery/recovery_eef_batch.R --seed $SEED --iter $ITERS_PER_JOB --output $OUT_FILE; echo 'Job $i Done. Press Enter to close.'; read"

    if [ "$i" -eq 1 ]; then
        tmux send-keys -t $SESSION_NAME "$CMD" C-m
    else
        tmux split-window -t $SESSION_NAME
        tmux select-layout -t $SESSION_NAME tiled
        tmux send-keys -t $SESSION_NAME "$CMD" C-m
    fi
    
    echo "Started Job $i (Seed $SEED)"
done

tmux attach-session -t $SESSION_NAME
