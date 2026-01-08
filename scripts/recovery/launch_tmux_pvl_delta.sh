#!/bin/bash

# PVL-Delta Parameter Recovery - Tmux Launcher
#
# This script launches multiple R processes in a tmux session to run
# parameter recovery in parallel. Each pane runs a batch with different
# random seeds so we can aggregate results afterwards.
#
# Usage from project root
#   ./scripts/recovery/launch_tmux_pvl_delta.sh

SESSION_NAME="PVL_Delta_Recovery"
NUM_JOBS=10
ITER_PER_JOB=10
BASE_SEED=3001
OUTPUT_DIR="outputs/recovery/pvl_delta"

mkdir -p "$OUTPUT_DIR"

# Check if session exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching..."
    tmux attach-session -t $SESSION_NAME
    exit 0
fi

echo "Creating new session $SESSION_NAME"
tmux new-session -d -s $SESSION_NAME -n "Recovery"

# Launch first job in the initial pane
SEED=$BASE_SEED
OUTPUT="$OUTPUT_DIR/batch_${SEED}.rds"
CMD="Rscript scripts/recovery/recovery_pvl_delta_batch.R --seed $SEED --iter $ITER_PER_JOB --output $OUTPUT; echo 'Job $SEED Done. Press Enter to close.'; read"
tmux send-keys -t $SESSION_NAME "$CMD" C-m

# Launch remaining jobs in new panes
for ((i=2; i<=NUM_JOBS; i++)); do
    SEED=$((BASE_SEED + i - 1))
    OUTPUT="$OUTPUT_DIR/batch_${SEED}.rds"
    CMD="Rscript scripts/recovery/recovery_pvl_delta_batch.R --seed $SEED --iter $ITER_PER_JOB --output $OUTPUT; echo 'Job $SEED Done. Press Enter to close.'; read"
    
    tmux split-window -t $SESSION_NAME
    tmux select-layout -t $SESSION_NAME tiled
    tmux send-keys -t $SESSION_NAME "$CMD" C-m
done

echo "Launched $NUM_JOBS parallel jobs ($ITER_PER_JOB iterations each = $((NUM_JOBS * ITER_PER_JOB)) total)"
echo "Attach with tmux attach -t $SESSION_NAME"
tmux attach-session -t $SESSION_NAME
