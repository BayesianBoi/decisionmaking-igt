#!/bin/bash

# ==============================================================================
# Model Fitting Tmux Launcher
# ==============================================================================
# This script launches 9 parallel model fitting jobs (3 Models x 3 Groups)
# in a tmux session.
#
# Usage: ./launch_tmux_fit.sh
# ==============================================================================

SESSION_NAME="Model_Fitting"
OUTPUT_DIR="outputs"

# Models and Groups to iterate over
MODELS=("pvl_delta" "orl" "eef")
GROUPS=("HC" "Amph" "Hero")
# Note: "Cannabis" group is not currently in the data loader map.

# Check if session exists
if tmux has-session -t $SESSION_NAME 2>/dev/null; then
    echo "Session $SESSION_NAME already exists. Attaching..."
    tmux attach-session -t $SESSION_NAME
    exit 0
fi

# Create new session
echo "Creating new session: $SESSION_NAME"
# Create first window but don't run command yet
tmux new-session -d -s $SESSION_NAME -n "Manager" "echo 'Starting 9 fitting jobs...'; read"

JOB_COUNT=0

for MODEL in "${MODELS[@]}"; do
    for GROUP in "${GROUPS[@]}"; do
        JOB_COUNT=$((JOB_COUNT + 1))
        
        # Script path
        SCRIPT="scripts/fitting/fit_${MODEL}.R"
        
        # Command to run
        CMD="Rscript $SCRIPT $GROUP; echo 'Job ${MODEL}-${GROUP} Done. Press Enter to close.'; read"
        
        # Determine pane logic
        # We'll create a new window for every 4 jobs to avoid overcrowding, 
        # or just tiled panes in one window. 9 panes is fine in one window but tight.
        # Let's split them into 3 windows (one per model) for cleanliness?
        # Actually, let's do one window per MODEL.
        
        WINDOW_NAME="${MODEL}"
        
        # Check if window exists
        tmux list-windows -t $SESSION_NAME | grep -q "${WINDOW_NAME}"
        if [ $? -ne 0 ]; then
            # Create new window for this model
            tmux new-window -t $SESSION_NAME -n "${WINDOW_NAME}"
            # First pane in this window
            tmux send-keys -t $SESSION_NAME:${WINDOW_NAME} "$CMD" C-m
        else
            # Split existing window
            tmux split-window -t $SESSION_NAME:${WINDOW_NAME} 
            tmux select-layout -t $SESSION_NAME:${WINDOW_NAME} tiled
            tmux send-keys -t $SESSION_NAME:${WINDOW_NAME} "$CMD" C-m
        fi
        
        echo "Launched: $MODEL - $GROUP"
    done
done

# Kill the initial "Manager" window (index 0) if we want, or keep it.
# Let's keep it as a landing page.

echo "All $JOB_COUNT jobs launched."
tmux attach-session -t $SESSION_NAME
