#!/bin/bash
set -e  # Exit on error

# KEGG Database Update Script
# This script fetches the latest KEGG modules data, regenerates graphs,
# and updates the package data files.

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Default values
UPDATE_DATE="${UPDATE_DATE:-$(date +'%B %d, %Y')}"
BRANCH_DATE="${BRANCH_DATE:-$(date +'%Y-%m-%d')}"
MAX_WORKERS="${MAX_WORKERS:-5}"
DELAY="${DELAY:-0.5}"
MAX_RETRIES="${MAX_RETRIES:-5}"

info "KEGG Database Update - $UPDATE_DATE"
info "================================================"

# Step 1: Extract old format files from modules_table.tsv
info "Step 1: Extracting old format files for change detection..."
if [ -f "kegg_pathways_completeness/pathways_data/modules_table.tsv" ]; then
    mkdir -p old_format

    # Extract definitions (module:definition)
    tail -n +2 kegg_pathways_completeness/pathways_data/modules_table.tsv | \
        awk -F'\t' '{print $1":"$2}' > old_format/all_pathways.txt

    # Extract names (module:name)
    tail -n +2 kegg_pathways_completeness/pathways_data/modules_table.tsv | \
        awk -F'\t' '{print $1":"$3}' > old_format/all_pathways_names.txt

    # Extract classes (module:class)
    tail -n +2 kegg_pathways_completeness/pathways_data/modules_table.tsv | \
        awk -F'\t' '{print $1":"$4}' > old_format/all_pathways_class.txt

    info "Extracted old format files for comparison"
else
    warn "No existing modules_table.tsv found. First-time setup."
fi

# Step 2: Fetch latest KEGG modules data
info "Step 2: Fetching latest KEGG modules data..."
if [ -d "old_format" ]; then
    fetch_modules_data \
        -o fetched_data \
        --old-definitions old_format/all_pathways.txt \
        --old-names old_format/all_pathways_names.txt \
        --old-classes old_format/all_pathways_class.txt \
        --max-workers "$MAX_WORKERS" \
        --delay "$DELAY" \
        --max-retries "$MAX_RETRIES"
else
    fetch_modules_data \
        -o fetched_data \
        --max-workers "$MAX_WORKERS" \
        --delay "$DELAY" \
        --max-retries "$MAX_RETRIES"
fi

# Check if changed.tsv was generated
if [ -f "fetched_data/changed.tsv" ]; then
    CHANGED_COUNT=$(tail -n +2 fetched_data/changed.tsv | wc -l | tr -d ' ')
    HAS_CHANGES="true"
    info "Found $CHANGED_COUNT changed modules"
else
    CHANGED_COUNT=0
    HAS_CHANGES="false"
    info "No changes detected"
fi

# Count total modules
TOTAL_MODULES=$(tail -n +2 fetched_data/modules_table.tsv | wc -l | tr -d ' ')
info "Total modules: $TOTAL_MODULES"

# Export variables for later use
export HAS_CHANGES
export CHANGED_COUNT
export TOTAL_MODULES

# Step 3: Update graphs
if [ "$HAS_CHANGES" = "true" ]; then
    info "Step 3: Performing incremental graph update for $CHANGED_COUNT changed modules..."
    make_graphs \
        -i fetched_data/modules_table.tsv \
        -o graphs_output \
        -e kegg_pathways_completeness/pathways_data/graphs.pkl \
        -c fetched_data/changed.tsv \
        -v
else
    info "Step 3: No changes detected, but verifying graphs are up to date..."
    make_graphs \
        -i fetched_data/modules_table.tsv \
        -o graphs_output \
        -e kegg_pathways_completeness/pathways_data/graphs.pkl \
        -v
fi

# Step 4: Update package data files
info "Step 4: Updating package data files..."
cp fetched_data/modules_table.tsv kegg_pathways_completeness/pathways_data/
cp graphs_output/graphs.pkl kegg_pathways_completeness/pathways_data/
cp fetched_data/mediators_list.txt kegg_pathways_completeness/pathways_data/
info "Updated modules_table.tsv, graphs.pkl and mediators_list.txt"

# Step 5: Update README with module count
info "Step 5: Updating README with module count..."
if [ -f "README.md" ]; then
    sed -i.bak "s/The current version includes [0-9]* KEGG modules (updated .*)\./The current version includes $TOTAL_MODULES KEGG modules (updated $UPDATE_DATE)./" README.md
    rm -f README.md.bak  # Remove backup file created by sed
    info "Updated README.md"
else
    warn "README.md not found, skipping update"
fi

# Step 6: Generate change summary
if [ "$HAS_CHANGES" = "true" ]; then
    info "Step 6: Generating change summary..."
    {
        echo "# KEGG Database Update Summary"
        echo ""
        echo "**Update Date**: $UPDATE_DATE"
        echo "**Total Modules**: $TOTAL_MODULES"
        echo "**Changed Modules**: $CHANGED_COUNT"
        echo ""
        echo "## Changed Modules"
        echo ""
        echo "| Module | Definition | Name | Class |"
        echo "|--------|------------|------|-------|"
        tail -n +2 fetched_data/changed.tsv | head -20 | \
            awk -F'\t' '{print "| " $1 " | " $2 " | " $3 " | " $4 " |"}'
    } > update_summary.md

    if [ "$CHANGED_COUNT" -gt 20 ]; then
        echo "" >> update_summary.md
        echo "*... and $((CHANGED_COUNT - 20)) more*" >> update_summary.md
    fi

    {
        echo ""
        echo "## Update Method"
        echo ""
        echo "- Used incremental update (only regenerated changed modules)"
        echo "- Reused existing graphs for unchanged modules"
        echo "- Source of truth: \`modules_table.tsv\`"
    } >> update_summary.md

    info "Generated update_summary.md"
else
    info "Step 6: No changes to summarize"
    # Create a minimal summary even if no changes
    {
        echo "# KEGG Database Update Summary"
        echo ""
        echo "**Update Date**: $UPDATE_DATE"
        echo "**Total Modules**: $TOTAL_MODULES"
        echo "**Changed Modules**: 0"
        echo ""
        echo "No changes detected in KEGG modules data."
    } > update_summary.md
fi

info "================================================"
info "Database update completed successfully!"
info "Changed modules: $CHANGED_COUNT"
info "Total modules: $TOTAL_MODULES"
