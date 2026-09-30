#!/bin/bash
# list_samples.sh <group_dir>
#   Print what the pipeline sees in a group: sample <TAB> R1 <TAB> R2.
#   Worth running before a job that takes hours - a sample missing here stays missing.
NO_STEP_LOG=1
LIST_ONLY=1        # a person looking: full detail, and nothing refused before it prints
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/lib/common.sh"

[ $# -eq 1 ] || { sed -n '2,4p' "$0"; exit 1; }
init_group "$1"

# Print the group first, then say what is wrong with it. A pre-flight look that passed a layout
# the pipeline is going to refuse would be worse than no pre-flight at all.
rc=0
out="$(list_samples)" || rc=1
[ -z "$out" ] || printf '%s\n' "$out"
check_sample_ids "$out" "$1" || rc=1
exit "$rc"
