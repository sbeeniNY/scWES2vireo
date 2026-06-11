#!/usr/bin/env bash
# Map an LSF job status to cluster-generic values: running / success / failed.
#
# Two hard rules:
#  1. Every LSF query is wrapped in `timeout`. On a busy cluster bjobs/bhist can
#     hang; without a timeout the status subprocess never returns and the whole
#     Snakemake driver freezes (it blocks waiting on this script). On timeout we
#     return empty and fall through to "running" so the next poll retries.
#  2. Only POSITIVE evidence leaves "running":
#       success <- bjobs DONE  or bhist "Done successfully"
#       failed  <- bjobs EXIT/ZOMBI  or bhist Exited/TERM_/exit code
#       running <- anything else, including not-found-yet / timed-out query.
#     Never default an unknown job to success (that restarts running jobs ->
#     duplicate LSF jobs + driver death). A finished job is reliably caught by
#     `bjobs -a` (kept for CLEAN_PERIOD, ~1h >> poll interval), so reporting
#     "running" on a transient gap cannot hang a genuinely completed job.
jobid="$1"

status=""
for _ in 1 2; do
    status=$(timeout 30 bjobs -a -noheader -o stat "$jobid" 2>/dev/null | tail -n1 | tr -d '[:space:]')
    [[ -n "$status" ]] && break
    sleep 2
done

case "$status" in
  PEND|PSUSP|SSUSP|USUSP|WAIT|RUN|PROV) echo running; exit 0 ;;
  DONE)                                 echo success; exit 0 ;;
  EXIT|ZOMBI)                           echo failed;  exit 0 ;;
esac

# Rare: not in bjobs (transient gap, or aged past CLEAN_PERIOD). bhist as a fast,
# timeout-guarded safety net (default window; do NOT use -n 0 -> scans all event
# logs and can hang).
hist=$(timeout 30 bhist -l "$jobid" 2>/dev/null)
if grep -q 'Done successfully' <<<"$hist"; then
    echo success
elif grep -Eq 'Exited|TERM_|exit code' <<<"$hist"; then
    echo failed
else
    echo running
fi
