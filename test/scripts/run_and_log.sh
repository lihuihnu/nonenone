#!/usr/bin/env bash
set -euo pipefail

# Run one test command, preserve its complete output, show a clean console
# stream, and return the tested command's real exit status.
#
# On the target cluster Slurm's PMI2 plugin may print malformed-key parser
# diagnostics although MPI_Init and the MPMC test itself succeed.  The key bytes
# can contain CR characters or even embedded newlines, so exact line-based
# regular expressions are not reliable.  We therefore identify only the two
# known PMI2 parser-message families by their fixed prefix and, when a malformed
# key is split across lines, suppress its continuation through the terminating
# "in req" fragment.
#
# Nothing is discarded from disk: <log>.raw always stores the byte stream from
# the launcher.  <log> and the console contain the filtered view.
#
# Disable filtering when debugging Slurm/PMI2 itself:
#   MPMC_FILTER_PMI2_WARNINGS=0 make run-full ...
#
# Usage:
#   run_and_log.sh <log-file> <command> [args ...]

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <log-file> <command> [args ...]" >&2
    exit 2
fi

log_file="$1"
shift
mkdir -p "$(dirname "${log_file}")"
raw_log="${log_file}.raw"
filter_pmi2="${MPMC_FILTER_PMI2_WARNINGS:-1}"

# Keep the filtering logic in one place so the console and normal log are
# generated identically.  LC_ALL=C makes awk operate byte-wise even when the
# malformed PMI2 key contains invalid/non-text bytes.
filter_stream() {
    MPMC_LOG_FILE="${log_file}" LC_ALL=C awk '
        BEGIN {
            filtered = 0
            continuation = 0
        }

        {
            line = $0
            sub(/\r$/, "", line)

            # A malformed key can contain an embedded newline.  Once a known
            # "no value for key" message starts without its terminating
            # "in req", suppress continuation lines until that terminator.
            if (continuation) {
                if (index(line, "in req") > 0) {
                    continuation = 0
                }
                next
            }

            if (index(line, "slurmstepd: error: mpi/pmi2: no value for key") == 1) {
                filtered++
                if (index(line, "in req") == 0) {
                    continuation = 1
                }
                next
            }

            if (index(line, "slurmstepd: error: mpi/pmi2: value not properly terminated in client request") == 1) {
                filtered++
                next
            }

            print line
            fflush()
        }

        END {
            if (filtered > 0) {
                printf("[INFO] filtered %d known Slurm PMI2 parser message(s); raw output: %s.raw\n",
                       filtered, ENVIRON["MPMC_LOG_FILE"])
            }
        }
    '
}

# Do not let `set -e` terminate before PIPESTATUS[0] is captured.  stdout and
# stderr are deliberately merged because yhrun/slurmstepd diagnostics are sent
# on stderr.
set +e
if [[ "${filter_pmi2}" == "1" ]]; then
    "$@" 2>&1 | tee "${raw_log}" | filter_stream | tee "${log_file}"
    rc=${PIPESTATUS[0]}
else
    "$@" 2>&1 | tee "${raw_log}" | tee "${log_file}"
    rc=${PIPESTATUS[0]}
fi
set -e

exit "${rc}"
