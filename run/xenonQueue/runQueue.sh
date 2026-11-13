#!/bin/bash
# Runs a queue of xenon cases one after the other in the background (runQueue.py); survives logout.
#   ./runQueue.sh N queueNextSteps.py [--only PATTERN ...] [--force]    N = number of cores
#   ./runQueue.sh stop                                                  stop the queue and its solver
# Progress: tail -f queue.log; state of each run: python3 runQueue.py queueNextSteps.py --list
# Starting again resumes a stopped run and skips the finished ones.
# OpenFOAM environment: the sourced one, else $OPENFOAM_BASHRC (default /usr/lib/openfoam/openfoam2606).

cd "${0%/*}" || exit
repo=$(cd ../.. && pwd)

running() { [ -f queue.pid ] && kill -0 "$(cat queue.pid)" 2>/dev/null; }

if [ "$1" = stop ]
then
    if running
    then
        kill -- -"$(cat queue.pid)" && echo "queue stopped (started again, it resumes the stopped run)"
    else
        echo "no queue running"
    fi
    rm -f queue.pid
    exit 0
fi

n=$1
case "$n" in
    ''|*[!0-9]*) echo "usage: $0 N <queue.py> [--only PATTERN ...] [--force]   or   $0 stop"; exit 1 ;;
esac
shift
[ -f "$1" ] || { echo "queue file '$1' not found"; exit 1; }
running && { echo "a queue is already running (pid $(cat queue.pid)); $0 stop first"; exit 1; }

if [ -z "$WM_PROJECT_DIR" ]
then
    . "${OPENFOAM_BASHRC:-/usr/lib/openfoam/openfoam2606/etc/bashrc}" \
        || { echo "OpenFOAM not found: source its etc/bashrc or set OPENFOAM_BASHRC"; exit 1; }
fi
. "$repo/etc/sourceme.sh"
export PATH=$FOAM_USER_APPBIN:$PATH
command -v xenonTransferInterFoam > /dev/null \
    || { echo "xenonTransferInterFoam not built: run ./Allwmake in $repo"; exit 1; }

# own process group, so that stop also ends mpirun and the solver
setsid nohup python3 runQueue.py "$@" --np "$n" >> queue.log 2>&1 < /dev/null &
echo $! > queue.pid
echo "queue started on $n cores (pid $!)"
echo "progress: tail -f queue.log    runs: python3 runQueue.py $1 --list    stop: $0 stop"
