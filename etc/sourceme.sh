#------------------------------------------------------------------------------
# Source this into your shell (do NOT execute it) after sourcing OpenFOAM's
# own etc/bashrc, every time you open a new shell that needs to build or run
# anything in this repository:
#
#   source /path/to/this/repo/etc/sourceme.sh
#
# Why this exists: Allwmake/Allwclean export WM_PROJECT_USER_DIR (and the
# FOAM_USER_APPBIN/FOAM_USER_LIBBIN paths derived from it) so the build
# itself works regardless of where this repo lives or what the local
# username is - but those exports only affect Allwmake's own subshell, not
# the interactive shell that ran it. Without sourcing this file too, your
# shell keeps whatever WM_PROJECT_USER_DIR/PATH it had before, which is why
# built solver names aren't found and a manual `wmake`/`wclean` on a single
# directory fails to find -llagrangianMyIntermediate/-llagrangianMyTurbulence
# even right after a successful Allwmake.
#------------------------------------------------------------------------------

# Respect an already-correct WM_PROJECT_USER_DIR (e.g. set by Allwmake,
# which knows its own location reliably) rather than recomputing it from
# ${BASH_SOURCE:-$0}: that fallback only works under bash - under /bin/sh
# (dash, no BASH_SOURCE) it resolves relative to $0, which breaks if the
# caller has already cd'd elsewhere before sourcing this file.
if [ -z "$WM_PROJECT_USER_DIR" ]; then
    projDir="$(cd "$(dirname "${BASH_SOURCE:-$0}")/.." && pwd)"
    export WM_PROJECT_USER_DIR="$projDir"
    unset projDir
fi
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"

export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"

#------------------------------------------------------------------------------
