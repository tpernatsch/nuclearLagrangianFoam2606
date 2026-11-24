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

# Under bash, BASH_SOURCE reliably gives this file's own location no matter
# what WM_PROJECT_USER_DIR already holds - and it must win: OpenFOAM's own
# etc/bashrc (sourced beforehand, as instructed above) sets a non-empty
# default ($HOME/$WM_PROJECT/$USER-$WM_PROJECT_VERSION) of its own, so "if
# unset" would never recompute it for an interactive `source etc/sourceme.sh`
# and silently keep pointing at that unrelated default instead of this repo.
#
# Under /bin/sh (dash, no BASH_SOURCE - i.e. Allwmake/Allwclean, which export
# WM_PROJECT_USER_DIR=$PWD themselves right before sourcing this file) there
# is no reliable own-location fallback once the caller has already cd'd, so
# trust whatever is already set instead.
if [ -n "$BASH_SOURCE" ]; then
    projDir="$(cd "$(dirname "$BASH_SOURCE")/.." && pwd)"
    export WM_PROJECT_USER_DIR="$projDir"
    unset projDir
elif [ -z "$WM_PROJECT_USER_DIR" ]; then
    projDir="$(cd "$(dirname "$0")/.." && pwd)"
    export WM_PROJECT_USER_DIR="$projDir"
    unset projDir
fi
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"

export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH"

#------------------------------------------------------------------------------
