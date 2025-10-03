#!/bin/bash 
#-----------------------------------------------------------------------------#
#           Group on Data Assimilation Development - GDAD/CPTEC/INPE          #
#-----------------------------------------------------------------------------#
#BOP
#
# !SCRIPT: config_smg.sh
#
# !DESCRIPTION:
#   Script used for SMG configuration and installation.
#
# !CALLING SEQUENCE:
#   ./config_smg.sh <option>
#   To learn more about the available options, execute:
#   ./config_smg.sh help
#
# !REVISION HISTORY:
#
#   20 Dec 2017 - J. G. de Mattos - Initial Version
#   Oct 2023 - J. A. Aravequia - Added HPC system identification
#
# !REMARKS:
#   - SMG paths are defined in etc/paths.sh
#   - Functions are contained in etc/smg_setup.sh
#
#EOP
#-----------------------------------------------------------------------------#
#BOC
# Resolve absolute path to this script (no symlink issues)
SMG_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
export SMG_ROOT

[[ -t 1 ]] && { C_ERR=$'\033[1;31m'; C_RST=$'\033[0m'; } || { C_INFO=; C_RST=; }
init="$SMG_ROOT/etc/__init__.sh"
[[ -r "$init" ]] || { printf '%s[ERROR]%s Missing %s\n' "$C_ERR" "$C_RST" "$init" >&2; exit 2; }
. "$init" || { printf '%s[ERROR]%s Failed to load %s\n' "$C_ERR" "$C_RST" "$init" >&2; exit 3; }

#BOP
# !FUNCTION: execute_function
# !DESCRIPTION:
#   Verifies if the function exists in `smg_setup.sh` and executes it.
#   Passes all arguments to the function, not just the name.
#EOP
execute_function() {
    local function_name=$1
    shift  # Remove function name to pass only actual arguments
    # Debug prints should go through the logging helpers; raw echo here is noisy.
    # _dump_cli "$@"   # uncomment if you want an argv snapshot here
    if declare -f "$function_name" > /dev/null; then
        _log_debug "Running: $function_name $*"
        "$function_name" "$@"  # Execute function with all remaining arguments
        local exit_code=$?
        if [[ $exit_code -ne 0 ]]; then
            _log_err "Function $function_name failed with exit code $exit_code."
            exit $exit_code
        fi
    else
        _log_err "Unknown option: $function_name"
        if declare -F show_help >/dev/null 2>&1; then
            show_help "$SMG_ROOT/etc/smg_setup.sh"
        else
            echo "Available commands:"
            # fallback simples
            declare -F | awk '{print $3}' | grep -vE '^(execute_function|main|_|show_help|get_doc|list_funcs|get_help_block|show_help_func)$' | sort
        fi
        exit 1
    fi

}
#BOP
# !FUNCTION: main
# !DESCRIPTION:
#   Main execution logic of the script, handling input arguments.
#EOP
main() {

    # Load functions from the external file
    export SMG_SETUP_AS_LIBRARY=1
    source "${SMG_ROOT}/etc/smg_setup.sh"
    unset SMG_SETUP_AS_LIBRARY
    
    # Default values if not set by the user
    export compgsi=${compgsi:-false}
    export compang=${compang:-false}
    export compbam=${compbam:-false}
    export compinctime=${compinctime:-false}
    # _dump_cli "$@"   # uncomment for a one-line argv dump at entry

    # Execute the requested function if it exists
    execute_function "$@"
}
# Call the main function
main "$@"

#EOC
#-----------------------------------------------------------------------------#

