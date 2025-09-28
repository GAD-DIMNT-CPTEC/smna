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
[[ -t 1 ]] && { C_INFO=$'\033[1;34m'; C_RST=$'\033[0m'; } || { C_INFO=; C_RST=; }
printf '%s[INFO]%s Installation path: SMG_ROOT=%s\n' "$C_INFO" "$C_RST" "$SMG_ROOT"

#BOP
# !FUNCTION: execute_function
# !DESCRIPTION:
#   Verifies if the function exists in `smg_setup.sh` and executes it.
#   Passes all arguments to the function, not just the name.
#EOP
execute_function() {
    local function_name=$1
    shift  # Remove function name to pass only actual arguments

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
    source "${SMG_ROOT}/etc/smg_setup.sh"

    if [[ $# -eq 0 ]]; then
        _log_warn "No arguments were passed!"
        show_help "$SMG_ROOT/etc/smg_setup.sh"
        exit 1
    fi

    local option=$1
    echo "[INFO] Selected option: $option"

    # Default values if not set by the user
    export compgsi=${compgsi:-false}
    export compang=${compang:-false}
    export compbam=${compbam:-false}
    export compinctime=${compinctime:-false}
    
    # Execute the requested function if it exists
    execute_function "$@"
}
# Call the main function
main "$@"

#EOC
#-----------------------------------------------------------------------------#

