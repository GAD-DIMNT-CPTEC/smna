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

RootDir=$(dirname "${BASH_SOURCE}")
export SMG_ROOT=${RootDir}

echo "[INFO] Installation path: SMG_ROOT=$SMG_ROOT"

#BOP
# !FUNCTION: detect_hpc_system
# !DESCRIPTION:
#   Identifies the HPC system and sets global variables accordingly.
#EOP
detect_hpc_system() {
    local sys_info=$(uname -a)

    if echo "$sys_info" | grep -q "cray_ari_s"; then
        export hpc_name="XC50"
        export SUB="cray"
        echo "[INFO] Detected: Cray XC50"
    
    elif echo "$sys_info" | grep -q "egeon-login1.cptec.inpe.br"; then
        export hpc_name="egeon"
        export SUB="egeon"
        echo "[INFO] Detected: EGEON Cluster"
    
    else
        echo "[ERROR] Unknown machine: $(hostname)"
        echo "[ACTION] 1) Add the machine to the defined systems in etc/mach/"
        echo "[ACTION] 2) Define an option for it in the function copy_fixed_files inside etc/smg_setup.sh"
        return 1
    fi
}

#BOP
# !FUNCTION: execute_function
# !DESCRIPTION:
#   Verifies if the function exists in `smg_setup.sh` and executes it.
#EOP
execute_function() {
    local function_name=$1

    if declare -f "$function_name" > /dev/null; then
        echo "[ OK ] Running: $function_name"
        vars_export
        "$function_name"
    else
        echo "[FAIL] Unknown option: $function_name"
        vars_export
        help
    fi
}

#BOP
# !FUNCTION: main
# !DESCRIPTION:
#   Main execution logic of the script, handling input arguments.
#EOP
main() {
    if [[ $# -eq 0 ]]; then
        echo "[WARNING] No arguments were passed!"
        help
        exit 1
    fi

    local option=$1
    echo "[INFO] Selected option: $option"

    # Default values if not set by the user
    export compgsi=${compgsi:-1}
    export compang=${compang:-1}
    export compbam=${compbam:-1}
    # Call HPC system detection function
    detect_hpc_system

    # Load functions from the external file
    source "${SMG_ROOT}/etc/smg_setup.sh"
    
    # Checks if Conda is active and deactivates it if necessary
    disable_conda
    
    # Execute the requested function if it exists
    execute_function "$option"
}

# Call the main function
main "$@"

#EOC
#-----------------------------------------------------------------------------#

