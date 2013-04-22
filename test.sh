#!/bin/bash -Eu
      # -E: ERR trap is inherited by shell functions.
      # -u: Treat unset variables as an error when substituting.
      #
      # Example script for handling bash errors.  Exit on error.  Trap exit.
      # This script is supposed to run in a subshell.
      # See also: http://fvue.nl/wiki/Bash:_Error_handling

      #  Trap non-normal exit signals: 1/HUP, 2/INT, 3/QUIT, 15/TERM, ERR
      trap onexit 1 2 3 15 ERR


      #--- onexit() -----------------------------------------------------
      #  @param $1 integer  (optional) Exit status.  If not set, use `$?'

      function onexit() {
          local exit_status=${1:-$?}
          echo Exiting $0 with $exit_status
          exit $exit_status
      }


      # myscript
      $1

      # Allways call `onexit' at end of script
      onexit