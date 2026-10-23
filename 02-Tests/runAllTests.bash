#!/bin/bash
#
# Run the standard test scripts in 02-Tests and print a summary.
#
# Usage, from the top directory of the repository after "source startup.bash":
#   bash 02-Tests/runAllTests.bash               # run all tests
#   bash 02-Tests/runAllTests.bash DriftTst ...  # run the named tests only
#
# Each test is judged to pass if it exits with status 0.  The output of
# each test is written to 99-Scratch/TestLogs/<test>.log.  Plots are made
# with the non-interactive "Agg" backend so that no windows are opened.
# A complete run takes about 15 minutes.

if [ -z "${LhARAOpticsPATH}" ]; then
  echo "ERROR: environment not set up; execute \"source startup.bash\" first."
  exit 1
fi
if [ "$(pwd -P)" != "$(cd "${LhARAOpticsPATH}" && pwd -P)" ]; then
  echo "ERROR: run from the top directory of the repository:"
  echo "       cd ${LhARAOpticsPATH}"
  exit 1
fi

export MPLBACKEND=Agg
LOGDIR="99-Scratch/TestLogs"
mkdir -p "${LOGDIR}"

#.. Tests that read files written by other tests run last, in this order;
#   cProfileStats is a utility that needs a profile written by cProfile:
LAST="RunDRACOsimulation ReadDRACOsimu RunLIONsimulation ReadLIONsimu
      plotTrajectory RunLhARAsimulation ReadLhARAsimu TestWrite TestRead
      RunPIONdecayTST decayCHAIN-vis"
SKIP="cProfileStats"

if [ $# -gt 0 ]; then
  TESTS="$*"
else
  TESTS=""
  for f in 02-Tests/*.py; do
    t=$(basename "$f" .py)
    case " $(echo ${LAST} ${SKIP}) " in
      *" $t "*) ;;
      *) TESTS="${TESTS} $t" ;;
    esac
  done
  TESTS="${TESTS} $(echo ${LAST})"
fi

nPass=0
nFail=0
Failed=""
for t in ${TESTS}; do
  Args=""
  if [ "$t" == "decayCHAIN-vis" ]; then
    Args="-i 99-Scratch/pionDECAYtst.dat -o 99-Scratch/pionDECAYvis.pdf -n 100"
  fi
  Start=$(date +%s)
  python "02-Tests/$t.py" ${Args} < /dev/null > "${LOGDIR}/$t.log" 2>&1
  rc=$?
  Time=$(( $(date +%s) - Start ))
  if [ ${rc} -eq 0 ]; then
    nPass=$(( nPass + 1 ))
    printf "  PASS  %-36s %5ds\n" "$t" "${Time}"
  else
    nFail=$(( nFail + 1 ))
    Failed="${Failed} $t"
    printf "  FAIL  %-36s %5ds   (see ${LOGDIR}/$t.log)\n" "$t" "${Time}"
  fi
done

echo
echo "  Tests passed: ${nPass}, failed: ${nFail}"
if [ ${nFail} -gt 0 ]; then
  echo "  Failed:${Failed}"
  exit 1
fi
exit 0
