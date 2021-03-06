if [ ! -f <RUN_DIRPATH>/job-is-complete.txt ]; then
  basedir=`pwd`
  cd <RUN_DIRPATH>

  if [ "$submit_cmd" = "" ]; then
    if <COMPILE_ONLY> ; then
	  ./<BATCH_FILENAME> --compile | tee submit.log
	else
      ./<BATCH_FILENAME> --compile --execute | tee submit.log
	fi
  else
    if [ ! -f "job-in-queue.txt" ] && [ ! -f "job-is-running.txt" ]; then
      if [ ! -f <BIN_FILEPATH> ]; then
        ./<BATCH_FILENAME> --compile
      fi
      if ! <COMPILE_ONLY> ; then
        export BATCH_EXECUTE_MODE=execute
        if ! eval "$submit_cmd" ./<BATCH_FILENAME> ; then
          echo "Submission failed for: <RUN_DIRPATH>/<BATCH_FILENAME>"
          exit 1
        fi
        touch "job-in-queue.txt"
        unset BATCH_EXECUTE_MODE
      fi
    fi
  fi
  cd "$basedir"
fi
