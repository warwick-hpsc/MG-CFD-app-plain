if [ ! -f <RUN_DIRPATH>/Times.csv ]; then
  basedir=`pwd`
  cd <RUN_DIRPATH>

  if [ "$submit_cmd" = "" ]; then
    ./<BATCH_FILENAME> > submit.log
  else
    if [ ! -f "job-in-queue.txt" ] && [ ! -f "job-is-running.txt" ]; then
      if [[ `hostname` == *"login"* ]]; then
        if [ ! -f <BIN_FILEPATH> ]; then
          ## Run batch script on login node to perform compilation, which will exit before app execution.
          ./<BATCH_FILENAME>
        fi
        if ! eval "$submit_cmd" ./<BATCH_FILENAME> ; then
          echo "Submission failed for: <RUN_DIRPATH>/<BATCH_FILENAME>"
          exit 1
        fi
        touch "job-in-queue.txt"
      else
        if ! eval "$submit_cmd" ./<BATCH_FILENAME> ; then
          echo "Submission failed for: <RUN_DIRPATH>/<BATCH_FILENAME>"
          exit 1
        fi
      fi
    fi
  fi
  cd "$basedir"
fi