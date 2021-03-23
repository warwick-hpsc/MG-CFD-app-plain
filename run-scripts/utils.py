from sets import Set
import copy, subprocess

def check_papi_batch_compatibility(papi_batch, event_type):
  compatible = False

  if event_type != "PRESET" and event_type != "NATIVE":
    print("ERROR: 'event_type' must be PRESET or NATIVE")
    quit()

  cmd = ["papi_event_chooser", event_type]
  for e in papi_batch:
    cmd.append(e)

  try:
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    output_lines = output.split('\n')
    for line in output_lines:
      if line.startswith("event_chooser.c"):
        if line.endswith("PASSED"):
          compatible = True
          break
      elif line.startswith("Total events reported:"):
        compatible = True
        break
  except OSError as error:
    ## papi_event_chooser not present. Assume events are compatible
    return True
  except subprocess.CalledProcessError as error:
    # exception_handled = False
    # output_lines = error.output.split('\n')
    # for e in papi_batch:
    #   s = "Event {0} can't be counted with others".format(e)
    #   for line in output_lines:
    #     if line.startswith(s):
    #       exception_handled = True
    #       break
    #   if exception_handled:
    #     break

    # for e in papi_batch:
    #   s = "Event {0} can't be found".format(e)
    #   for line in output_lines:
    #     if line.startswith(s):
    #       event_exists = False
    #       exception_handled = True

    # if not exception_handled:
    #   print "Execution failure detected:"
    #   print "  STDERR: '" + error.output + "'"
    #   sys.exit(-1)

    ## Update: PAPI 5.5.1 does not give a reason for incompatibility between 
    ## native events:
    exception_handled = True

  return compatible

def batch_papi_events(papi_preset_events, papi_native_events, batch_cache_filepath=None):
  papi_event_sets_to_batches = {}
  if not batch_cache_filepath is None and os.path.isfile(batch_cache_filepath):
    current_key = frozenset([])
    with open(batch_cache_filepath) as infile:
      for line in infile:
        line = line.rstrip('\n')
        if line.startswith("key;"):
          current_key = line.split(";")[1:]
          current_key = ';'.join(current_key)
          current_key = current_key.split(",")
          current_key = frozenset(current_key)
          papi_event_sets_to_batches[current_key] = []
        elif line.startswith("batch;"):
          batch = line.split(";")[1:]
          batch = ';'.join(batch)
          batch = batch.split(",")
          papi_event_sets_to_batches[current_key].append(batch)

  papi_preset_event_batches = []
  if len(papi_preset_events) > 0:
    if frozenset(papi_preset_events) in papi_event_sets_to_batches.keys():
      papi_preset_event_batches = papi_event_sets_to_batches[frozenset(papi_preset_events)]
    else:
      print("  batching preset events")
      for e in papi_preset_events:
        e_processed = False
        L = len(papi_preset_event_batches)
        for i in range(L):
          b_copy = copy.deepcopy(papi_preset_event_batches[i])
          b_copy.add(e)
          if check_papi_batch_compatibility(b_copy, "PRESET"):
            papi_preset_event_batches[i].add(e)
            e_processed = True
            break
        if not e_processed:
          papi_preset_event_batches += [Set([e])]
      papi_event_sets_to_batches[frozenset(papi_preset_events)] = papi_preset_event_batches

  papi_native_event_batches = []
  if len(papi_native_events) > 0:
    if frozenset(papi_native_events) in papi_event_sets_to_batches.keys():
      papi_native_event_batches = papi_event_sets_to_batches[frozenset(papi_native_events)]
    else:
      for e in papi_native_events:
        e_processed = False
        L = len(papi_native_event_batches)
        for i in range(L):
          b_copy = copy.deepcopy(papi_native_event_batches[i])
          b_copy.add(e)
          if check_papi_batch_compatibility(b_copy, "NATIVE"):
            papi_native_event_batches[i].add(e)
            e_processed = True
            break
        if not e_processed:
          papi_native_event_batches += [Set([e])]
      papi_event_sets_to_batches[frozenset(papi_native_events)] = papi_native_event_batches

  if not batch_cache_filepath is None:
    with open(batch_cache_filepath, "w") as outfile:
      for k in papi_event_sets_to_batches.keys():
        outfile.write("key;" + ",".join(k) + '\n')
        batches = papi_event_sets_to_batches[k]
        for b in batches:
          outfile.write("batch;" + ",".join(b) + '\n')
      outfile.flush()

  batched_papi_events = papi_preset_event_batches + papi_native_event_batches
  return batched_papi_events
