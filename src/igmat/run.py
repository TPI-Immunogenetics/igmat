import os
import sys
import traceback
import multiprocessing
import signal
import shutil
import time
import argparse
import logging

from . import logger
from . import fasta
from . import helpers
from . import counter
from . import igmat

from igmat.hmm.manager import Manager

# Install custom logger
logger.install()

class ServiceExit(Exception):

  """
  Custom exception which is used to trigger the clean exit
  of all running processes and the main program.
  """
  pass

def serviceShutdown(signum, frame):
  logging.warning('Caught signal {signum}'.format(signum=signum))
  raise ServiceExit

class outputWorker(multiprocessing.Process):

  def __init__(self, outputQueue, kwargs):
    super(outputWorker, self).__init__()
    self.kwargs = kwargs
    self._stop = multiprocessing.Event()
    self.outputQueue = outputQueue

  def run(self):

    # Open output handles
    bedHandle = open(self.kwargs['annotation'] + '.bed', 'w+') if self.kwargs['annotation'] else None
    fastaHandle = open(self.kwargs['annotation'] + '.fa', 'w+') if self.kwargs['annotation'] else None
    alignHandle = open(self.kwargs['annotation'] + '.align.fa', 'w+') if self.kwargs['annotation'] else None

    # Open the log file
    logHandle = open(self.kwargs['log'], 'w+') if self.kwargs['log'] else None

    while True:

      # Stop thread
      if self.isStopped():
        break

      if self.outputQueue.empty():
        time.sleep(0.1)
        continue

      # Get data
      data = self.outputQueue.get()

      # Store log file
      if logHandle:
        logHandle.write('{name}\t{status}\n'.format(name=data['name'], status=data['status']))

      # The sequence has been annotated
      if data['status'] == 'annotated':

        # Verbose
        if self.kwargs['verbose']:
          print('Query \'{name}\':\n'.format(name=data['name']))
          print(data['data'])

        # Store bed file
        if bedHandle:
          for feature in data['data'].annotations():
            bedHandle.write('{0}\t{1}\t{2}\t{3}\n'.format(
              data['name'].replace(' ', '_'), 
              feature['start']-data['data'].start, 
              feature['stop']+1-data['data'].start, 
              feature['type']
            ))

        # Store fasta file
        if fastaHandle:
          sequence = data['sequence'][ data['data'].start:(data['data'].end+1)]
          fastaHandle.write('>{name}\n{sequence}\n'.format(name=data['name'].replace(' ', '_'), sequence=sequence))

        # Store fasta alignment file
        if alignHandle:
          sequence = data['data'].sequence
          alignHandle.write('>{name}\n{sequence}\n'.format(name=data['name'].replace(' ', '_'), sequence=sequence))

      # Job done
      self.outputQueue.task_done()

    # Thread terminated
    logging.info('Output Process terminated')

    # Close handle
    if bedHandle:
      bedHandle.close()

    if fastaHandle:
      fastaHandle.close()

    if alignHandle:
      alignHandle.close()

    if logHandle:
      logHandle.close()

  def stop(self):
    self._stop.set()

  def isStopped(self):
    return self._stop.is_set()

class processWorker(multiprocessing.Process):

  def __init__(self, index, inputQueue, outputQueue, kwargs):
    super(processWorker, self).__init__()
    self.index = index
    self.kwargs = kwargs
    self._stop = multiprocessing.Event()
    self.inputQueue = inputQueue
    self.outputQueue = outputQueue
    
  def process_nucleotide(self, sequence, model, threshold=0, restrict=[], hmmerpath=None):
    
    frame_list = fasta.translate(sequence)
    
    # Keep the result with the smaller evalue
    result = None
    for frame, sequence in frame_list.items():
      try:
        dataList = igmat.annotate(sequence, 
          model=model,
          restrict=restrict,
          threshold=threshold)
        for idx, value in enumerate(dataList):
          if result is None or value.evalue < result.evalue:
            result = dataList
      except:
        pass
      
    if result is None:
      raise Exception('Unable to annotate sequence')
    
    return result

  def run(self):

    while True:

      # Stop thread
      if self.isStopped():
        break

      sequence = None
      try:
        sequence = self.inputQueue.pop()
      except Exception as e:
        break   # List is empty

      self.kwargs['result'].increment('total')

      # Processing sequence
      logging.info("Processing sequence \'{name}\' with worker #{index}".format(index=self.index, name=sequence.getName()))
      try:

        # Process sequence
        if sequence.getType() == 'protein':
          dataList = igmat.annotate(sequence, 
            self.kwargs['dataset'],
            restrict=self.kwargs['restrict'],
            threshold=self.kwargs['bit_score_threshold'])
          if dataList is None:
            raise Exception('No annotation found!')
        elif sequence.getType() == 'nucleotide':
          dataList = self.process_nucleotide(
            sequence=sequence.getSequence(), 
            model=self.kwargs['dataset'], 
            restrict=self.kwargs['restrict'], 
            threshold=self.kwargs['bit_score_threshold'])
        else:
          raise Exception('Invalid sequence type')
        
        for idx, value in enumerate(dataList):
          sequence_name = sequence.getName() if len(dataList) == 1 else '{0}_{1}'.format(sequence.getName(),idx+1)
          self.outputQueue.put({
            'data': value,
            'status': 'annotated',
            'sequence': sequence.getSequence(),
            'name': sequence_name
          })

        # All done
        self.kwargs['result'].increment('success')

      except Exception as e:
        logging.error('Error for sequence {name}: {error}'.format(name=sequence.getName(), error=str(e)))

        # Append a result anyway
        self.outputQueue.put({
          'data': None,
          'status': 'failed',
          'name': sequence.getName(),
          'sequence': sequence.getSequence()
        })

        # Increment the number of failed jobs
        self.kwargs['result'].increment('failed')

    logging.info('Worker #{index} stopped'.format(index=self.index))

  def stop(self):
    self._stop.set()

  def isStopped(self):
    return self._stop.is_set()

def run(input, model, restrict=[], logPath=None, annotationPath=None, ncpu=1, bit_score_threshold=80, verbose=False, hmmerpath=None):

  if not input:
    raise Exception('No input provided')

  # Check if there should be some restriction as to which chain types should be numbered.
  # If it is not the imgt scheme they want then restrict to only igs (otherwise you'll hit assertion errors)
  types_to_chains = {"ig":["H","K","L"], "tr":["A", "B","G","D"], "heavy":["H"], "light":["K","L"] }
  restrictList = []
  for r in restrict:
   restrictList += types_to_chains[r] if r in types_to_chains else [r]

  # TODO: code cleanup
  restrict = restrictList

  # # Get process start time
  start_time = time.time()

  # Create the manager object
  manager = Manager(hmmerpath)

  # Load the HMMR model
  dataset = manager.load(model)

  # Initialize the result counter
  result = counter.Counter(['total', 'success', 'failed'])

  manager = multiprocessing.Manager()
  workQueue = manager.list()

  # Parse the input sequence or fasta file.
  if os.path.isfile(input):

    # Read the sequences. All are read into memory currently...
    for sequence in fasta.parse(input):
      workQueue.append(sequence)

  elif isinstance(input, str): # Single sequence
    workQueue.append(fasta.sequence('Input sequence', input))

  # Register the signal handlers
  signal.signal(signal.SIGTERM, serviceShutdown)
  signal.signal(signal.SIGINT, serviceShutdown)

  # Generate the output queue
  outputQueue = multiprocessing.JoinableQueue()

  # Initialize the output process
  outputProcess = outputWorker(outputQueue, kwargs={
    'log': logPath,
    'annotation': annotationPath,
    'verbose': verbose
  })
  outputProcess.start()

  # Initialize the processes
  processList = []
  try:

    # Start the processes
    logging.info('Starting {size} worker processes'.format(size=ncpu))
    for i in range(ncpu):

      # Create a new process worker
      logging.info('Starting worker #{index}'.format(index=i))
      worker = processWorker(i, workQueue, outputQueue, kwargs={
        'restrict': restrict,
        'dataset':dataset,
        'bit_score_threshold': bit_score_threshold,
        'result':result
      })
      worker.start()

      # Append to processes list
      processList.append(worker)

    # Keep the main process running
    while True:

      # Check if process are running
      terminate = True
      for worker in processList:
        if worker.is_alive():
          terminate = False
          break

      # All workers are terminated
      if terminate:
        break

      # Wait
      time.sleep(0.1)

    # Wait for the output worker to finish
    outputQueue.join()
    outputProcess.stop()
    outputProcess.join()

  except ServiceExit:
      
    # Send stop signal
    logging.info('Stopping {size} processes'.format(size=len(processList)))
    for worker in processList:
      worker.stop()

  # Record the end time
  time_taken = time.time() - start_time

  # Convert to minutes and hours
  seconds = time_taken
  minutes = seconds / 60
  hours = minutes / 60
  print('Total execution time: {}:{}:{:.2f}'.format(int(hours), int(minutes), seconds))
  print('Total: {}; Success: {}; Failed: {}'.format(result.value('total'), result.value('success'), result.value('failed')))
  
  # sys.exit(0)
