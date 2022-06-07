
import copy
import math

class Mapper():

  def __init__(self, sequence):
    self._map = []
    self._solutionList = []

    self._regionList = [
      {'name': 'FR1', 'start': 0, 'stop': 26, 'index': 0},
      {'name': 'CDR1', 'start': 26, 'stop': 38, 'index': 1},
      {'name': 'FR2', 'start': 38, 'stop': 55, 'index': 2},
      {'name': 'CDR2', 'start': 55, 'stop': 65, 'index': 3},
      {'name': 'FR3', 'start': 65, 'stop': 104, 'index': 4},
      {'name': 'CDR3', 'start': 104, 'stop': 117, 'index': 5},
      {'name': 'FR4', 'start': 117, 'stop': 127, 'index': 6}
    ]
    self._regionMap = {}
    self._sequence = sequence

    for i in range(len(self._regionList)):
      name = self._regionList[i]['name']
      self._map.append({
        'index': i,
        'name': name,
        'list': []
        })
      self._regionMap[name] = i

  def append(self, region, item):
    # Get the name index
    index = self._regionMap[region] if region in self._regionMap else None
    if index == None:
      raise Exception('Invalid region name: {0}'.format(region))

    self._map[index]['list'].append(item)

  def _createRegion(self, name):
    '''
    Create an empty region as placeholder
    '''

    region_index = self._regionMap[name]
    region = self._regionList[region_index]
    stateRegion = {
      'score': 0,
      'start': -1,
      'stop': -1,
      'size': 0,
      'vector': []
    }
    for i in range(region['start']+1, region['stop']+1):
      stateRegion['vector'].append({'hmm': i, 'type': '-', 'idx': None, 'score': 0})

    return stateRegion

  def _debugState(self, state_vector):

    for region_index in range(len(self._regionList)):

      for state in state_vector[region_index]['vector']:

        residue = self._sequence[ state['idx'] ] if state['idx'] != None else '-'
        print('{0}\t{1}\t{2}\t{3}\t{4}'.format(
          state['hmm'],
          state['type'],
          state['idx'],
          residue,
          self._regionList[region_index]['name']
          ))

  def _recursiveProcess(self, solution, index):

    # We reached the end of the map
    region_size = len(self._regionList)
    if index >= region_size:

      # Check if all regions are continuous
      isValid = True
      prev_item = None
      for i in range(region_size):
        # prev_item = solution[(i-1)] if i > 0 else None
        # next_item = solution[(i+1)] if i+1 < region_size else None
        # curr_item = solution[i]

        curr_item = solution[i]
        next_item = None
        for j in range(i+1, region_size):
          next_item = solution[j]
          if next_item['size'] > 0:
            break

        # Check if previous region is contiguous
        if (prev_item and prev_item['size'] > 0 and curr_item['size'] > 0) and prev_item['stop'] > curr_item['start']:
          isValid = False
          break

        # Check if next region is contiguous
        if (next_item and next_item['size'] > 0 and curr_item['size'] > 0) and next_item['start'] < curr_item['stop']:
          isValid = False
          break

        prev_item = curr_item

      # Append solution
      if isValid:
        self._solutionList.append(solution)

      return True

    # Iterate current block solutions
    for i in range(len(self._map[index]['list'])):

      # Get current item
      curr_item = self._map[index]['list'][i]

      # Append current solution
      item = copy.deepcopy(solution)
      # item = solution
      item.append(copy.deepcopy(curr_item))
      self._recursiveProcess(item, index+1)

    return False

  def process(self):

    # Add empty region where missing
    for i in range(len(self._map)):
      
      # if len(self._map[i]['list']) == 0:
      # print('Missing region: {0}'.format(self._map[i]['name']))
      self._map[i]['list'].append(self._createRegion(self._map[i]['name']))
      # print(i, self._map[i]['name'],len(self._map[i]['list']))

    # Generate all combination of regions
    for i in range(len(self._map[0]['list'])):

      solution = []
      solution.append(self._map[0]['list'][i])

      # Create a new solution
      self._recursiveProcess(solution, 1)

    solutionScore = []
    for k in range(len(self._solutionList)):

      # The current solution
      state_vector = self._solutionList[k]

      # Check for non contiguous regions
      # This fixes problems when aligning two or more hmms and something is overlapping
      prev_valid = None
      prev_region = None
      for j in range(len(self._regionList)):

        region_name = self._regionList[j]['name']
        
        first_valid = None
        last_valid = None
        for state in state_vector[j]['vector']:
          first_valid = first_valid if first_valid and first_valid['idx'] else state
          first_valid = first_valid if first_valid and first_valid['idx'] else None
          last_valid = state if state['idx'] else last_valid

        # This region is not contiguous, it needs to be removed
        delta_idx = first_valid['idx']-prev_valid['idx'] if (prev_valid and first_valid) else None
        if prev_valid and first_valid and (delta_idx != 1):
          
          # region_remove = region_name if region_name in ['CDR1', 'CDR2', 'CDR3'] else prev_region
          
          # Select which region needs to be removed
          region_remove = None
          if region_name and prev_region:
            region_remove = region_name if region_name in ['CDR1', 'CDR2', 'CDR3'] else prev_region
            if region_name[0] == prev_region[0]:
              prev_region_index = self._regionMap[prev_region]
              region_remove = region_name if state_vector[j]['size'] < state_vector[prev_region_index]['size'] else prev_region

          
          if region_remove:
            region_remove_index = self._regionMap[region_remove]
            state_vector[region_remove_index]['vector'] = []
            state_vector[region_remove_index]['score'] = 0
            state_vector[region_remove_index]['size'] = 0
            for i in range(self._regionList[region_remove_index]['start'], self._regionList[region_remove_index]['stop']):
              state_vector[ region_remove_index ]['vector'].append({
                'hmm': i+1,
                'type': '-',
                'idx': None,
                'score': 0
                })

        # Update for next
        # prev_region = region_name
        prev_region = region_name if last_valid else prev_region
        prev_valid = last_valid if last_valid else prev_valid
  
      # Apply IMGT annotation numbering
      for region_name in ['CDR1', 'CDR2', 'CDR3']:

        # Get all valid residues
        region_index = self._regionMap[region_name]
        region_size = len(state_vector[region_index]['vector'])
        residueList = [ x for x in state_vector[region_index]['vector'] if (x['type'] != 'd' and x['idx']) ]
        if len(residueList) == 0 or (len(residueList) == len(state_vector[region_index]['vector'])):
          continue

        leftSize = math.ceil(len(residueList)/2)
        rightSize = len(residueList)-leftSize
        gapSize = region_size-len(residueList)

        orderedList = []
        # hmmCount = 0
        hmmCount = min(x['hmm'] for x in residueList)-1
        for i in range(region_size):

          if i >= leftSize and i < (gapSize+leftSize):

            if item['type'] != '-':
              hmmCount += 1

            item = {
              'hmm': hmmCount,
              'type': 'd',
              'score': 0,
              'idx': None
            }
          else:
            item = residueList.pop(0)
            if item['type'] == 'm':
              hmmCount += 1
            
            item['hmm'] = hmmCount
            

          orderedList.append(item)
          # hmmCount += 1

        # Update state vector
        state_vector[region_index]['vector'] = orderedList
      
      # Fill in the missing bits
      state_start=-1
      state_end=0
      prev = None
      region_index = 0
      while region_index < len(self._regionList):
        region_pos = 0

        while (region_index < len(self._regionList) and region_pos < len(state_vector[region_index]['vector'])):

          # Update start and end position on the query sequence
          state_start = state_start if (state_start != None and state_start >= 0) else state_vector[region_index]['vector'][region_pos]['idx']
          state_end = state_vector[region_index]['vector'][region_pos]['idx'] if state_vector[region_index]['vector'][region_pos]['idx'] else state_end
          
          # Set current position
          curr = {
            'region': region_index,
            'pos': region_pos,
            'hmm': state_vector[region_index]['vector'][region_pos]['hmm'],
            'idx': state_vector[region_index]['vector'][region_pos]['idx'],
            'type': state_vector[region_index]['vector'][region_pos]['type']
          }

          # There is nothing before or this is a valid position
          if (curr['type'] != '-') or (prev is None or prev['idx'] is None):
            # position += 1
            region_pos += 1
            prev = curr
            continue

          # This is an invalid position, try to find the next valid position
          last = None
          while region_index < len(self._regionList):
            while region_index < len(self._regionList) and region_pos < len(state_vector[region_index]['vector']):
              last = { 
                'region': region_index, 
                'pos': region_pos, 
                'hmm': state_vector[region_index]['vector'][region_pos]['hmm'], 
                'idx': state_vector[region_index]['vector'][region_pos]['idx'], 
                'type': state_vector[region_index]['vector'][region_pos]['type']
              }
              if last['type'] != '-':
                break

              last = None
              region_pos += 1


            # We found a valid position, we are breaking also this loop
            if last and last['type'] != '-':
              break

            region_index += 1
            region_pos = 0

          if last is None:
            prev = curr
            region_pos += 1
            continue

          if not last['idx']:
            continue

          idx = prev['idx']+1
          state_region = curr['region']
          state_pos = curr['pos']
          hmm_pos = curr['hmm']
          while idx < last['idx']:

            if state_pos >= len(state_vector[state_region]['vector']):
              state_pos = 0;
              state_region += 1
              continue

            if state_vector[state_region]['vector'][state_pos]['type'] != '-':
              state_vector[state_region]['size'] += 1
              state_vector[state_region]['vector'].insert(state_pos,{
                'hmm': hmm_pos,
                'idx': idx,
                'type': 'i',
                'score': 0
              })
              state_pos += 1
            else:
              hmm_pos = state_vector[state_region]['vector'][state_pos]['hmm']
              residue_score = state_vector[state_region]['vector'][state_pos]['score'] if state_vector[state_region]['vector'][state_pos]['idx'] == idx else 5
              state_vector[state_region]['size'] += 1
              state_vector[state_region]['vector'][state_pos] = {
                'hmm': hmm_pos,
                'type': 'm',
                'idx': idx,
                'score': residue_score
              }
              
              state_pos += 1
              # hmm_pos += 1

            idx += 1
            if state_pos >= len(state_vector[state_region]['vector']):
              state_pos = 0;
              state_region += 1

        # Increment region index
        region_index += 1

      # Update scores
      for j in range(len(self._regionList)):
        score = sum(state['score'] for state in state_vector[j]['vector'])
        state_vector[j]['score'] = score

      # Store result
      solutionScore.append({
        'start': state_start,
        'end': state_end
      })

    # Validate regions
    k = 0
    resultList = []
    while k < len(self._solutionList):

      score = 0
      solutionMask = ''
      for h in range(len(self._regionList)):
        score += self._solutionList[k][h]['score']
        solutionMask += '1' if self._solutionList[k][h]['size'] > 0 else '0'
        
      # Check if the solution is valid
      solutionMask.strip('0')
      resultList.append({
        'score': score,
        'mask': solutionMask,
        'start': solutionScore[k]['start'],
        'end': solutionScore[k]['end'],
        'index': k
      })

      # Valid solution
      k += 1

    # Sort by score
    resultList.sort(key=lambda x: x['score'], reverse=True)
    if not resultList:
      return None

    # Generate the state list
    state_vector = []
    solution_index = resultList[0]['index']
    for i in range(len(self._solutionList[solution_index])):
      for state in self._solutionList[solution_index][i]['vector']:
        state_vector.append(((state['hmm'], state['type']), state['idx']))
        
    return {
      'state': state_vector,
      'start': resultList[0]['start'],
      'end': resultList[0]['end'],
      'score': resultList[0]['score']
    }
