from src.igmat import igmat

class TestAnnotation:

  def test_input(self):
    
    # The sequence that needs to be annotated
    sequence = "DVQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWYRQAPGKEREWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCGVRSIREYWGQGTQVTVSSHHHHHH"

    # Generate the list of results
    resultList = igmat.annotate(sequence, 'IMGT')
    if not resultList:
      raise Exception('No result found')

    # Iterate the list of results
    for result in resultList:

      # Print the resulting sequence
      print(result)

      # Print some details
      print(result.sequence)
      print(result.start)
      print(result.end)
      print(result.type)

      # Print the annotations
      for feature in result.annotations():
        print('Annotation {type}: {start}-{stop}'.format(
          type=feature['type'],
          start=feature['start'],
          stop=feature['stop']
          ))