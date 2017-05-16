import csv
class Ontology(object):
  """An Ontology has two member variables:

  Attributes:
    names (dict): 
      An (ID,name) dictionary 
    hierarchy (dict):
      a (parentID,[childID]) dictionary 
  """
  def __init__(self,ontology_file_name):
      names = dict()
      hierarchy = dict()
      with open(ontology_file_name) as f:
        rdr = csv.reader(f)
        header = next(rdr)
        for row in rdr:
          ID = int(row[0])
          names[ID] = row[2]
          if row[3]:
            parentID = int(row[3])
            if hierarchy.get(parentID):
                hierarchy[parentID].append(ID)
            else:
              hierarchy[parentID] = [ID]

      self.names = names
      self.hierarchy = hierarchy

  def get_all_regionIDs(self,regionID):
    IDs = [regionID]
    all_IDs = []
    while IDs:
      # print 'IDs is',IDs
      ID = IDs.pop()
      all_IDs.append(ID)
      child_IDs = self.hierarchy.get(ID)
      if child_IDs:
        IDs += child_IDs
    return all_IDs
