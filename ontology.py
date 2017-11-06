import os
import csv
class Ontology(object):
  """An Ontology has two member variables:

  Attributes:
    names (dict): 
      An (ID,name) dictionary 
    hierarchy (dict):
      a (parentID,[childID]) dictionary 
    parents
      a child to parent dictionary
  """
  def __init__(self,ontology_file_name):
      names = dict()
      parents = dict()
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
            
            parents[ID] = parentID

      self.names = names
      self.hierarchy = hierarchy
      self.parents = parents

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

  def get_enclosing_regions(self,regionID):
    IDs = [regionID]
    all_IDs = []
    while IDs:
      # print 'IDs is',IDs
      ID = IDs.pop()
      all_IDs.append(ID)
      parent_ID = self.parents.get(ID)
      if parent_ID is not None:
        IDs.append(parent_ID)
    return all_IDs


def main():
    o = Ontology(os.path.join(config.scriptLocation, "data",  "Ontology.csv"))
    print(o.names[4008])
    
    print(o.hierarchy.get(4008))
    print(o.get_enclosing_regions(4008))
    print(o.get_enclosing_regions(4896))
    print([o.names[x] for x in o.get_enclosing_regions(4896)])
    cortex_divisions = [str(x) for x in o.hierarchy.get(4008)]
    enclosing_regions = [str(x) for x in o.get_enclosing_regions(4896)]
    cortex_subdivision = set(enclosing_regions).intersection(cortex_divisions)
    print(cortex_subdivision)
    cortex_subdivision= next(iter(cortex_subdivision))
    print(cortex_subdivision)
    cortex_subdivision = o.names[int(cortex_subdivision)]
    print(cortex_subdivision)




if __name__ == '__main__':
    import config
    main()
