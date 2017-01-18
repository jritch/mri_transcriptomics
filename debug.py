import read_data
import scipy
import pickle

def save_pickle(item):
  file_name = 'data.pickle'
  with open('data.pickle', 'wb') as f:
    pickle.dump({'item':item}, f)
  return file_name

def load_pickle(file_name):
    item = pickle.load(file_name)['item']
    return item

def test_pickle(item):
    f = save_pickle(item)
    item2 = load_pickle(item)
    print item2
    return

def print_out_hierarchy(coords,coord_to_region_map,ontology,flat_mri_data):
    for ID,name in ontology.names.iteritems():
        coords_list = read_data.get_coords_from_region_id(ID,coords,coord_to_region_map,ontology)
        #print name, ID, coords_list
        #print coords_list

    '''
      data = (coords,coord_to_region_map,o,flat_t1t2_ratio_data)

     
    '''
    return

def debug_region_id(coords,coord_to_region_map,ontology,flat_mri_data):
    
    # the length of the region IDs list from the entire brain is the same as the total number of region IDs
    print ontology.get_all_regionIDs(12959)
    return

def debug_ranking(coords,coord_to_region_map,o,flat_t1t2_ratio_data):
  #delete the file if it already exists
  os.unlink('C:/Users/Jacob/Google Drive/4th Year/Thesis/regions_ranked_by_t1overt2_ratio.txt')
  with open('C:/Users/Jacob/Google Drive/4th Year/Thesis/regions_ranked_by_t1overt2_ratio.txt','w') as f:    
    #debug.print_out_hierarchy(data[0],data[1],data[2],data[3])
    debug.debug_region_id(coords,coord_to_region_map,o,flat_t1t2_ratio_data)
    t2 = time.clock()
    print t2-t1
    f.write(pprint.pformat(rank_regions_by_intensity(coords,coord_to_region_map,o,flat_t1t2_ratio_data)))

    pass

import math


def main():
    test_pickle([1,2,3])
    return

if __name__ == '__main__':
    main()