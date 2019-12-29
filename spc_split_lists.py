"""Split glaciers into lists to run on separate nodes on the supercomputer"""

# Built-in libraries
import argparse
import os
# External libraries
import numpy as np
import pickle
# Local libraries
import globaldebris_input as input


def getparser():
    """
    Use argparse to add arguments from the command line
    
    Parameters
    ----------
    n_batches (optional) : int
        number of nodes being used on the supercomputer
    ignore_regionname (optional) : int
        switch to ignore region name or not (1 ignore it, 0 use region)
    add_cal : int
        switch to add "Cal" to the batch filenames such that calibration and simulation can be run at same time
    option_ordered : int
        option to keep glaciers ordered or to grab every n value for the batch
        (the latter helps make sure run times on each core are similar as it removes any timing differences caused by 
         regional variations)
        
    Returns
    -------
    Object containing arguments and their respective values.
    """
    parser = argparse.ArgumentParser(description="run calibration in parallel")
    # add arguments
    parser.add_argument('-n_batches', action='store', type=int, default=1,
                        help='number of nodes to split the glaciers amongst')
    parser.add_argument('-ignore_regionname', action='store', type=int, default=0,
                        help='switch to include the region name or not in the batch filenames')
    parser.add_argument('-option_ordered', action='store', type=int, default=1,
                        help='switch to keep lists ordered or not')
    return parser


def split_list(lst, n=1, option_ordered=1):
    """
    Split list into batches for the supercomputer.
    
    Parameters
    ----------
    lst : list
        List that you want to split into separate batches
    n : int
        Number of batches to split glaciers into.
    
    Returns
    -------
    lst_batches : list
        list of n lists that have sequential values in each list
    """
    # If batches is more than list, then there will be one glacier in each batch
    if option_ordered == 1:
        if n > len(lst):
            n = len(lst)
        n_perlist_low = int(len(lst)/n)
        n_perlist_high = int(np.ceil(len(lst)/n))
        lst_copy = lst.copy()
        count = 0
        lst_batches = []
        for x in np.arange(n):
            count += 1
            if count <= len(lst) % n:
                lst_subset = lst_copy[0:n_perlist_high]
                lst_batches.append(lst_subset)
                [lst_copy.remove(i) for i in lst_subset]
            else:
                lst_subset = lst_copy[0:n_perlist_low]
                lst_batches.append(lst_subset)
                [lst_copy.remove(i) for i in lst_subset]
    
    else:
        if n > len(lst):
            n = len(lst)
    
        lst_batches = [[] for x in np.arange(n)]
        nbatch = 0
        for count, x in enumerate(lst):
            if count%n == 0:
                nbatch = 0
    
            lst_batches[nbatch].append(x)
            
            nbatch += 1
    
            
    return lst_batches   


if __name__ == '__main__':
    parser = getparser()
    args = parser.parse_args()   
        
    # Count glaciers in existing batch
    batch_list = []
    count_latlons = 0
    batch_str = 'latlon_batch_'

    # check files
    for i in os.listdir():
            
        if args.ignore_regionname == 0:
            check_str = input.roi + '_' + batch_str
        elif args.ignore_regionname == 1:
            check_str = batch_str
        
        # List batch fns and count total number of glaciers
        if i.startswith(check_str) and i.endswith('.pkl'):
            with open(i, 'rb') as f:
                latlon_list = pickle.load(f)
                batch_list.append(i)
            
            count_latlons += len(latlon_list)
    
    #%%    
    # Check if need to update old batch files or not
    #  (different number of glaciers or batches)
    if count_latlons != len(input.latlon_list) or args.n_batches != len(batch_list):
        # Delete old files
        for i in batch_list:
            os.remove(i)
            
        # Split list of of lat/lons
        # Lat/lon lists to pass for parallel processing
        latlon_lsts = split_list(input.latlon_list, n=args.n_batches, option_ordered=args.option_ordered)
    
        # Export new lists
        for n in range(len(latlon_lsts)):
            if args.ignore_regionname == 0:
                batch_fn = input.roi + '_' + batch_str + str(n) + '.pkl'
            elif args.ignore_regionname == 1:
                batch_fn = batch_str + str(n) + '.pkl'
                
            print('Batch', n, ':\n', batch_fn, '\n')
            with open(batch_fn, 'wb') as f:
                pickle.dump(latlon_lsts[n], f)