from __future__ import print_function

import argparse
import logging
import os
import shutil
import sys
import time

def main():
    """
    This will try to run transgene from system arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=False, default=1)
    # parser.add_argument('--n_iter', dest='n_iter', type=int, help='Total number of iterations.', required=False, default=10)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    params = parser.parse_args()

    # Now that all the processes have completed, we need to process the worker specific vcfs
    worker_files = {i: '.worker_%s.txt' % i for i in range(0, params.cores)}
    with open(params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(worker_files[filename]) as w_file:
                first_line = True
                for line in w_file:
                    # Only print the first line header once
                    if !first_line || filename != ".worker_0.txt"
                        print(line, file=outfile, end='')
                    if first_line:
                        first_line = False

    for filename in worker_files:
        # Delete the temp files now that we're done with them.
        os.remove(worker_files[filename])

if __name__ == '__main__':
    main()
