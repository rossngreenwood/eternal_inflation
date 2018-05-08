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
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=True)
    # parser.add_argument('--n_iter', dest='n_iter', type=int, help='Total number of iterations.', required=False, default=10)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    parser.add_argument('--output_dir', dest='output_dir', type=str, help='Output directory.', required=False, default='')
    params = parser.parse_args()

    # Now that all the processes have completed, we need to process the worker specific vcfs
    worker_files = {i: '.worker_%s.txt' % i for i in range(0, params.cores)}
    with open(params.output_dir + params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(params.output_dir + worker_files[filename]) as w_file:
                header_line = w_file.readline()
		if filename == 0:
		    print(header_line, file=outfile, end='')
                for line in w_file:
                    if line[0] == "3" and line[1] == ",":
                        print(line, file=outfile, end='')

    #for filename in worker_files:
    #    # Delete the temp files now that we're done with them.
    #    os.remove(params.output_dir + worker_files[filename])

if __name__ == '__main__':
    main()
