from __future__ import print_function

import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=True)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    parser.add_argument('--output_dir', dest='output_dir', type=str, help='Output directory.', required=False, default='')
    params = parser.parse_args()

    worker_files = {i: '.worker_%s.txt' % i for i in range(0, params.cores)}

    worker_seeds = ""
    for filename in worker_files:
        with open(params.output_dir + worker_files[filename]) as w_file:
            w_header = w_file.readline()
            worker_seeds = worker_seeds + w_header.split(',')[-2] + ','
    worker_seeds = worker_seeds[:-1] # Remove last comma

    with open(params.output_dir + params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(params.output_dir + worker_files[filename]) as w_file:
                header_line = w_file.readline()
        		if filename == 0:
                    # Multiply the number of iterations for each worker by number of workers
                    (n_iter,header_line) = header_line.split(',',1)
                    n_iter = str(float(n_iter)*params.cores)[:-2]
                    # Write header line only once, followed by a line of seeds
        		    print(n_iter + ',' + header_line, file=outfile, end='')
                    print(worker_seeds, file=outfile, end='')
                for line in w_file:
                    print(line, file=outfile, end='')

    #for filename in worker_files:
    #    # Delete the temp files now that we're done with them.
    #    os.remove(params.output_dir + worker_files[filename])

if __name__ == '__main__':
    main()
