from __future__ import print_function

import argparse
import os
import sys

def main():

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=False,default=24)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    parser.add_argument('--output_dir', dest='output_dir', type=str, help='Output directory.', required=False, default='')
    parser.add_argument('--truncate', dest='flag_truncate', type=int, help='Truncate output file?', required=False, default=0)
    params = parser.parse_args()
    
    with open(params.output_dir + '.worker_0.txt') as w_file:
        header_line = w_file.readline()
        header_list = header_line.split(',')
        cores = int(float(header_list[1]))

    worker_files = {i: '.worker_%s.txt' % i for i in range(0, cores)}

    worker_seeds = ""
    for filename in worker_files:
        with open(params.output_dir + worker_files[filename]) as w_file:
            w_header = w_file.readline()
            worker_seeds = worker_seeds + w_header.split(',')[-2] + ','
    worker_seeds = worker_seeds[:-1] + '\r\n' # Remove last comma

    with open(params.output_dir + params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(params.output_dir + worker_files[filename]) as w_file:
                header_line = w_file.readline()
                if filename == 0:
                    # Multiply the number of iterations for each worker by number of workers
                    header_list = header_line.split(',')
                    n_iter = str(float(header_list[0])*cores)[:-2]
                    # Write header line only once, followed by a line of seeds
                    print(n_iter + ',' + ','.join(header_list[1:13]+[header_list[13]]), file=outfile, end='')
                for line in w_file:
                    if params.flag_truncate == 0 or line[:2] == '3,':
                        print(line, file=outfile, end='')

    #for filename in worker_files:
    #    # Delete the temp files now that we're done with them.
    #    os.remove(params.output_dir + worker_files[filename])

if __name__ == '__main__':
    main()
