from __future__ import print_function

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from functools import partial
from multiprocessing import Manager, Pool


def rngenerator(worker_id, worker_iter, mh, mv, kmax):
    """
    Generate a random ross n g

    :param int worker_id: The id for this worker
    :param int mh:
    :param int mv:
    :param int kmax:
    """
    worker_id = str(worker_id)
    logging.info('worker %s is up and running.' % worker_id)
    worker_outfile = '.worker_%s.txt' % worker_id

    #call = ['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-nojvm', '-minimize', '-r', '\"test(\'', worker_outfile, '\'),exit\"']
    call = ['matlab', '-nodisplay','-nosplash','-nodesktop','-nojvm','-minimize','-r',
        'eis_wrapper(\'' + worker_outfile + '\',' + worker_iter + ',' + mv + ',' + mh + ',' + kmax + ')']
    subprocess.check_call(call)

    logging.info('worker %s received signal to go down.' % worker_id)


def main():
    """
    This will try to run transgene from system arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--input_file', dest='input_file', type=str, help='Input filename.', required=True)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=False, default=1)
    parser.add_argument('--n_iter', dest='n_iter', type=int, help='Total number of iterations.', required=False, default=10)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    params = parser.parse_args()

    # Gather parameters from input file
    file_object = open(params.input_file)
    file_contents = file_object.readline().split(",")
    if len(file_contents) != 3:
        print("Failed to read file")
    else:
        mv = float(file_contents[0])
        mh = float(file_contents[1])
        kmax = int(file_contents[2])
        print(mv,mh,kmax)

    pool = Pool(processes=params.cores)
    rngenerator_partial = partial(rngenerator,
                                  worker_iter=str(params.n_iter/params.cores),
                                  mh=str(mh),
                                  mv=str(mv),
                                  kmax=str(kmax))
    pool.map(rngenerator_partial, range(0, params.cores))
    pool.close()
    pool.join()

    os.system("echo \"Waiting 10 seconds...\"")
    time.sleep(10)

    # Now that all the processes have completed, we need to process the worker specific vcfs
    worker_files = {i: '.worker_%s.txt' % i for i in range(0, params.cores)}
    with open(params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(worker_files[filename]) as w_file:
                for line in w_file:
                    print(line, file=outfile, end='')
    # for filename in worker_files:
    #     # Delete the temp files now that we're done with them.
    #     os.remove(worker_files[filename])

if __name__ == '__main__':
    main()
