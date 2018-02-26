from __future__ import print_function

import argparse
import logging
import os
import shutil
import subprocess
import sys
import time
from random import randint
from functools import partial
from multiprocessing import Manager, Pool


def rngenerator(worker_id, worker_iter, mh, mv, kmax, gamma, measure, n_tunnel_max, lambdascreen, rho_Lambda_thres, fixQ, Nafter, seed, n_recycle, output_dir):
    """
    Generate a random ross n g

    :param int worker_id: The id for this worker
    :param int worker_iter:
    :param float mh:
    :param float mv:
    :param int kmax:
    :param float gamma:
    :param str measure:
    :param int n_tunnel_max:
    :param bool lambdascreen:
    :param float rho_Lambda_thres:
    :param bool fixQ:
    :param float Nafter:
    """
    worker_id = str(worker_id)
    logging.info('worker %s is up and running.' % worker_id)
    worker_outfile = output_dir + ('.worker_%s.txt' % worker_id)

    #call = ['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-nojvm', '-minimize', '-r', '\"test(\'', worker_outfile, '\'),exit\"']
    call = ['sh','/hb/software/apps/matlab/bin/matlab', '-nodisplay','-nosplash','-wait','-nodesktop','-nojvm','-minimize','-r',
        'warning(\'off\');eis_wrapper(\'' +
        worker_outfile + '\',' +
        worker_iter + ',' +
        mv + ',' +
        mh + ',' +
        kmax + ',' +
        gamma + ',' +
        '\'' + measure + '\',' +
        n_tunnel_max + ',' +
        lambdascreen + ',' +
        rho_Lambda_thres + ',' +
        fixQ + ',' +
        Nafter + ',' +
        str(randint(1,1000)) + ',' +
        n_recycle +
        ');exit']
    subprocess.check_call(call)

    logging.info('worker %s received signal to go down.' % worker_id)


def main():
    """
    This will try to run transgene from system arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--input_file', dest='input_file', type=str, help='Input filename.', required=True)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=False, default=1)
    # parser.add_argument('--n_iter', dest='n_iter', type=int, help='Total number of iterations.', required=False, default=10)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    parser.add_argument('--output_dir', dest='output_dir', type=str, help='Output directory.', required=False, default='')
    params = parser.parse_args()

    # Default parameter values
    n_iter       = 1000
    mv           = 1.0
    mh           = 1.0
    kmax         = 30
    gamma        = 0.0
    measure      = 'B'
    n_tunnel_max = 3
    lambdascreen = True
    rho_Lambda_thres = 0.0000001
    fixQ         = False
    Nafter       = 60
    seed         = randint(1,1000);
    n_recycle    = 4;

    # Gather parameters from input file
    file_object = open(params.input_file)
    file_contents = file_object.readline().split(",")
    if len(file_contents) >= 1:
        n_iter       = int(  file_contents[0])
    if len(file_contents) >= 2:
        mv           = float(file_contents[1])
    if len(file_contents) >= 3:
        mh           = float(file_contents[2])
    if len(file_contents) >= 4:
        kmax         = int(  file_contents[3])
    if len(file_contents) >= 5:
        gamma        = float(file_contents[4])
    if len(file_contents) >= 6:
        measure      = str(  file_contents[5])
    if len(file_contents) >= 7:
        n_tunnel_max = int(  file_contents[6])
    if len(file_contents) >= 8:
        lambdascreen = int( file_contents[7])
    if len(file_contents) >= 9:
        rho_Lambda_thres = float( file_contents[8])
    if len(file_contents) >= 10:
        fixQ         = int( file_contents[9])
    if len(file_contents) >= 11:
        Nafter       = float(file_contents[10])
    if len(file_contents) >= 12:
        seed         = int(  file_contents[11])
    if len(file_contents) >= 13:
        n_recycle    = int(  file_contents[12])
    print(n_iter,mv,mh,kmax,gamma,measure,n_tunnel_max,
        lambdascreen,rho_Lambda_thres,fixQ,Nafter,seed,n_recycle)

    pool = Pool(processes=params.cores)
    rngenerator_partial = partial(rngenerator,
        worker_iter=str(n_iter/params.cores),
        mh=str(mh),
        mv=str(mv),
        kmax=str(kmax),
        gamma=str(gamma),
        measure=str(measure),
        n_tunnel_max=str(n_tunnel_max),
        lambdascreen=str(lambdascreen),
        rho_Lambda_thres=str(rho_Lambda_thres),
        fixQ=str(fixQ),
        Nafter=str(Nafter),
        seed=str(seed),
        n_recycle=str(n_recycle),
        output_dir=params.output_dir)
    pool.map(rngenerator_partial, range(0, params.cores))
    pool.close()
    pool.join()

    os.system("echo \"Combining output...\"")
    # time.sleep(10)

    # Now that all the processes have completed, we need to process the worker specific vcfs
    worker_files = {i: '.worker_%s.txt' % i for i in range(0, params.cores)}
    with open(params.output_dir + params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(params.output_dir + worker_files[filename]) as w_file:
                header_line = w_file.readline()
                if filename == 0:
                    print(header_line, file=outfile, end='')
                for line in w_file:
                    print(line, file=outfile, end='')
    for filename in worker_files:
        # Delete the temp files now that we're done with them.
        os.remove(params.output_dir + worker_files[filename])

if __name__ == '__main__':
    main()
