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
import platform


def rngenerator(worker_id, worker_iter, mh, mv, kmax, gamma, measure, n_tunnel_max, lambdascreen, rho_Lambda_thres, fixQ, Nafter, seed, n_recycle):
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
    worker_outfile = '.worker_%s.txt' % worker_id
# /hb/software/apps/matlab/bin/
    #call = ['matlab', '-nodisplay', '-nosplash', '-nodesktop', '-nojvm', '-minimize', '-r', '\"test(\'', worker_outfile, '\'),exit\"']
    if platform.system() != "Windows":
        call = ['sh','/hb/software/apps/matlab/bin/matlab', '-nodisplay','-nosplash','-wait','-nodesktop','-nojvm','-minimize','-r',
            '\"eis_wrapper(\'' +
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
            ');exit\"']
    else:
        call = ['matlab', '-nodisplay','-nosplash','-wait','-nodesktop','-nojvm','-minimize','-r',
            '\"eis_wrapper(\'' +
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
            ');exit\"']
    subprocess.check_call(call)

def main():
    """
    This will try to run transgene from system arguments
    """
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--task_id',dest='task_id',type=int,help='Task ID',required=True)
    parser.add_argument('--input_file', dest='input_file', type=str, help='Input filename.', required=True)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=False, default=1)
    # parser.add_argument('--n_iter', dest='n_iter', type=int, help='Total number of iterations.', required=False, default=10)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
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
    print(n_iter,mv,mh,kmax,gamma,measure,n_tunnel_max,lambdascreen,rho_Lambda_thres,fixQ,Nafter,seed,n_recycle)

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
                                  n_recycle=str(n_recycle))
    rngenerator_partial(params.task_id)

if __name__ == '__main__':
    main()
