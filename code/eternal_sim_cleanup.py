from __future__ import print_function

import argparse
import os
import sys
import numpy as np
from pyparsing import nestedExpr
from math import isnan

def main():

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--cores', dest='cores', type=int, help='Number of cores to use.', required=False,default=24)
    parser.add_argument('--output_file', dest='output_file', type=str, help='Output filename.', required=False, default='outfile.txt')
    parser.add_argument('--output_dir', dest='output_dir', type=str, help='Output directory.', required=False, default='')
    parser.add_argument('--truncate', dest='flag_truncate', type=int, help='Truncate output file?', required=False, default=0)
    parser.add_argument('--cut_params', dest='cutParams', type=str, help='File describing how to truncate data', required=False, default='')
    parser.add_argument('--cut_bounds', dest='cutBounds', type=str, help='File describing how to truncate data', required=False, default='')
    params = parser.parse_args()
    
    nested_braces = nestedExpr('[',']')
    cutParams = nested_braces.parseString(params.cutParams).asList()[0]
    cutBounds = nested_braces.parseString(params.cutBounds).asList()[0]
    cutBounds = [[float(j) for j in i] for i in cutBounds]
    
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
    
    lncount = 0
    with open(params.output_dir + params.output_file, 'w') as outfile:
        for filename in worker_files:
            with open(params.output_dir + worker_files[filename]) as w_file:
                header_line = w_file.readline()
                if filename == 0:
                    # Multiply the number of iterations for each worker by number of workers
                    header_list = header_line.split(',')
                    n_iter = str(float(header_list[0])*float(header_list[1]))[:-2]
                    # Write header line only once, followed by a line of seeds
                    print(n_iter + ',' + ','.join(header_list[1:13]+[header_list[13]]) + '\r\n', file=outfile, end='')
                for line in w_file:
                    # Enforce various means of truncating the data
                    if 'count' in params.cutParams:
                        ind = [i for i,x in enumerate(cutParams) if x == 'count'][0]
                        maxCount = cutBounds[ind][0]
                        if lncount+1 > maxCount:
                            break
                    data_list = line.split(',')
                    if 'record_flag' in cutParams:
                        ind = [i for i,x in enumerate(cutParams) if x == 'record_flag'][0]
                        if not float(data_list[0]) in cutBounds[ind]:
                            continue
                    if 'status' in cutParams:
                        ind = [i for i,x in enumerate(cutParams) if x == 'status'][0]
                        if not float(data_list[1]) in cutBounds[ind]:
                            continue
                    if data_list[0] in ['2','3']: 
                        if 'log_tunnel_rate' in cutParams:
                            log_tunnel_rate = float(data_list[5])
                            ind = [i for i,x in enumerate(cutParams) if x == 'log_tunnel_rate'][0]
                            if isnan(log_tunnel_rate):
                                continue
                        if 'flag_hm' in cutParams:
                            flag_hm = float(data_list[6])
                            ind = [i for i,x in enumerate(cutParams) if x == 'flag_hm'][0]
                            if not flag_hm in cutBounds[ind]:
                                continue
                        if 'flag_fv' in cutParams:
                            flag_fv = float(data_list[4])
                            ind = [i for i,x in enumerate(cutParams) if x == 'flag_fv'][0]
                            if not flag_fv in cutBounds[ind]:
                                continue
                    if data_list[0] == '3':
                        if 'Q' in cutParams:
                            Q = float(data_list[11])
                            ind = [i for i,x in enumerate(cutParams) if x == 'Q'][0]
                            bounds = cutBounds[ind]
                            #if not type(bounds[0]) is int:
                            #    Amin = 3.089-bounds[0]*0.036
                            #    Amax = 3.089+bounds[1]*0.036
                            #    if Q < np.sqrt(np.exp(Amin)/(1e10)):
                            #        continue
                            #    if Q > np.sqrt(np.exp(Amax)/(1e10)):
                            #        continue
                            #else:
                            if Q < bounds[0] or Q > bounds[1]:
                                continue
                        if 'n_s' in cutParams:
                            n_s = float(data_list[13])
                            ind = [i for i,x in enumerate(cutParams) if x == 'n_s'][0]
                            bounds = cutBounds[ind]
                            if not type(bounds[0]) is int:
                                nmin = 0.9655-bounds[0]*0.0062
                                nmax = 0.9655+bounds[0]*0.0062
                                if n_s < nmin or n_s > nmax:
                                    continue
                            else:
                                if n_s < bounds[0] or n_s > bounds[1]:
                                    continue
                        if 'r' in cutParams:
                            r = float(data_list[12])
                            ind = [i for i,x in enumerate(cutParams) if x == 'r'][0]
                            bounds = cutBounds[ind]
                            if r < bounds[0] or r > bounds[1]:
                                continue
                        if 'Lambda' in cutParams:
                            Lambda = float(data_list[18])
                            ind = [i for i,x in enumerate(cutParams) if x == 'Lambda'][0]
                            bounds = cutBounds[ind]
                            if Lambda < bounds[0] or Lambda > bounds[1]:
                                continue
                        if 'flag_top' in cutParams:
                            ind = [i for i,x in enumerate(cutParams) if x == 'flag_top'][0]
                            if not float(data_list[10]) in cutBounds[ind]:
                                continue
                    lncount += 1
                    print(line, file=outfile, end='')

    #for filename in worker_files:
    #    # Delete the temp files now that we're done with them.
    #    os.remove(params.output_dir + worker_files[filename])

if __name__ == '__main__':
    main()
