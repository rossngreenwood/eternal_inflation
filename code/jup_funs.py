import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sys import stdout
from os import name as osname

def getBinnedFractions(run_ids,binParam,bins):

    if osname == 'nt':
        data_dir = 'C:/Users/Ross/Documents/data/'
    else:
        data_dir = '~/eternal_inflation/data/'

    data_names = ['record_flag','status','N','rho_offset','flag_fv_eternal', \
                'log_tunnel_rate', 'flag_hawking_moss','rho_false', \
                'numStochEpochs', 'NSinceStoch','numTopolEpochs','Q','r', \
                'n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure', \
                'n_tunnel_max','lambdascreen', 'rho_Lambda_thres','fixQ', \
                'Nafter','seed','n_recycle'];
    frac_names = ['mv','mh','len','n_iter','success','Q','r','ns','alpha', \
                'rho_Lambda','lgOk', 'stoch1','stoch2','stochAtExit','fv', \
                'fv_hm','top','wildcard']

    n_run = len(run_ids)

    print('Loading data...         |')
    count = 0

    fractions = pd.DataFrame(np.zeros([len(bins)+1,len(frac_names)]), columns=frac_names)

    for i in range(0,n_run):

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_t_' + ('%04d' % run_ids[i]) + '.txt'

        try:
            meta = pd.read_csv(fname,header=None,names=meta_names,nrows=1)
        except:
            continue

        if meta.shape[0] == 0:
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            stdout.write('#')

        data = pd.read_csv(fname, skiprows=2,header=None,names=data_names)

        if data.shape[0] == 0:
            continue

        if not binParam in ['mv','mh']:
            ibin = np.digitize(data[binParam],bins)
        else:
            ibin = np.digitize(meta[binParam][0],bins)*np.ones((len(data),1))

        for ind in range(0,len(bins)+1):

            ilog = (ibin == ind)
            if not any(ilog):
                continue

            fractions['len'][ind]         += data[ilog].shape[0]

            # Fractions for observables
            fractions['Q'][ind]           += sum(np.logical_and(data[ilog]['Q'] > np.sqrt(np.exp(3.089-0.036)/(1e10)), \
                                              data[ilog]['Q'] < np.sqrt(np.exp(3.089+0.036)/(1e10))))
            fractions['r'][ind]           += sum(data[ilog]['r'] < 0.114)
            fractions['ns'][ind]          += sum(np.logical_and(data[ilog]['n_s'] > 0.9655-0.0062, \
                                              data[ilog]['n_s'] < 0.9655+0.0062))
            fractions['alpha'][ind]       += sum(np.logical_and(data[ilog]['alpha'] > -0.0057-0.0071, \
                                                                data[ilog]['alpha'] <  0.0084+0.0082))
            fractions['rho_Lambda'][ind]  += sum(data[ilog]['rho_Lambda'] == 0)
            fractions['lgOk'][ind]        += sum(data[ilog]['lgOk'] < -2)

            # Fractions for eternal inflation
            fractions['stoch1'][ind]      += sum(data[ilog]['numStochEpochs'] == 1)
            fractions['stoch2'][ind]      += sum(data[ilog]['numStochEpochs'] == 2)
            fractions['stochAtExit'][ind] += sum(data[ilog]['NSinceStoch'] == 0)
            fractions['fv'][ind]          += sum(data[ilog]['flag_fv_eternal'] > 0)
            fractions['top'][ind]         += sum(data[ilog]['numTopolEpochs'] > 0)
            fractions['fv_hm'][ind]       += sum(data[ilog]['flag_hawking_moss'] > 0)

    # Divide by total counts to obtain fractions
    for j in fractions.index:
        if fractions['len'][j] > 0:
            fractions['stoch1'][j]      /= fractions['len'][j]
            fractions['stoch2'][j]      /= fractions['len'][j]
            fractions['stochAtExit'][j] /= fractions['len'][j]
            fractions['top'][j]         /= fractions['len'][j]
            fractions['fv'][j]          /= fractions['len'][j]
            fractions['fv_hm'][j]       /= fractions['len'][j]
            fractions['wildcard'][j]    /= fractions['len'][j]
            fractions['Q'][j]           /= fractions['len'][j]
            fractions['r'][j]           /= fractions['len'][j]
            fractions['ns'][j]          /= fractions['len'][j]
            fractions['alpha'][j]       /= fractions['len'][j]
            fractions['rho_Lambda'][j]  /= fractions['len'][j]
            fractions['lgOk'][j]        /= fractions['len'][j]

    print(' Done')

    return fractions

def getBinnedFractions2D(run_ids,binParam,bins,weights=None,merit=True,massBounds=((-np.inf,np.inf),(-np.inf,np.inf)),cutParams=None,cutBounds=None):

    if osname == 'nt':
        data_dir = 'C:/Users/Ross/Documents/data/'
    else:
        data_dir = '~/eternal_inflation/data/'

    data_names = ['record_flag','status','N','rho_offset','flag_fv_eternal', \
                'log_tunnel_rate', 'flag_hawking_moss','rho_false', \
                'numStochEpochs', 'NSinceStoch','numTopolEpochs','Q','r', \
                'n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda','weight_E'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure', \
                'n_tunnel_max','lambdascreen', 'rho_Lambda_thres','fixQ', \
                'Nafter','seed','n_recycle'];
    frac_names = ['mv','mh','len','counts','n_iter','success','Q','r','ns','alpha', \
                'rho_Lambda','lgOk', 'stoch1','stoch2','stochAtExit','fv', \
                'fv_hm','top','wildcard','NStoch','stochNearMax','nsalpha']

    n_run = len(run_ids)

    #if merit:
    sumfun = np.sum
    #else:
    #    sumfun = np.mean
    
    if weights is None:
        w = [1]*n_run
    else:
        w = weights

    print('Loading data...         |')
    count = 0

    fractions = pd.DataFrame(np.zeros([np.prod([len(b)+1 for b in bins]),len(frac_names)]), columns=frac_names)

    maxN = np.zeros([fractions.shape[0],1])
    NStoch = np.zeros([fractions.shape[0],1])
    valid_runs = np.zeros([fractions.shape[0],1])

    for i in range(0,n_run):

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_t_' + ('%04d' % run_ids[i]) + '.txt'

        try:
            meta = pd.read_csv(fname,header=None,names=meta_names,nrows=1)
        except:
            continue

        if meta.shape[0] == 0:
            continue

        #print('%.2f' % meta['mh'][0])
        if meta['mv'][0] < massBounds[0][0] or meta['mv'][0] > massBounds[0][1]:
            print(meta['mv'][0])
            continue
        if meta['mh'][0] < massBounds[1][0] or meta['mh'][0] > massBounds[1][1]:
            print(meta['mh'][0])
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            if osname == 'nt':
                print('#',end='')
            else:
                print('%d' % ((count-1)*4))

        data = pd.read_csv(fname, skiprows=2,header=None,names=data_names)

        if data.shape[0] == 0:
            continue
        if not merit:
            w[i] *= meta['n_iter'][0]/data.shape[0]

        ibin = [0 for x in binParam]
        for iparam in range(0,len(binParam)):
            if not binParam[iparam] in ['mv','mh']:
                ibin[iparam] = np.digitize(data[binParam[iparam]],bins[iparam])
                ibin[iparam][np.isnan(data[binParam[iparam]])] = -1
            else:
                ibin[iparam] = np.digitize(meta[binParam[iparam]][0],bins[iparam])*np.ones((len(data),1)[0])

        icut = np.greater(data['record_flag'],-1)
        if cutParams is not None:
            for ic in range(len(cutParams)):
                p = data[cutParams[ic]]
                icut = np.logical_and.reduce((icut,p >= cutBounds[ic][0], p <= cutBounds[ic][1]))

        for ind0 in range(0,len(bins[0])+1):
            for ind1 in range(0,len(bins[1])+1):

                ilog = np.logical_and.reduce((ibin[0] == ind0, ibin[1] == ind1, icut))

                ind = ind0*(len(bins[1])+1)+ind1

                if not any(ilog):
                    continue

                fractions['counts'][ind]      += data[ilog].shape[0]
                fractions['len'][ind]         += w[i]*data[ilog].shape[0]

                # Fractions for observables
                fractions['Q'][ind]           += w[i]*sumfun(np.logical_and(data[ilog]['Q'] > np.sqrt(np.exp(3.089-0.036)/(1e10)), \
                                                  data[ilog]['Q'] < np.sqrt(np.exp(3.089+0.036)/(1e10))))
                #fractions['Q'][ind]           += w[i]*sumfun(np.logical_and(data[ilog]['Q'] > 3.16E-3,data[ilog]['Q'] < 7.94E-3))
                fractions['r'][ind]           += w[i]*sumfun(data[ilog]['r'] < 0.114)
                #fractions['nsalpha'][ind]          += w[i]*sumfun(np.logical_and.reduce((data[ilog]['n_s'] > 0.9655-0.0062, \
                #                                  data[ilog]['n_s'] < 0.9655+0.0062, \
                #                                  data[ilog]['alpha'] > -0.0057-0.0071, data[ilog]['alpha'] < 0.0084+0.0082)))
                fractions['ns'][ind]          += w[i]*sumfun(np.logical_and(data[ilog]['n_s'] > 0.9655-0.0062, \
                                                  data[ilog]['n_s'] < 0.9655+0.0062))
                fractions['alpha'][ind]       += w[i]*sumfun(np.logical_and(data[ilog]['alpha'] > -0.0057-0.0071, \
                                                                    data[ilog]['alpha'] <  0.0084+0.0082))
                fractions['rho_Lambda'][ind]  += w[i]*sumfun(data[ilog]['rho_Lambda'] == 0)
                fractions['lgOk'][ind]        += w[i]*sumfun(data[ilog]['lgOk'] < -2)

                # Fractions for eternal inflation
                fractions['stoch1'][ind]      += w[i]*sumfun(data[ilog]['numStochEpochs'] == 1)
                fractions['stoch2'][ind]      += w[i]*sumfun(data[ilog]['numStochEpochs'] == 2)
                fractions['stochAtExit'][ind] += w[i]*sumfun(data[ilog]['NSinceStoch'] == 0)
                fractions['fv'][ind]          += w[i]*sumfun(data[ilog]['flag_fv_eternal'] > 0)
                fractions['top'][ind]         += w[i]*sumfun(data[ilog]['numTopolEpochs'] > 0)
                fractions['fv_hm'][ind]       += w[i]*sumfun(data[ilog]['flag_hawking_moss'] > 0)
                fractions['wildcard'][ind]       += w[i]*sumfun(np.logical_or(data[ilog]['numTopolEpochs'] > 0,data[ilog]['numStochEpochs']>0))
                #ilog = np.logical_and(ilog,data['N']>0)
                ilog = np.logical_and(ilog,np.logical_not(np.isnan(data['NSinceStoch'])))
                newMax = np.max(data[ilog]['N'])
                if newMax > maxN[ind]:
                    maxN[ind] = newMax
                fractions['NStoch'][ind]      += w[i]*sumfun(data[ilog]['N'])
                NStoch[ind]      += w[i]*np.sum(ilog)
                fractions['stochNearMax'][ind] += w[i]*np.sum(data[ilog]['N']-data[ilog]['NSinceStoch']-55 > 1)

    # Divide by total counts to obtain fractions
    for j in fractions.index:
        nsamp = fractions['len'][j]
        if nsamp > 0:
            fractions['stoch1'][j]      /= nsamp
            fractions['stoch2'][j]      /= nsamp
            fractions['stochAtExit'][j] /= nsamp
            fractions['top'][j]         /= nsamp
            fractions['fv'][j]          /= nsamp
            fractions['fv_hm'][j]       /= nsamp
            fractions['wildcard'][j]    /= nsamp
            fractions['Q'][j]           /= nsamp
            fractions['r'][j]           /= nsamp
            fractions['ns'][j]          /= nsamp
            fractions['alpha'][j]       /= nsamp
            fractions['rho_Lambda'][j]  /= nsamp
            fractions['lgOk'][j]        /= nsamp
            fractions['nsalpha'][j]     /= nsamp
            fractions['NStoch'][j]      -= maxN[j]
            #fractions['NStoch'][j]      /= np.max([1,nsamp])
            fractions['NStoch'][j]      /= np.max([1,(NStoch[j]-1)])
            fractions['stochNearMax'][j]/= NStoch[j]

    print(' Done')

    return fractions

def getBin2D(run_ids,binParam,bins):

    if osname == 'nt':
        data_dir = 'C:/Users/Ross/Documents/data/'
    else:
        data_dir = '~/eternal_inflation/data/'

    data_names = ['record_flag','status','N','rho_offset','flag_fv_eternal', \
                'log_tunnel_rate', 'flag_hawking_moss','rho_false', \
                'numStochEpochs', 'NSinceStoch','numTopolEpochs','Q','r', \
                'n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure', \
                'n_tunnel_max','lambdascreen', 'rho_Lambda_thres','fixQ', \
                'Nafter','seed','n_recycle'];

    n_run = len(run_ids)

    print('Loading data...         |')
    count = 0

    dataBin = pd.DataFrame(np.zeros([0,len(data_names)]), columns=data_names)

    for i in range(0,n_run):

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_t_' + ('%04d' % run_ids[i]) + '.txt'

        try:
            meta = pd.read_csv(fname,header=None,names=meta_names,nrows=1)
        except:
            continue

        if meta.shape[0] == 0:
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            stdout.write('#')

        data = pd.read_csv(fname, skiprows=2,header=None,names=data_names)

        if data.shape[0] == 0:
            continue

        ibin = [0 for x in binParam]
        for iparam in range(0,len(binParam)):
            if not binParam[iparam] in ['mv','mh']:
                ibin[iparam] = np.digitize(data[binParam[iparam]],bins[iparam])
                ibin[iparam][np.isnan(data[binParam[iparam]])] = -1
            else:
                ibin[iparam] = np.digitize(meta[binParam[iparam]][0],bins[iparam])*np.ones((len(data),1)[0])
        print(np.sum(ibin[0]==1))
        for ind0 in range(1,2):
            for ind1 in range(1,2):

                ilog = np.logical_and(ibin[0] == ind0, ibin[1] == ind1)

                ind = ind0*(len(bins[1])+1)+ind1

                if not any(ilog):
                    continue

                dataBin = dataBin.append(data[ilog])

    print(' Done')

    return dataBin

def massBin2D(run_ids,cutParams=None,cutBounds=None):

    if osname == 'nt':
        data_dir = 'C:/Users/Ross/Documents/data/'
    else:
        data_dir = '~/eternal_inflation/data/'

    data_names = ['record_flag','status','N','rho_offset','flag_fv_eternal', \
                'log_tunnel_rate', 'flag_hawking_moss','rho_false', \
                'numStochEpochs', 'NSinceStoch','numTopolEpochs','Q','r', \
                'n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda','weight_E'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure', \
                'n_tunnel_max','lambdascreen', 'rho_Lambda_thres','fixQ', \
                'Nafter','seed','n_recycle'];
    frac_names = ['mv','mh','len','n_iter','success','Q','r','ns','alpha', \
                'rho_Lambda','lgOk', 'stoch1','stoch2','stochAtExit','fv', \
                'fv_hm','top','wildcard']

    n_run = len(run_ids)

    print('Loading data...         |')
    count = 0

    fractions = pd.DataFrame(np.zeros([n_run,len(frac_names)]), columns=frac_names)

    for i in range(0,n_run):

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_t_' + ('%04d' % run_ids[i]) + '.txt'

        try:
            meta = pd.read_csv(fname,header=None,names=meta_names,nrows=1)
        except:
            continue

        if meta.shape[0] == 0:
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            stdout.write('#')

        data = pd.read_csv(fname, skiprows=2,header=None,names=data_names)

        icut = np.greater(data['record_flag'],-1)
        if cutParams is not None:
            for ic in range(len(cutParams)):
                p = data[cutParams[ic]]
                icut = np.logical_and.reduce((icut,p >= cutBounds[ic][0], p <= cutBounds[ic][1]))
        
        ind = i
        fractions['mv'][ind] = meta['mv'][0]
        fractions['mh'][ind] = meta['mh'][0]

        fractions['n_iter'][ind]      += meta['n_iter'][0]

        if data[icut].shape[0] == 0:
            continue

        fractions['len'][ind]         += data[icut].shape[0]
        fractions['success'][ind]     += data[icut].shape[0]

        # Fractions for observables
        fractions['Q'][ind]           += sum(np.logical_and(data[icut]['Q'] > np.sqrt(np.exp(3.089-0.036)/(1e10)), \
                                          data[icut]['Q'] < np.sqrt(np.exp(3.089+0.036)/(1e10))))
        fractions['r'][ind]           += sum(data[icut]['r'] < 0.064)
        fractions['ns'][ind]          += sum(np.logical_and(data[icut]['n_s'] > 0.9655-0.0062, \
                                          data[icut]['n_s'] < 0.9655+0.0062))
        fractions['alpha'][ind]       += sum(np.logical_and(data[icut]['alpha'] > -0.0057-0.0071, \
                                                            data[icut]['alpha'] <  0.0084+0.0082))
        fractions['rho_Lambda'][ind]  += sum(data[icut]['rho_Lambda'] == 0)
        fractions['lgOk'][ind]        += sum(data[icut]['lgOk'] < -2)

        # Fractions for eternal inflation
        fractions['stoch1'][ind]      += sum(data[icut]['numStochEpochs'] == 1)
        fractions['stoch2'][ind]      += sum(data[icut]['numStochEpochs'] == 2)
        fractions['stochAtExit'][ind] += sum(data[icut]['NSinceStoch'] == 0)
        fractions['fv'][ind]          += sum(data[icut]['flag_fv_eternal'] > 0)
        fractions['top'][ind]         += sum(data[icut]['numTopolEpochs'] > 0)
        fractions['fv_hm'][ind]       += sum(data[icut]['flag_hawking_moss'] > 0)
        fractions['wildcard'][ind]    += sum(np.logical_or(data[icut]['numTopolEpochs'] > 0,data[icut]['numStochEpochs']>0))
    
    # Divide by total counts to obtain fractions
    print(fractions.shape)
    for j in fractions.index:
        if fractions['n_iter'][j] > 0:
            fractions['success'][j]     /= fractions['n_iter'][j]
        if fractions['len'][j] > 0:
            fractions['stoch1'][j]      /= fractions['len'][j]
            fractions['stoch2'][j]      /= fractions['len'][j]
            fractions['stochAtExit'][j] /= fractions['len'][j]
            fractions['fv'][j]          /= fractions['len'][j]
            fractions['fv_hm'][j]       /= fractions['len'][j]
            fractions['wildcard'][j]    /= fractions['len'][j]
            fractions['Q'][j]           /= fractions['len'][j]
            fractions['r'][j]           /= fractions['len'][j]
            fractions['ns'][j]          /= fractions['len'][j]
            fractions['top'][j]         /= fractions['len'][j]
            fractions['alpha'][j]       /= fractions['len'][j]
            fractions['rho_Lambda'][j]  /= fractions['len'][j]
            fractions['lgOk'][j]        /= fractions['len'][j]

    print(' Done')

    return fractions

def getMarginalFractions(run_ids,model=False,rho_thres=0.02):

    if osname == 'nt':
        data_dir = 'C:/Users/Ross/Documents/data/'
    else:
        data_dir = '~/eternal_inflation/data/'

    data_names = ['record_flag','status','N','rho_offset','flag_fv_eternal', \
                'log_tunnel_rate', 'flag_hawking_moss','rho_false', \
                'numStochEpochs', 'NSinceStoch','numTopolEpochs','Q','r', \
                'n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure', \
                'n_tunnel_max','lambdascreen', 'rho_Lambda_thres','fixQ', \
                'Nafter','seed','n_recycle'];
    frac_names = ['mv','mh','len','n_iter','success','Q','r','ns','alpha', \
                'rho_Lambda','lgOk', 'stoch1','stoch2','stochAtExit','fv', \
                'fv_hm','wildcard']

    n_run = len(run_ids)

    print('Loading data...         |')
    count = 0

    fractions = pd.DataFrame(np.zeros([n_run,len(frac_names)]), columns=frac_names)

    for i in range(0,n_run):

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_t_' + ('%04d' % run_ids[i]) + '.txt'

        try:
            meta = pd.read_csv(fname,header=None,names=meta_names,nrows=1)
        except:
            continue

        if meta.shape[0] == 0:
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            stdout.write('#')

        data = pd.read_csv(fname, skiprows=2,header=None,names=data_names)

        scale_match = fractions.query(('mv == %f' % meta['mv'][0]) + (' & mh == %f' % meta['mh'][0]));

        if scale_match.shape[0] > 0:
            ind = scale_match.index[0]
        else:
            ind = i
            fractions['mv'][ind] = np.round(meta['mv'][0],decimals=5)
            fractions['mh'][ind] = np.round(meta['mh'][0],decimals=5)

        if model:
            data = data.query('not(rho_Lambda > %G)' % (np.power(float(meta['mv'][0]),4)*rho_thres))
            data = data.query('not(rho_offset > %G)' % (np.power(float(meta['mv'][0]),4)*rho_thres))
            data = data.query('not(rho_offset <= 0)')
            data = data.query('not(N < 0)')
        else:
            data = data.query('Q > %f' % np.sqrt(np.exp(3.089-0.036)/(1e10)))
            data = data.query('Q < %f' % np.sqrt(np.exp(3.089+0.036)/(1e10)))
            data = data.query('r < 0.114')
            data = data.query('n_s > %f' % (0.9655-0.0062))
            data = data.query('n_s < %f' % (0.9655+0.0062))

        fractions['n_iter'][ind]      += meta['n_iter'][0]

        if data.shape[0] == 0:
            continue

        fractions['len'][ind]         += data.shape[0]
        fractions['success'][ind]     += data.shape[0]

        # Fractions for observables
        fractions['Q'][ind]           += sum(np.logical_and(data['Q'] > np.sqrt(np.exp(3.089-0.036)/(1e10)), \
                                          data['Q'] < np.sqrt(np.exp(3.089+0.036)/(1e10))))
        fractions['r'][ind]           += sum(data['r'] < 0.114)
        fractions['ns'][ind]          += sum(np.logical_and(data['n_s'] > 0.9655-0.0062, \
                                          data['n_s'] < 0.9655+0.0062))
        fractions['alpha'][ind]       += sum(np.logical_and(data['alpha'] > -0.0057-0.0071, \
                                                            data['alpha'] <  0.0084+0.0082))
        fractions['rho_Lambda'][ind]  += sum(data['rho_Lambda'] == 0)
        fractions['lgOk'][ind]        += sum(data['lgOk'] < -2)

        # Fractions for eternal inflation
        fractions['stoch1'][ind]      += sum(data['numStochEpochs'] == 1)
        fractions['stoch2'][ind]      += sum(data['numStochEpochs'] == 2)
        fractions['stochAtExit'][ind] += sum(data['NSinceStoch'] == 0)
        fractions['fv'][ind]          += sum(data['flag_fv_eternal'] > 0)
        fractions['top'][ind]         += sum(data['numTopolEpochs'] > 1)
        fractions['fv_hm'][ind]       += sum(data['flag_hawking_moss'] > 0)

    df = fractions[fractions.mv != 0] # Remove would-be duplicate entries

    # Divide by total counts to obtain fractions
    print(df.shape)
    for j in df.index:
        if df['n_iter'][j] > 0:
            df['success'][j]     /= df['n_iter'][j]
        if df['len'][j] > 0:
            df['stoch1'][j]      /= df['len'][j]
            df['stoch2'][j]      /= df['len'][j]
            df['stochAtExit'][j] /= df['len'][j]
            df['fv'][j]          /= df['len'][j]
            df['fv_hm'][j]       /= df['len'][j]
            df['wildcard'][j]    /= df['len'][j]
            df['Q'][j]           /= df['len'][j]
            df['r'][j]           /= df['len'][j]
            df['ns'][j]          /= df['len'][j]
            df['alpha'][j]       /= df['len'][j]
            df['rho_Lambda'][j]  /= df['len'][j]
            df['lgOk'][j]        /= df['len'][j]

    print(' Done')

    return df

def getModelFractions_old(run_ids,new_output=True,rho_thres=0.02):

    if osname == 'nt':
        data_dir = 'C:/Users/Ross/Documents/data/'
    else:
        data_dir = '~/eternal_inflation/data/'

    col_names = ['record_flag','status','N','rho_offset','flag_fv_eternal','log_tunnel_rate', \
                 'flag_hawking_moss','rho_false','numStochEpochs', \
                 'NSinceStoch','numTopolEpochs','Q','r','n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure','n_tunnel_max','lambdascreen', \
                  'rho_Lambda_thres','fixQ','Nafter','seed','n_recycle'];

    n_run = len(run_ids)

    print('Loading data...         |')
    count = 0

    fractions = pd.DataFrame(np.zeros([n_run,16]), \
                             columns=['mv','mh','success','Q','r','ns','alpha','rho_Lambda','lgOk', \
                                      'stoch1','stoch2','stochAtExit','fv','fv_hm','top','wildcard'])

    for i in range(0,n_run):

        # print('outfile_t_' + ('%04d' % run_ids[i]) + '.txt')

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_t_' + ('%04d' % run_ids[i]) + '.txt'

        try:
            meta = pd.read_csv(fname,header=None,names=meta_names,nrows=1)
        except:
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            stdout.write('#')

        if meta.shape[0] == 0:
            continue

        # Create grid of mass scales

        data = pd.read_csv(fname, skiprows=2,header=None,names=data_names)

        fractions['mv'][i] = meta['mv'][0]
        fractions['mh'][i] = meta['mh'][0]

        data = data.query('not(rho_Lambda > %G)' % (np.power(float(meta['mv'][0]),4)*rho_thres))
        data = data.query('not(rho_offset > %G)' % (np.power(float(meta['mv'][0]),4)*rho_thres))
        data = data.query('not(rho_offset <= 0)')
        data = data.query('not(N < 0)')

        if data.shape[0] == 0:
            continue

        fractions['success'][i]     = data.shape[0]/meta['n_iter'][0]

        fractions['Q'][i]           = sum(np.logical_and(data['Q'] > np.sqrt(np.exp(3.089-0.036)/(1e10)), \
                                          data['Q'] < np.sqrt(np.exp(3.089+0.036)/(1e10))))/data.shape[0]
        fractions['r'][i]           = sum(data['r'] < 0.114)/data.shape[0]
        fractions['ns'][i]          = sum(np.logical_and(data['n_s'] > 0.9655-0.0062, \
                                          data['n_s'] < 0.9655+0.0062))/data.shape[0]
        fractions['alpha'][i]       = sum(np.logical_and(data['alpha'] > -0.0057-0.0071, \
                                          data['alpha'] <  0.0084+0.0082))/data.shape[0]
        fractions['rho_Lambda'][i]  = sum(data['rho_Lambda'] == 0)/data.shape[0]
        fractions['lgOk'][i]        = sum(data['lgOk'] < -2)/data.shape[0]

        # Fractions for eternal inflation
        fractions['stoch1'][i]      = sum(np.logical_and(data['numStochEpochs'] == 1, \
                                           data['NSinceStoch'] != 0.314))/data.shape[0]
        fractions['stoch2'][i]      = sum(np.logical_and(data['numStochEpochs'] == 2, \
                                           data['NSinceStoch'] != 0.314))/data.shape[0]
        fractions['stochAtExit'][i] = sum(data['NSinceStoch'] == 0)/data.shape[0]
        # fractions['fv'][i]          = sum(data['flag_fv_eternal'] > 0)/data.shape[0]
        fractions['fv'][i]          = sum(np.logical_and(data['log_tunnel_rate'] > -10000000, \
                                        data['log_tunnel_rate'] < 9/4/3.14159))/data.shape[0]
        fractions['top'][i]         = sum(data['numTopolEpochs'] > 0)/data.shape[0]
        fractions['fv_hm'][i]       = sum(data['flag_hawking_moss'] > 0)/data.shape[0]
        fractions['wildcard'][i]    = sum(np.logical_and(data['numStochEpochs'] > 0, \
                                        data['flag_hawking_moss'] == 0))/data.shape[0]
    print(' Done')

    return fractions

def fracPlotObservables(fractions):

    dfrac = fractions.query('mv != 0')

    fig, axs = plt.subplots(nrows=1,ncols=2,sharey=False,figsize=(14,6))

    # null_inds = (fractions['mv'] != 0)

    plt.sca(axs[0])
    mv2_mh = np.power(dfrac['mv'],2)/np.power(dfrac['mh'],1)
    plt.xlabel('$m_v^2/m_h$')
    plt.plot(mv2_mh,dfrac['Q'],         'o',alpha=0.3)
    plt.plot(mv2_mh,dfrac['r'],         'o',alpha=0.1)
    plt.plot(mv2_mh,dfrac['ns'],        'o',alpha=0.1)
    plt.plot(mv2_mh,dfrac['success'],   'd',alpha=0.1)
    # lt.plot(mv2_mh,dfrac['Q']*dfrac['r']*dfrac['success'],'kD',markerfacecolor='None')
    plt.legend(['$Q$','$r$','$n_s$','$N_e > 55$'],loc='lower left')
    plt.loglog()

    plt.ylabel('Fraction')

    plt.sca(axs[1])
    mv2_mh = dfrac['mh']
    plt.xlabel('$m_h$')
    plt.plot(mv2_mh,dfrac['Q'],         'o',alpha=0.1)
    plt.plot(mv2_mh,dfrac['r'],         'o',alpha=0.3)
    plt.plot(mv2_mh,dfrac['ns'],         'o',alpha=0.3)
    plt.plot(mv2_mh,dfrac['success'],   'd',alpha=0.3)
    # plt.plot(mv2_mh,dfrac['Q']*dfrac['r']*dfrac['success'],'kD',markerfacecolor='None')
    plt.legend(['$Q$','$r$','$n_s$','$N_e > 55$'],loc='lower left')
    plt.loglog()

    # plt.savefig(data_dir + '../pyplotfig.png')

def fracPlotStochastic(fractions):

    fig, axs = plt.subplots(nrows=1,ncols=1,sharey=False,figsize=(6,6))

    dfrac = fractions.query('mv != 0')

    plt.sca(axs)
    mv2_mh = np.power(dfrac['mv'],1)/np.power(dfrac['mh'],-1)
    plt.plot(mv2_mh,dfrac['stoch1']+dfrac['stoch2'],     's',alpha=0.4)
    plt.plot(mv2_mh,dfrac['top'],     'd',alpha=0.4)
    plt.plot(mv2_mh,dfrac['fv'],'o',alpha=0.4)
    plt.xlabel(u"$m_v m_h$")
    plt.legend(['Stochastic','Topological','False Vacuum'],loc='lower right')
    plt.loglog()

    plt.ylabel('Fraction')

    # plt.sca(axs[1])
    # mv2_mh = np.power(dfrac['mv'],1)/np.power(dfrac['mh'],-1)
    # # plt.plot(mv2_mh,dfrac['fv'],'x',alpha=0.4)
    # # plt.plot(mv2_mh,dfrac['fv_hm'],'+',alpha=0.4)
    # plt.plot(mv2_mh,dfrac['wildcard'],'+',alpha=0.4)
    # plt.xlabel(u"$m_v m_h$")
    # # plt.legend(['False vacuum','Hawking-Moss'],loc='upper left')
    # plt.loglog()

def fracPlotStochastics(fractions1,fractions2):

    fig, axs = plt.subplots(nrows=1,ncols=1,sharey=False,figsize=(6,6))

    dfrac1 = fractions1.query('mv != 0')
    dfrac2 = fractions2.query('mv != 0')

    plt.sca(axs)

    mv2_mh = np.power(dfrac1['mv'],2)/np.power(dfrac1['mh'],-1)
    plt.plot(mv2_mh,dfrac1['stoch1']+dfrac1['stoch2'],     'bo',alpha=0.4)
    # plt.plot(mv2_mh,dfrac1['stochAtExit'],     'bd',alpha=0.4)

    mv2_mh = np.power(dfrac2['mv'],2)/np.power(dfrac2['mh'],-1)
    plt.plot(mv2_mh,dfrac2['stoch1']+dfrac2['stoch2'],     'rs',alpha=0.4)
    # plt.plot(mv2_mh,dfrac2['stochAtExit'],     'rs',alpha=0.4)

    plt.xlabel(u"$m_v^2 m_h$")
    plt.legend([u"0.1 $m_v^4$ Allowance", u"0.1 $m_v^4$ Allowance", 'No Allowance', 'No Allowance'],loc='lower left')
    plt.loglog()

def loadData(run_ids,rho_thres=0.02):

    data_dir = 'C:/Users/Ross/Documents/data/';

    col_names = ['record_flag','status','N','rho_offset','flag_fv_eternal','log_tunnel_rate', \
                 'flag_hawking_moss','rho_false','numStochEpochs', \
                 'NSinceStoch','numTopolEpochs','Q','r','n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda'];
    meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure','n_tunnel_max','lambdascreen', \
                  'rho_Lambda_thres','fixQ','Nafter','seed','n_recycle'];

    n_run = len(run_ids)

    print('Loading data...         |')
    count = 0

    for i in range(0,n_run):

        try:
            meta = pd.read_csv(data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt', \
                                 header=None,names=meta_names,nrows=1)
        except:
            continue

        if i == 0 and meta.shape[0] > 0:
            print('Measure: %s' % meta['measure'][0])

        if np.floor(i/(n_run/25)) >= count:
            count = count + 1
            stdout.write('#')

        if meta.shape[0] == 0:
            continue

        data = pd.read_csv(data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt', \
                         skiprows=2,header=None,names=col_names)

        if data.shape[0] == 0:
            continue

        if i == 0:
            dataAll = data
        else:
            dataAll = dataAll.append(data)

    print(' Done')

    return dataAll


# col_names = ['record_flag','status','N','rho_offset','flag_fv_eternal','log_tunnel_rate','flag_hm','V_false','numStochEpochs', \
#              'NSinceStoch','numTopolEpochs','Q','r','n_s','alpha','n_t','dlgrho','lgOk','rho_Lambda'];
# meta_names = ['n_iter','cores','mv','mh','m_Pl','kmax','gamma','measure','n_tunnel_max','lambdascreen', \
#               'rho_Lambda_thres','Nafter','seed','n_recycle'];
