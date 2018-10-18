import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sys import stdout
from os import name as osname

def getBinnedFractions(run_ids,binParam,bins,rho_thres=0.02):

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
            fractions['top'][ind]         += sum(data[ilog]['numTopolEpochs'] > 1)
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

def massBin2D(run_ids):

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

    fractions = pd.DataFrame(np.zeros([n_run,len(frac_names)]), columns=frac_names)

    for i in range(0,n_run):

        if osname == 'nt':
            fname = data_dir + 'outfile_t_' + ('%04d' % run_ids[i]) + '.txt'
        else:
            fname = data_dir + 'out_' + ('%04d' % run_ids[i]) + '/outfile_Q_' + ('%04d' % run_ids[i]) + '.txt'

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

        ind = i
        fractions['mv'][ind] = np.round(meta['mv'][0],decimals=5)
        fractions['mh'][ind] = np.round(meta['mh'][0],decimals=5)

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
        fractions['top'][ind]         += sum(data['numTopolEpochs'] > 0)
        fractions['fv_hm'][ind]       += sum(data['flag_hawking_moss'] > 0)

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
