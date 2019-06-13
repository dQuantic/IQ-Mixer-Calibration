#sys.path.append("/Applications/Labber/Script")
#sys.path.append("C:\\Program Files (x86)\\Labber\\Script")

#Try to charge a Measurement template

# set output path and filenames
#now = datetime.datetime.now()
#sScriptsPath = os.path.dirname(os.path.abspath(__file__))
#sParentFolderPath = os.path.split(os.path.split(sScriptsPath)[0])[0] # ugh, this is ugly
#sOutPath = sParentFolderPath+'\\{:d}\\{:02d}\\Data_{:02d}{:02d}\\'.format(now.year,now.month,now.month,now.day)

#Meas = ScriptTools.MeasurementObject(os.path.join(sPath, 'Mixer_Char_Template.hdf5'),
                                     #os.path.join(sPath, 'IQ_Out.hdf5'))

                                     
import Labber
import time, numpy as np
from Labber import ScriptTools

#We connect to the local host
client=Labber.connectToServer('localhost')

#We read which instruments are currently connected to the server
instruments = client.getListOfInstrumentsString()
for inst in instruments:
    print (inst)

#We define the instruments we will be using for the IQ Mixer calibration
#We can add options such as set congig/get config at startup
Keysight_PXI_AWG=client.connectToInstrument('Keysight PXI AWG', dict(interface='PXI', address= 3,startup='Set config', lock=True))
SingleQubitPulse=client.connectToInstrument('Single-Qubit Pulse Generator',dict(name='Single-Qubit Pulse Generator' ))
RF_CW=client.connectToInstrument('Rohde&Schwarz RF Source', dict(interface='TCPIP', address='169.254.2.20'))
RF_IQ=client.connectToInstrument('Rohde&Schwarz RF Source', dict(interface='TCPIP', address='169.254.111.73'))
RF_Spectrometer=client.connectToInstrument('Rohde&Schwarz Spectrum Analyzer', dict(interface='TCPIP', address='169.254.2.21'))

#We start the defined instruments!!
Keysight_PXI_AWG.startInstrument()
RF_CW.startInstrument()
RF_IQ.startInstrument()
RF_Spectrometer.startInstrument()
SingleQubitPulse.startInstrument()


# set output path and filenames
#now = datetime.datetime.now()
#sScriptsPath = os.path.dirname(os.path.abspath(__file__))
#sParentFolderPath = os.path.split(os.path.split(sScriptsPath)[0])[0] # ugh, this is ugly
#sOutPath = sParentFolderPath+'\\{:d}\\{:02d}\\Data_{:02d}{:02d}\\'.format(now.year,now.month,now.month,now.day)

#Meas = ScriptTools.MeasurementObject(os.path.join(sPath, 'Mixer_Char_Template.hdf5'),
                                     #os.path.join(sPath, 'IQ_Out.hdf5'))


#C:\Users\penlabs.icn2\Labber\Data\2019\06\Data_0606\IQ_Mixer_1749(2)_6_6_LO_template.hdf5



#Define the Optimizer method
def nelder_mead(fun, x0,
                initial_step=0.1,
                no_improve_thr=10e-6, no_improv_break=10,
                maxiter=0, min_thr=-np.inf,
                alpha=1., gamma=2., rho=-0.5, sigma=0.5,
                verbose=False,args = ()):
    '''
    parameters:
        fun (function): function to optimize, must return a scalar score
            and operate over a numpy array of the same dimensions as x0
        x0 (numpy array): initial position
        initial_step (float/np array): determines the stepsize to construct
            the initial simplex. If a float is specified it uses the same
            value for all parameters, if an array is specified it uses
            the specified step for each parameter.
        no_improv_thr,  no_improv_break (float, int): break after
            no_improv_break iterations with an improvement lower than
            no_improv_thr
        maxiter (int): always break after this number of iterations.
            Set it to 0 to loop indefinitely.
        alpha (float): reflection coefficient
        gamma (float): expansion coefficient
        rho (float): contraction coefficient
        sigma (float): shrink coefficient
            For details on these parameters see Wikipedia page
    return: tuple (best parameter array, best score)
    Pure Python/Numpy implementation of the Nelder-Mead algorithm.
    Reference: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
    '''
    # init
    x0 = np.array(x0)  # ensures algorithm also accepts lists
    dim = len(x0)
    prev_best = fun(np.append(x0,args))
    res_list = []
    no_improv = 0
    res = [[x0, prev_best]]
    if type(initial_step) is float:
        initial_step_matrix = np.eye(dim)*initial_step
    elif (type(initial_step) is list) or (type(initial_step) is np.ndarray):
        if len(initial_step) != dim:
            raise ValueError('initial_step array must be same lenght as x0')
        initial_step_matrix = np.diag(initial_step)
    else:
        raise TypeError('initial_step ({})must be list or np.array'.format(
                        type(initial_step)))

    for i in range(dim):
        x = copy.copy(x0)
        x = x + initial_step_matrix[i]
        score = fun(np.append(x,args))
        res.append([x, score])

    # simplex iter
    iters = 0
    while 1:
        # order
        res.sort(key=lambda x: x[1])
        best = res[0][1]
        res_list.append(res[0])

        if best < min_thr:
            # Conclude success, break the loop
            if verbose:
                print('The minimum threshold reached after {} rounds,'.format(
                      iters) + 'concluding succesful convergence')
            break

        # break after maxiter
        if maxiter and iters >= maxiter:
            # Conclude failure break the loop
            if verbose:
                print('max iterations exceeded, optimization failed')
            break
        iters += 1

        if best < prev_best - no_improve_thr:
            no_improv = 0
            prev_best = best
        else:
            no_improv += 1

        if no_improv >= no_improv_break:
            # Conclude success, break the loop
            if verbose:
                print('No improvement registered for {} rounds, '.format(
                      no_improv_break) + 'concluding succesful convergence')
            break


        # centroid
        x0 = [0.] * dim
        for tup in res[:-1]:
            for i, c in enumerate(tup[0]):
                x0[i] += c / (len(res)-1)

        # reflection
        xr = x0 + alpha*(x0 - res[-1][0])
        rscore = fun(np.append(xr,args))
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]
            res.append([xr, rscore])
            continue

        # expansion
        if rscore < res[0][1]:
            xe = x0 + gamma*(x0 - res[-1][0])
            escore = fun(np.append(xe,args))
            if escore < rscore:
                del res[-1]
                res.append([xe, escore])
                continue
            else:
                del res[-1]
                res.append([xr, rscore])
                continue

        # contraction
        xc = x0 + rho*(x0 - res[-1][0])
        cscore = fun(np.append(xc,args))
        if cscore < res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # reduction
        x1 = res[0][0]
        nres = []
        for tup in res:
            redx = x1 + sigma*(tup[0] - x1)
            score = fun(np.append(redx,args))
            nres.append([redx, score])
        res = nres

    # once the loop is broken evaluate the final value one more time as
    # verification
    fun(np.append(res[0][0],args))
    return res[0], res_list