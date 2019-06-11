#sys.path.append("/Applications/Labber/Script")
#sys.path.append("C:\\Program Files (x86)\\Labber\\Script")

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
