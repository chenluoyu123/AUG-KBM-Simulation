# setup for the input.tgyro
inputtgyro=OMFIT['INPUTS']['input.tgyro']
inputtgyro['TGYRO_RMIN']=0.2
inputtgyro['TGYRO_RMAX']=0.8
OMFIT['SCRIPTS']['main.py'].run()
#OMFIT['WholeInfo']['SCRIPTS']['inputgyro2cgyro.py']=OMFITpythonTask('/home/jianx/mymodule/WholeInfo/SCRIPTS/inputgyro2cgyro.py')
#OMFIT['WholeInfo']['SCRIPTS']['inputgyro2cgyro.py'].run()