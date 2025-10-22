# this script is used to update the scripts
scripts=root['SCRIPTS']
plots=root['PLOTS']
dirscripts='/home/users/xiangjian/mymodule/CGYRO_SCAN/SCRIPTS/'
dirplots='/home/users/xiangjian/mymodule/CGYRO_SCAN/PLOTS/'
for item in scripts.keys():
    scripts[item]=OMFITpythonTask(dirscripts+item)
for item in plots['CGYROalone'].keys():
    if item.split('.')[-1]=='py':
        plots['CGYROalone'][item]=OMFITpythonTask(dirplots+'CGYROalone/'+item)
for item in plots['CGYROalone']['assist'].keys():
    plots['CGYROalone']['assist'][item]=OMFITpythonTask(dirplots+'CGYROalone/assist/'+item)
for item in plots['CGYROscan'].keys():
    if item.split('.')[-1]=='py':
        plots['CGYROscan'][item]=OMFITpythonTask(dirplots+'CGYROscan/'+item)
for item in plots['CGYROscan']['eigen'].keys():
    plots['CGYROscan']['eigen'][item]=OMFITpythonTask(dirplots+'CGYROscan/eigen/'+item)
for item in plots['CGYROscan']['assist'].keys():
    plots['CGYROscan']['assist'][item]=OMFITpythonTask(dirplots+'CGYROscan/assist/'+item)
