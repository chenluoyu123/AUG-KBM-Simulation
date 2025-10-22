# this script is used to update the scripts
scripts=root['SCRIPTS']
plots=root['PLOTS']
wholeinfordir='/home/users/xiangjian/mymodule/WholeInfo/'
dirscripts=wholeinfordir+'SCRIPTS/'
dirplots=wholeinfordir+'PLOTS/'
for item in scripts.keys():
    scripts[item]=OMFITpythonTask(dirscripts+item)
for item in plots.keys():
    plots[item]=OMFITpythonTask(dirplots+item)

