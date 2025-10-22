tgyroout=OMFIT['WholeInfo']['OUTPUTS']['TGYRO']
for item in tgyroout.keys():
    if item[0:4]=='out.':
        tgyroout[item].deploy('/fusion/projects/xpsi/east-d3d_experiments/Jian/collegues/matthias/nsh-kefit-xj/tgyroout')