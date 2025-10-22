inputcgyro=OMFIT['INPUTS']['input.cgyro']
alpha=np.sum([inputcgyro['DENS_'+str(m)]*inputcgyro['TEMP_'+str(m)]*\
              (inputcgyro['DLNNDR_'+str(m)]+inputcgyro['DLNTDR_'+str(m)]) for m in np.arange(1,inputcgyro['N_SPECIES']+1) ])
alpha=alpha*inputcgyro['BETAE_UNIT']*inputcgyro['Q']**2*inputcgyro['RMAJ']
print(alpha)