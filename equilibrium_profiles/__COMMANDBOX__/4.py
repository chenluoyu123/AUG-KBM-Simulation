for key in OMFIT['OUTPUTS']['Profiles_gen'].keys():
    OMFIT['OUTPUTS']['Profiles_gen'][key].deploy(dirpath+'Profiles_gen/'+key)
OMFIT['OUTPUTS']['TGYRO'].deploy(dirpath+'TGYRO')