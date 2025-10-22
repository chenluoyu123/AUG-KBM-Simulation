#OMFIT['WholeInfo']=OMFIT.load('/home/users/xiangjian/mymodule/WholeInfo/OMFITsave.txt')
dirpath='/home/users/xiangjian/augwork/baseinfo/'
OMFIT['INPUTS']['input.gacode']=OMFITgacode(dirpath+'input.gacode')
inputtgyro=OMFIT['INPUTS']['input.tgyro']
inputtgyro['LOC_N_ION']=1