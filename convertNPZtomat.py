import numpy as np
import scipy.io
fpath = "E:/2021-09-02-aged-test/testsession211015/2021-10-15-30min/cluster_rejection_mask.npz"
data = np.load(fpath)
#scipy.io.savemat("path/to/matfile.mat", mdict=file)
scipy.io.savemat("E:/2021-09-02-aged-test/testsession211015/2021-10-15-30min/cluster_rejection_mask.mat", mdict=data)
