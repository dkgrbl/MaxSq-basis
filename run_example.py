# author: Denis Kopylov


import numpy as np
import matplotlib.pyplot as plt

import msq_funcs
import plot_funcs


import os
FolderPath = os.path.dirname(os.path.realpath(__file__)) + "/" 



def run_example():
 
    ####
    #### Load data
    ####

    PDC =  np.load(os.path.join(FolderPath,"data_63.npz"))
    cov_matrix_initial = PDC["cov_matrix"]

    msq_funcs.check_covariance_matrix(cov_matrix_initial)


    # ####
    # #### Compute first N_modes from MSq basis 
    # ####

    N_modes = 200
    U_msq = msq_funcs.build_MSq_basis_swap(cov_matrix_initial, N=N_modes)
    covariance_matrix_msq = msq_funcs.transfrom_covariance_with_U(cov_matrix_initial, U_msq)
    
    ####
    #### Plot covariance matrix, squeezing and 3 first modes
    ####

    plot_funcs.plot_covariance_matrix_xxpp_log(cov_matrix_initial, title="Covariance, initial")
    plot_funcs.plot_covariance_matrix_xxpp_log(covariance_matrix_msq, title="Covariance, Msq-basis")
    plot_funcs.plot_squeezing(covariance_matrix_msq)
    
    plot_funcs.plot_mode_profile(U_msq, n=0)
    plot_funcs.plot_mode_profile(U_msq, n=1)
    plot_funcs.plot_mode_profile(U_msq, n=2)
    plot_funcs.plot_mode_profile(U_msq, n=43)
    # plot_funcs.plot_mode_profile(U_msq, n=42)

 
if __name__ == '__main__':
    run_example()
    plt.show()
