#Metal-Mott insulator interfaces


Published:

    "Metal-Mott insulator interfaces"
    http://arxiv.org/abs/1606.00936

PREINSTALL:

    CTQMC (continuous-time quantum Monte Carlo) for DMFT calculation
        -"http://hauleweb.rutgers.edu/downloads/"
        -See also additional preinstallations such as python and C++
        -Download "dmft_w2k.tgz" (version 2015)
        -Extract:  "tar -zxvf dmft_w2k.tgz"
        -Install according to "http://hauleweb.rutgers.edu/tutorials/Installation.html"
        -"export WIEN_DMFT_ROOT=where_your_bin_folder_is"



Execution:

"phase_diagram/cubic_DMFT.py"

    - Calculate phase diagram of 3D-cubic Hubbard model
    - options:
        U   : onsite Hubbard interaction
        mu  : chemical potential
        beta    : inverse temperature
        Mstep   : number of Monte-Carlo steps

"transition/mit.py"
    
    - Perform the DMFT calculation with fixed U/mu/beta parameter in the coexistence regime across given number of sites
    - There should exist boundary conditions at the ends of the sites at given U/mu/beta:
        "Sig_i.out": insulating self-energy, boundary condition of the rightmost site
        "Sig_m.out": metallic self-energy, boundary condition of the leftmost site
    - options:
        U : onsite interaction
        mu : chemical potential
        beta : inverse temperature
        Nitt : number of iteration
        Nsite : number of total sites across x-axis
        Mstep : number of Monte-Carlo steps
        
