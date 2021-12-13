## N-1_SCUC
N-1 security-constrained unit commitment (N-1 SCUC) that considers all contingency-case network constraints.

This AMPL program solves an N-1 SCUC problem that considers all contingency-case (N-1) network constraints. Reserve requirment is still included in the Master problem, which may not be necessary since N-1 constraints are enforced but may help reduce the number of iterations. 

This N-1 SCUC is solved by Bender decomposition, implemented in the code SCUC.run. The code file, N-1_SCUC.mod, defines the models for the Master problem and slave sub-problems, as well as the cut. For this work, N-1_SCUC.mod is not the main runnable code while the code file, SCUC.run, is the script main runnable program.

This program takes 1,325 seconds (~22 minutes) on a laptop: Intel(R) Core(TM) i7-8850H CPU @ 2.60GHz, 32 GB RAM, Windows 10.

This program lays the foundation for the following paper: 

[Arun Venkatesh Ramesh, Xingpeng Li and Kory Hedman, “An Accelerated-Decomposition Approach for Security-Constrained Unit Commitment with Corrective Network Reconfiguration,” IEEE Transactions on Power Systems, online early access, Jul. 2021. (DOI: 10.1109/TPWRS.2021.3098771)](https://ieeexplore.ieee.org/document/9492752)

## Contact:
Dr. Xingpeng Li

University of Houston

Email: xli83@central.uh.edu

Website: https://rpglab.github.io/


## License:
This work is licensed under the terms of the Creative Commons Attribution 4.0 (CC BY 4.0) license. 
https://creativecommons.org/licenses/by/4.0/


## Disclaimer:
The author doesn’t make any warranty for the accuracy, completeness, or usefulness of any information disclosed; and the author assumes no liability or responsibility for any errors or omissions for the information (data/code/results etc) disclosed.
