## CHANGE LOG:
    -   Created the script for redocking of Astex Diverse Dataset (approx. 01/23/2026)
    -   Adjusting of GNINA's arguments for the minimization for Iteration_1.sh & Iteration_2.sh (approx. 01/28/2026)
    -   Implemented a pose selection criteria for Iteration_3.sh (02/02/2026)
    -   Recoding of the algorithm in Iteration_4.sh for PopOS! Since the custom kernel drivers for CachyOS caused issues in booting of the laptop (02/06/2026) 
    -   Fixed a minimization problem where RMSD values would be evaluated as "inf" by implementing a self-dependent pose selection criteria for Iteration_5.sh   (02/06/2026)

## BENCHMARK LOG (CONVERGENCE OF MEAN RMSD, MEDIAN RMSD, SUCCESS RATE):
Iteration_1.sh - MEAN: 2.6183101377049, MEDIAN: 2.8717171291866, SUCCESS RATE: 100%
Iteration_2.sh - 
Iteration_3.sh - 
Iteration_4.sh - 
Iteration_5.sh - 