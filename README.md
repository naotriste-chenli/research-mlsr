## CHANGE LOG:
    -   Created the script for redocking of Astex Diverse Dataset (approx. 01/23/2026)
    -   Adjusting of GNINA's arguments for the minimization for Iteration_1.sh & Iteration_2.sh (approx. 01/28/2026)
    -   Implemented a pose selection criteria for Iteration_3.sh (02/02/2026)
    -   Recoding of the algorithm in Iteration_4.sh for PopOS! Since the custom kernel drivers for CachyOS caused issues in booting of the laptop (02/06/2026) 
    -   Fixed a minimization problem where RMSD values would be evaluated as "inf" when using --cnn_scoring refinement and --cnn_scoring all by implementing a self-dependent pose selection criteria for Iteration_4.sh (02/06/2026)

## RMSD OF REDOCKING BENCHMARK LOG:
    - Iteration_1.sh - MEAN: 2.8946391114754, MEDIAN: 1.8659                   , SUCCESS RATE: 98.82% (84/85 proteins were docked with high confidence)
    - Iteration_2.sh - MEAN: inf            , MEDIAN: 2.1491813233766 (invalid), SUCCESS RATE: 0%
    - Iteration_3.sh - MEAN: 1.8643846007287, MEDIAN: 1.45537                  , SUCCESS RATE: 84.70% (72/85 proteins were docked with high confidence)
    - Iteration_4.sh - MEAN: 1.0283359858491, MEDIAN: 0.951621                 , SUCCESS RATE: 95.29% (81/85 proteins were docked with high confidence)