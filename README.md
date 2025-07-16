# Optimized S-Box Modeling Techniques 

This repository contains several S-box modeling techniques for evaluating and improving the security of S-boxes used in symmetric cryptography algorithms. The main focus is generating, reducing, and augmenting inequalities for various S-boxes, including 4-bit and 5-bit configurations.

## Directory Structure

- **Direct Inequality Generation**: This section contains scripts and results related to the direct generation of inequalities for different S-boxes.
  - `4-bit_sboxes_256_256/`: Contains inequality generation results for 4-bit S-boxes with 256 as the range for coefficients.
  - `4-bit_sboxes_500_500/`: Contains inequality generation results for 4-bit S-boxes with 500 as the range for coefficients.
  - `direct_inequality_generation.py`: The main script for generating inequalities for S-boxes.

- **Greedy Generation and Reduction**: Implements a greedy approach to generate and reduce inequalities.
  - `Results/`: Contains the greedy generation and reduction process results for various S-boxes.
  - `greedy_generation_and_reduction.py`: Script implementing the greedy algorithm for inequality generation and reduction.

- **Iterative Inequality Augmentation**: Focuses on iteratively improving the generated inequalities for different S-boxes.
  - `4-bit_sboxes/`: Contains results for 4-bit S-boxes.
  - `5-bit_sboxes/`: Contains results for 5-bit S-boxes.
  - `iterative_inequality_augmentation.py`: Script for augmenting inequalities iteratively.

- **Modified Greedy Approach**: Implements a modified version of the greedy algorithm for generating and reducing inequalities.
  - `Results/`: Contains results of the modified greedy approach.
  - `modified_greedy_approach.py`: Script implementing the modified greedy approach.

> All results in the above directories were computed using the Difference Distribution Table (DDT).

---
<!-- 
- **SBOXES**: Contains JSON files for different S-boxes.
  - `4_bit_sboxes.json`: List of 4-bit S-boxes.
  - `5_bit_sboxes.json`: List of 5-bit S-boxes.

### Additional Cryptanalysis Metrics
These directories contain results based on different S-box analysis metrics:

- **BCT**: Uses the Boomerang Connectivity Table.
- **DPT**: Uses the Division Property Table.
- **LAT**: Uses the Linear Approximation Table.

Each directory includes outputs from all three modeling techniques:
- Direct Inequality Generation
- Iterative Inequality Augmentation
- Modified Greedy Approach -->