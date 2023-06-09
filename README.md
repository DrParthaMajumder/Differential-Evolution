# Differential-Evolution (Modified)

# Developer: Dr. Partha Majumder

# Contact Details: parthamajpk@gmail.com

# Short note on differential evolution:
Differential Evolution (DE) is a popular and powerful evolutionary algorithm used for global optimization of functions. It is known for its simplicity, robustness, and effectiveness in finding optimal solutions in various domains. DE operates by maintaining a population of candidate solutions, typically represented as vectors in a search space. The algorithm iteratively improves the population by applying mutation, crossover, and selection operators. In the mutation stage, DE generates new candidate solutions by combining existing solutions. It does this by adding a scaled difference vector between two randomly selected individuals to a third individual in the population. This process promotes exploration of the search space. The crossover stage involves combining the mutated individuals with the original individuals to create trial solutions. A binomial crossover scheme is commonly used, where each component of the trial solution is selected from either the mutated or original solution based on a crossover probability. In the selection stage, the trial solutions are evaluated using an objective function, and the individuals with better fitness replace their corresponding original solutions in the population. This selection process ensures that the population evolves towards better solutions over successive iterations. DE has proven to be highly effective in solving various optimization problems, including parameter estimation, function optimization, and engineering design. It is known for its versatility, as it does not rely on derivative information and can handle both continuous and discrete search spaces. 
Overall, Differential Evolution offers a simple yet powerful approach to global optimization. By combining mutation, crossover, and selection operators, DE efficiently explores the search space and converges towards optimal solutions.

# License
Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
