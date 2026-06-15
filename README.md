# EE 3006 Transformer Design - Bonus Optimization

This repository contains the MATLAB optimization code for the design of a single-phase 1 kVA 220/110 V, 50 Hz shell-type transformer. The script utilizes a **Differential Evolution (DE)** algorithm to find the Pareto-optimal design that minimizes both material cost and total analytical losses.

## Project Overview

In power transformer design, reducing losses (operational cost) and material usage (capital cost) are often conflicting objectives. This code employs a multi-objective optimization approach converted into a single weighted objective function. The optimization ensures that the final design is physically and electromagnetically feasible by applying strict penalty constraints.

### Design Variables
The Differential Evolution algorithm explores the design space by modifying the following 4 variables:
1. `x(1)`: Primary turns ($N_1$)
2. `x(2)`: Center limb width ($W_c$) in mm
3. `x(3)`: Primary conductor gauge (Standard AWG index)
4. `x(4)`: Secondary conductor gauge (Standard AWG index)

*Note: The secondary number of turns ($N_2$) is deterministically calculated based on the required 2:1 voltage ratio.*

### Optimization Constraints
To guarantee a practical and safe transformer design, the algorithm applies heavy penalty functions if any of the following limits are exceeded:
* **Magnetic Saturation:** Maximum flux density $B_{max} \le 1.50 \text{ T}$
* **Thermal Limit:** Current density $J \le 3.0 \text{ A/mm}^2$
* **Mechanical Feasibility:** Window fill factor $\le 0.45$ (Ensures the windings physically fit into the core window)
* **Voltage Ratio:** Allowed turns-ratio error $\le 2\%$
* **Cost Limit:** The optimized material cost must not exceed the baseline design cost.

## Requirements
* **MATLAB** (R2018a or newer is recommended). No additional toolboxes are required as the Differential Evolution algorithm is built from scratch within the script.

## How to Run
1. Open `transformer_optimization.m` in MATLAB.
2. Run the script.
3. The script will automatically execute the Differential Evolution loop for 140 generations with a population size of 70.
4. Progress and final comparisons will be printed in the **Command Window**.
5. The generated figures will be automatically saved to a folder named `Transformer_Report_Figures_Bonus` on your **Desktop**.

## Outputs

### 1. Command Window Summary
The script prints a detailed side-by-side comparison of the **Baseline Design** and the **Optimized Design**. It includes electrical parameters (current densities, resistances, losses), magnetic parameters ($B_{max}$), mechanical parameters (fill factor, dimensions, mass), and total material cost.

### 2. Exported Figures
The algorithm generates and saves four high-resolution (300 DPI) plots for report integration:
* `bonus_DE_convergence.png`: Shows the convergence of the best objective value ($F$) over 140 generations.
* `bonus_cost_loss.png`: A scatter plot demonstrating the trade-off between Material Cost and Total Analytical Loss among all feasible designs.
* `bonus_efficiency_cost.png`: Illustrates the relationship between analytical efficiency and material cost.
* `bonus_bmax_map.png`: A 2D map of the design space ($N_1$ vs. $W_c$) colored by the maximum magnetic flux density ($B_{max}$).

## Authors
* **Berk Kaan** - Electrical and Electronics Engineering, AGU
* **Berkay Mete** - Electrical and Electronics Engineering, AGU

*Developed for the EE 3006 Electromechanical Energy Conversion Laboratory Course.*
