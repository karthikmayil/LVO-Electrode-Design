# Design Principles to Govern Electrode Fabrication for the Lithium Trivanadate Cathode
Simulation results and code used for analysis and interactive plotting. Accompanies paper found at: https://iopscience.iop.org/article/10.1149/1945-7111/ab91c8

## Installation

**Important:** As of right now, these installation instructions work for those using a Mac/Linux. However, if you have jupyter notebooks set up on your Windows system, you can install the necessary packages in requirements.txt and open the notebook Design Principles to Govern Electrode Fabrication for the Lithium Trivanadate Cathode [Accompanying Code].ipynb in this repository. 

To use the interactive tool, follow the instructions below:

1) Make sure python3 is installed. The below command should work.

<code>python3 --version</code>

2) Clone this repository. Either download it from github, or type the following into your terminal:
 
<code>git clone https://github.com/karthikmayil/LVO-Electrode-Design</code>

Change into this directory:

<code>cd LVO-Electrode-Design</code>

3) The following commands need only be done once in terminal. Note: The last step, which installs necessary packages, might take a bit longer than the first two which will be very quick.

<code>chmod +x install.sh</code>

<code>chmod +x run.sh</code>

<code>./install.sh</code>

4) Open the tool and interact! This same line can be executed to open the interactive tool again.

<code>./run.sh</code>

## What can I do with this tool?

### Optimal Design

This tool can be used to create plots that show the optimal design of a cathode to optimize a selected performance metric (energy, power, capacity for the electrode and for a cell). These plots show how optimal design parameter values vary with discharge rate, and include a sensitivity analysis. For example, a plot like below can be generated, which shows the optimal design as well as a shaded region showing sensitivity of performance to the design parameters.

![Example 1](/example_images/Figure1.png "Optimal Design")

### Rate Capability

Rate capability of electrodes prepared differently may be of interest as well. Plots like the one shown below can be created for cell and electrode power and energy densities, as well as electrode capacty delivered over a range of C-rates.

![Example 2](/example_images/Figure2.png "Rate Capability")

### Detailed Design Guide

Looking at slices of the simulated data via contour maps like the one below can be very helpful in understanding the effects of different design paramters at different discharge rates. Choose the performance metric, and pick a design parameter to fix at a given discharge rate and explore a 2-D parameter space.

![Example 3](/example_images/Figure3.png "Detailed Design Guide")
