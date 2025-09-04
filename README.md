# mdpi-liver
This code accompanies our paper "Modeling competitive dynamics of healthy and cancerous liver cells with YAP hyperactivation."

- - - 

Description of files:

genFig1_7.m - this file generates Figures 1-7 from the paper (including all of their subfigures).

genFig8.m - this file generates Figure 8 from the paper.

runner.m - run this file to perform data-fitting.

liverMultiStart.m - this is called by the "runner.m" file during the data-fitting process; it runs fmincon multiple times with random initial guesses in the parameter space.

liverFunctional.m - this is calld by the "runner.m" file during the data-fitting process; it calculates the value of the functional we are minimizing during the data-fitting process.

parameter_control.m - this file contains parameter values used in many simulations when generating figures.

dimensional_odes.m - this file describes Model 5 in our paper, used for finding parameter estimates when fitting to cancer data presented by Anton et al.

odes_after_fitting.m - this file describes the ode system we developed in our paper and includes a third equation tracking the total volume of the liver (d(L+C)/dt).

twoD_odes.m - this file describes the ode system (Model 3) developed in our paper.

- - -
