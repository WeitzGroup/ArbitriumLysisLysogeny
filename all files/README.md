# ArbitriumLysisLysogeny

Time-scales modulate optimal lysis-lysogeny decision switches and near-term phage fitness

Code to find optimal switching strategy from lysis to lysogeny to maximize fitness. The probability of lysogeny is modelled as a function of the optimization paramters and the state. Based on the cost function (here the fitness of the phage), the optimal parameters are computed through gradient descent using Armijo step-size. 

Navigating the code:

main_all.m calls helper function files to generate optimal switching and shaping parameters for different time horizons/resource conditions/arbitrium production rate etc. All optimal parameters have been generated and saved as mat files for quick generation of results. Once the parameters are saved, the plot_all.m file needs to be run to generate plots using this data. 

There are 3 options for running main_all.m (line 3, option = X). Option 1 corresponds to iterating over different values of R0 (initial resource concentration) and J(resource influx); option 2 corresponds to iterating over different production rates of arbitrium - from 10^5 to 10^8 (relavant to Fig. S1); option 3 is for lower production rate of 10^7 - for Fig. 3(b) in manuscript.

Option 1 generates all mat files of the format 'optimal_params_r0_X_j_Y.mat' (where R0 = X and J = Y) 
Option 2 generates mat file 'optimal_params_r0_40_j_0_sweep.mat' (sweeping over different production rates of arbitrium - from 10^5 to 10^8)
Option 3 generates mat file 'optimal_params_r0_40_j_0_low_production.mat' (for production rate of arbitrium of 10^7)


Helper functions and description:

monod.m - returns the monod function value (used in system dynamics)

probability.m - returns probability of lysogeny based on current state and optimization parameter values

cost_function.m - returns the cost based on the cost function used in the optimization algorithm

forward_euler.m - computes and returns the state trajectory based on system dynamics

compute_derivative_system.m - returns the derivative of the state at time t based on state variables using the system dynamics

backward_euler_costate.m - returns the costate trajectory computed backwards from final time based on costate dynamics

compute_derivative_costate.m - returns the costate derivative at current time

grad_J.m = returns the gradient of the cost w.r.t. the optimization parameters

compute_dfdtheta.m - returns the derivative of the state dynamics w.r.t. the parameters

compute_dfdx_system.m - returns the derivative of the state dynamics w.r.t. the state parameters

minimum_step.m - takes Armijo step based on gradient w.r.t. optimization parameters, cost function and armijo conditions
 

Generating optimal parameters:
Run file named 'main_all.m'; generates a mat file with saved optimal parameters. 

Parameters to edit (present in file main_all.m): 

Line 5 changes the production rate (from 10^5 to 10^8 for Fig. S1, currently set to single value of 5 * 10^7); Fig. 3 used production rates of 10^7 and 5 * 10^7; other results only use base param of 5 * 10^7

Data files:

File name 'optimal_params_r0_X_j_Y.mat' stores the optimal parameters for initial conditions of R0 = X and pars.J = Y; these optimal params are used for generating the figures. 

File name 'optimal_params_r0_40_j_0_low_production.mat' stores the optimal params for a production rate of 10^7 (instead of 5 * 10^7); used for Fig. 3b

File name 'optimal_params_r0_40_j_0_sweep.mat' stores optimal params for different production rates varying from 10^5 to 10^8. Used for Fig. S1

 
Figures:

Run file named 'plot_all.m' to generate all figures in paper - no changes needed to file (all relevant data to generate figures is already included in repo as .mat files). These mat files can be generated via main_all.m by changing the option (line 3 of 'main_all.m'). 

Figure 2:
Temperate phage-bacteria infection dynamics for different fixed probabilities of lysogeny (P= 0, P= 0.5 or P= 1 wherePis the probability of lysogeny) for 48 hours with an MOI of 0.01. 

Figure 3:
Comparison  of  the  optimal  lysis-lysogeny  decision  response  functions  (function  of  arbitrium  moleculeconcentration)  given  variations  in  time  horizon  (from  12  hours  to  48  hours).

Figure 4:
Temperate  phage-bacterial  infection  dynamics  for  an  optimal  switching  strategy  in  the  lysis-lysogenydecision as a function of the final time horizon. 

Figure 5:
Phage  fitness  comparison  for  optimal  and  fixed  probabilities  of  lysogeny. 

Figure 6:
Effect  of  resource  level  (change  in  initial  concentration  and  influx)  on  optimal  lysis-lysogeny  switching point. 



