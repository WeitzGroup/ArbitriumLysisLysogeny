# ArbitriumLysisLysogeny

Time-scales modulate optimal lysis-lysogeny decision switches and near-term phage fitness

Code to find optimal switching strategy from lysis to lysogeny to maximize fitness. The probability of lysogeny is modelled as a function of the optimization paramters and the state. Based on the cost function (here the fitness of the phage), the optimal parameters are computed through gradient descent using Armijo stepp-size.

Navigating the code:

main_all.m calls helper function files to generate optimal switching and shaping parameters for different time horizons/resource conditions/arbitrium production rate etc. All optimal parameters have been generated and saved as mat files for quick generation of results. Once the parameters are saved, the plot_all.m file needs to be run to generate plots using this data. 

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

Line 40 changes the production rate from 10^5 to 10^8 (Relevant for Fig. S1); Fig. 3 used production rates of 10^7 and 5*10^7; other results only use base param of 5*10^7

Line 46 changes time horizon (currently fixed at 24 hours, can be done from 12 to 48 hours)

 
Figures:

Run file named 'plot_all.m' to generate all figures in paper - no changes needed to file (all relevant data to generate figures is already included in repo as .mat files).

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



