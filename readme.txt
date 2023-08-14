The Dynare models can be found in the files "Hansen2020_NotLinearized_Ramsey.mod" and "Corr_Hansen2020_Linearized_OSR.mod",
where the former contains the linear model used for the results on Ramsey optimal policy, and the latter contains the
corrected linearized model used for the results on Taylor rules.

To run the models, the files "CB_Ramsey_calc.m" and "CB_OSR_calc.m" can simply be executed in the directory, which produces 
all non looping results used in the thesis. The results requiring looping over different parameter values are outsourced into
the files "Loop_Ramsey.m" and "Loop_OSR.m", due to significantly longer run times. The possible configurations for each file are
commented in the respective header. The file "Plot_weights.m" simply produces the non-dynamic results for the weights on the quadratic
gaps equations.

The above files are able to produce all graphical output used in the thesis, as well as save all results from the simulations in the
Matlab workspace. The outputs in the workspace are saved via structures and are explained below:

CB_regimes : saves the output from all optimal policy simulations
policy : saves the output from all taylor rule simulations
L : saves the loss on gaps for each regime/taylor rule
CEL : saves the consumption equivalent loss on gaps for each regime/taylor rule
To : saves the T_0 values for each regime/taylor rule
L_inf : saves the loss on the inflation gap for each regime/taylor rule
L_out : saves the loss on the output gap for each regime/taylor rule
L_cons : saves the loss on the consumption gap for each regime/taylor rule
CEL_inf : saves the consumption equivalent loss on the inflation gap for each regime/taylor rule
CEL_out : saves the consumption equivalent loss on the inflation gap for each regime/taylor rule
CEL_cons : saves the consumption equivalent loss on the inflation gap for each regime/taylor rule
results_cell : this cell is saved within the structures CB_regimes and policy and contains the results from looping over parameters

The structures hold all results, so they can simply be called by specifying the respective regime. So for example L.ric will show the
loss for the Ricardian regime and CEL.key the consumption equivalent loss for the Keynesian regime. The possible regimes/taylor rule
configurations are:

Ramsey: "ric", "key", "avg", "optim"
Taylor: "ric", "key", "ineq", "zero"