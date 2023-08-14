//Endogenous variables
var 
Y Yn A Cr Ck N R Pi w d tr tk 
YS
log_lin_Ck log_lin_Cr log_lin_Y log_lin_Pi log_lin_YS log_lin_A
@#if CB_objective == "ric"
    YS_ric log_lin_YS_ric
@#elseif CB_objective == "key"
    YS_key log_lin_YS_key
@#endif
;

//Exogenous variables
varexo  
eps;

//Parameters
parameters 
lam tau delt gam alph bet chi phi thetaP psiP rho phiPi phiY phiC Tp AB
YB CrB CkB
CrSB CkSB lamS Way Wy Wpi Wdc
@#if CB_objective == "ric"
    Way_ric Wy_ric
@#elseif CB_objective == "key"
    Way_key Wy_key 
@#endif
;

//Inialization of parameter values
AB  = 1;    //steady_state of exogenous process A

lam  = 0.4;  //Share of Keynesian agents
tau  = 0.93;  //Redistribution of Profits (1- share of dividends taxed)
delt  = 1;  //Redistribution of Profits (share of dividends taxed)
gam  = 1.67;  //Degree of Skill Bias
alph  = 0.25;  //Profits Share
bet  = 0.9925;  //Discount Factor
chi  = 1;  //Labor Disutility
phi  = 1;   //Inverse Frisch Elasticity
thetaP  = 9;   //CES Elasticity
psiP  = 372.8;   // Price Adjustment
rho  = 0.9;   // Persistence of Shock

phiPi = 1.5;   //Monetary Policy Inflation 
phiY = 0.125;   //Monetary Policy Output gap
phiC = 0;   //Monetary Policy Consumption gap

Tp  = 1/(thetaP - 1); //from CES Elasticity (monopolistic distortions)

//Helper Parameters
YB=  		 0.897735;
CrB= 		 1.10646;
CkB= 		 0.58465;

CrSB = CrB/YB;
CkSB = CkB/YB;
lamS = lam*CkB/YB;

Way = (1+phi)/(1-alph) + (lam-lamS)/(1-lamS)^2 * ((1-(1-tau)*delt)*(1-alph)/CkSB)^2 * ((1+phi)/(1-alph))^2 - (lam-lamS)/(1-lamS)*gam*(1-alph)/CkSB * (1-((1-(1-tau)*delt)*(1-alph)/CkSB)) * (1+phi)/(1-alph);
Wy = (1+phi)/(1-alph) + (lam-lamS)/(1-lamS)^2 * ((1-alph)*(1-(1-tau)*delt)/CkSB)^2 * ((1+phi)/(1-alph))^2;
Wpi = psiP * (1 - (lam-lamS)/(1-lamS)*(1-(1-tau)*delt/CkSB));
Wdc = lamS*(1-lamS);

@#if CB_objective == "ric"
    Way_ric = (1+phi)/(1-alph);
    Wy_ric = (1+phi)/(1-alph);
@#elseif CB_objective == "key"
    Way_key = (1+phi)/(1-alph) + (1-CkB/YB)/(1-CkB/YB)^2 * ((1-(1-tau)*delt)*(1-alph)/CkSB)^2 * ((1+phi)/(1-alph))^2 - (1-CkB/YB)/(1-CkB/YB)*gam*(1-alph)/CkSB * (1-((1-(1-tau)*delt)*(1-alph)/CkSB)) * (1+phi)/(1-alph);
    Wy_key = (1+phi)/(1-alph) + (1-CkB/YB)/(1-CkB/YB)^2 * ((1-alph)*(1-(1-tau)*delt)/CkSB)^2 * ((1+phi)/(1-alph))^2;
@#endif

model;
    //Euler Ricardian
    Cr(+1) = bet*R*Cr/Pi(+1);
                                 
    //"Euler" Keynesian
    Ck = (A/AB)^(-gam)*w*N - Tp*Y + tk;

    //Labour Supply
    w = chi*N^phi*Y;

    //Production Function
    Y = A*N^(1-alph);

    //dividend yield
    d = (1+Tp)*Y - w*N - psiP/2*Y*(Pi-1)^2;

    //Phillips Curve
    Pi*(Pi-1) = bet*Cr/Cr(+1)*Y(+1)/Y*Pi(+1)*(Pi(+1)-1) + thetaP/psiP*(1/(1-alph)*w/A^(1/(1-alph))*Y^(alph/(1-alph)) - (1+Tp)*(thetaP-1)/thetaP);

    //Fiscal Policy (redistribution)
    tr = (delt*d - lam*tk)/(1-lam);

    //Fiscal Policy (tax rate)
    tk = (1-tau)*delt*d;

    //Goods Market Clearing
    Y = lam*Ck + (1-lam)*Cr + psiP/2*Y*(Pi-1)^2;

    //Natural Output
    Yn = A*((1+Tp)*(1-alph)*(thetaP-1)/thetaP*(1/chi))^((1-alph)/(1+phi));

    //Productivity shock 
    A = AB*exp(eps);


    //Helper variables
    YS = A^(Way/Wy);

    log_lin_Ck = log(Ck/steady_state(Ck));
    log_lin_Cr = log(Cr/steady_state(Cr));
    log_lin_Pi = log(Pi/steady_state(Pi));
    log_lin_Y = log(Y/steady_state(Y));
    log_lin_YS = log(YS/steady_state(YS));
    log_lin_A = log(A/steady_state(A));

    //CB objectives
    @#if CB_objective == "ric"
        YS_ric = A^(Way_ric/Wy_ric);
        log_lin_YS_ric = log(YS_ric/steady_state(YS_ric));
    @#elseif CB_objective == "key"
        YS_key = A^(Way_key/Wy_key);
        log_lin_YS_key = log(YS_key/steady_state(YS_key));
    @#endif
end;

steady_state_model;
    A = AB;
    Yn = A*((1+Tp)*(1-alph)*(thetaP-1)/thetaP*(1/chi))^((1-alph)/(1+phi));
    Pi = 1;

    R = 1/bet;
    Y = A*((1+Tp)*(1-alph)*(thetaP-1)/thetaP*(1/chi))^((1-alph)/(1+phi));
    N = (Y/A)^(1/(1-alph));
    w = chi*N^phi*Y;
    d = (1+Tp)*Y - w*N - psiP/2*Y*(Pi-1)^2;
    tk = (1-tau)*delt*d;
    tr = (delt*d-lam*tk)/(1-lam);
    Ck = w*N - Tp*Y + tk;
    Cr = (Y-lam*Ck-psiP/2*Y*(Pi-1)^2)/(1-lam);

    YS = A^(Way/Wy);

    log_lin_Ck = 0;
    log_lin_Cr = 0;
    log_lin_Pi = 0;
    log_lin_Y = 0;
    log_lin_YS = 0;
    log_lin_A = 0;
    
    @#if CB_objective == "ric"
        YS_ric = A^(Way_ric/Wy_ric);
        log_lin_YS_ric = 0;
    @#elseif CB_objective == "key"
        YS_key = A^(Way_key/Wy_key);
        log_lin_YS_key = 0;
    @#endif
end;

@#if CB_objective == "ric"
    planner_objective(log(Cr) - N^(1+phi)/(1+phi));
@#elseif CB_objective == "key"
    planner_objective(log(Ck) - N^(1+phi)/(1+phi));
@#elseif CB_objective == "avg"
    planner_objective(log(lam*Ck+(1-lam)*Cr) - N^(1+phi)/(1+phi));
@#elseif CB_objective == "optim"
    planner_objective(lam*log(Ck)+(1-lam)*log(Cr) - N^(1+phi)/(1+phi));
@#endif

ramsey_model(planner_discount=bet);

//check and calc. steady state
steady;

//check residuals
resid;

//check eigenvalues
check;

//produce shocks
epsV = 0.01; %change for different shock

shockVs = [epsV];
for i = 1:99
    shockVs(i+1) = shockVs(i)*rho;
end

shocks;
    var eps;
    periods 1:100;
    values (shockVs);
end;

//simulation
perfect_foresight_setup(periods=100);
perfect_foresight_solver(no_homotopy);