//All variables expressed in percent deviations from steady state


//Endogenous variables
var 
YV YnV AV CrV CkV NV RV PiV wV dV trV tkV
YVS dCV dYV;

//Exogenous variables
varexo  
eps;

//Parameters
parameters 
lam tau delt gam alph bet chi phi thetaP psiP rho phiPi phiY phiC Tp AB
YB YnB CrB CkB NB RB PiB wB dB trB tkB
CrSB CkSB lamS Way Wy Wpi Wdc
@#if policy == "ric" && specific_OSR_weights
    Wy_spec Wpi_spec Wdc_spec
@#elseif policy == "key" && specific_OSR_weights
    Wy_spec Wpi_spec Wdc_spec
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
phiC = 0;    //Monetary Policy consumption gap

Tp  = 1/(thetaP - 1); //from CES Elasticity (monopolistic distortions)

//Steady-State Parameters:
YB=  		 0.897735;
YnB= 		 0.897735;
CrB= 		 1.10646;
CkB= 		 0.58465;
NB=  		 0.866025;
RB=  		 1.00756;
PiB= 		 1;
wB=  		 0.777461;
dB=  		 0.33665;
trB= 		 0.545374;
tkB= 		 0.0235655;

//Helper Parameters:
CrSB = CrB/YB;
CkSB = CkB/YB;
lamS = lam*CkB/YB;

Way = (1+phi)/(1-alph) + (lam-lamS)/(1-lamS)^2 * ((1-(1-tau)*delt)*(1-alph)/CkSB)^2 * ((1+phi)/(1-alph))^2 - (lam-lamS)/(1-lamS)*gam*(1-alph)/CkSB * (1-((1-(1-tau)*delt)*(1-alph)/CkSB)) * (1+phi)/(1-alph);
Wy = (1+phi)/(1-alph) + (lam-lamS)/(1-lamS)^2 * ((1-alph)*(1-(1-tau)*delt)/CkSB)^2 * ((1+phi)/(1-alph))^2;
Wpi = psiP * (1 - (lam-lamS)/(1-lamS)*(1-(1-tau)*delt/CkSB));
Wdc = lamS*(1-lamS);

@#if policy == "ric" && specific_OSR_weights
    Wy_spec = (1+phi)/(1-alph);
    Wpi_spec = psiP;
    Wdc_spec = 0;
@#elseif policy == "key" && specific_OSR_weights
    Wy_spec = (1+phi)/(1-alph) + 1/(1-CkSB) * ((1-alph)*(1-(1-tau)*delt)/CkSB)^2 * ((1+phi)/(1-alph))^2;
    Wpi_spec = psiP * (1 - (1-CkSB)/(1-CkSB)*(1-(1-tau)*delt/CkSB));
    Wdc_spec = CkSB*(1-CkSB);
@#endif

model(linear);
    //Euler Ricardian
    CrV = CrV(+1) - RV + PiV(+1);
                                 
    //"Euler" Keynesian
    CkB*CkV/(wB*NB) + Tp*YB*YV/(wB*NB) - tkB*tkV/(wB*NB) = -gam*AV + wV + NV;

    //Labour Supply
    wV = phi*NV + YV;

    //Production Function
    YV = AV + (1-alph)*NV;

    //dividend yield
    dB*dV = (1+Tp)*YB*YV - wB*NB*(wV+NV) - psiP/2*YB*(PiB^2*(2*PiV+YV) - 2*PiB*(PiV+YV) + YV);

    //Phillips Curve
    (2*PiB^2-PiB)*PiV = bet*(PiB^2-PiB)*(YV(+1)-YV-CrV(+1)+CrV) + bet*(2*PiB^2-PiB)*PiV(+1) + thetaP/psiP*wB/AB^(1/(1-alph))*YB^(alph/(1-alph))*(wV + alph/(1-alph)*YV - AV/(1-alph));

    //Fiscal Policy (redistribution)
    dV = (1-lam)*trB*trV/(delt*dB) + lam*tkB*tkV/(delt*dB);

    //Fiscal Policy (tax rate)
    tkV = dV;

    //Goods Market Clearing
    YV = lam*CkB*CkV/YB + (1-lam)*CrB*CrV/YB + psiP/2*(PiB^2*(2*PiV+YV) - 2*PiB*(PiV+YV) + YV);

    //Natural Output
    YnV = AV;

    //Monetary Policy
    @#if policy == "ric"
        RV = phiPi*PiV + phiY*YV + phiC*CrV;
    @#elseif policy == "key"
        RV = phiPi*PiV + phiY*YV + phiC*CkV;
    @#elseif policy == "ineq"
        RV = phiPi*PiV + phiY*YV + phiC*(CrV-CkV);
    @#elseif policy == "zero"
        RV = phiPi*PiV + phiY*YV + 0*(CrV-CkV);
    @#endif

    //Productivity shock
    AV = rho*AV(-1) + eps;

    //Helper variables
    YVS = Way/Wy*AV;
    
    dYV = YV - YVS;
    @#if policy == "ric"
        dCV = CrV;
    @#elseif policy == "key"
        dCV = CkV;
    @#elseif policy == "ineq"
        dCV = CrV - CkV;
    @#elseif policy == "zero"
        dCV = CrV - CkV;
    @#endif
end;

//calc. and check steady state
steady;

//check residuals
resid;

//check eigenvalues
check;

//produce shocks
shocks;
var eps;
stderr 0.01;
end;

//sensitivity analysis
@#if sense
    estimated_params;
    phiPi, 1.5, -2, 2;
    phiY, 0.125, -2, 2;
    phiC, -0.1, -2, 2;
    end;

    varobs
    YV YnV AV CrV CkV NV RV PiV wV dV trV tkV;

    dynare_sensitivity;
@#endif

//simulation

@#if policy == "zero"
   stoch_simul(irf=100, nograph);
@#elseif (policy == "ric" || policy == "key") && specific_OSR_weights
    optim_weights;
        PiV Wpi_spec;
        dYV Wy_spec;
        dCV Wdc_spec;
    end;
    
    osr_params phiC;
    
    osr_params_bounds;
        phiC, -2, 2;
    end;
    
    osr(opt_algo=9, irf=100, nograph);
    oo_.osr.optim_params;
@#else
    optim_weights;
        PiV Wpi;
        dYV Wy;
        dCV Wdc;
    end;
    
    osr_params phiC;
    
    osr_params_bounds;
        phiC, -2, 2;
    end;
    
    osr(opt_algo=9, irf=100, nograph);
    oo_.osr.optim_params;
@#endif