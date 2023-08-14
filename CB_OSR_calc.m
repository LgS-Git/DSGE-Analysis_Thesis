
% Define what policies to run
policy_types = ["ric", "key", "ineq", "zero"];
%policy_types = ["ric", "key", "ineq"];

% Toggle on for policy specific OSR weights
specific_OSR_weights = 1;

% Toggle on to start sensitivity analysis
sense = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calc Taylor policies
for k=1:length(policy_types)
    dynare_command = 'dynare Corr_Hansen2020_Linearized_OSR.mod -Dpolicy="'+policy_types(k)+'" -Dspecific_OSR_weights='+specific_OSR_weights + ' -Dsense='+sense;
    eval(dynare_command);
    fname = "OSR_" + policy_types(k);
    save(fname);
    clearvars -except policy_types specific_OSR_weights sense;
end

%Load Taylor policies
policy = struct();
for k=1:length(policy_types)
    fname = "OSR_" + policy_types(k); 
    policy.(policy_types(k)) = load(fname);
    del_command = "delete " + fname + ".mat";
    eval(del_command);
end

%calc loss
L = struct();
L_inf = struct();
L_out = struct();
L_cons = struct();
To = struct();
for k=1:length(policy_types)
    %create beta vector
    betV = zeros(length(policy.(policy_types(k)).YV_eps),1);
    for m = 1:length(policy.(policy_types(k)).YV_eps)
        betV(m) = policy.(policy_types(k)).bet^(m-1);
    end
    %load variables
    lam = policy.(policy_types(k)).lam;
    lamS = policy.(policy_types(k)).lamS;
    tau = policy.(policy_types(k)).tau;
    delt = policy.(policy_types(k)).delt;
    alph = policy.(policy_types(k)).alph;
    CkSB = policy.(policy_types(k)).CkSB;
    phi = policy.(policy_types(k)).phi;
    gam = policy.(policy_types(k)).gam;

    Wpi = policy.(policy_types(k)).Wpi;
    Wy = policy.(policy_types(k)).Wy;
    Wdc = policy.(policy_types(k)).Wdc;

    PiV = policy.(policy_types(k)).PiV_eps;
    YV = policy.(policy_types(k)).YV_eps;
    YVS = policy.(policy_types(k)).YVS_eps;
    CrV = policy.(policy_types(k)).CrV_eps;
    CkV = policy.(policy_types(k)).CkV_eps;
    AV = policy.(policy_types(k)).AV_eps;

    %calc quadratic loss
    L.(policy_types(k)) = sum(betV*(-0.5).*(Wpi*PiV.^2 + Wy*(YV - YVS).^2 + Wdc*(CrV - CkV).^2));
    L_inf.(policy_types(k)) =  sum(betV*(-0.5).*(Wpi*PiV.^2));
    L_out.(policy_types(k)) =  sum(betV*(-0.5).*(Wy*(YV - YVS).^2));
    L_cons.(policy_types(k)) =  sum(betV*(-0.5).*(Wdc*(CrV - CkV).^2));

    %calc T_0
    To.(policy_types(k)) = (lam-lamS)/(1-lamS) * (1-(1-tau)*delt)*(1-alph)/CkSB * betV.*(((1+phi)/(1-alph).*(YV-AV) + (0.5 + lamS/(1-lamS)*(1-(1-tau)*delt)*(1-alph)/CkSB) * ((1+phi)/(1-alph).*(YV-AV)).^2 - lamS/(1-lamS)*gam*(1-alph)/CkSB.*AV*(1+phi)/(1-alph).*(YV-AV)));
end


%calc consumption equivalent loss compared to phi_C=0
if ismember("zero", policy_types)
    CEL = struct();
    CEL_inf = struct();
    CEL_out = struct();
    CEL_cons = struct();
    CEW = struct();
    CEW_notip = struct();
    for k=1:length(policy_types)
        %load variables
        bet = policy.(policy_types(k)).bet;
        
        L_x = L.(policy_types(k));
        L_x_inf = L_inf.(policy_types(k));
        L_x_out = L_out.(policy_types(k));
        L_x_cons = L_cons.(policy_types(k));

        L_zero = L.zero;
        L_zero_inf = L_inf.zero;
        L_zero_out = L_out.zero;
        L_zero_cons = L_cons.zero;
        
        %calc consumption equivalence
        CEL.(policy_types(k)) =  exp((1-bet)*(L_zero-L_x))-1;
        CEL_inf.(policy_types(k)) = exp((1-bet)*(L_zero_inf-L_x_inf))-1;
        CEL_out.(policy_types(k)) = exp((1-bet)*(L_zero_out-L_x_out))-1;
        CEL_cons.(policy_types(k)) = exp((1-bet)*(L_zero_cons-L_x_cons))-1;
    end
end



%Plots

%%plot config
x0=10;
y0=10;
width = 400;
height = 400;
x = -1:length(policy.(policy_types(1)).PiV_eps)+1-2;
xlimV = [-2 10];

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
for k=1:length(policy_types)
    plot(x, vertcat(0, policy.(policy_types(k)).PiV_eps)*100);
end
hold off
xlim(xlimV);
xlabel('Quarters');
ylabel('\pi');
%title('Inflation Gap');
legend(policy_types,'FontSize',10);

figure
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
for k=1:length(policy_types)
    plot(x, vertcat(0, policy.(policy_types(k)).dYV_eps)*100);
end
hold off
xlim(xlimV);
xlabel('Quarters');
ylabel('y-y^{*}');
%title('Output Gap');
legend(policy_types,'FontSize',10);

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
for k=1:length(policy_types)
    plot(x, vertcat(0, (policy.(policy_types(k)).CrV_eps - policy.(policy_types(k)).CkV_eps)*100));
end
hold off
xlim(xlimV);
xlabel('Quarters');
ylabel('c_{r}-c_{k}');
%title('Consumption Gap');
legend(policy_types,'FontSize',10);

%clear vars
clearvars -except policy_types policy L L_inf L_out L_cons CEL CEL_inf CEL_out CEL_cons To