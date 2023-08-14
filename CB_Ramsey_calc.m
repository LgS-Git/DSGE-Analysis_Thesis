
% Define what regimes to run
CB_regimes_types = ["ric", "key", "avg", "optim"];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc and save CB_regimes
for k=1:length(CB_regimes_types)
    dynare_command = sprintf('dynare Hansen2020_NotLinearized_Ramsey.mod -DCB_objective="%s"', CB_regimes_types(k));
    eval(dynare_command);
    fname = "Ramsey_" + CB_regimes_types(k);
    save(fname);
    clearvars -except CB_regimes_types;
end

%Load CB_regimes
CB_regimes = struct();
for k=1:length(CB_regimes_types)
    fname = "Ramsey_" + CB_regimes_types(k);
    CB_regimes.(CB_regimes_types(k)) = load(fname);
    del_command = "delete " + fname + ".mat";
    eval(del_command);
end

%calculate loss function for regime types
L = struct();
L_inf = struct();
L_out = struct();
L_cons = struct();
To = struct();
for k=1:length(CB_regimes_types)
    %create beta vector
    betV = zeros(length(CB_regimes.(CB_regimes_types(k)).Y)-2,1);
    for m = 1:length(CB_regimes.(CB_regimes_types(k)).Y)-2
        betV(m) = CB_regimes.(CB_regimes_types(k)).bet^(m-1);
    end
    %load variables
    lam = CB_regimes.(CB_regimes_types(k)).lam;
    lamS = CB_regimes.(CB_regimes_types(k)).lamS;
    tau = CB_regimes.(CB_regimes_types(k)).tau;
    delt = CB_regimes.(CB_regimes_types(k)).delt;
    alph = CB_regimes.(CB_regimes_types(k)).alph;
    CkSB = CB_regimes.(CB_regimes_types(k)).CkSB;
    phi = CB_regimes.(CB_regimes_types(k)).phi;
    gam = CB_regimes.(CB_regimes_types(k)).gam;

    Wpi = CB_regimes.(CB_regimes_types(k)).Wpi;
    Wy = CB_regimes.(CB_regimes_types(k)).Wy;
    Wdc = CB_regimes.(CB_regimes_types(k)).Wdc;

    PiV = CB_regimes.(CB_regimes_types(k)).log_lin_Pi(2:end-1);
    YV = CB_regimes.(CB_regimes_types(k)).log_lin_Y(2:end-1);
    YVS = CB_regimes.(CB_regimes_types(k)).log_lin_YS(2:end-1);
    CrV = CB_regimes.(CB_regimes_types(k)).log_lin_Cr(2:end-1);
    CkV = CB_regimes.(CB_regimes_types(k)).log_lin_Ck(2:end-1);
    AV = CB_regimes.(CB_regimes_types(k)).log_lin_A(2:end-1);

    %calc quadratic loss
    L.(CB_regimes_types(k)) = sum(betV*(-0.5).*(Wpi*PiV.^2 + Wy*(YV - YVS).^2 + Wdc*(CrV - CkV).^2));
    L_inf.(CB_regimes_types(k)) =  sum(betV*(-0.5).*(Wpi*PiV.^2));
    L_out.(CB_regimes_types(k)) =  sum(betV*(-0.5).*(Wy*(YV - YVS).^2));
    L_cons.(CB_regimes_types(k)) =  sum(betV*(-0.5).*(Wdc*(CrV - CkV).^2));

    %calc T_0
    To.(CB_regimes_types(k)) = (lam-lamS)/(1-lamS) * (1-(1-tau)*delt)*(1-alph)/CkSB * betV.*(((1+phi)/(1-alph).*(YV-AV) + (0.5 + lamS/(1-lamS)*(1-(1-tau)*delt)*(1-alph)/CkSB) * ((1+phi)/(1-alph).*(YV-AV)).^2 - lamS/(1-lamS)*gam*(1-alph)/CkSB.*AV*(1+phi)/(1-alph).*(YV-AV)));
end

%calc consumption equivalent loss, compared to optim
if ismember("optim", CB_regimes_types)
    CEL = struct();
    CEL_inf = struct();
    CEL_out = struct();
    CEL_cons = struct();
    for k=1:length(CB_regimes_types)
        %load variables
        bet = CB_regimes.(CB_regimes_types(k)).bet;

        L_x = L.(CB_regimes_types(k));
        L_x_inf = L_inf.(CB_regimes_types(k));
        L_x_out = L_out.(CB_regimes_types(k));
        L_x_cons = L_cons.(CB_regimes_types(k));

        L_optim = L.optim;
        L_optim_inf = L_inf.optim;
        L_optim_out = L_out.optim;
        L_optim_cons = L_cons.optim;

        %calc CEL
        CEL.(CB_regimes_types(k)) =  exp((1-bet)*(L_optim-L_x))-1;
        CEL_inf.(CB_regimes_types(k)) = exp((1-bet)*(L_optim_inf-L_x_inf))-1;
        CEL_out.(CB_regimes_types(k)) = exp((1-bet)*(L_optim_out-L_x_out))-1;
        CEL_cons.(CB_regimes_types(k)) = exp((1-bet)*(L_optim_cons-L_x_cons))-1;
    end
end


%Plots

%%plot config
x0=10;
y0=10;
width = 400;
height = 400;
x = -1:length(CB_regimes.(CB_regimes_types(1)).log_lin_Pi*100)-2;

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
for k=1:length(CB_regimes_types)
    plot(x, CB_regimes.(CB_regimes_types(k)).log_lin_Pi*100);
end
hold off
xlim([-2 10]);
xlabel('Quarters');
ylabel('\pi');
%title('Inflation Gap');
legend(CB_regimes_types,'FontSize',10);

figure
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
for k=1:length(CB_regimes_types)
    plot(x, CB_regimes.(CB_regimes_types(k)).log_lin_Y*100 - CB_regimes.(CB_regimes_types(k)).log_lin_YS*100);
end
hold off
xlim([-2 10]);
xlabel('Quarters');
ylabel('y-y^{*}');
%title('Output Gap');
legend(CB_regimes_types,'FontSize',10);

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
for k=1:length(CB_regimes_types)
    plot(x, CB_regimes.(CB_regimes_types(k)).log_lin_Cr*100 - CB_regimes.(CB_regimes_types(k)).log_lin_Ck*100);
end
hold off
xlim([-2 10]);
xlabel('Quarters');
ylabel('c_{r}-c_{k}');
%title('Consumption Gap');
legend(CB_regimes_types,'FontSize',10);

%clear variables
clearvars -except CB_regimes CB_regimes_types L L_inf L_out L_cons CEL CEL_inf CEL_out CEL_cons To