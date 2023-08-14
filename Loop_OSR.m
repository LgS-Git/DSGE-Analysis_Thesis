
% Define different taylor rule types
policy_types = ["ric", "key", "ineq", "zero"];
%policy_types = ["ric", "key", "ineq"];

% Toggle on for policy specific OSR weights
specific_OSR_weights = 1;

% Toogle on which parameter to loop (only one at a time)
gam_On = 0;
lam_On = 0;
tau_On = 1;

% Adjust range of param loop
gammas = 1:0.025:2;
lambdas = 0:0.025:1;
taus = 0:0.025:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc and save policy
if gam_On
        paramS = gammas;
        param = 'gam';
elseif lam_On %#ok<*UNRCH> 
        paramS = lambdas;
        param = 'lam';
elseif tau_On
        paramS = taus;
        param = 'tau';
end
for k=1:length(policy_types)
first_time = 1;
results_cell = cell(1,length(paramS));
    for i = 0:length(paramS)
        if first_time
            dynare_command = 'dynare Corr_Hansen2020_Linearized_OSR.mod -Dpolicy="'+policy_types(k)+'" -Dspecific_OSR_weights='+specific_OSR_weights + ' -Dsense=0';
            eval(dynare_command);
            first_time=0;
        else
            set_param_value(param,paramS(i));
            try
                if policy_types(k) == "zero"
                    [info, oo_] = stoch_simul(M_, options_, oo_, var_list_);
                else
                    osr(M_.endo_names,M_.osr.param_names,M_.osr.variable_indices,M_.osr.variable_weights);
                end
                results_cell{i} = oo_;
            catch
                results_cell{i} = NaN;
                disp("Computation fails for " + param + " = " + string(paramS(i)));
            end
        end
    end
    fname = "OSR_" + policy_types(k);
    save(fname);
    clearvars -except policy_types gammas lambdas taus paramS param gam_On tau_On lam_On specific_OSR_weights;
end



%load policy
policy = struct();
for k=1:length(policy_types)
    fname = "OSR_" + policy_types(k);
    policy.(policy_types(k)) = load(fname);
    del_command = "delete " + fname + ".mat";
    eval(del_command);
end


%calculate loss function
L = struct();
L_inf = struct();
L_out = struct();
L_cons = struct();
for k=1:length(policy_types)
    %create beta vector
    betV = zeros(length(policy.(policy_types(k)).YV_eps),1);
    for m = 1:length(policy.(policy_types(k)).YV_eps)
        betV(m) = policy.(policy_types(k)).bet^(m-1);
    end

    %load variables
    Wpi = policy.(policy_types(k)).Wpi;
    Wy = policy.(policy_types(k)).Wy;
    Wdc = policy.(policy_types(k)).Wdc;
    
    L_perParam = zeros(1,length(policy.(policy_types(k)).results_cell));
    L_perParam_inf = zeros(1,length(policy.(policy_types(k)).results_cell));
    L_perParam_out = zeros(1,length(policy.(policy_types(k)).results_cell));
    L_perParam_cons = zeros(1,length(policy.(policy_types(k)).results_cell));
    
    for i=1:length(policy.(policy_types(k)).results_cell)
        %load variables per gamma
        if isa(policy.(policy_types(k)).results_cell{1,i}, 'struct')
            CkV = policy.(policy_types(k)).results_cell{1,i}.irfs.CkV_eps';
            CrV = policy.(policy_types(k)).results_cell{1,i}.irfs.CrV_eps';
            YV = policy.(policy_types(k)).results_cell{1,i}.irfs.YV_eps';
            PiV = policy.(policy_types(k)).results_cell{1,i}.irfs.PiV_eps';
            YVS = policy.(policy_types(k)).results_cell{1,i}.irfs.YVS_eps';
        else
            CkV = nan(length(policy.(policy_types(k)).YV_eps),1);
            CrV = nan(length(policy.(policy_types(k)).YV_eps),1);
            YV = nan(length(policy.(policy_types(k)).YV_eps),1);
            PiV = nan(length(policy.(policy_types(k)).YV_eps),1);
            YVS = nan(length(policy.(policy_types(k)).YV_eps),1);
        end
    
        %calc quadratic loss for each gamma
        L_perParam(1,i) =  sum(betV*(-0.5).*(Wpi*PiV.^2 + Wy*(YV - YVS).^2 + Wdc*(CrV - CkV).^2));
        L_perParam_inf(1,i) =  sum(betV*(-0.5).*(Wpi*PiV.^2));
        L_perParam_out(1,i) =  sum(betV*(-0.5).*(Wy*(YV - YVS).^2));
        L_perParam_cons(1,i) =  sum(betV*(-0.5).*(Wdc*(CrV - CkV).^2));
    end
    L.(policy_types(k)) = L_perParam;
    L_inf.(policy_types(k)) = L_perParam_inf;
    L_out.(policy_types(k)) = L_perParam_out;
    L_cons.(policy_types(k)) = L_perParam_cons;
end

if ismember("zero", policy_types)
    CEL = struct();
    CEL_inf = struct();
    CEL_out = struct();
    CEL_cons = struct();
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

%%Plot config
x0 = 10;
y0 = 10;
width = 400;
height = 400;
x = linspace(paramS(1), paramS(end), length(paramS));

if ismember("zero", policy_types) && lam_On
    %exclude "zero"
    idx = policy_types ~= "zero";
    policy_types = policy_types(idx);
    
    %creates a reference line due to small differences
    reference_line = zeros(size(CEL.(policy_types(1))));
    for k = 1:length(policy_types)
        reference_line = reference_line + CEL.(policy_types(k));
    end
    reference_line = reference_line / length(policy_types);
    
    % Calculate the differences between each line and the reference line
    for k = 1:length(policy_types)
        CEL_diff.(policy_types(k)) = CEL.(policy_types(k)) - reference_line;
    end

    %Plot
    figure;
    set(gcf, 'position', [x0, y0, width, height]);
    hold on;
    for k = 1:length(policy_types)
        %figure;
        %set(gcf, 'position', [x0, y0, width, height]);
        plot(x, CEL_diff.(policy_types(k)));
        %legend(policy_types(k));
    end
    hold off
    %set x-axis
    xticks(lambdas(1):0.1:lambdas(end));
    xlabel('\lambda')
    legend(policy_types);

elseif ismember("zero", policy_types) && tau_On
    %exclude "zero"
    idx = policy_types ~= "zero";
    policy_types = policy_types(idx);
    
    %creates a reference line due to small differences
    reference_line = zeros(size(CEL.(policy_types(1))));
    for k = 1:length(policy_types)
        reference_line = reference_line + CEL.(policy_types(k));
    end
    reference_line = reference_line / length(policy_types);
    
    % Calculate the differences between each line and the reference line
    for k = 1:length(policy_types)
        CEL_diff.(policy_types(k)) = CEL.(policy_types(k)) - reference_line;
    end

    %Plot
    %plot
    figure;
    set(gcf, 'position', [x0, y0, width, height]);
    hold on;
    for k = 1:length(policy_types)
        plot(x, CEL_diff.(policy_types(k)));
    end
    hold off
    %set x-axis
    xticks(taus(1):0.1:taus(end));
    xlabel('\tau')
    legend(policy_types);

elseif ismember("zero", policy_types) && gam_On
    %exclude "zero"
    idx = policy_types ~= "zero";
    policy_types = policy_types(idx);
    
    %creates a reference line due to small differences
    reference_line = zeros(size(CEL.(policy_types(1))));
    for k = 1:length(policy_types)
        reference_line = reference_line + CEL.(policy_types(k));
    end
    reference_line = reference_line / length(policy_types);
    
    % Calculate the differences between each line and the reference line
    for k = 1:length(policy_types)
        CEL_diff.(policy_types(k)) = CEL.(policy_types(k)) - reference_line;
    end

    %Plot
    %plot
    figure;
    set(gcf, 'position', [x0, y0, width, height]);
    hold on;
    for k = 1:length(policy_types)
        plot(x, CEL_diff.(policy_types(k)));
    end
    hold off
    %set x-axis
    xticks(gammas(1):0.1:gammas(end));
    xlabel('\gamma')
    legend(policy_types);
end

%clearvars
clearvars -except policy L L_inf L_out L_cons CEL CEL_inf CEL_out CEL_cons