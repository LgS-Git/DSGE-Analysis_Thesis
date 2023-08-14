
% Define which CB_regimes to run
CB_regimes_types = ["ric", "key", "avg", "optim"];

% Toogle on which parameter to loop (only one at a time)
gam_On = 1;
lam_On = 0;
tau_On = 0;

% Adjust range of parameter loops
gammas = 0:0.05:2.5;
lambdas = 0:0.025:1;
taus = 0.5:0.01:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calc and save CB_regimes
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
for k=1:length(CB_regimes_types)
    %loop over param
    first_time = 1;
    results_cell = cell(1,length(paramS));
    for i = 0:length(paramS)
        if first_time
            dynare_command = 'dynare Hansen2020_NotLinearized_Ramsey.mod -DCB_objective="' + CB_regimes_types(k) + '"';
            eval(dynare_command);
            first_time=0;
        else
            set_param_value(param,paramS(i));
            try
                perfect_foresight_setup;
                perfect_foresight_solver;
                if oo_.deterministic_simulation.status
                    results_cell{i} = oo_;
                else
                    results_cell{i} = NaN;
                    disp("Computation fails for " + param + " = " + string(paramS(i)));
                end
            catch
                results_cell{i} = NaN;
                disp("Computation fails for " + param + " = " + string(paramS(i)));
            end
        end
    end
    fname = "Ramsey_" + CB_regimes_types(k);
    save(fname);
    clearvars -except CB_regimes_types gammas lambdas taus paramS param gam_On lam_On tau_On;
end

%load CB_regimes
CB_regimes = struct();
for k=1:length(CB_regimes_types)
    fname = "Ramsey_" + CB_regimes_types(k);
    CB_regimes.(CB_regimes_types(k)) = load(fname);
    del_command = "delete " + fname + ".mat";
    eval(del_command);
end

%calculate loss function
L = struct();
L_inf = struct();
L_out = struct();
L_cons = struct();
for k=1:length(CB_regimes_types)
    %create beta vector
    betV = zeros(length(CB_regimes.(CB_regimes_types(k)).Y)-2,1);
    for m = 1:length(CB_regimes.(CB_regimes_types(k)).Y)-2
        betV(m) = CB_regimes.(CB_regimes_types(k)).bet^(m-1);
    end

    %load variables
    Wpi = CB_regimes.(CB_regimes_types(k)).Wpi;
    Wy = CB_regimes.(CB_regimes_types(k)).Wy;
    Wdc = CB_regimes.(CB_regimes_types(k)).Wdc;
    
    L_perParam = zeros(1,length(CB_regimes.(CB_regimes_types(k)).results_cell));
    L_perParam_inf = zeros(1,length(CB_regimes.(CB_regimes_types(k)).results_cell));
    L_perParam_out = zeros(1,length(CB_regimes.(CB_regimes_types(k)).results_cell));
    L_perParam_cons = zeros(1,length(CB_regimes.(CB_regimes_types(k)).results_cell));
    
    for i=1:length(CB_regimes.(CB_regimes_types(k)).results_cell)
        %load variables per param
        if isa(CB_regimes.(CB_regimes_types(k)).results_cell{1,i}, 'struct')
            CkV = CB_regimes.(CB_regimes_types(k)).results_cell{1,i}.endo_simul(14,2:end-1).';
            CrV = CB_regimes.(CB_regimes_types(k)).results_cell{1,i}.endo_simul(15,2:end-1).';
            YV = CB_regimes.(CB_regimes_types(k)).results_cell{1,i}.endo_simul(16,2:end-1).';
            PiV = CB_regimes.(CB_regimes_types(k)).results_cell{1,i}.endo_simul(17,2:end-1).';
            YVS = CB_regimes.(CB_regimes_types(k)).results_cell{1,i}.endo_simul(18,2:end-1).';
        else
            CkV = nan(length(CB_regimes.(CB_regimes_types(k)).Y(2:end-1)),1);
            CrV = nan(length(CB_regimes.(CB_regimes_types(k)).Y(2:end-1)),1);
            YV = nan(length(CB_regimes.(CB_regimes_types(k)).Y(2:end-1)),1);
            PiV = nan(length(CB_regimes.(CB_regimes_types(k)).Y(2:end-1)),1);
            YVS = nan(length(CB_regimes.(CB_regimes_types(k)).Y(2:end-1)),1);
        end
    
        %calc quadratic loss for each param
        L_perParam(1,i) =  sum(betV*(-0.5).*(Wpi*PiV.^2 + Wy*(YV - YVS).^2 + Wdc*(CrV - CkV).^2));
        L_perParam_inf(1,i) =  sum(betV*(-0.5).*(Wpi*PiV.^2));
        L_perParam_out(1,i) =  sum(betV*(-0.5).*(Wy*(YV - YVS).^2));
        L_perParam_cons(1,i) =  sum(betV*(-0.5).*(Wdc*(CrV - CkV).^2));
    end
    L.(CB_regimes_types(k)) = L_perParam;
    L_inf.(CB_regimes_types(k)) = L_perParam_inf;
    L_out.(CB_regimes_types(k)) = L_perParam_out;
    L_cons.(CB_regimes_types(k)) = L_perParam_cons;
end

%calc consumption equivalent loss per param, compared to optim
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

        %calc CEL's
        CEL.(CB_regimes_types(k)) =  exp((1-bet)*(L_optim-L_x))-1;
        CEL_inf.(CB_regimes_types(k)) = exp((1-bet)*(L_optim_inf-L_x_inf))-1;
        CEL_out.(CB_regimes_types(k)) = exp((1-bet)*(L_optim_out-L_x_out))-1;
        CEL_cons.(CB_regimes_types(k)) = exp((1-bet)*(L_optim_cons-L_x_cons))-1;
    end
end
    
%Plots

%%Plot config
x0 = 10;
y0 = 10;
width = 400;
height = 400;
x = linspace(paramS(1), paramS(end), length(paramS));

if ismember("optim", CB_regimes_types) && gam_On
    %exclude "optim"
    idx = CB_regimes_types ~= "optim";
    CB_regimes_types = CB_regimes_types(idx);
    
    for k = 1:length(CB_regimes_types)
        % construct Fill space variables
        colors = {'r','g','b'};
        xFill = [x, fliplr(x)];
        
        Y = [CEL_inf.(CB_regimes_types(k)); CEL_out.(CB_regimes_types(k)); CEL_cons.(CB_regimes_types(k))];
        [csort,idx] = sort(Y(:,40),'descend');
        y = Y(idx,:);
        % Create Y coordinates for fill:
        yFillPos = y(csort > 0,:);
        yFillNeg = y(csort <= 0,:);
        yFill = [y, fliplr([yFillPos(2:end,:); zeros(~isempty(yFillNeg) + ~isempty(yFillPos),length(x)); yFillNeg(1:end-1,:)])];
    
        %draw figure
        figure;
        set(gcf, 'position', [x0, y0, width, height]);
        hold on;
        cel_handle = plot(x, CEL.(CB_regimes_types(k)), 'k', 'Linewidth', 1.5);
        fill_handles = gobjects(1, length(CB_regimes_types));
        for i = 1:length(CB_regimes_types)
            plot(x, y(i,:), colors{idx(i)});
            % Fill between curves:
            fill_handles(idx(i)) = fill(xFill, yFill(i,:), colors{idx(i)}, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
        end
        %set x-axis
        xticks(gammas(1):0.5:gammas(end));
        xlabel('\gamma')
        
        %lines
        xline(1.67,"--");
        yline(0);
    
        %create legend
        legend_labels = {'CEL', 'Inflation Gap', 'Output Gap', 'Inequality Gap'};
        patch_handles = gobjects(1, 3);
        for i = 1:length(CB_regimes_types)
            patch_handles(i) = patch(NaN, NaN, colors{i}, 'FaceAlpha', 0.1);
        end
        legend([cel_handle, patch_handles], legend_labels);
    end

elseif ismember("optim", CB_regimes_types) && lam_On
    %exclude "optim"
    idx = CB_regimes_types ~= "optim";
    CB_regimes_types = CB_regimes_types(idx);
    
    %creates a reference line due to small differences
    reference_line = zeros(size(CEL.(CB_regimes_types(1))));
    for k = 1:length(CB_regimes_types)
        reference_line = reference_line + CEL.(CB_regimes_types(k));
    end
    reference_line = reference_line / length(CB_regimes_types);
    
    %calculate the differences between each line and the reference line
    for k = 1:length(CB_regimes_types)
        CEL_diff.(CB_regimes_types(k)) = CEL.(CB_regimes_types(k)) - reference_line;
    end

    %plot
    figure;
    set(gcf, 'position', [x0, y0, width, height]);
    hold on;
    for k = 1:length(CB_regimes_types)
        plot(x, CEL_diff.(CB_regimes_types(k)));
    end
    hold off
    %set x-axis
    xticks(lambdas(1):0.1:lambdas(end));
    xlabel('\lambda')
    legend(CB_regimes_types);
    
    %plot the differences between key and avg
    if ismember("key", CB_regimes_types) && ismember("avg", CB_regimes_types)    

        close_line_reference = (CEL.key + CEL.avg) / 2;

        %calculate the differences between each close line and the close_line_reference
        CEL_diff_close_key = CEL.key - close_line_reference;
        CEL_diff_close_avg = CEL.avg - close_line_reference;
        
        figure;
        set(gcf, 'position', [x0, y0, width, height]);
        hold on;
        plot(x, CEL_diff_close_key, 'color', [0.8500 0.3250 0.0980]);
        plot(x, CEL_diff_close_avg, 'color', [0.9290 0.6940 0.1250]);
        hold off;
        %set x-axis
        xticks(lambdas(1):0.1:lambdas(end));
        xlabel('\lambda');
        legend('key', 'avg');
    end

elseif ismember("optim", CB_regimes_types) && tau_On
    %exclude "optim"
    idx = CB_regimes_types ~= "optim";
    CB_regimes_types = CB_regimes_types(idx);
    
    reference_line = zeros(size(CEL.(CB_regimes_types(1))));
    for k = 1:length(CB_regimes_types)
        reference_line = reference_line + CEL.(CB_regimes_types(k));
    end
    reference_line = reference_line / length(CB_regimes_types);
    
    %calculate the differences between each line and the reference line
    for k = 1:length(CB_regimes_types)
        CEL_diff.(CB_regimes_types(k)) = CEL.(CB_regimes_types(k)) - reference_line;
    end

    %plot
    figure;
    set(gcf, 'position', [x0, y0, width, height]);
    hold on;
    for k = 1:length(CB_regimes_types)
        plot(x, CEL_diff.(CB_regimes_types(k)));
    end
    hold off
    %set x-axis
    xticks(taus(1):0.1:taus(end));
    xlabel('\tau')
    legend(CB_regimes_types);

    %plot the differences between ric and avg
    if ismember("ric", CB_regimes_types) && ismember("avg", CB_regimes_types)    

        close_line_reference = (CEL.ric + CEL.avg) / 2;

        %calculate the differences between each close line and the close_line_reference
        CEL_diff_close_ric = CEL.ric - close_line_reference;
        CEL_diff_close_avg = CEL.avg - close_line_reference;
        
        figure;
        set(gcf, 'position', [x0, y0, width, height]);
        hold on;
        plot(x, CEL_diff_close_ric);
        plot(x, CEL_diff_close_avg, 'color', [0.9290 0.6940 0.1250]);
        hold off;
        %set x-axis
        xticks(lambdas(1):0.1:lambdas(end));
        xlabel('\tau');
        legend('ric', 'avg');
    end
end

%clearvars
clearvars -except CB_regimes L L_inf L_out L_cons CEL CEL_inf CEL_out CEL_cons