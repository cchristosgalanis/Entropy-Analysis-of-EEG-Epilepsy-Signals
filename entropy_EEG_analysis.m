% Christos Galanis AEM:2111 

%main function for project
function entropy_EEG_analysis1
    close all;
    clc;

    %trying to load data from file
    try
        %loading file 
        load BonnEEGdata.mat
        disp("Data has been loaded correctly") %display file 
    catch
        error("Error while loading BonnEGGdata.mat file. Proccess has been terminated");
    end


    fprintf("\n--- Compute Entropy ---")
    %using Shannon's equation to determine the entropy 
    try
        %calculating entropy for category Z, category:-> healthy, open eyes
        entropyZ = EEG_Tools.ShEn(Z_EEG_data);
        
        %calculating entropy for category O, category:-> healthy, cloded
        %eyes
        entropyO = EEG_Tools.ShEn(O_EEG_data);

        %calculating entropy for category N, category:-> epileptic, healthy
        %hemishpere
        entropyN = EEG_Tools.ShEn(N_EEG_data);

        %calculating entropy for category F, category:->epileptic, epileptic
        %focus
        entropyF = EEG_Tools.ShEn(F_EEG_data);

        %calculating entropy for category S, category:->epileptic,seizure
        entropyS = EEG_Tools.ShEn(S_EEG_data);

    catch
        error("An occured while calculating entropy for categories. Process has been terminated.")
    end

    %call EnPlot for entropies
    EEG_Tools.EnPlot(entropyZ,entropyO,entropyN,entropyF,entropyS);

% ------ Non Linear Analysis -------
    
    fprintf('\n--- Starting Non Linear Analysis  ---\n');
    
    % max_tau 
    max_tau = 50;
    max_dim = 15;
    fs = 173.61; %change to correct frequency
    maxiter = 100; %steps forward to calculate divergence
    %meanperiod = 50; %using optimal tau (Theiler window)
    tlinear = 3:8; %or 2:8 , being hardcored but will be hybrid
    
    categories = {Z_EEG_data, O_EEG_data, N_EEG_data, F_EEG_data, S_EEG_data};
    names = {'Z (Healthy Open)', 'O (Healthy Closed)', 'N (Inter-Hemi)', 'F (Inter-Focus)', 'S (Ictal)'};
    colors = {'b', 'c', 'g', 'm', 'r'};
    
    %figure for AMI curve
    ami_fig = figure('Name', 'Average Mutual Information Analysis','Color','w');
    hold on; grid on;
    title('Average Mutual Information (AMI)');
    xlabel('Time Delay (\tau)');
    ylabel('AMI (bits)');

    %figure for FNN curve
    fnn_fig = figure('Name',"FNN Analysis",'Color','w');
    title('FNN Fractiono vs Embedding Dimension');
    xlabel("Embedding Dimension");
    ylabel("FNN (%)");
    grid on; hold on;

    %figure for LLE_divergence curve
    lle_fig = figure('Name',"LLE Divergence Analysis",'Color','w');
    hold on; grid on;
    title('Average Logarithmic Divergence (LLE)');
    xlabel('Time Steps (k)');
    ylabel('Divergence d(k)');
    
    % Loop through each category
    for i = 1:length(categories)
        current_data = categories{i};
        current_name = names{i};
        current_color = colors{i};
        
        fprintf('Processing Category: %s ...\n', current_name);
        
        % ---- Calculation/Plot for AMI/optimal tau ----
        %calculate ami matrix
        ami_matrix = EEG_Tools.Compute_AMI_Matrix(current_data, max_tau);
        
        %finding optimal tau for each category
        tau_vector = EEG_Tools.Find_Tau_From_Matrix(ami_matrix);
        
        fprintf(' Average Optimal Tau: %.2f\n', mean(tau_vector));
        
        %plot for curve
        figure(ami_fig);
        EEG_Tools.Plot_AMI_Curve(ami_matrix, current_name, current_color);

        % ---- Calculation/Plot for Embedding Dimension ----
        %finding optimal embedding dimension m
        [num_signal,~] = size(current_data);
        m_vector = zeros(num_signal,1);
        fnn_matrix = zeros(num_signal,max_dim);
        lle_vector = zeros(num_signal,1);
        d_matrix = zeros(num_signal,maxiter);

        for s = 1:num_signal
            signal = current_data(s,:);
            tau = tau_vector(s);

            [opt_m,fnn_curve] = EEG_Tools.embeddingDimensionFnn(signal,tau,max_dim);
            m_vector(s) = opt_m;
            fnn_matrix(s,:) = fnn_curve;

            d_curve = EEG_Tools.lyarosenstein(signal,opt_m,tau,tau,maxiter);
            d_matrix(s,:) = d_curve; %savine curve

            %compute LLE(slope)
            F = polyfit(tlinear,d_curve(tlinear),1);
            lle_vector(s) = F(1)*fs;
            
        end
        fprintf(" Average Embedding Dimension: %.2f\n", mean(m_vector));
        fprintf(" Average LLE: %.4f\n",mean(lle_vector));

        figure(fnn_fig);
        EEG_Tools.plotFNN(fnn_matrix,current_name,current_color);

        figure(lle_fig);
        EEG_Tools.plotLLE_Divergence(d_matrix,current_name,current_color);
    end
    
    %final for AMI
    figure(ami_fig);
    legend('Location', 'northeast');
    hold off;
    %saveas(gcf, 'AMI_Analysis_Plot.png');

    %final for FNN
    figure(fnn_fig);
    legend('Location','northeast');
    hold off;
    %saveas(gcf, 'FNN_Analysis_Plot.png');
    disp('Non Linear Analysis has been completed');

    %final for LLE
    figure(lle_fig);
    legend('Location', 'northeast');
    hold off;
    %saveas(gcf, 'LLE_DIvergence_Plot.png');
end

