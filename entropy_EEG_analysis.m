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
        entropyZ = ShEn(Z_EEG_data);
        
        %calculating entropy for category O, category:-> healthy, cloded
        %eyes
        entropyO = ShEn(O_EEG_data);

        %calculating entropy for category N, category:-> epileptic, healthy
        %hemishpere
        entropyN = ShEn(N_EEG_data);

        %calculating entropy for category F, category:->epileptic, epileptic
        %focus
        entropyF = ShEn(F_EEG_data);

        %calculating entropy for category S, category:->epileptic,seizure
        entropyS = ShEn(S_EEG_data);

    catch
        error("An occured while calculating entropy for categories. Process has been terminated.")
    end

    %call EnPlot for entropies
    %EnPlot(entropyZ,entropyO,entropyN,entropyF,entropyS);

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
        ami_matrix = Compute_AMI_Matrix(current_data, max_tau);
        
        %finding optimal tau for each category
        tau_vector = Find_Tau_From_Matrix(ami_matrix);
        
        fprintf(' Average Optimal Tau: %.2f\n', mean(tau_vector));
        
        %plot for curve
        figure(ami_fig);
        Plot_AMI_Curve(ami_matrix, current_name, current_color);

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

            [opt_m,fnn_curve] = embeddingDimensionFnn(signal,tau,max_dim);
            m_vector(s) = opt_m;
            fnn_matrix(s,:) = fnn_curve;

            d_curve = lyarosenstein(signal,opt_m,tau,tau,maxiter);
            d_matrix(s,:) = d_curve; %savine curve

            %compute LLE(slope)
            F = polyfit(tlinear,d_curve(tlinear),1);
            lle_vector(s) = F(1)*fs;
            
        end
        fprintf(" Average Embedding Dimension: %.2f\n", mean(m_vector));
        fprintf(" Average LLE: %.4f\n",mean(lle_vector));

        figure(fnn_fig);
        plotFNN(fnn_matrix,current_name,current_color);

        figure(lle_fig);
        plotLLE_Divergence(d_matrix,current_name,current_color);
    end
    
    %final for AMI
    % figure(ami_fig);
    % legend('Location', 'northeast');
    % hold off;
    %saveas(gcf, 'AMI_Analysis_Plot.png');

    %final for FNN
    % figure(fnn_fig);
    % legend('Location','northeast');
    % hold off;
    %saveas(gcf, 'FNN_Analysis_Plot.png');
    disp('Non Linear Analysis has been completed');

    %final for LLE
    %final for AMI
    figure(lle_fig);
    legend('Location', 'northeast');
    hold off;
    %saveas(gcf, 'LLE_DIvergence_Plot.png');
end

%function for calculate entropy for each category
function en_res = ShEn(EGG_data)
    %firstly we have to quoantize datas in order to use shannon's equation
    %find numbers of signals in each category (make our code more hybrid for future use)
    samples = size(EGG_data,1);
    
    %fullfill an 1D array with zeros
    en_res = zeros(samples,1);
    
    %calculating histcount and each entropy
    for signal = 1:samples
        temp_signal = EGG_data(signal,:); %getting each signal from datas

        %calculate histcount for signal
        temp = histcounts(temp_signal,'Normalization','probability');
        % keep values that aren't 0, because this gonna lead to inf values
        % to Shannon's calculation
        temp = temp(temp > 0);
           
        %using Shannon's equation to calculate entropy for each signal
        en_res(signal) = -sum(temp.*log2(temp));
    end
end

%EnPlot function for creating diagram
function EnPlot(enZ,enO,enN,enF,enS)
    category_label = {enZ,enO,enN,enF,enS};
    %makes code more stable for future use
    categories = length(category_label);

    category_name = {'Z','O','N','F','S'}; %using each category for plot names
    colours = {'b', 'w', 'g', 'm', 'r'}; %colors for each category
    
    %creating figure for Shannon Entropy
    figure('Name','EEG Entropy Distribution');
    hold on;
    grid on;
    
    title("Entropy by Shannon's Equation");
    xlabel("Entropy's Value");
    ylabel("Probabiliy Density");

    %fullfill 1D arrays with zeros
    mean_values = zeros(categories,1);
    std_values = zeros(categories,1);
    
    %calculate mean and std values for each entropy category
    for i = 1:categories
        temp_mean_values = mean(category_label{i});
        temp_std_values = std(category_label{i});

        %saving temp values to arrays
        mean_values(i) = temp_mean_values;
        std_values(i) = temp_std_values;
        
        %creating range values for x
        x = linspace(temp_mean_values - 4*temp_std_values, temp_mean_values + 4*temp_std_values, 100);

        %calculate values for y
        y = normpdf(x,temp_mean_values,temp_std_values);
        
        %Plot each entropy value
        plot(x,y,'Color',colours{i},'LineWidth',2,'DisplayName',category_name{i});
    end

    for i = 1:categories
        category = category_name{i};
        mean_value = mean_values(i);
        std_value = std_values(i);
        fprintf("\n Mean Entropy Value for category %s is: %4f", category,mean_value)
        fprintf("\n Std Entropy Value for category %s is: %4f",category,std_value)
    end
    legend('Location','northeast')
    hold off;
    %saveas(gcf, 'Shannon_Entropy.png');
end

% ----- Functions for AMI calculation -----

% function to calculate AMI matrix
function ami_matrix = Compute_AMI_Matrix(signal_data, max_tau)
    % Get dimensions
    [num_signals, ~] = size(signal_data);
    n_bins = 32;
    
    ami_matrix = zeros(num_signals, max_tau);
    
    % Loop through each signal
    for s = 1:num_signals
        signal = signal_data(s, :);
        
        
        for t = 1:max_tau
            %shifted vectors
            x = signal(1:end-t);
            y = signal(t+1:end);
            
            % 2D Histogram
            N = histcounts2(x, y, n_bins);
            p_xy = N / sum(N(:));
            
            % Marginals
            p_x = sum(p_xy, 2);
            p_y = sum(p_xy, 1);
            
            %ami's calculation
            px_py = p_x * p_y;
            nzs = p_xy > 0;
            ami_matrix(s, t) = sum(p_xy(nzs) .* log2(p_xy(nzs) ./ px_py(nzs)));
        end
    end
end

% function to calculate optimal tau for each signal category
function optimal_tau_vec = Find_Tau_From_Matrix(ami_matrix)
    [num_signals, max_tau] = size(ami_matrix);
    optimal_tau_vec = zeros(num_signals, 1);
    
    for s = 1:num_signals
        % get the AMI curve for this specific signal
        curve = ami_matrix(s, :);
        
        % Find first local minimum
        found_tau = 1; 
        for t = 2:(max_tau - 1)
            if curve(t) < curve(t-1) && curve(t) < curve(t+1)
                found_tau = t;
                break; 
            end
        end
        optimal_tau_vec(s) = found_tau;
    end
end

% function to plot ami curve
function Plot_AMI_Curve(ami_matrix, name, color)
    % mean value for each category
    mean_curve = mean(ami_matrix);
    
    x_axis = 1:length(mean_curve);
    
    plot(x_axis, mean_curve, 'Color', color, 'LineWidth', 2, 'DisplayName', name);
end


%function to calculate embedding (m)
function [optimal_m, fnn_ratio] = embeddingDimensionFnn(signal, tau, max_dim)
    % Kennel's implementation parameters
    r_tol = 15;  % distance tolerance
    a_tol = 2;   % attractor tolerance

    N = length(signal);
    fnn_ratio = zeros(1, max_dim);
    sigma = std(signal); % standard deviation for normalization

    for m = 1:max_dim
        m_curr = N - m*tau;
        
        Y = zeros(m_curr, m);
        for i = 1:m
            
            Y(:, i) = signal(1 + (i-1)*tau : m_curr + (i-1)*tau);
        end

        % using knnsearch
        [idx, dist_matrix] = knnsearch(Y, Y, 'K', 2);

        nearest_idx = idx(:, 2);       % index of neighbor
        dist_in_m = dist_matrix(:, 2); % distance in dimension m

        false_count = 0;
        next_tau = m * tau;

        for k = 1:m_curr
            % value on next dimension (m+1)
            val_point_next = signal(k + next_tau);
            
            val_neighbor_next = signal(nearest_idx(k) + next_tau);

            % distance increase in next dimension
            dist_increase = abs(val_point_next - val_neighbor_next);

            % distance drops down of machine number
            if dist_in_m(k) < 1e-10
                r_i = 0; % Avoid division by zero
            else
                r_i = dist_increase / dist_in_m(k);
            end

            % get to attractor's tolerance size
            s_i = dist_increase / sigma;

            if (r_i > r_tol) || (s_i > a_tol)
                false_count = false_count + 1;
            end
        end
        
        % matrix of ratios
        fnn_ratio(m) = (false_count / m_curr) * 100; 
    end
    
    % find first m drops down to 1% 
    threshold = 1; % 1% threshold
    opt_m = find(fnn_ratio < threshold);
    
    if isempty(opt_m)
        % if never drops down the threshhold
        [~, optimal_m] = min(fnn_ratio);
    else
        optimal_m = opt_m(1);
    end
end

%function to plot fnn
function plotFNN(fnn_matrix,name,color)
    mean_curve = mean(fnn_matrix);

    x_axis = 1:length(mean_curve);

    plot(x_axis,mean_curve,'Color',color,'LineWidth',2,'DisplayName',name);
end

%Built in Rosenstein algorithm for calculating LLE(largest lyapunov exponent)
function d = lyarosenstein(x,m,tao,meanperiod,maxiter) 
% d:divergence of nearest trajectoires
% x:signal
% tao:time delay
% m:embedding dimension

    N=length(x);
    M=N-(m-1)*tao;
    Y=psr_deneme(x,m,tao);

    for i=1:M
        x0=ones(M,1)*Y(i,:);
        distance=sqrt(sum((Y-x0).^2,2));
        for j=1:M
            if abs(j-i)<=meanperiod
                distance(j)=1e10;
            end
        end
    [neardis(i) nearpos(i)]=min(distance);
    end

    for k=1:maxiter
        maxind=M-k;
        evolve=0;
        pnt=0;
        for j=1:M
            if j<=maxind && nearpos(j)<=maxind
                dist_k=sqrt(sum((Y(j+k,:)-Y(nearpos(j)+k,:)).^2,2));
                if dist_k~=0
                    evolve=evolve+log(dist_k);
                    pnt=pnt+1;
                end
            end
        end
        if pnt > 0
            d(k)=evolve/pnt;
        else
            d(k)=0;
        end
    end
end

function Y=psr_deneme(x,m,tao,npoint)
    N=length(x);
    if nargin == 4
    M=npoint;
    else
    M=N-(m-1)*tao;
    end

    Y=zeros(M,m); 

    for i=1:m
        Y(:,i)=x((1:M)+(i-1)*tao)';
    end
end

%function to plot divergence (Lyapunov Exponent curve)
function plotLLE_Divergence(d_matrix,current_name,current_color)
    mean_curve = mean(d_matrix,1);
    maxiter = length(mean_curve);
    plot(1:maxiter,mean_curve,'Color',current_color,'LineWidth',1.5,'DisplayName',current_name);
end