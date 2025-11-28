% Christos Galanis AEM:2111 

%main function for project
function entropy_EEG_analysis
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

    %displaying calculated values for each category
    fprintf("Entropy Values:\n")
    fprintf("Calculated entropies for category Z are:")
    fprintf("%4f ",entropyZ)
    fprintf("\n")
    fprintf("\n Calculated entropies for category O are:")
    fprintf("%4f ",entropyO)
    fprintf("\n")
    fprintf("\n Calculated entropies for category N are:")
    fprintf("%4f ",entropyN)
    fprintf("\n")
    fprintf("\n Calculated entropies for category F are:")
    fprintf("%4f ", entropyF)
    fprintf("\n")
    fprintf("\n Calculated entropies for category S are:")
    fprintf("%4f ", entropyS)
    fprintf("\n")

    %call EnPlot for entropies
    EnPlot(entropyZ,entropyO,entropyN,entropyF,entropyS);

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
    colours = {'b', 'k', 'g', 'm', 'r'}; %colors for each category
    
    %creating figure
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
end