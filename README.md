# Entropy Analysis of EEG Epilepsy Signals

### 💡 Project Overview

* **Course:** Digital Signal Processing
* **Author:** Christos Galanis
* **Software/Version:** MatLab 2025b

The main objective of this project is the **quantitative distinction of five physiological brain states** (Z, O, N, S, F) using the value of **entropy** as a characteristic feature. The results from various entropy calculation approaches must confirm that entropy can be used as a **biomarker** for the diagnosis and monitoring of epilepsy. This specific work utilizes the **Shannon Equation** to calculate the entropy.

---

### 🧠 Theoretical Background

Entropy is a measure of **disorder**, **complexity**, or the **rate of information** within a signal. It allows us to quantify how unpredictable or chaotic an electrical signal recorded from the brain is.

* **High Entropy:** Associated with healthy, complex, and uncoordinated brain function, where independent neuronal groups are utilized.
* **Low Entropy:** Associated with signals that are organized, predictable, and synchronized.

#### Expected Results (Based on Neuroscience)

Based on Neurology, the value of entropy is expected to be **reduced during an epileptic seizure (S)** compared to a healthy state (Z or O).

* **Reason:** The extreme coordination (hypersynchronization) of neuronal activity that occurs during a seizure. This hypersynchronization makes the EEG signal highly organized and predictable, reducing complexity and leading to a smaller entropy value.
* **Normal State (Z or O):** Neuronal activity is asynchronous and chaotic, exhibiting high entropy.

#### Shannon Entropy Equation

Shannon Entropy ($H$) is the foundation of Information Theory, representing the **average uncertainty** contained in an information source (like EEG).

$$H(x)=-\sum_{i=1}^{n}p_{i}\log_{2}(p_{i})$$

---

### 💾 Data

The **Bonn EEG Dataset** is used for this work, a database widely used for comparing EEG analysis algorithms.

* **Source:** Department of Epileptology at the University of Bonn.
* **Signals:** 500 single-channel EEG signals, 100 from each category (state): **Z, O, N, F, and S**.
* **Duration/Sampling:** Each signal lasts 23.6 seconds, consisting of 4097 samples. The sampling rate was $173.61 \text{ Hz}$.

---

### 💻 Code Implementation Flow (MatLab 2025b)

The code is generally based on three main functions:
1.  `entropy_EEG_analysis` (Main function)
2.  `ShEn` (Shannon Entropy calculation)
3.  `EnPlot` (Graph visualization)

#### General Code Flow

1.  Load the database.
2.  Call the **`ShEn`** function, which returns the entropy values for each signal category.
3.  Call the **`EnPlot`** function to visualize the data.

#### `ShEn` Function Analysis (Entropy Calculation)

The `ShEn` function calculates the entropy for every signal in each category. Since the Shannon Equation requires a **Probability Distribution** for discrete values, and EEG signals are continuous time series, the following process is necessary:

1.  **Quantization/Normalization:** The built-in MatLab function **`histcounts`** is used for quantization. It is called with the arguments `('Normalization','Probability')` to directly return a probability vector $p=\{p_{1},p_{2},...,p_{n}\}$.
2.  **Logarithm Check:** A critical step is executed to prevent the entropy from taking an infinite value (due to $\log_2(0)$):
    ```matlab
    temp = temp(temp > 0);
    ```
3.  **Final Calculation:** The Shannon Equation is applied to the resulting probability vector.

#### `EnPlot` Function Analysis (Plotting)

The `EnPlot` function plots the **Normal Distribution** curve for the entropy values of each category.

1.  **Statistical Calculation:** For each category, the **Mean ($\mu$)** and **Standard Deviation ($\sigma$)** of the entropy values are calculated.
2.  **Axis Generation:**
    * The x-axis values (entropy range) are generated using `linspace` based on $\mu \pm 4\sigma$.
    * The y-axis values (Probability Density Function) are calculated using the MatLab function `normpdf`.

---

### 📊 Results and Conclusions

#### Mean and Standard Deviation Table

| Category | Description | Mean Entropy ($\mu$) | Standard Deviation ($\sigma$) |
| :--- | :--- | :--- | :--- |
| **Z** | Open eyes | 4.2805 | 0.2365 |
| **O** | Closed eyes | 4.3683 | 0.2630 |
| **N** | Healthy hemisphere | 4.2980 | 0.3166 |
| **F** | Epileptic focus | 4.2657 | 0.3278 |
| **S** | Seizure crisis | 4.1745 | 0.2583 |

#### Key Findings

* **Category O ($\mu=4.3683$):** Has the **highest mean** (highest entropy).
* **Category S ($\mu=4.1745$):** Has the **lowest mean** (lowest entropy).
* **Category Z ($\sigma=0.2365$):** Has the **smallest standard deviation** (highest, narrowest curve).
* **Category F ($\sigma=0.3278$):** Has the **largest standard deviation** (flattest curve, highest uncertainty/instability).

The results for healthy states (Z, O) showing high entropy and the seizure state (S) showing low entropy **validate the neuroscience theory**.

#### Conclusion on Entropy as a Biomarker

Entropy can be used as a **biomarker** because it:
1.  **Correlates with Neuroscience:** High entropy $\rightarrow$ healthy brain; Low entropy $\rightarrow$ pathological condition.
2.  **Reflects Dynamics:** Intermediate states (N, F) are positioned between the two extremes, indicating an ability to sense the transition from stability to instability, characterized by their large standard deviations ($\sigma$).
3.  **Shows Discriminative Ability:** The difference in entropy between the two extreme states (healthy vs. seizure) is **statistically significant**.

#### ⚠️ Limitation of Shannon Equation

The quantitative separation between the categories is small, attributed to the **Shannon Equation's** method of calculation:
* **Logarithmic Compression:** The use of the logarithm compresses the information, making minimal changes in entropy difficult to detect.
* **Averaging Effect:** It calculates the **average uncertainty** of the entire signal, overlooking crucial instantaneous events that occur in intermediate states (N, F).

---

### 📚 Bibliography

1.  [https://www.researchgate.net/publication/396448622\_Collapse\_of\_Complexity\_in\_Epileptic\_Seizures\_A\_Neurotopological\_Analysis\_with\_Entropy\_Homology\_and\_Ricci\_Curvature](https://www.researchgate.net/publication/396448622_Collapse_of_Complexity_in_Epileptic_Seizures_A_Neurotopological\_Analysis\_with\_Entropy\_Homology\_and\_Ricci\_Curvature)
2.  [https://journals.ametsoc.org/view/journals/atot/35/5/jtech-d-17-0056.1.xml](https://journals.ametsoc.org/view/journals/atot/35/5/jtech-d-17-0056.1.xml)
3.  [https://www.researchgate.net/publication/2173230\_The\_Indefinite\_Logarithm\_Logarithmic\_Units\_and\_the\_Nature\_of\_Entropy](https://www.researchgate.net/publication/2173230_The\_Indefinite\_Logarithm\_Logarithmic\_Units\_and\_the\_Nature\_of\_Entropy)
