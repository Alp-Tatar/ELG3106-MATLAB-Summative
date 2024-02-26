% Part 1: Reflection Spectrum Calculation
n_air = 1;
n_layer1 = 1.4;
n_layer2 = 2.62;
n_layer3 = 3.5;
RI_1 = (n_air - n_layer1) / (n_air + n_layer1);
RI_2 = (n_layer1 - n_layer2) / (n_layer1 + n_layer2);
RI_3 = (n_layer2 - n_layer3) / (n_layer2 + n_layer3);
wavelength_c = 650;
WL_1 = wavelength_c / n_layer1;
WL_2 = wavelength_c / n_layer2;
d_1 = WL_1 / 4;
d_2 = WL_2 / 4;

n_points = 1001;
data = zeros(1, n_points);
L_val = 400;

for k_val = 1:n_points
    delta_1 = (2 * pi * n_layer1 * d_1) / L_val;
    delta_2 = (2 * pi * n_layer2 * d_2) / L_val;
    
    M_11 = exp(1i * delta_1);
    M_12 = exp(1i * delta_2);
    M_21 = exp(-1i * delta_1);
    M_22 = exp(-1i * delta_2);
    
    tao_21 = M_12 * (M_11 * RI_1 + M_21 * RI_2) + M_22 * RI_3 * (M_11 * RI_1 * RI_2 + M_21);
    tao_11 = M_12 * (M_11 + M_21 * RI_1 * RI_2) + M_22 * RI_3 * (M_11 * RI_2 + M_21 * RI_1);
    gamma_sys = tao_21 / tao_11;
    R = 100 * (abs(gamma_sys) * abs(gamma_sys));
    data(1, k_val) = R;
    L_val = L_val + 1;
end

% Plotting the Reflection Spectrum
figure;
hold on;
title('Reflection Spectrum as a Function of Wavelength');
xlabel('Wavelength (nm)');
ylabel('Reflection (%)');
grid on;
plot(400:1:1400, data,'k');

% Part 2: Power Calculation and Plotting
wavelength_array = 400:10:1400;
reflectivity_array = zeros(1, length(wavelength_array));
transmissivity_array = zeros(1, length(wavelength_array));
power_array = zeros(1, length(wavelength_array));
total_power = 0;

% Material properties
ref_coeff_air_layer1 = ref_coefficient(n_air, n_layer1);
ref_coeff_layer1_layer2 = ref_coefficient(n_layer1, n_layer2);
ref_coeff_layer2_cell = ref_coefficient(n_layer2, n_cell);

trans_coeff_air_layer1 = trans_coefficient(n_air, n_layer1);
trans_coeff_layer1_layer2 = trans_coefficient(n_layer1, n_layer2);
trans_coeff_layer2_cell = trans_coefficient(n_layer2, n_cell);

q_air_layer1 = (1/trans_coeff_air_layer1) * [1, ref_coeff_air_layer1; ref_coeff_air_layer1, 1];
q_layer1_layer2 = (1/trans_coeff_layer1_layer2) * [1, ref_coeff_layer1_layer2; ref_coeff_layer1_layer2, 1];
q_layer2_cell = (1/trans_coeff_layer2_cell) * [1, ref_coeff_layer2_cell; ref_coeff_layer2_cell, 1];

wavelength_air = 650;
d_layer1 = 0.25 * wavelength_air / n_layer1;
d_layer2 = 0.25 * wavelength_air / n_layer2;

i_val = 1;
for wavelength_val = wavelength_array
    p_layer1 = prop_matrix(n_layer1, d_layer1, wavelength_val);
    p_layer2 = prop_matrix(n_layer2, d_layer2, wavelength_val);
    t_matrix = q_air_layer1 * p_layer1 * q_layer1_layer2 * p_layer2 * q_layer2_cell;
    reflection_coefficient = abs(t_matrix(2, 1) / t_matrix(1, 1))^2;
    reflectivity_array(i_val) = reflection_coefficient * 100;
    transmission_coefficient = 1 - reflection_coefficient;
    transmissivity_array(i_val) = transmission_coefficient;
    intensity_val = calc_intensity(wavelength_val);
    power_val = transmission_coefficient * intensity_val;
    power_array(i_val) = power_val;
    total_power = total_power + power_val * 10;
    i_val = i_val + 1;
end

% Plotting Power vs. Wavelength
figure;
plot(wavelength_array, power_array, 'k');
ylabel('Instantaneous Power');
xlabel('Wavelength (nm)');
title('Power Spectrum from 400nm to 1400nm');
disp(['The total found power from 400nm to 1400nm is ' num2str(total_power)]);
grid on;
xlim([400, 1400]);
ylim([0.2, 1.4]);

wavelength_array_new = 200:10:2200;

reflectivity_array_new = zeros(1, length(wavelength_array_new));
transmissivity_array_new = zeros(1, length(wavelength_array_new));
power_array_new = zeros(1, length(wavelength_array_new));
total_power_new = 0;

ref_coeff_air_layer1 = ref_coefficient(n_air, n_layer1);
ref_coeff_layer1_layer2 = ref_coefficient(n_layer1, n_layer2);
ref_coeff_layer2_cell = ref_coefficient(n_layer2, n_cell);

trans_coeff_air_layer1 = trans_coefficient(n_air, n_layer1);
trans_coeff_layer1_layer2 = trans_coefficient(n_layer1, n_layer2);
trans_coeff_layer2_cell = trans_coefficient(n_layer2, n_cell);

q_air_layer1 = (1/trans_coeff_air_layer1) * [1, ref_coeff_air_layer1; ref_coeff_air_layer1, 1];
q_layer1_layer2 = (1/trans_coeff_layer1_layer2) * [1, ref_coeff_layer1_layer2; ref_coeff_layer1_layer2, 1];
q_layer2_cell = (1/trans_coeff_layer2_cell) * [1, ref_coeff_layer2_cell; ref_coeff_layer2_cell, 1];

wavelength_air = 650;
d_layer1 = 0.25 * wavelength_air / n_layer1;
d_layer2 = 0.25 * wavelength_air / n_layer2;

i_val = 1;
for wavelength_val = wavelength_array_new
    p_layer1 = prop_matrix(n_layer1, d_layer1, wavelength_val);
    p_layer2 = prop_matrix(n_layer2, d_layer2, wavelength_val);
    t_matrix = q_air_layer1 * p_layer1 * q_layer1_layer2 * p_layer2 * q_layer2_cell;
    reflection_coefficient = abs(t_matrix(2, 1) / t_matrix(1, 1))^2;
    reflectivity_array_new(i_val) = reflection_coefficient * 100;
    transmission_coefficient = 1 - reflection_coefficient;
    transmissivity_array_new(i_val) = transmission_coefficient;
    intensity_val = calc_intensity(wavelength_val);
    power_val = transmission_coefficient * intensity_val;
    power_array_new(i_val) = power_val;
    total_power_new = total_power_new + power_val * 10;
    i_val = i_val + 1;
end

% Plot
figure;
plot(wavelength_array_new, power_array_new, 'k');
ylabel('Instantaneous Power');
xlabel('Wavelength (nm)');
title('Power Distribution from 200nm to 2200nm');
disp(['The total found power from 200nm to 2200nm is ' num2str(total_power_new)]);
grid on;
xlim([200, 2200]);
ylim([0.04, 1.4]);

% Part 3: Refractive Index Variation and Power Calculation
lambda_c = 650;
wavelength_vals = 200:1:2199;
g01 = (n_air - n_layer1) / (n_air + n_layer1);  % reflection coefficient between air and first layer
t01 = 2 * n_air / (n_air + n_layer1); % transmission coefficient between air and first layer
q01 = (1 / t01) * [1, g01; g01, 1];

Power_vals = zeros(1, length(wavelength_vals));
PowerC_vals = zeros(1, 200);
indice_array_vals = 1.3:0.01:3.29;
n2_val = 1.3;

for n_val = 1:length(indice_array_vals)
    t12 = 2 * n_layer1 / (n_layer1 + n2_val);
    t23 = 2 * n2_val / (n2_val + n_layer3);
    g12 = (n_layer1 - n2_val) / (n_layer1 + n2_val);
    g23 = (n2_val - n_layer3) / (n2_val + n_layer3);
    q12 = (1 / t12) * [1, g12; g12, 1];
    q23 = (1 / t23) * [1, g23; g23, 1];
    L_val = 200;
    for k_val = 1:length(wavelength_vals)
        delta_m = (pi / 2) * (lambda_c / L_val);
        p1 = [exp(1i * delta_m), 0; 0, exp(-1i * delta_m)];
        T = q01 * p1 * q12 * p1 * q23;
        g_val = T(2, 1) / T(1, 1);
        r_val = abs(g_val)^2;
        pow_val = (((1 - r_val) * (6.16 * 10^15)) / ((L_val^5) * (exp(2484 / L_val) - 1)));
        Power_vals(1, k_val) = pow_val;
        L_val = L_val + 1;
    end
    PowerC_vals(1, n_val) = sum(Power_vals);
    n2_val = n2_val + 0.01;
end

% Plotting Power vs. Refractive Indices
figure;
subplot(2, 1, 1);
plot(indice_array_vals, PowerC_vals, 'k');
title('Power vs. Refractive Indices');
ylabel('Power (W/m^2)');
xlabel('Refractive Index');
xlim([1.3, 3.3]);

% Comparing Reflectivity for Different Refractive Indices
Reflectivity_n23 = zeros(1, length(wavelength_vals));
Reflectivity_n262 = zeros(1, length(wavelength_vals));
n2_val = 2.3;

for k_val = 1:length(wavelength_vals)
    delta_m = (pi / 2) * (lambda_c / L_val);
    p1 = [exp(1i * delta_m), 0; 0, exp(-1i * delta_m)];
    T = q01 * p1 * q12 * p1 * q23;
    g_val = T(2, 1) / T(1, 1);
    r_val = abs(g_val)^2;
    Reflectivity_n23(1, k_val) = r_val;
    L_val = L_val + 1;
end

n2_val = 2.62;
for k_val = 1:length(wavelength_vals)
    delta_m = (pi / 2) * (lambda_c / L_val);
    p1 = [exp(1i * delta_m), 0; 0, exp(-1i * delta_m)];
    T = q01 * p1 * q12 * p1 * q23;
    g_val = T(2, 1) / T(1, 1);
    r_val = abs(g_val)^2;
    Reflectivity_n262(1, k_val) = r_val;
    L_val = L_val + 1;
end

% Plotting Reflectivity vs. Wavelength for Different Refractive Indices
subplot(2, 1, 2);
plot(wavelength_vals, Reflectivity_n23 * 100, 'k');
hold on;
plot(wavelength_vals, Reflectivity_n262 * 100, 'r');
title('Reflectivity vs. Wavelength for n=2.3 & 2.62');
xlabel('Wavelength (nm)');
ylabel('Reflectivity (%)');
xlim([200, 2200]);

function calc_ref_coefficient = ref_coefficient(n1, n2)
    calc_ref_coefficient = (n1 - n2) / (n2 + n1);
end

function calc_trans_coefficient = trans_coefficient(n1, n2)
    calc_trans_coefficient = (2 * n1) / (n1 + n2);
end

function calc_propagation_matrix = prop_matrix(n, d, wavelength)
    calc_propagation_matrix = [exp(1i * 2 * pi * n * d / wavelength), 0; 0, exp(-(1i * 2 * pi * n * d / wavelength))];
end

function calc_intensity = calc_intensity(wavelength)
    calc_intensity = (6.16 * 10^15) / ((wavelength^5) * (exp(2484 / wavelength) - 1));
end
