% Constants and initial values
refractive_index_air = 1;
refractive_index_1 = 1.4;
refractive_index_2 = 1.4;
refractive_index_3 = 3.15;
refractive_index_4 = 3.5;
wavelength_c = 650;
transmission_1 = 2 * refractive_index_air / (refractive_index_air + refractive_index_1);
transmission_4 = 2 * refractive_index_3 / (refractive_index_3 + refractive_index_4);
Q_1 = calculateQ(refractive_index_air, refractive_index_1, transmission_1);
Q_4 = calculateQ(refractive_index_3, refractive_index_4, transmission_4);
wavelength_1 = wavelength_c / refractive_index_1;
thickness_1 = wavelength_1 / 4;
wavelength_3 = wavelength_c / refractive_index_3;
thickness_3 = wavelength_3 / 4;

% Initialize variables
num_steps_1 = 2001;
Range_power_1 = 200:1:2200;
Range_n_2 = 1.4:0.01:3.0;
data_power_1 = zeros(1, num_steps_1);
data_n_2_2k = zeros(1, length(Range_n_2));

% Loop for 200-2200 nm
data_n_2_2k = calculatePowerTransmission(Range_n_2, refractive_index_1, refractive_index_2, refractive_index_3, num_steps_1, data_power_1, data_n_2_2k, Q_1, Q_4, wavelength_c);

% Initialize variables for the second loop
num_steps_2 = 1001;
Range_power_2 = 400:1:1400;
Range_n_2 = 1.4:0.01:3.0;
data_power_2 = zeros(1, num_steps_2);
data_n_2_1k = zeros(1, length(Range_n_2));

% Loop for 400-1400 nm
data_n_2_1k = calculatePowerTransmission(Range_n_2, refractive_index_1, refractive_index_2, refractive_index_3, num_steps_2, data_power_2, data_n_2_1k, Q_1, Q_4, wavelength_c);

% Plotting
plotPowerTransmission(Range_n_2, data_n_2_1k, data_n_2_2k);

% Reflectivity calculations and plotting
plotReflectivity(wavelength_c, refractive_index_1, refractive_index_2, refractive_index_3, refractive_index_4, Q_1);

function Q = calculateQ(n1, n2, t)
    % Calculate the Q matrix for a given layer
    reflection_coefficient = (n1 - n2) / (n1 + n2);
    Q = (1 / t) * [1 reflection_coefficient; reflection_coefficient 1];
end

function data_n_2 = calculatePowerTransmission(Range_n_2, n_1, n_2, n_3, t, data_power, data_n_2, Q_1, Q_4, wavelength_c)
    % Calculate power transmission for a range of refractive indices
    for x = 1:length(Range_n_2)
        L = 200;
        t_2 = 2 * n_1 / (n_1 + n_2);
        t_3 = 2 * n_2 / (n_3 + n_2);
        reflection_2 = (n_1 - n_2) / (n_1 + n_2);
        reflection_3 = (n_2 - n_3) / (n_2 + n_3);
        Q_2 = calculateQ(n_1, n_2, t_2);
        Q_3 = calculateQ(n_2, n_3, t_3);

        % Rest of the loop
        for k = 1:t
            WL_2 = wavelength_c / n_2;
            d_2 = WL_2 / 4;
            delta_m = (pi / 2) * (wavelength_c / L);

            Prop_M = [exp(1i * delta_m) 0; 0 exp(-1i * delta_m)];
            T_Mat = (Q_1 * Prop_M * Q_2 * Prop_M * Q_3 * Prop_M * Q_4);
            gamma_sys = T_Mat(2, 1) / T_Mat(1, 1);
            R = abs(gamma_sys) * abs(gamma_sys);
            T = 1 - R;

            Power = T * (6.16e15 / (power(L, 5) * (exp(2484 / L) - 1)));
            data_power(1, k) = Power;
            L = L + 1;
        end

        data_n_2(1, x) = sum(data_power);
        n_2 = n_2 + 0.01;
    end
end

function plotPowerTransmission(Range_n_2, data_n_2_1k, data_n_2_2k)
    % Plot power transmission for different wavelength ranges
    figure
    hold on
    title('Power Transmission at Varying Refractive Indices n2 ');
    xlabel('n2');
    ylabel('Power [W/m^2]');
    plot(Range_n_2, data_n_2_1k, 'r');
    plot(Range_n_2, data_n_2_2k, 'k');
    legend('Wavelength (200-2200 nm)', 'Wavelength (400-1400 nm)');
    grid on;
end

function plotReflectivity(wavelength_c, n1, n2, n3, n4, Q_1)
    % Plot reflectivity for different refractive indices
    wavelength_range = 400:1399;
    Reflectivity_n2_19 = zeros(1, length(wavelength_range));
    L = 200;
    n2 = 2.3;
    t12 = 2 * n1 / (n1 + n2);
    t23 = 2 * n2 / (n2 + n3);
    t34 = 2 * n4 / (n3 + n4);
    g12 = (n1 - n2) / (n1 + n2);
    g23 = (n2 - n3) / (n2 + n3);
    g34 = (n3 - n4) / (n3 + n4);
    q12 = (1 / t12) * [1 g12; g12 1];
    q23 = (1 / t23) * [1 g23; g23 1];
    q34 = (1 / t34) * [1 g34; g34 1];

    for k = 1:length(wavelength_range)
        delta_m = (pi / 2) * (wavelength_c / L);
        p1 = [exp(1i * delta_m) 0; 0 exp(-1i * delta_m)];
        T = (Q_1 * p1 * q12 * p1 * q23 * p1 * q34);
        g = T(2, 1) / T(1, 1);
        r = abs(g)^2;
        Reflectivity_n2_19(1, k) = r;
        L = L + 1;
    end

    % Plot Reflectivity
    figure
    hold on
    title('Reflectivity (%) for Varying Wavelength Contrasting the Reflectivity of n=2.3 and n=2.62');
    xlabel('Wavelength (nm)');
    ylabel('Reflectivity (%)');
    plot(wavelength_range, Reflectivity_n2_19 * 100, 'r');
    
    % Reflectivity calculations and plotting for n=2.62
    Reflectivity_n2_62 = zeros(1, length(wavelength_range));
    L = 200;
    n2 = 2.62;
    t12 = 2 * n1 / (n1 + n2);
    t23 = 2 * n2 / (n2 + n3);
    g12 = (n1 - n2) / (n1 + n2);
    g23 = (n2 - n3) / (n2 + n3);
    q12 = (1 / t12) * [1 g12; g12 1];
    q23 = (1 / t23) * [1 g23; g23 1];

    for k = 1:length(wavelength_range)
        delta_m = (pi / 2) * (wavelength_c / L);
        p1 = [exp(1i * delta_m) 0; 0 exp(-1i * delta_m)];
        T = (Q_1 * p1 * q12 * p1 * q23);
        g = T(2, 1) / T(1, 1);
        r = abs(g)^2;
        Reflectivity_n2_62(1, k) = r;
        L = L + 1;
    end
    
    % Plot Reflectivity
    plot(wavelength_range, Reflectivity_n2_62 * 100, 'k');
    legend('Reflectivity for n=2.3', 'Reflectivity for n=2.62');
    grid on;
    xlim([400, 1400]);
end
