clear, clc;

signal_size = 10000000;    % size of signal
upsample_rate = 8;         % upsampling rate
M = 64;                    % M-ary QAM
k = log2(M);               % number of bits per subcarrier
x = 0:0.5:10;              % signal_noise_ratio
y = [];                    % error rates
epochs = 100;              % number of epochs per loop

for signal_noise_ratio = 0:0.5:10
    
    average_error_rate = 0;  % average error rate
    
    for epoch = 1:epochs
        
        % generate random binary signal
        signal = randi([0, 1], signal_size, 1);

        % increase the length of signal to a multiple of k
        if mod(signal_size, k) ~= 0
            completion = zeros(k - mod(signal_size, k), 1);
            signal = [signal; completion];
        end

        % constellation of 64-QAM
        qam_list = [-7 + 7i, -7 + 5i, -7 + 1i, -7 + 3i, -7 - 7i, -7 - 5i,...
            -7 - 1i, -7 - 3i, -5 + 7i, -5 + 5i, -5 + 1i, -5 + 3i, -5 - 7i,...
            -5 - 5i, -5 - 1i, -5 - 3i, -1 + 7i, -1 + 5i, -1 + 1i, -1 + 3i,...
            -1 - 7i, -1 - 5i, -1 - 1i, -1 - 3i, -3 + 7i, -3 + 5i, -3 + 1i,...
            -3 + 3i, -3 - 7i, -3 - 5i, -3 - 1i, -3 - 3i, 7 + 7i, 7 + 5i, 7 + 1i,...
            7 + 3i, 7 - 7i, 7 - 5i, 7 - 1i, 7 - 3i, 5 + 7i, 5 + 5i, 5 + 1i,...
            5 + 3i, 5 - 7i, 5 - 5i, 5 - 1i, 5 - 3i, 1 + 7i, 1 + 5i, 1 + 1i,...
            1 + 3i, 1 - 7i, 1 - 5i, 1 - 1i, 1 - 3i, 3 + 7i, 3 + 5i, 3 + 1i,...
            3 + 3i, 3 - 7i, 3 - 5i, 3 - 1i, 3 - 3i];

        % a series of powers of 2
        power_series = ones(k, 1);
        for i = 2:k
            power_series(i) = 2 * power_series(i - 1);
        end

        % modulate signal with 64-QAM and upsample
        signal_matrix = reshape(signal, length(signal) / k, k);
        signal_symbols = signal_matrix * power_series;
        modulated_signal = zeros(length(signal) / k, 1);
        upsampled_signal = zeros(length(modulated_signal) * upsample_rate, 1);
        modulated_signal = qam_list(signal_symbols + 1).';
        upsampled_signal(1:upsample_rate:length(upsampled_signal)) =...
            modulated_signal;

        % compute coefficients of raised cosine filter
        beta = 0.5;           % rolloff factor
        span = 11;            % number of symbols
        sps = upsample_rate;  % number of samples per symbol
        delay = span * sps / 2;
        t = (-delay:delay) / sps;
        filter_coefficients = zeros(1, length(t));
        mid_point = find(t == 0);
        filter_coefficients(mid_point) = -(pi * (beta - 1) - 4 * beta) / (pi *...
            sps);
        zero_denom_indices = find(abs(abs(4 * beta * t) - 1) < sqrt(eps));
        filter_coefficients(zero_denom_indices) = 1 / (2 * pi * sps) * (pi *...
            (beta + 1) * sin(pi * (beta + 1) / (4 * beta)) - 4 * beta * sin(pi *...
            (beta - 1) / (4 * beta)) + pi * (beta - 1) * cos(pi * (beta - 1) /...
            (4 * beta)));
        indices = 1:length(t);
        indices([mid_point, zero_denom_indices]) = [];
        t = t(indices);
        filter_coefficients(indices) = -4 * beta / sps * (cos((1 + beta) * pi...
            * t) + sin((1 - beta) * pi * t) ./ (4 * beta * t)) ./ (pi * ((4 *...
            beta * t) .^ 2 - 1));
        filter_coefficients = filter_coefficients / sqrt(sum(...
            filter_coefficients .^ 2));

        % filter signal
        filtered_signal = conv(upsampled_signal, filter_coefficients);

        % add noise
        noise_power = 1 / (10 ^ (signal_noise_ratio / 10));
        noisy_signal = filtered_signal + sqrt(noise_power / 2) * (randn(size(...
            filtered_signal)) + 1i * randn(size(filtered_signal)));

        % filter received signal
        filtered_noisy_signal = conv(noisy_signal, filter_coefficients);

        % downsample signal
        downsampled_signal = filtered_noisy_signal(1:upsample_rate:length(...
            filtered_noisy_signal));
        downsampled_signal = downsampled_signal((length(downsampled_signal)-...
            length(modulated_signal))/2+1:(length(downsampled_signal)+...
            length(modulated_signal))/2);

        % demodulate signal
        distance = abs(downsampled_signal - qam_list);
        [min_distances, received_symbols] = min(distance, [], 2);
        received_symbols = received_symbols - 1;
        received_matrix = zeros(length(received_symbols), k);
        for i = k:-1:1
            received_matrix(:, i) = floor(received_symbols / power_series(i));
            received_symbols = mod(received_symbols, power_series(i));
        end
        received_signal = received_matrix(:);

        % calculate error rate
        total_errors = sum((signal ~= received_signal) * 1);
        error_rate = total_errors / length(signal);
        average_error_rate = average_error_rate + error_rate;
    end
    
    average_error_rate = average_error_rate / epochs;
    y = [y, average_error_rate];
end

plot(x, y);
xlabel('Signal noise ratio');
ylabel('Error rate');

set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, ...
    'PaperPosition', [0 0 screenposition(3:4)], ...
    'PaperSize', screenposition(3:4));
print -dpdf -painters error_rate_wrt_snr