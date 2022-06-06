clear, clc;
signal_size = 10000000;
upsample_rate = 8;
signal_noise_ratio = 10;
M = 64;
k = log2(M);
signal = randi([0, 1], signal_size, 1);
if mod(signal_size, k) ~= 0
    completion = zeros(k - mod(signal_size, k), 1);
    signal = [signal; completion];
end
signal_matrix = reshape(signal, length(signal) / k, k);
signal_symbols = bi2de(signal_matrix);
modulated_signal = qammod(signal_symbols, M);
upsampled_signal = upsample(modulated_signal, upsample_rate);
filter_coefficients = rcosdesign(0.5, 11, upsample_rate, 'sqrt').';
filtered_signal = conv(upsampled_signal, filter_coefficients);
noisy_signal = awgn(filtered_signal, signal_noise_ratio);
filtered_noisy_signal = conv(noisy_signal, filter_coefficients);
downsampled_signal = downsample(filtered_noisy_signal, upsample_rate);
downsampled_signal = downsampled_signal((length(downsampled_signal)-...
    length(modulated_signal))/2+1:(length(downsampled_signal)+...
    length(modulated_signal))/2);
received_symbols = qamdemod(downsampled_signal, M);
received_matrix = de2bi(received_symbols, k);
received_signal = received_matrix(:);
error_calculator = comm.ErrorRate;
error_rate = error_calculator(signal, received_signal)
