clear, clc;
signal_size = 10000000;
upsample_rate = 8;
signal_noise_ratio = 10;
signal = randi([0, 1], signal_size, 1);
BPSK_modulator = comm.BPSKModulator;
modulated_signal = real(BPSK_modulator(signal));
upsampled_signal = upsample(modulated_signal, upsample_rate);
filter_coefficients = rcosdesign(0.5, 11, upsample_rate, 'sqrt').';
filtered_signal = conv(upsampled_signal, filter_coefficients);
noisy_signal = awgn(filtered_signal, signal_noise_ratio);
filtered_noisy_signal = conv(noisy_signal, filter_coefficients);
downsampled_signal = downsample(filtered_noisy_signal, upsample_rate);
BPSK_demodulator = comm.BPSKDemodulator;
demodulated_signal = BPSK_demodulator(downsampled_signal);
received_signal = (demodulated_signal >= 0.5) * 1;
received_signal = received_signal((length(received_signal)-length(signal...
    ))/2+1:(length(received_signal)+length(signal))/2);
error_calculator = comm.ErrorRate;
error_rate = error_calculator(signal, received_signal)