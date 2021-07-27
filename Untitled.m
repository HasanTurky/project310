
clear all
close all
clc

%% channel parameters
L =10;                 %length in m
a = 6 ;                %transmitter height
b = 11 ;                %receiver height
h = 16;                %surface height
o = 100;                 %order
r =1.33;                 %surface reflection coefficient
%% Simulation Parameters
%Moduluation method: QPSK, 8PSK, 16QAM
modulation_method = 1 ; %1, 2, or 3

%fft bin Size
fftSize = 64;
nFFT = 64;

%Size of cycle prefix extension

cyclic_prefix_extension = 16; %

%number of channel taps.

channel_tap = 2; % for time delay in channel

% mod_order = find(ismember(mod_methods,mod_method));

%% Read the image and convert it into binary format.
im1 = imread('eagle2.jpg');
im2 = imread('im2.jpg')
im_bin1 = dec2bin(im1(:))';
im_bin2 = dec2bin(im2(:))'
im_bin = [im_bin1 im_bin2];
im_bin = im_bin(:);
im_bin_stored = im_bin;

%% interleaving and encoding
im_bin = im_bin-'0';
im_bin = encode(im_bin,7,4,'hamming/binary');
im_bin = randintrlv(im_bin,3107);
im_bin = char(im_bin+'0');

%% Binary stream to symbols
% 1. parse binary stream into mod_order bit symbols
% 2. pads input signal to appropriate length

if modulation_method ==1
    mod_order = 2;
elseif modulation_method == 2
    mod_order = 3;
else
    mod_order = 4;
end
symbol_remainder = mod(mod_order-mod(length(im_bin),mod_order),mod_order);
zero_padding = repmat('0',symbol_remainder,1);
im_bin_zero_padded = [im_bin;zero_padding];
constellation_data = reshape(im_bin_zero_padded,mod_order,length(im_bin_zero_padded)/mod_order)';
constellation_symbol_id = bin2dec(constellation_data);

%% modulation technique

% Phase shift keying about unit circle 
if modulation_method == 1 
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quadrature = sin(n+pi/4);
    modulation_ID = (in_phase + quadrature*1i);
end

if modulation_method == 2
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quadrature = sin(n+pi/4);
    modulation_ID = (in_phase + quadrature*1i);
end

if modulation_method == 3
    mod_ind = sqrt(2^mod_order);
    %n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
    quadrature = repmat(linspace(-1,1,mod_ind)',1,mod_ind);
    modulation_ID = (in_phase(:) + quadrature(:)*1i);
end


X = modulation_ID(constellation_symbol_id+1);



%% move to time domain
dft_remainder = mod(nFFT-mod(length(X),nFFT),nFFT);

if modulation_method == 3
    X_zero_padded = [X;zeros(dft_remainder,1)];
else
    X_zero_padded = [X';zeros(dft_remainder,1)];
end
X_parallel = reshape(X_zero_padded,fftSize,length(X_zero_padded)/fftSize);
x = ifft(X_parallel);

x_cyclic_prfx_ext = [x(end-cyclic_prefix_extension+1:end,:);x];
x_s   = x_cyclic_prfx_ext(:);

%% Add Channel Noises
% Calculate data power
data_pwr = mean(abs(x_s.^2));

% Add noise to the channel
u = 2; %velocity of wind in m/s
Nw = 20.5 + 22.4*log10(u); % wind noise in dB re 1 pW
Nwp = 10^(Nw/10)*1e-12;
f = 0.528;
Ns=Nw+20.7-15.9*log(f);
Nsp = 10^(Ns/10)*1e-12;
accoustic_freq = 1e5 %MHz
% Nt=-15+20*log10(accoustic_freq);
% Ntp = 10^(Nt/10)*1e-12;
a = 5 + 5.7 * (5-u);
b = 50 + 2.4 * (5-u);
rr = 0.4;
Nr = b + a * log10(rr); % raining noise
Nrp = 10^(Nr/10)*1e-12;


noise_pwr = Nwp + Nsp + Nrp %+Ntp;
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;
x_s_noise = x_s + noise;

% Measure SNR
snr = 10*log10(mean(abs(x_s.^2))/mean(abs(noise.^2)));

%% Apply multipath fading channel 
g = exp(-(0:channel_tap-1));
g = g/norm(g);
a = mpm(L,a,b,h,o,r);
x_s_noise_fading = conv(x_s_noise,a,'same');
x_s_noise_fading = conv(x_s_noise_fading,g,'same');

%% Use FFT to move to frequency domain
% Remove cyclic prefix extension and shift from serial to parallel
x_p = reshape(x_s_noise_fading,fftSize+cyclic_prefix_extension,length(x_s_noise_fading)/(fftSize+cyclic_prefix_extension));
x_p_cpr = x_p(cyclic_prefix_extension+1:end,:);

% Move to frequency domain
X_hat_blocks = fft(x_p_cpr);


%% channels estimation 

    if channel_tap > 1
            G = X_hat_blocks(:,1)./X_parallel(:,1);
            X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));
    end


%% Symbol demodulation
% remove fft padding 
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-dft_remainder);

%Recover data from modulated symbols
A=[real(modulation_ID) imag(modulation_ID)];
if (size(A,2)>2)
    A=[real(modulation_ID)' imag(modulation_ID)'];
end
reconstructed_symbols = knnsearch(A,[real(X_hat) imag(X_hat)])-1;

%Parse to binary stream to remove symbol padding
reconstructed_symbol_constellation = dec2bin(reconstructed_symbols);
reconstructed_binary_image = reshape(reconstructed_symbol_constellation',numel(reconstructed_symbol_constellation),1);
reconstructed_binary_image = reconstructed_binary_image(1:end-symbol_remainder);
%% deinterleaving
reconstructed_binary_image = reconstructed_binary_image-'0';
reconstructed_binary_image = randdeintrlv(reconstructed_binary_image,3107);
reconstructed_binary_image = decode(reconstructed_binary_image,7,4,'hamming/binary');
reconstructed_binary_image = char(reconstructed_binary_image+'0');
%% bit error rate
ber = sum(abs(reconstructed_binary_image-im_bin_stored))/length(im_bin_stored);

%% recover image
reconstructed_image = reshape(reconstructed_binary_image,8,numel(reconstructed_binary_image)/8);
reconstructed_image = uint8(bin2dec(reconstructed_image'));
reconstructed_image = reshape(reconstructed_image,size(im1));

%% result observation

% Original image
subplot(1,2,1);
imshow(im1);
title('Transmitted Image');

% Recovered image
subplot(1,2,2);
imshow(reconstructed_image);
title(sprintf('\\Received Image\n \\rmBit Error Rate: %.2g',ber));
