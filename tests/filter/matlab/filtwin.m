function [ph_out] = filtwin(ph, alpha, beta, n_win, n_pad, low_pass_wavelength)

n_fft = n_win + n_pad

freq0=1/low_pass_wavelength;
freq_i=-(n_fft)/n_fft/2:1/n_fft:(n_fft-2)/n_fft/2;
butter_i=1./(1+(freq_i/freq0).^(2*5));
low_pass=butter_i'*butter_i;
low_pass=fftshift(low_pass);

B = gausswin(7) * gausswin(7)';
ph_bit = zeros(n_fft);

ph_bit(1:n_win, 1:n_win) = ph(1:n_win, 1:n_win);
ph_fft = fft2(ph_bit);
H = abs(ph_fft);
H = ifftshift( filter2(B, fftshift(H)));
meanH = median(H(:));
if meanH~=0
    H = H/meanH;
end
H = H.^alpha;
H = H-1;
H(H<0) = 0;

G = H * beta + low_pass;

%M = ph_fft.*G;
%M(1:5,1:5)
%M(14:18,14:18)
%M(28:32,28:32)

ph_filt = ifft2(ph_fft.*G);
ph_out = ph_filt(1:n_win, 1:n_win);

