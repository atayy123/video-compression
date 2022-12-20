% get first 50 frames of the video
Y = yuv_import_y('foreman_qcif/foreman_qcif.yuv',[176 144],50);

% calculate distortion and rate for different stepsizes
d = zeros(4,1);
en = zeros(4,1);

for i = 3:6
    [d(i-2),en(i-2)] = intraframe_coder(Y, 2^i);
end

% the result is bits/pel, we need kbits/s
[n,m] = size(Y{1});
sz = n*m;
entro = en*sz*30/1000;

% plot distortion per rate
figure
plot(d, entro)
grid on
xlabel("Distortion")
ylabel("Rate (kbit/s)")

% calculate PSNR
psnr = 10*log10(255^2./d);
% plot rate-PSNR curve
figure
plot(entro, psnr)
grid on
ylabel("PSNR")
xlabel("Rate (kbit/s)")
