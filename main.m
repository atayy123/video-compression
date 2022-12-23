% get first 50 frames of the video
Y = yuv_import_y('foreman_qcif/foreman_qcif.yuv',[176 144],50);
for i=1:20
    figure
    image(Y{i});
end
%% Intra-Frame Video Coder
% calculate distortion and rate for different stepsizes
d_intra = zeros(4,1);
en_intra = zeros(4,1);

for i = 3:6
    [d_intra(i-2),en_intra(i-2)] = intraframe_coder(Y, 2^i);
end

% the result is bits/pel, we need kbits/s
[n,m] = size(Y{1});
sz = n*m;
entro_intra = en_intra*sz*30/1000;

% % plot distortion per rate
% figure
% plot(d, entro_intra)
% grid on
% xlabel("Distortion")
% ylabel("Rate (kbit/s)")

% calculate PSNR
psnr_intra = 10*log10(255^2./d_intra);
% plot rate-PSNR curve
figure
plot(entro_intra, psnr_intra)
grid on
hold on
ylabel("PSNR")
xlabel("Rate (kbit/s)")

%% Conditional Replenishment Video Coder
% calculate distortion and rate for different stepsizes
d_mc = zeros(4,1);
en_mc = zeros(4,1);

for i = 3:6
    [d_mc(i-2),en_mc(i-2)] = conditional_rep(Y, 2^i, en_intra);
end

entro_cond = en_mc*30/1000;
psnr_cond = 10*log10(255^2./d_mc);
plot(entro_cond, psnr_cond)
legend('Intra-Frame Video Coder','Conditional Replenishment Video Coder')

%% Video coder with motion compensation
d_mc = zeros(4,1);
en_mc = zeros(4,1);

for i = 3:6
    stepsize = 2^i;
    [d_mc(i-2),en_mc(i-2)] = coder_mc(Y, stepsize, en_intra(log2(stepsize)-2));
end
% bit rate in kbit and PSNR
entro_mc = en_mc*30/1000;
psnr_mc = 10*log10(255^2./d_mc);
plot(entro_mc, psnr_mc)


