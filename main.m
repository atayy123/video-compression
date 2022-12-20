% get first 50 frames of the video
Y = yuv_import_y('foreman_qcif/foreman_qcif.yuv',[176 144],50);
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
d_cond = zeros(4,1);
en_cond = zeros(4,1);

for i = 3:6
    [d_cond(i-2),en_cond(i-2)] = conditional_rep(Y, 2^i, en_intra);
end

entro_cond = en_cond*sz*30/1000;
psnr_cond = 10*log10(255^2./d_cond);
plot(entro_cond,psnr_cond)