function beta = extract_feature(impath)
% Extracting exponent feature of ROI
% 
% Author: Haifeng Xu
% March07/2018

% Read reference image and ROI-detection
im = imread(impath);        %read image
im = im2double(im);         %convert to double
im = imresize(im, .25);     %re-size for speed up
isMLO = true;               %true for MLO view
isFFDM = true;              %true only for FFDM images

% Segment breast boundary and chest wall:
[mask, contour, cwall] = segBreast(im, isMLO, isFFDM);

% Detect maximum squared ROI, and resize it to square
mask_SQ = sqmax(mask);
imc = im (mask_SQ);
L = sqrt(length(imc(:)));
imc = reshape(imc, [L, L]);
% figure(1);
% imshow(imc);

% Applying Hanning window function:
w = hanning (L);
w = w.* w';
imc = imc.*w;

% Applying 2-D Fourier Transform:
mask_SQ_DFT = fft2(imc);
spectral_density = abs(mask_SQ_DFT).^2;
%imagesc(spectral_density)
close
surf(spectral_density)
spectral_density = fftshift(spectral_density);
%imagesc(log(spectral_density))
%figure, imshow(spectral_density);

% Creating frequency domain matrix as [u, v]
fvec = linspace (-0.5*L, 0.5*L, L);
[u, v]= meshgrid(fvec, fvec);
f = sqrt(u.^2 + v.^2);
%figure, imagesc(f)

% Creating empty vector for meanfrequency and mean power spectural
mean_freq_vector = ones(1, 8);
mean_power_spe = ones(1, 8);

% Average the power spectrum along a radial slice
R = 0.5*int64(L);
for i = 1:1:8
    fmax = f(R, R-(i+1)*0.1*R);
    fmin = f(R, R-i*0.1*R);
    select = (f<fmax) & (f>fmin);
    %figure, imshow(select)
    mean_frequency = mean(mean(select));
    mean_freq_vector(i) = mean_frequency;
    
    %selected_power_sp = spectral_density * select;
    %figure, imshow(selected_power_sp)
    %put a logic indexing then it will be selected directly
    meanvalue = mean(spectral_density(select));
    mean_power_spe(i) = meanvalue;
end

% Estimate the exponent beta from the least squares fit of 
% the power law spectrum up to the Nyquist frequency
figure, plot(mean_freq_vector, mean_power_spe,'+');
beta = polyfit(log(mean_freq_vector), log(mean_power_spe), 1);

end

