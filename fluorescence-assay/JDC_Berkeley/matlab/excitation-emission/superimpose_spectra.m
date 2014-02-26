function spectrum2_scaled = superimpose_spectra(spectrum1_wavelengths, spectrum1, spectrum2_wavelengths, spectrum2)

% find shared wavelengths
[common_wavelengths, indices1, indices2] = intersect(spectrum1_wavelengths, spectrum2_wavelengths);

% align by best fit
spectrum2_scaled = spectrum2 * mean(spectrum1(indices1)) / mean(spectrum2(indices2));

return
