% Plot excitation, emission, and absorbance spectra.

clear;

% Load data.
filename = '../../2012-05-15/2012-05-15 plate B3 spectrum after titration.txt';

nrows = 4;
ncols = 6; 

%titles = {'200 nM Src', 'buffer', '200 nM Src + 5 uM bosutinib', '5 uM bosutinib', '200 nM Src + 5 uM vandetanib', '5 uM vandetanib', '200 nM Src + 5 uM erlotinib', '5 uM erlotinib'};

column_titles = {'2 uM bosutinib + 375 nM Src', '2 uM bosutinib', '2 uM bosutinib + 150 nM Src', '2 uM bosutinib', '2 uM bosutinib + 60 nM Src', '2 uM bosutinib'};
row_titles = {'replicate 1', 'replicate 2', 'replicate 3', 'replicate 4'};

% Read spectra.

infile = fopen(filename, 'r');
[absorbance_wavelengths, absorbance_spectrum] = get_spectrum(infile, '^absorbance:Spectrum$', nrows, ncols);
[excitation_wavelengths_top, excitation_spectrum_top] = get_spectrum(infile, '^excitation spectra \d+ emission:EX Spectrum$', nrows, ncols);
[excitation_wavelengths_bottom, excitation_spectrum_bottom] = get_spectrum(infile, '^excitation spectra emission \d+ bottom:EX Spectrum$', nrows, ncols);
[emission_wavelengths_top, emission_spectrum_top] = get_spectrum(infile, '^emission spectrum excite \d+:EM Spectrum$', nrows, ncols);
[emission_wavelengths_bottom, emission_spectrum_bottom] = get_spectrum(infile, '^emission spectrum excite \d+ bottom:EM Spectrum$', nrows, ncols);
fclose(infile);


% Plot absorbance and excitation spectra.
clf;

ny = nrows;
nx = ncols;

for row = 1:nrows
  for col = 1:ncols
    subplot(ny, nx, ncols*(row-1) + col);
    hold on;
    abs = squeeze(absorbance_spectrum(row, col, :));

    ex_top = squeeze(excitation_spectrum_top(row, col, :));
    ex_bottom = squeeze(excitation_spectrum_bottom(row, col, :));

    em_top = squeeze(emission_spectrum_top(row, col, :));
    em_bottom = squeeze(emission_spectrum_bottom(row, col, :));

    ex_top = superimpose_spectra(absorbance_wavelengths, abs, excitation_wavelengths_top, ex_top);
    ex_bottom = superimpose_spectra(absorbance_wavelengths, abs, excitation_wavelengths_bottom, ex_bottom);

    em_top = em_top * sum(ex_top) / sum(em_top);
    em_bottom = em_bottom * sum(ex_bottom) / sum(em_bottom);

    plot(absorbance_wavelengths, abs, 'k-');

    plot(excitation_wavelengths_top, ex_top, 'b-');
    plot(excitation_wavelengths_bottom, ex_bottom, 'b--');    

    plot(emission_wavelengths_top, em_top, 'r-');
    plot(emission_wavelengths_bottom, em_bottom, 'r--');    

    xlabel('wavelength / nm');    

    if (row == 1) && (col == 1)
      legend('abs', 'ex top', 'ex bottom', 'em top', 'em bottom', 'location', 'northeast');
    end

    title(sprintf('%s; %s', char(column_titles{col}), char(row_titles{row})));
    axis([200 700 0 0.5]);

  end
end

filename = 'excitation-emission-absorbance.eps';
exportfig(gcf, filename, 'width', 22, 'height', 17, 'color', 'cmyk');
system(sprintf('epstopdf "%s"', filename));
