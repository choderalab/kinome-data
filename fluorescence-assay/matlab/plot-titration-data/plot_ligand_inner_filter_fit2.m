% Check whether primary inner filter effects could be contributing to observed fluorecence behavior for ligands in buffer.

clear;

% Load data.
filename = '../../2012-05-09/2012-05-08 1-2 dilution series Src 400 uM plate A4 after spin down.txt';

infile = fopen(filename, 'r');

%regexp_string = ':280,340';
%regexp_string = ':280,440';
%regexp_string = ':280,460';
regexp_string = ':280,480';
%regexp_string = ':350,440';
%regexp_string = ':350,460';
%regexp_string = ':350,480';
%regexp_string = 'read 3 gain 125:280,480';
%regexp_string = 'gain 125:350,480';

markersize = 15;

ligand_concentrations = 0.5 * 10e-6 * 2.^(-(0:11)); % ligand concentrations in cell (M)
protein_concentration = 100e-9; % protein concentration (M)

nrows = 8; % number of rows
ncolumns = 12; % number of columns
data = zeros(nrows, ncolumns);

while ~feof(infile)
  line = fgetl(infile);

  S = regexp(line, regexp_string);  
  if length(S) > 0
    % Match found.

    % Read subsequent data
    line = fgetl(infile); % skip header line
    for row = 1:nrows
      line = fgetl(infile); 
      row_data = sscanf(line, '%*c\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
      
      if length(row_data) < 8
	break
      end

      % Store.
      data(row, :) = row_data;
    end      

  end
end
fclose(infile);

% normalize data
%data = data / max(max(data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

% Fit using inner filter effect model.
well_volume = 200e-6; % L
well_area = 0.3165 / 100 / 100; % well area (m^2)
path_length = well_volume / 1000 / well_area; % path length (m)
observed_ligand_fluorescence = data(4,:);
observed_buffer_fluorescence = data(2,:);

[computed_ligand_fluorescence_nonlinear, computed_buffer_fluorescence_nonlinear, extinction_coefficient_nonlinear] = fit_ligand_inner_filter(well_volume, path_length, ligand_concentrations, observed_ligand_fluorescence, observed_buffer_fluorescence, 'nonlinear');
[computed_ligand_fluorescence, computed_buffer_fluorescence, extinction_coefficient] = fit_ligand_inner_filter(well_volume, path_length, ligand_concentrations, observed_ligand_fluorescence, observed_buffer_fluorescence, 'linear');

plot(ligand_concentrations, observed_ligand_fluorescence, 'b.', 'markersize', markersize);
hold on;
plot(ligand_concentrations, computed_ligand_fluorescence, 'b-');
plot(ligand_concentrations, computed_ligand_fluorescence_nonlinear, 'b--');

plot(ligand_concentrations, observed_buffer_fluorescence, 'g.', 'markersize', markersize);
plot(ligand_concentrations, computed_buffer_fluorescence, 'g-');
set(gca, 'xscale', 'log');

%legend('observed ligand', 'linear', sprintf('nonlinear (\\epsilon = %.0f cm^{-1} M^{-1})', extinction_coefficient_nonlinear), 'observed buffer', 'fit buffer', 'location', 'northwest');
legend('observed ligand', 'I_{obs} = I_0 \Phi_f \epsilon b c', sprintf('I_{obs} = I_0 \\Phi_f (1 - e^{-\\epsilon b c}) (\\epsilon = %.0f cm^{-1} M^{-1})', extinction_coefficient_nonlinear), 'observed buffer', 'fit buffer', 'location', 'northwest');

title(sprintf('bosutinib%s', regexp_string));
xlabel('[bosutinib] / M');
ylabel('fluorescence (AU)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('primary inner filter effect fit - bosutinib in buffer%s.eps', regexp_string);
exportfig(gcf, filename, 'width', 10, 'height', 7.5, 'color', 'cmyk');
system(sprintf('epstopdf "%s"', filename));