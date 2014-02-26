% plot serial dilution data

clear;

% Load data.
filename = '../../2012-05-09/2012-05-08 1-2 dilution series Src 400 uM plate A4 before spin down.txt';
filename = '../../2012-05-09/2012-05-08 1-2 dilution series Src 400 uM plate A4 after spin down.txt';

infile = fopen(filename, 'r');

%regexp_string = ':280,340';
%regexp_string = ':280,440';
%regexp_string = ':280,460';
%regexp_string = ':280,480';
%regexp_string = ':350,440';
%regexp_string = ':350,460';
regexp_string = ':350,480';

ligand_concentrations = 0.5 * 10e-6 * 2.^(-(0:11)); % ligand concentrations in cell (M)
protein_concentration = 200e-9; % protein concentration (M)

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

      % Store.
      data(row, :) = row_data;
    end      

  end
end
fclose(infile);

% normalize data
data = data / max(max(data));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

markersize = 10;
protein_name = 'Src';
well_volume = 200e-6; % L

observed_protein_fluorescence = data(1,:);
observed_buffer_fluorescence = data(2,:);

% Fix pipetting error.
observed_protein_fluorescence(11) = observed_protein_fluorescence(12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bosutinib

subplot(3,2,1);
ligand_name = 'bosutinib';
observed_mixture_fluorescence = data(3,:);
observed_ligand_fluorescence = data(4,:);
semilogx(ligand_concentrations, observed_protein_fluorescence, 'g.', ligand_concentrations, observed_ligand_fluorescence, 'k.', ligand_concentrations, observed_buffer_fluorescence, 'b.', ligand_concentrations, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
xlabel('ligand concentration (M)');
title(ligand_name);
ylabel('fluorescence (AU)');
legend(protein_name, ligand_name, 'buffer', sprintf('%s:%s', protein_name, ligand_name), 'location', 'west');
axis([1e-9 1e-5 0 1]);
set(gca, 'XTick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]);

[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_buffer_fluorescence, computed_mixture_fluorescence] = fit_dilution_data(well_volume, protein_concentration, ligand_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_buffer_fluorescence, observed_mixture_fluorescence)
hold on;
plot(ligand_concentrations, computed_protein_fluorescence, 'g-', ligand_concentrations, computed_ligand_fluorescence, 'k-', ligand_concentrations, computed_buffer_fluorescence, 'b-', ligand_concentrations, computed_mixture_fluorescence, 'r-');
ligand_name
Kd
title(sprintf('%s%s : Kd = %.3f nM, [P] = %.3f nM', ligand_name, regexp_string, Kd / 1e-9, protein_concentration / 1e-9));

subplot(3,2,2);
semilogx(ligand_concentrations, observed_ligand_fluorescence, 'k.', 'markersize', markersize); 
xlabel('ligand concentration (M)');
ylabel('fluorescence (AU)');
title(sprintf('%s%s', ligand_name, regexp_string))
axis([1e-9 1e-5 0 max(observed_ligand_fluorescence)]);
set(gca, 'XTick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]);;

% Fit using inner filter effect model.
well_area = 0.3165 / 100 / 100; % well area (m^2)
path_length = well_volume / 1000 / well_area; % path length (m)
%[computed_ligand_fluorescence, computed_buffer_fluorescence, extinction_coefficient] = fit_ligand_inner_filter(well_volume, path_length, ligand_concentrations, observed_ligand_fluorescence, observed_buffer_fluorescence);
hold on;
plot(ligand_concentrations, computed_ligand_fluorescence, 'k-');
%legend('observed', sprintf('fit (\\epsilon = %.0f cm^{-1} M^{-1})', extinction_coefficient), 'location', 'northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vandetanib

subplot(3,2,3);
ligand_name = 'vandetanib';
observed_mixture_fluorescence = data(5,:);
observed_ligand_fluorescence = data(6,:);
semilogx(ligand_concentrations, observed_protein_fluorescence, 'g.', ligand_concentrations, observed_ligand_fluorescence, 'k.', ligand_concentrations, observed_buffer_fluorescence, 'b.', ligand_concentrations, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
xlabel('ligand concentration (M)');
title(ligand_name);
ylabel('fluorescence (AU)');
legend(protein_name, ligand_name, 'buffer', sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
axis([1e-9 1e-5 0 1]);
set(gca, 'XTick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]);

[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_buffer_fluorescence, computed_mixture_fluorescence] = fit_dilution_data(well_volume, protein_concentration, ligand_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_buffer_fluorescence, observed_mixture_fluorescence)
hold on;
plot(ligand_concentrations, computed_protein_fluorescence, 'g-', ligand_concentrations, computed_ligand_fluorescence, 'k-', ligand_concentrations, computed_buffer_fluorescence, 'b-', ligand_concentrations, computed_mixture_fluorescence, 'r-');
ligand_name
Kd
title(sprintf('%s%s : Kd = %.3f nM, [P] = %.3f nM', ligand_name, regexp_string, Kd / 1e-9, protein_concentration / 1e-9));

subplot(3,2,4);
semilogx(ligand_concentrations, observed_ligand_fluorescence, 'k.', 'markersize', markersize); 
xlabel('ligand concentration (M)');
ylabel('fluorescence (AU)');
title(sprintf('%s%s', ligand_name, regexp_string));
axis([1e-9 1e-5 0 max(observed_ligand_fluorescence)]);
set(gca, 'XTick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]);

% Fit using inner filter effect model.
well_area = 0.3165 / 100 / 100; % well area (m^2)
path_length = well_volume / 1000 / well_area; % path length (m)
%[computed_ligand_fluorescence, computed_buffer_fluorescence, extinction_coefficient] = fit_ligand_inner_filter(well_volume, path_length, ligand_concentrations, observed_ligand_fluorescence, observed_buffer_fluorescence);
hold on;
plot(ligand_concentrations, computed_ligand_fluorescence, 'k-');
%legend('observed', sprintf('fit (\\epsilon = %.0f cm^{-1} M^{-1})', extinction_coefficient), 'location', 'northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% erlotinib

subplot(3,2,5);
ligand_name = 'erlotinib';
observed_mixture_fluorescence = data(7,:);
observed_ligand_fluorescence = data(8,:);

semilogx(ligand_concentrations, observed_protein_fluorescence, 'g.', ligand_concentrations, observed_ligand_fluorescence, 'k.', ligand_concentrations, observed_buffer_fluorescence, 'b.', ligand_concentrations, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
xlabel('ligand concentration (M)');
title(ligand_name);
ylabel('fluorescence (AU)');
legend(protein_name, ligand_name, 'buffer', sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
axis([1e-9 1e-5 0 1]);
set(gca, 'XTick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]);

[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_buffer_fluorescence, computed_mixture_fluorescence] = fit_dilution_data(well_volume, protein_concentration, ligand_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_buffer_fluorescence, observed_mixture_fluorescence)
hold on;
plot(ligand_concentrations, computed_protein_fluorescence, 'g-', ligand_concentrations, computed_ligand_fluorescence, 'k-', ligand_concentrations, computed_buffer_fluorescence, 'b-', ligand_concentrations, computed_mixture_fluorescence, 'r-');
ligand_name
Kd
title(sprintf('%s%s : Kd = %.3f nM, [P] = %.3f nM', ligand_name, regexp_string, Kd / 1e-9, protein_concentration / 1e-9));

subplot(3,2,6);
semilogx(ligand_concentrations, observed_ligand_fluorescence, 'k.', 'markersize', markersize); 
xlabel('ligand concentration (M)');
ylabel('fluorescence (AU)');
title(sprintf('%s%s', ligand_name, regexp_string));
axis([1e-9 1e-5 0 max(observed_ligand_fluorescence)]);
set(gca, 'XTick', [1e-9, 1e-8, 1e-7, 1e-6, 1e-5]);

% Fit using inner filter effect model.
well_area = 0.3165 / 100 / 100; % well area (m^2)
path_length = well_volume / 1000 / well_area; % path length (m)
%[computed_ligand_fluorescence, computed_buffer_fluorescence, extinction_coefficient] = fit_ligand_inner_filter(well_volume, path_length, ligand_concentrations, observed_ligand_fluorescence, observed_buffer_fluorescence);
hold on;
plot(ligand_concentrations, computed_ligand_fluorescence, 'k-');
%legend('observed', sprintf('fit (\\epsilon = %.0f cm^{-1} M^{-1})', extinction_coefficient), 'location', 'northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('serial-dilution%s.eps', regexp_string);
exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));

%filename = sprintf('serial-dilution%s.png', regexp_string);
%print('-dpng', '-r300', filename);

