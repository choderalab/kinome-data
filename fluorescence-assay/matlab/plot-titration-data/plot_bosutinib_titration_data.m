% plot excitation and emission spectrum

clear;

% Load data.
filename = '../../2012-05-15/2012-05-15 plate B3 titration.txt';

infile = fopen(filename, 'r');

fluorescence_wavelengths = '280,480';  % top
%fluorescence_wavelengths = '280,480\[2\]';  % bottom

ninjections = 12; % number of injections
nrows = 4; % number of rows
ncolumns = 6; % number of columns
data = zeros(ninjections+1, nrows, ncolumns);

i = 1;
j = 1;
injection = 0;
while ~feof(infile)
  line = fgetl(infile);

  regexp_string = sprintf('^injection (\\d+) fluorescence:%s$', fluorescence_wavelengths);
  S = regexp(line, regexp_string);  
  if length(S) > 0
    % Match found.

    % Read subsequent data
    line = fgetl(infile); % skip header line

    for row = 1:nrows
      line = fgetl(infile); 
      disp(sprintf('%1d > %s', row, line));
      tokens = regexp(line, '\W', 'split');

      for column = 1:ncolumns
	value_string = char(tokens{column+1});
	if strcmp(value_string, 'OVRFLW')==1
	  value = NaN;
	else
	  value = sscanf(value_string, '%f');
	end
	
	% Store.
	data(injection+1, row, column) = value;
      end
    end      
    
    injection = injection + 1;
  end
end
fclose(infile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

buffer_columns = [2, 4, 6];

ligand_name = 'bosutinib';
nreplicates = 4;

markersize = 10;

well_area = 0.1586 / 100 / 100; % m^2

initial_volume = 75e-6; % L
injection_volumes = 10e-6 * ones(1,12); % L
injection_concentrations = 20e-6 * 2.^(-(11:-1:0)); % M
well_volumes = cumsum([initial_volume injection_volumes]);
path_lengths = (well_volumes * 1000 / (100^3)) / well_area ; % m
total_ligand_quantities = cumsum([0 injection_volumes.*injection_concentrations]);
total_ligand_concentrations = total_ligand_quantities ./ well_volumes;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_protein_concentration_estimate = 1e-6; % M

subplot(3,3,1);
hold on;
protein_column = 1; 
buffer_column = 2;
protein_name = 'Src';
for replicate = 1:nreplicates
  observed_mixture_fluorescence = squeeze(data(:,replicate,protein_column))'; 
  observed_ligand_fluorescence = squeeze(data(:,replicate,buffer_column))'; 
  indices = find(~isnan(observed_mixture_fluorescence));
  ninjections_fit = length(indices)-1;

  plot(0:ninjections, observed_ligand_fluorescence, 'k.', 0:ninjections, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
  
  % Fit data.
  %[Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data3(initial_volume, initial_protein_concentration_estimate, injection_volumes(1,1:(length(indices)-1)), injection_concentrations(1,1:(length(indices)-1)), observed_ligand_fluorescence(1,indices), observed_mixture_fluorescence(1,indices));
  %[Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data4(initial_volume, initial_protein_concentration_estimate, injection_volumes(1,1:(length(indices)-1)), injection_concentrations(1,1:(length(indices)-1)), observed_ligand_fluorescence(1,indices), observed_mixture_fluorescence(1,indices));
  [Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data5(initial_volume, initial_protein_concentration_estimate, injection_volumes(1,1:(length(indices)-1)), injection_concentrations(1,1:(length(indices)-1)), observed_ligand_fluorescence(1,indices), observed_mixture_fluorescence(1,indices));
  disp(sprintf('replicate %d / %d : Kd = %f nM, [P]_0 = %f nM', replicate, nreplicates, Kd / 1e-9, protein_concentration / 1e-9));
  % Plot fit.
  plot(0:ninjections_fit, computed_mixture_fluorescence, 'r:');
  plot(0:ninjections_fit, computed_ligand_fluorescence, 'k:');
  title(sprintf('1 uM : Kd = %.1e M', Kd));  
end
legend('bosutinib into buffer', 'bosutinib into protein', 'location', 'northwest');
title('1 uM');
set(gca, 'XTick', 0:ninjections);
ylabel('fluorescence');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);

subplot(3,3,2);
hold on;
for replicate = 1:nreplicates
  observed_mixture_fluorescence = squeeze(data(:,replicate,protein_column))'; 
  plot(0:ninjections, observed_mixture_fluorescence, 'r.-', 'markersize', markersize); 
end
title('mixture zoomed');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);
set(gca, 'XTick', 0:ninjections);
subplot(3,3,3);
hold on;
for replicate = 1:nreplicates
  observed_ligand_fluorescence = squeeze(data(:,replicate,buffer_column))'; 
  plot(0:ninjections, observed_ligand_fluorescence, 'k.-', 'markersize', markersize);
end
title('ligand zoomed');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);
set(gca, 'XTick', 0:ninjections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_protein_concentration_estimate = 400e-9; % M

% DEBUG: Fix one NaN
data(5,1,4) = 0.5 * (data(4,1,4) + data(6,1,4));

subplot(3,3,4);
hold on;
protein_column = 3; 
buffer_column = 4;
protein_name = 'Src';
for replicate = 1:nreplicates
  observed_mixture_fluorescence = squeeze(data(:,replicate,protein_column))';
  observed_ligand_fluorescence = squeeze(data(:,replicate,buffer_column))';
  indices = find(~isnan(observed_mixture_fluorescence));
  ninjections_fit = length(indices)-1;

  plot(0:ninjections, observed_ligand_fluorescence, 'k.', 0:ninjections, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
  
  % Fit data.
  [Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data5(initial_volume, initial_protein_concentration_estimate, injection_volumes(1,1:(length(indices)-1)), injection_concentrations(1,1:(length(indices)-1)), observed_ligand_fluorescence(1,indices), observed_mixture_fluorescence(1,indices));
  disp(sprintf('replicate %d / %d : Kd = %f nM, [P]_0 = %f nM', replicate, nreplicates, Kd / 1e-9, protein_concentration / 1e-9));
  % Plot fit.
  plot(0:ninjections_fit, computed_mixture_fluorescence, 'r:');
  plot(0:ninjections_fit, computed_ligand_fluorescence, 'k:');
  title(sprintf('1 uM : Kd = %.1e M', Kd));
end
title('400 nM');
set(gca, 'XTick', 0:ninjections);
ylabel('fluorescence');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);

subplot(3,3,5);
hold on;
for replicate = 1:nreplicates
  observed_mixture_fluorescence = squeeze(data(:,replicate,protein_column))'; 
  plot(0:ninjections, observed_mixture_fluorescence, 'r.-', 'markersize', markersize); 
end
title('mixture');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);
set(gca, 'XTick', 0:ninjections);
subplot(3,3,6);
hold on;
for replicate = 1:nreplicates
  observed_ligand_fluorescence = squeeze(data(:,replicate,buffer_column))'; 
  plot(0:ninjections, observed_ligand_fluorescence, 'k.-', 'markersize', markersize);
end
title('ligand');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);
set(gca, 'XTick', 0:ninjections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initial_protein_concentration_estimate = 160e-9; % M

subplot(3,3,7);
hold on;
protein_column = 5; 
buffer_column = 6;
protein_name = 'Src';
for replicate = 1:nreplicates
  observed_mixture_fluorescence = squeeze(data(:,replicate,protein_column))'; 
  observed_ligand_fluorescence = squeeze(data(:,replicate,buffer_column))'; 
  indices = find(~isnan(observed_mixture_fluorescence));
  ninjections_fit = length(indices)-1;

  plot(0:ninjections, observed_ligand_fluorescence, 'k.', 0:ninjections, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
  
  % Fit data.
  [Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data5(initial_volume, initial_protein_concentration_estimate, injection_volumes(1,1:(length(indices)-1)), injection_concentrations(1,1:(length(indices)-1)), observed_ligand_fluorescence(1,indices), observed_mixture_fluorescence(1,indices));
  disp(sprintf('replicate %d / %d : Kd = %f nM, [P]_0 = %f nM', replicate, nreplicates, Kd / 1e-9, protein_concentration / 1e-9));
  % Plot fit.
  plot(0:ninjections_fit, computed_mixture_fluorescence, 'r:');
  plot(0:ninjections_fit, computed_ligand_fluorescence, 'k:');
  title(sprintf('1 uM : Kd = %.1e M', Kd));
end
title('160 nM');
set(gca, 'XTick', 0:ninjections);
ylabel('fluorescence');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);

subplot(3,3,8);
hold on;
for replicate = 1:nreplicates
  observed_mixture_fluorescence = squeeze(data(:,replicate,protein_column))'; 
  plot(0:ninjections, observed_mixture_fluorescence, 'r.-', 'markersize', markersize); 
end
title('mixture');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);
set(gca, 'XTick', 0:ninjections);
subplot(3,3,9);
hold on;
for replicate = 1:nreplicates
  observed_ligand_fluorescence = squeeze(data(:,replicate,buffer_column))'; 
  plot(0:ninjections, observed_ligand_fluorescence, 'k.-', 'markersize', markersize);
end
title('ligand');
oldaxis = axis; axis([0 ninjections 0 oldaxis(4)]);
set(gca, 'XTick', 0:ninjections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('2012-05-15 bosutinib titration-top.eps');
exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
system(sprintf('epstopdf "%s"', filename));
