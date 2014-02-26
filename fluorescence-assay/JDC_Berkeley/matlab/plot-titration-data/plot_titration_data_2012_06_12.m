% Plot titration data from 2012-06-12 experiment.

clear;

% Load data.

titrant_dilution_indices = [6 5 4 3 2 1 1];
ninjections = length(titrant_dilution_indices); % number of injections
titrant_concentrations = 10e-6 * 2.^(-(titrant_dilution_indices - 1)); % concentrations of titrant (M)
titrant_volumes = 10e-6 * ones(1, ninjections); % titrant volumes (L)
protein_volume = 100e-6; % initial protein volume (L)
protein_concentration = 250e-9; % initial protein concentration (M)

%ligand_row = 2; ligand_name = 'quinazoline';
%ligand_row = 3; ligand_name = 'staurosporine';
ligand_row = 4; ligand_name = 'bosutinib';
%ligand_row = 5; ligand_name = 'vandetinib';
%ligand_row = 6; ligand_name = 'erlotinib';
%ligand_row = 7; ligand_name = 'gefitinib';
%ligand_row = 8; ligand_name = 'imatinib';

nrows = 8; % number of rows
ncolumns = 12; % number of columns
data = zeros(ninjections+1, nrows, ncolumns);
for injection_index = 0:ninjections
  injection_index
  
  %filename = sprintf('../../2012-06-12/2012-06-12 plate C2 titration %d.txt', injection_index);
  filename = sprintf('../../2012-06-12/2012-06-12 plate C1 titration %d.txt', injection_index);
  
  infile = fopen(filename, 'r');
  
  % Skip to second part of file (weird concatenation)

  for section = 1:2
    regexp_string = 'Procedure Summary';
    while ~feof(infile)
      line = fgetl(infile);    
      S = regexp(line, regexp_string); 
      if length(S) > 0
	% match found
	break
      end
    end
  end

  %regexp_string = 'F280/480 100 top:280,480'; % data to look for
  regexp_string = 'F280/480 130 top:280,480'; % data to look for

  i = 1;
  j = 1;
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
	% TODO: Break up line to handle OVRFLW

	% Store.
	data(injection_index+1, row, :) = row_data;
      end      

    end
  end
  fclose(infile);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

buffer_row = 1;

markersize = 10;

injection_volumes = 10e-6 * ones(1,10); % L
injection_concentrations = [1e-9, 10e-9, 100e-9, 1e-6, 10e-6, 100e-6, 100e-6, 100e-6, 100e-6, 100e-6]; % M
initial_volume = 150e-6; % L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buffer_column = 1;
initial_protein_concentration_estimate = 1e-6; % M

subplot(3,3,1);
protein_column = 2; 
protein_name = 'Src';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1
observed_ligand_fluorescence = (ligand ./ buffer)' - 1
observed_mixture_fluorescence = (signal ./ buffer)' - 1
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
title('1 uM');
set(gca, 'XTick', 0:10);
ylabel('relative fluorescence');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('1 uM : Kd = %.1e M', Kd));

subplot(3,3,2);
protein_column = 3; 
protein_name = 'p38';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
set(gca, 'XTick', 0:10);
title('1 uM');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.1e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('1 uM : Kd = %.1e M', Kd));

subplot(3,3,3);
protein_column = 4; 
protein_name = 'hCAII';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
set(gca, 'XTick', 0:10);
title('1 uM');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.1e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('1 uM : Kd = %.1e M', Kd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buffer_column = 5;
initial_protein_concentration_estimate = 1e-6 / 2.5; % M

subplot(3,3,4);
protein_column = 6; 
protein_name = 'Src';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
title('400 nM');
set(gca, 'XTick', 0:10);
ylabel('relative fluorescence');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %.1e M, [P]_0 = %.1e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('400 nM : Kd = %.1e M', Kd));

subplot(3,3,5);
protein_column = 7; 
protein_name = 'p38';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
set(gca, 'XTick', 0:10);
title('400 nM');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.3e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('400 nM : Kd = %.1e M', Kd));

subplot(3,3,6);
protein_column = 8; 
protein_name = 'hCAII';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
set(gca, 'XTick', 0:10);
title('400 nM');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.3e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('400 nM : Kd = %.1e M', Kd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buffer_column = 9;
initial_protein_concentration_estimate = 1e-6 / 2.5^2; % M

subplot(3,3,7);
protein_column = 10; 
protein_name = 'Src';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
title('160 nM');
xlabel('injection');
set(gca, 'XTick', 0:10);
ylabel('relative fluorescence');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.3e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('160 nM : Kd = %.1e M', Kd));

subplot(3,3,8);
protein_column = 11; 
protein_name = 'p38';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
title('160 nM');
set(gca, 'XTick', 0:10);
xlabel('injection');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.3e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('160 nM : Kd = %.1e M', Kd));

subplot(3,3,9);
protein_column = 12; 
protein_name = 'hCAII';
signal = squeeze(data(:,ligand_row,protein_column)); buffer = squeeze(data(:,buffer_row,buffer_column)); protein = squeeze(data(:,buffer_row,protein_column)); ligand = squeeze(data(:,ligand_row,buffer_column)); 
observed_protein_fluorescence = (protein ./ buffer)' - 1;
observed_ligand_fluorescence = (ligand ./ buffer)' - 1;
observed_mixture_fluorescence = (signal ./ buffer)' - 1;
plot(0:10, observed_protein_fluorescence, 'g.', 0:10, observed_ligand_fluorescence, 'k.', 0:10, observed_mixture_fluorescence, 'r.', 'markersize', markersize); 
legend(protein_name, ligand_name, sprintf('%s:%s', protein_name, ligand_name), 'location', 'best');
title('160 nM');
set(gca, 'XTick', 0:10);
xlabel('injection');

% Fit data.
[Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence);
disp(sprintf('Kd = %e M, [P]_0 = %.3e M', Kd, protein_concentration));
% Plot fit.
hold on;
plot(0:10, computed_protein_fluorescence, 'g-');
plot(0:10, computed_ligand_fluorescence, 'k-');
plot(0:10, computed_mixture_fluorescence, 'r-');
title(sprintf('160 nM : Kd = %.1e M', Kd));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('%s-%s.eps', ligand_name, regexp_string);
exportfig(gcf, filename, 'width', 12, 'height', 12, 'color', 'cmyk');
system(sprintf('epstopdf %s', filename));


