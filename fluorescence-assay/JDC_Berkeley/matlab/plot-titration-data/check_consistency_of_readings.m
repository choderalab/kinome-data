% plot serial dilution data

clear;

% Load data.
filename = '../../2012-05-09/2012-05-09 re-read plate A1 from 2012-05-07.txt';

infile = fopen(filename, 'r');

%regexp_string = ':297,328';
regexp_string = ':280,340';
%regexp_string = ':280,440';
%regexp_string = ':280,460';
%regexp_string = ':280,480';
%regexp_string = ':350,440';
%regexp_string = ':350,460';
%regexp_string = ':350,480';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

imagesc(data);
colorbar;
