% plot dilution series check 

clear;

% Load data.
filename = '../../2012-05-09/2012-05-09 Phoenix 1-10 dilution series check.txt';

infile = fopen(filename, 'r');

%wavelengths = '350,440';
%wavelengths = '350,460';
wavelengths = '350,480';

nplates = 6; % number of plates, labeled D1 - D6
nrows = 8; % number of rows
ncolumns = 12; % number of columns
data = zeros(nplates, nrows, ncolumns);

while ~feof(infile)
  line = fgetl(infile);

  regexp_string = sprintf('D(\\d+):%s', wavelengths);
  S = regexp(line, regexp_string);  
  if length(S) > 0
    % Match found.
    S = regexp(line, regexp_string, 'tokens');
    plate_index = sscanf(char(S{1}), '%d');

    % Read subsequent data
    line = fgetl(infile); % skip header line

    for row = 1:nrows
      line = fgetl(infile); 
      row_data = sscanf(line, '%*c\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');

      % Store.
      data(plate_index, row, :) = row_data;
    end      

  end
end
fclose(infile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;

markersize = 10;


% set colormap
nmap = max(max(max(data)));
cmap = jet(nmap);
colormap(cmap);

concentrations = [10e-6, 1e-6, 100e-9, 10e-9, 1e-9, 100e-10];

ligand_names = {'buffer', 'quinazoline', 'staurosporine', 'bosutinib', 'vandetanib', 'erlotinib', 'gefitinib', 'imatinib'};

ny = 4;
nx = 2;
for row = 1:8
  subplot(ny, nx, row);
  
  %x = log10(squeeze(data(:,row,:)));
  a = squeeze(data(:,row,:));
  b = squeeze(data(:,1,:));
  x = a;
  imagesc(1:12, 1:6, x);
  colorbar;
  if floor((row+1)/2) == 4
    xlabel('well (all should be same)');
  end
  ylabel('dilution index');
  title(char(ligand_names{row}));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filename = sprintf('emission-spectra-280nm-excitation-readheights.eps');
%exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
%system(sprintf('epstopdf %s', filename));

filename = sprintf('plate-ligand-dilution-series-check:%s.png', wavelengths);
print('-dpng', '-r450', filename);

