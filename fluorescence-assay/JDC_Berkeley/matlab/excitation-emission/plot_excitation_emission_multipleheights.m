% plot excitation and emission spectrum at multiple heights

clear;

% Load data.
%filename = 'Emission spectrum - multiple excitation wavelengths_110610_170412.txt';
%filename = '../../2012-05-09/2012-05-09 plate A4 spectra - adjusted gain.txt';
filename = '../../2012-05-09/2012-05-09 plate A4 spectra - varying read height.txt';
infile = fopen(filename, 'r');

%excitation = [270 280 290 300 350];
%emission = 350:5:600;

titles = {'200 nM Src', 'buffer', '200 nM Src + 5 uM bosutinib', '5 uM bosutinib', '200 nM Src + 5 uM vandetanib', '5 uM vandetanib', '200 nM Src + 5 uM erlotinib', '5 uM erlotinib'};

% TODO: Discover excitation and emission wavelengths from file automatically.
excitation_wavelengths = [280];
emission_wavelengths = 310:5:550;
heights = [4 6 8 10]; % read heights in mm

nexcitation = length(excitation_wavelengths);
nemission = length(emission_wavelengths);
nrows = 8;
nheights = length(heights);
data = zeros(nheights,nrows,nexcitation,nemission);

height_index = 0;
i = 1;
j = 1;
while ~feof(infile)
  line = fgetl(infile);

  regexp_string = 'excite (\d+) nm bandwidth (\d+) nm height (\d+) mm:EM Spectrum';
  S = regexp(line, regexp_string);
  if length(S) > 0
    S = regexp(line, regexp_string, 'tokens');
    excitation_wavelength = sscanf(char(S{1}(1)), '%d');
    i = find(excitation_wavelengths == excitation_wavelength);
    bandwidth = sscanf(char(S{1}(2)), '%d');
    height = sscanf(char(S{1}(3)), '%d');
    height_index = find(heights == height); 
    continue
  end

  regexp_string = 'Wavelength \d+ \((\d+) nm\)';
  S = regexp(line, regexp_string);  
  if length(S) > 0
    S = regexp(line, regexp_string, 'tokens');
    emission_wavelength = sscanf(char(S{1}), '%d');

    j = find(emission_wavelengths == emission_wavelength);

    % Read subsequent data
    line = fgetl(infile); % skip

    for row = 1:nrows
      line = fgetl(infile); 
      value_string = line(3:8);
      if value_string == 'OVRFLW'
	value = NaN;
      else
	value = sscanf(value_string, '%f');
      end

      % Store.
      data(height_index,row,i,j) = value;
    end      
  end
end
fclose(infile);

% set colormap
nmap = max(max(max(data)));
cmap = zeros(nmap+1,3);
cmap(2:(nmap+1),:) = jet(nmap);
colormap(cmap);

ny = 4;
nx = 3;
clf;
for i = 1:4
  subplot(ny, nx, nx*(i-1)+1);
  row = 2*(i-1)+1;
  imagesc(emission_wavelengths, heights, squeeze(data(:,row,1,:)));
  title(titles(2*(i-1)+1));
  set(gca, 'ytick', heights);
  if i == 4
    xlabel('emission / nm');
  end
  ylabel('read height / mm');
  colorbar;

  subplot(ny, nx, nx*(i-1)+2);
  row = 2*(i-1)+2;
  imagesc(emission_wavelengths, heights, squeeze(data(:,row,1,:)));
  title(titles(2*(i-1)+2));
  set(gca, 'ytick', heights);
  if i == 4
    xlabel('emission / nm');
  end
  colorbar;

  subplot(ny, nx, nx*(i-1)+3);
  imagesc(emission_wavelengths, excitation_wavelengths, squeeze(data(:,row,1,:)) - squeeze(data(:,row-1,1,:)));
  title(sprintf('%s - %s', char(titles(2*(i-1)+2)), char(titles(2*(i-1)+1))));
  set(gca, 'ytick', heights);
  if i == 4
    xlabel('emission / nm');
  end
  colorbar;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filename = sprintf('emission-spectra-280nm-excitation-readheights.eps');
%exportfig(gcf, filename, 'width', 10, 'height', 10, 'color', 'cmyk');
%system(sprintf('epstopdf %s', filename));

filename = sprintf('excitation-emission-readheights.png');
print('-dpng', '-r300', filename);

