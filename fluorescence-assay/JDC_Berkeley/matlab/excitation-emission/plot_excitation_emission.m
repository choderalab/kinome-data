% plot excitation and emission spectrum

clear;

% Load data.
filename = 'Emission spectrum - multiple excitation wavelengths_110610_170412.txt';
filename = '../../2012-05-09/2012-05-09 plate A4 spectra - adjusted gain.txt';
infile = fopen(filename, 'r');

%excitation = [270 280 290 300 350];
% emission = 350:5:600;

titles = {'200 nM Src', 'buffer', '200 nM Src + 5 uM bosutinib', '5 uM bosutinib', '200 nM Src + 5 uM vandetanib', '5 uM vandetanib', '200 nM Src + 5 uM erlotinib', '5 uM erlotinib'};

% TODO: Discover excitation and emission wavelengths from file automatically.
excitation_wavelengths = [260 270 280 290 300 310];
emission_wavelengths = 290:5:550;

nexcitation = length(excitation_wavelengths);
nemission = length(emission_wavelengths);
nrows = 8;
data = zeros(nrows,nexcitation,nemission);

i = 1;
j = 1;
while ~feof(infile)
  line = fgetl(infile);

  regexp_string = 'excite (\d+) nm bandwidth (\d+)\s*\w*:EM Spectrum';
  S = regexp(line, regexp_string);
  if length(S) > 0
    S = regexp(line, regexp_string, 'tokens');
    excitation_wavelength = sscanf(char(S{1}(1)), '%d');
    i = find(excitation_wavelengths == excitation_wavelength);
    bandwidth = sscanf(char(S{1}(2)), '%d');
    continue
  end

  regexp_string = 'Wavelength \d+ \((\d+) nm\)';
  S = regexp(line, regexp_string);  
  if length(S) > 0
    S = regexp(line, regexp_string, 'tokens');
    emission_wavelength = sscanf(char(S{1}), '%d');

    if bandwidth ~= 9 
      continue
    end

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
      data(row,i,j) = value;
    end      
  end
end
fclose(infile);

% set colormap
nmap = max(max(max(data)));
cmap = zeros(nmap+1,3);
cmap(2:(nmap+1),:) = jet(nmap);
colormap(cmap);

%data = data(:,1:4,:);
%excitation = excitation(1:4);

% normalize
%for i = 1:4
%  data(:,i,:) = data(:,i,:) / max(max(squeeze(data(:,i,:))));
%end

%ny = 4;
%nx = 3;
%clf;
%for row = 1:nrows
%  subplot(ny,nx,row);
%  %imagesc(excitation_wavelengths, emission_wavelengths, squeeze(data(row,:,:)));
%  imagesc(emission_wavelengths, excitation_wavelengths, squeeze(data(row,:,:)));
%  title(titles{row});
%end

hfigure1 = figure(1)
buffer_row = 2;
ny = 4;
nx = 4;
clf;
for i = 1:4
  mixture_row = 2*(i-1)+1;
  ligand_row = 2*(i-1)+2;

  subplot(ny, nx, nx*(i-1)+1);
  imagesc(emission_wavelengths, excitation_wavelengths, squeeze(data(mixture_row,:,:)));
  title(titles(2*(i-1)+1));
  set(gca, 'ytick', excitation_wavelengths);
  if i == 4
    xlabel('emission / nm');
  end
  ylabel('excitation / nm');
  colorbar;

  subplot(ny, nx, nx*(i-1)+2);
  imagesc(emission_wavelengths, excitation_wavelengths, squeeze(data(ligand_row,:,:)));
  title(titles(2*(i-1)+2));
  set(gca, 'ytick', excitation_wavelengths);
  if i == 4
    xlabel('emission / nm');
  end
  colorbar;

  subplot(ny, nx, nx*(i-1)+3);
  imagesc(emission_wavelengths, excitation_wavelengths, squeeze(data(mixture_row,:,:)) - squeeze(data(ligand_row,:,:)));
  title(sprintf('(%s) - (%s)', char(titles(mixture_row)), char(titles(ligand_row))));
  set(gca, 'ytick', excitation_wavelengths);
  if i == 4
    xlabel('emission / nm');
  end
  colorbar;

  subplot(ny, nx, nx*(i-1)+4);
  imagesc(emission_wavelengths, excitation_wavelengths, squeeze(data(ligand_row,:,:)) - squeeze(data(buffer_row,:,:)));
  title(sprintf('(%s) - (%s)', char(titles(ligand_row)), char(titles(buffer_row))));
  set(gca, 'ytick', excitation_wavelengths);
  if i == 4
    xlabel('emission / nm');
  end
  colorbar;
end

filename = 'excitation-emission-images.eps';
exportfig(hfigure1, filename, 'width', 22, 'height', 17, 'color', 'cmyk');
system(sprintf('epstopdf "%s"', filename));


hfigure2 = figure(2)
buffer_row = 2;
ny = 4;
nx = 4;
clf;

legend_text = cell(1,nexcitation);
for i = 1:nexcitation
  legend_text{i} = sprintf('%d nm', excitation_wavelengths(i));
end

colors = {'r', 'g', 'm', 'c', 'y', 'k'};
ncolors = length(colors);
linestyles = {'-', ':', '--', '-.'};
nlinestyles = length(linestyles);

for i = 1:4
  mixture_row = 2*(i-1)+1;
  ligand_row = 2*(i-1)+2;

  subplot(ny, nx, nx*(i-1)+1);
  hold on;
  spectra = squeeze(data(mixture_row,:,:));
  for j = 1:nexcitation
    spectrum = spectra(j,:);
    indices = find(spectrum ~= 0);
    plot(emission_wavelengths(indices), spectrum(indices), sprintf('%s%s', colors{mod(j, ncolors)+1}, linestyles{mod(j, nlinestyles)+1}));  
  end    
  axis([emission_wavelengths(1) emission_wavelengths(end) min(min(spectra)) max(max(spectra))]);
  title(titles(2*(i-1)+1));
  if i == 4
    xlabel('emission / nm');
  end
  ylabel('fluorescence (AU)');
  set(gca, 'ytick', []);

  subplot(ny, nx, nx*(i-1)+2);
  hold on;
  spectra = squeeze(data(ligand_row,:,:));
  for j = 1:nexcitation
    spectrum = spectra(j,:);
    indices = find(spectrum ~= 0);
    plot(emission_wavelengths(indices), spectrum(indices), sprintf('%s%s', colors{mod(j, ncolors)+1}, linestyles{mod(j, nlinestyles)+1}));  
  end    
  axis([emission_wavelengths(1) emission_wavelengths(end) min(min(spectra)) max(max(spectra))]);
  title(titles(2*(i-1)+2));
  if i == 4
    xlabel('emission / nm');
  end
  set(gca, 'ytick', []);

  subplot(ny, nx, nx*(i-1)+3);
  hold on;
  spectra = squeeze(data(mixture_row,:,:)) - squeeze(data(ligand_row,:,:));
  for j = 1:nexcitation
    spectrum = spectra(j,:);
    indices = find(spectrum ~= 0);
    plot(emission_wavelengths(indices), spectrum(indices), sprintf('%s%s', colors{mod(j, ncolors)+1}, linestyles{mod(j, nlinestyles)+1}));  
  end    
  title(sprintf('(%s) - (%s)', char(titles(mixture_row)), char(titles(ligand_row))));
  if i == 4
    xlabel('emission / nm');
  end
  axis([emission_wavelengths(1) emission_wavelengths(end) min(min(spectra)) max(max(spectra))]);
  set(gca, 'ytick', []);
  if i == 1
    legend(legend_text, 'location', 'northeast');
  end

  if ligand_row ~= buffer_row
    subplot(ny, nx, nx*(i-1)+4);
    hold on;
    spectra = squeeze(data(ligand_row,:,:)) - squeeze(data(buffer_row,:,:));
    for j = 1:nexcitation
      spectrum = spectra(j,:);
      indices = find(spectrum ~= 0);
      plot(emission_wavelengths(indices), spectrum(indices), sprintf('%s%s', colors{mod(j, ncolors)+1}, linestyles{mod(j, nlinestyles)+1}));  
    end    
    axis([emission_wavelengths(1) emission_wavelengths(end) min(min(spectra)) max(max(spectra))]);
    title(sprintf('(%s) - (%s)', char(titles(ligand_row)), char(titles(buffer_row))));
    if i == 4
      xlabel('emission / nm');
    end
    set(gca, 'ytick', []);
  end
end

filename = 'excitation-emission-lines.eps';
exportfig(hfigure2, filename, 'width', 22, 'height', 17, 'color', 'cmyk');
system(sprintf('epstopdf "%s"', filename));
%print -dpng -r300 excitation-emission-lines.png

