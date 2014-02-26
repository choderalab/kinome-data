% plot excitation and emission spectrum

clear;

% Load data.
filename = 'Emission spectrum - multiple excitation wavelengths_110610_170412.txt';
infile = fopen(filename, 'r');

excitation = [270 280 290 300 350];
emission = 350:5:600;
nexcitation = length(excitation);
nemission = length(emission);
nrows = 8;
data = zeros(nrows,nexcitation,nemission);

i = 1;
j = 1;
while ~feof(infile)
  line = fgetl(infile);

  regexp_string = 'excite (\d+):EM Spectrum';
  S = regexp(line, regexp_string);
  if length(S) > 0
    S = regexp(line, regexp_string, 'tokens');
    excitation_wavelength = sscanf(char(S{1}), '%f');
    continue
  end

  regexp_string = 'Wavelength \d+ \((\d+) nm\)';
  S = regexp(line, regexp_string);  
  if length(S) > 0
    S = regexp(line, regexp_string, 'tokens');
    emission_wavelength = sscanf(char(S{1}), '%f');

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

    j = j + 1;
    if j > nemission
      j = 1;
      i = i + 1;
    end
    continue
  end
end
fclose(infile);

data = data(:,1:4,:);
excitation = excitation(1:4);

% normalize
%for i = 1:4
%  data(:,i,:) = data(:,i,:) / max(max(squeeze(data(:,i,:))));
%end
