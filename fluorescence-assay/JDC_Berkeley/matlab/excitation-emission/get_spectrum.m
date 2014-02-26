function [wavelengths, spectrum] = get_spectrum(infile, regexp_string, nrows, ncols)

% Read wavelengths
wavelengths = [];
frewind(infile);
while ~feof(infile)
  line = fgetl(infile);
  S = regexp(line, regexp_string);
  if length(S) > 0
    % Read wavelength line.
    line = fgetl(infile);
    S = regexp(line, '^Wavelength \d+ \((\d+) nm\)', 'tokens');
    wavelength = sscanf(char(S{1}(1)), '%d');
    wavelengths = [wavelengths wavelength];        
  end
end
nwavelengths = length(wavelengths);

% Allocate spectrum data.
spectrum = zeros(nrows, ncols, nwavelengths);
frewind(infile);
while ~feof(infile)
  line = fgetl(infile);
  S = regexp(line, regexp_string);
  if length(S) > 0
    % Read wavelength line.
    line = fgetl(infile);
    S = regexp(line, '^Wavelength (\d+) \((\d+) nm\)', 'tokens');
    wavelength_index = sscanf(char(S{1}(1)), '%d');
    % Skip header.
    line = fgetl(infile);
    % Read data.
    for row = 1:nrows
      line = fgetl(infile);
      tokens = regexp(line, '\t', 'split');
      for col = 1:ncols
	token = tokens{col+1};
	value = sscanf(token, '%f');
	% TODO: Handle overflow.
	spectrum(row, col, wavelength_index) = value;
      end
    end
  end
end

return
