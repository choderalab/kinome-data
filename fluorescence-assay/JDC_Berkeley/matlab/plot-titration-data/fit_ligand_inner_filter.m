function [computed_ligand_fluorescence, computed_buffer_fluorescence, extinction_coefficient] = fit_ligand_inner_filter(well_volume, path_length, ligand_concentrations, observed_ligand_fluorescence, observed_buffer_fluorescence, mode)
% Fit fluorescence of various fluorophore concentrations accounting for "primary inner filter effect", or nonlinearity of Beer's law.
%
% ARGUMENTS
%
% 

% Determine number of wells.
nwells = length(ligand_concentrations);

% Determine ligand quantities.
ligand_quantities = ligand_concentrations * well_volume;

% Determine well volumes.
well_volumes = well_volume * ones(1,nwells);

if strcmp(mode, 'nonlinear')==1
  % Define nonlinear equation for observed fluorescence accounting for primary inner filter effect.
  % 'IQ' includes intensity of incident light and quantum efficiency
  % From Eq. 1 of [1].
  % [1] Kao S, Asanov AN, and Oldham PB. A comparison of fluorescence inner-filter effects for different cell configurations. Instrumentation Science & Technology 26;375, 1998.
  observed_fluorescence = @(concentration, IQ, extinction_coefficient, path_length) IQ * (1 - exp(-extinction_coefficient * path_length * concentration));
elseif strcmp(mode, 'linear')==1
  % Linearized model.
  observed_fluorescence = @(concentration, IQ, extinction_coefficient, path_length) IQ * extinction_coefficient * path_length * concentration;
else
  error(sprintf('mode %s unknown', mode));
end

% Define an objective function on observed fluorescence data, including inner filter effect.
lsq_fit_error = @(ligand_fluorescence, buffer_fluorescence, extinction_coefficient) sum((observed_buffer_fluorescence - buffer_fluorescence).^2) + sum((observed_ligand_fluorescence - (observed_fluorescence(ligand_concentrations, ligand_fluorescence, extinction_coefficient, path_length) + buffer_fluorescence)).^2);

objective = @(x) lsq_fit_error(exp(x(1)), exp(x(2)), exp(x(3)));

% Create initial guess.
buffer_fluorescence_estimate = mean(observed_buffer_fluorescence);
ligand_fluorescence_estimate = abs(sum(observed_ligand_fluorescence - buffer_fluorescence_estimate) / sum(ligand_quantities));
extinction_coefficient_estimate = 75000/(1e-2);

% Solve for unknown parameters.
options = optimset('maxfunevals', 1e4, 'maxiter', 1e4);
x = log([ligand_fluorescence_estimate buffer_fluorescence_estimate extinction_coefficient_estimate]);
x = fminsearch(objective, x, options);

% Extract solution.
ligand_fluorescence = exp(x(1));
buffer_fluorescence = exp(x(2));
extinction_coefficient = exp(x(3));

% Compute model-based fluorecence to compare with observed fluorescence.
computed_ligand_fluorescence = observed_fluorescence(ligand_concentrations, ligand_fluorescence, extinction_coefficient, path_length) + buffer_fluorescence;
computed_buffer_fluorescence = buffer_fluorescence * ones(1,nwells);

% Convert extinction coefficient to cm^-1 M^-1
extinction_coefficient = extinction_coefficient * 1e-2;

return

