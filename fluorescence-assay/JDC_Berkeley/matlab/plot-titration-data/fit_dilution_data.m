function [Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_buffer_fluorescence, computed_mixture_fluorescence] = fit_dilution_data(well_volume, protein_concentration_estimate, ligand_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_buffer_fluorescence, observed_mixture_fluorescence)
% Fit the observed titration data with a 1:1 binding model.
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

% Solve for complex concentration given total ligand and protein concentrations.
% Kd = P*L / PL
% P + PL = Pt
% L + PL = Lt
% Kd = (Pt - PL)*(Lt - PL) / PL
% PL Kd = Pt Lt - PL Pt - PL Lt + PL^2
% PL^2 - (Pt + Lt + Kd) PL + Pt Lt = 0
% PL = ((Pt + Lt + Kd) +- sqrt((Pt + Lt + Kd)^2 - 4 Pt Lt)) / 2

% Pt and Lt are in M, result is in M.
complex_concentration = @(Pt,Lt,Kd) ((Pt + Lt + Kd) - sqrt((Pt + Lt + Kd).^2 - 4*Pt.*Lt)) / 2; % complex
ligand_concentration = @(Pt,Lt,Kd) Lt - complex_concentration(Pt,Lt,Kd); % free ligand
protein_concentration = @(Pt,Lt,Kd) Pt - complex_concentration(Pt,Lt,Kd); % free protein


% Define an objective function on observed fluorescence data.
lsq_fit_error = @(protein_fluorescence, ligand_fluorescence, buffer_fluorescence, complex_fluorescence, protein_quantity, Kd) sum((observed_protein_fluorescence - (protein_fluorescence*protein_quantity + buffer_fluorescence)).^2) + sum((observed_ligand_fluorescence - (ligand_quantities*ligand_fluorescence + buffer_fluorescence)).^2) + sum((observed_mixture_fluorescence - (complex_fluorescence.*well_volumes.*complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + ligand_fluorescence.*well_volumes.*ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + protein_fluorescence.*well_volumes.*protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + buffer_fluorescence)).^2);

objective = @(x) lsq_fit_error(exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), exp(x(6)));

% Create initial guess.
protein_quantity_estimate = protein_concentration_estimate * well_volume;
buffer_fluorescence_estimate = abs(mean(observed_buffer_fluorescence));
protein_fluorescence_estimate = abs((mean(observed_protein_fluorescence) - buffer_fluorescence_estimate) / protein_quantity_estimate);
ligand_fluorescence_estimate = abs(sum(observed_ligand_fluorescence - buffer_fluorescence_estimate) / sum(ligand_quantities));
complex_fluorescence_estimate = abs((observed_mixture_fluorescence(end) - buffer_fluorescence_estimate) / min(protein_quantity_estimate, ligand_quantities(end)));
Kd_estimate = 1e-6; 

% Solve for unknown parameters.
options = optimset('maxfunevals', 1e4, 'maxiter', 1e4);
x = log([protein_fluorescence_estimate ligand_fluorescence_estimate buffer_fluorescence_estimate complex_fluorescence_estimate protein_quantity_estimate Kd_estimate]);
x = fminsearch(objective, x, options);

% Extract solution.
protein_fluorescence = exp(x(1));
ligand_fluorescence = exp(x(2));
buffer_fluorescence = exp(x(3));
complex_fluorescence = exp(x(4));
protein_quantity = exp(x(5));
Kd = exp(x(6)); 

% Compute model-based fluorecence to compare with observed fluorescence.
computed_protein_fluorescence = protein_fluorescence * protein_quantity * ones(1,nwells) + buffer_fluorescence;
computed_ligand_fluorescence = ligand_fluorescence * ligand_quantities + buffer_fluorescence;
computed_mixture_fluorescence = complex_fluorescence * well_volumes .* complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + protein_fluorescence * well_volumes .* protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + ligand_fluorescence * well_volumes .* ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + buffer_fluorescence;
computed_buffer_fluorescence = buffer_fluorescence * ones(1,nwells);

protein_concentration = protein_quantity / well_volume;

return
