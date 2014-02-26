function [Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data5(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_ligand_fluorescence, observed_mixture_fluorescence)
% Fit the observed titration data with a 1:1 binding model.
%
% ARGUMENTS
%
% 

% Determine number of injections.
ninjections = length(injection_volumes);

% Determine ligand quantities.
ligand_quantities = cumsum([0 injection_volumes.*injection_concentrations]);

% Determine well volumes.
well_volumes = cumsum([initial_volume injection_volumes]);

% Ligand total concentrations.
ligand_total_concentrations = ligand_quantities ./ well_volumes;

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

protein_quantity = initial_protein_concentration_estimate * initial_volume;

% Define an objective function on observed fluorescence data.
error_fun = @(protein_fluorescence, ligand_fluorescence, complex_fluorescence, buffer_fluorescence, protein_quantity, Kd) sum((observed_ligand_fluorescence - buffer_fluorescence - ligand_quantities*ligand_fluorescence).^2) + sum((observed_mixture_fluorescence - buffer_fluorescence - complex_fluorescence.*well_volumes.*complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - ligand_fluorescence.*well_volumes.*ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - protein_fluorescence.*well_volumes.*protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd)).^2);

objective = @(x) error_fun(exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), exp(x(6)));

% Create initial guess.
protein_quantity_estimate = initial_protein_concentration_estimate * initial_volume;
protein_fluorescence_estimate = (observed_mixture_fluorescence(1) - observed_ligand_fluorescence(1)) / protein_quantity_estimate;
complex_fluorescence_estimate = (observed_mixture_fluorescence(end) - observed_ligand_fluorescence(end)) / min(protein_quantity_estimate, ligand_quantities(end));
buffer_fluorescence_estimate = observed_ligand_fluorescence(1);
ligand_fluorescence_estimate = (observed_ligand_fluorescence(end) - buffer_fluorescence_estimate) / ligand_quantities(end);
Kd_estimate = 1e-9; 

% Solve for unknown parameters.
x = log([protein_fluorescence_estimate, ligand_fluorescence_estimate, complex_fluorescence_estimate, buffer_fluorescence_estimate, protein_quantity_estimate, Kd_estimate]);
options = optimset('maxfunevals', 1e4, 'maxiter', 1e4);
x = fminsearch(objective, x, options);

% Extract solution.
protein_fluorescence = exp(x(1));
ligand_fluorescence = exp(x(2));
complex_fluorescence = exp(x(3));
buffer_fluorescence = exp(x(4));
protein_quantity = exp(x(5));
Kd = exp(x(6));

% Compute model-based fluorecence to compare with observed fluorescence.
computed_ligand_fluorescence = ligand_fluorescence * ligand_quantities + buffer_fluorescence;
computed_mixture_fluorescence = complex_fluorescence * well_volumes .* complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + protein_fluorescence * well_volumes .* protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + ligand_fluorescence * well_volumes .* ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + buffer_fluorescence;

protein_concentration = protein_quantity / well_volumes(1);

return
