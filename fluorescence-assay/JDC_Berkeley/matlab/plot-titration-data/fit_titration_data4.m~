function [Kd, protein_concentration, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data4(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_ligand_fluorescence, observed_mixture_fluorescence)
% Fit the observed titration data with a 1:1 binding model.
%
% ARGUMENTS
%
% 


% Determine number of injections.
ninjections = length(injection_volumes);

% Determine ligand quantities.
ligand_quantities = [0 cumsum(injection_volumes.*injection_concentrations)];

% Determine well volumes.
well_volumes = [initial_volume (initial_volume + cumsum(injection_volumes))];

total_ligand_concentrations = ligand_quantities ./ well_volumes;

well_area = 0.1586 / 100 / 100; % m^2
path_lengths = (well_volumes * 1000 / (100^3)) / well_area ; % m

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

errorfun = @(protein_fluorescence, ligand_fluorescence, complex_fluorescence, buffer_fluorescence, protein_quantity, Kd, epsilonb) sum((observed_ligand_fluorescence - buffer_fluorescence*well_volumes - ligand_fluorescence*(1-exp(-epsilonb.*path_lengths.*total_ligand_concentrations))).^2) + sum((observed_mixture_fluorescence - buffer_fluorescence*well_volumes - complex_fluorescence.*well_volumes.*complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - ligand_fluorescence*(1-exp(-epsilonb.*path_lengths.*ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd))) - protein_fluorescence.*well_volumes.*protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd)).^2);

objective = @(x) errorfun(exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), exp(x(6)), exp(x(7)));

% Create initial guess.
protein_quantity_estimate = initial_protein_concentration_estimate * initial_volume;
buffer_fluorescence_estimate = observed_ligand_fluorescence(1) / well_volumes(1);
protein_fluorescence_estimate = (observed_mixture_fluorescence(1) - observed_ligand_fluorescence(1)) / protein_quantity_estimate;
ligand_fluorescence_estimate = abs(observed_ligand_fluorescence(end) - buffer_fluorescence_estimate*well_volumes(end)) / ligand_quantities(end);
complex_fluorescence_estimate = (observed_mixture_fluorescence(end) - observed_ligand_fluorescence(end)) / min(protein_quantity_estimate, ligand_quantities(end));
Kd_estimate = 1e-6; 
epsilonb = 1.0 * 10e-3; % epsilon (M/cm) * b

% Solve for unknown parameters.
x = log([protein_fluorescence_estimate ligand_fluorescence_estimate complex_fluorescence_estimate buffer_fluorescence_estimate protein_quantity_estimate Kd_estimate epsilonb]);
options = optimset('maxfunevals', 1e4, 'maxiter', 1e4);
x = fminsearch(objective, x, options)

% Extract solution.
protein_fluorescence = exp(x(1));
ligand_fluorescence = exp(x(2));
complex_fluorescence = exp(x(3));
buffer_fluorescence = exp(x(4));
protein_quantity = exp(x(5));
Kd = exp(x(6));
epsilonb = exp(x(7));

% Compute model-based fluorecence to compare with observed fluorescence.
computed_ligand_fluorescence = ligand_fluorescence * (1 - exp(-epsilonb * total_ligand_concentrations .* path_lengths)) + buffer_fluorescence * well_volumes;
computed_mixture_fluorescence = complex_fluorescence * well_volumes .* complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + protein_fluorescence * well_volumes .* protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + ligand_fluorescence * (1 - exp(-epsilonb * ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) .* path_lengths)) + buffer_fluorescence * well_volumes;

protein_concentration = protein_quantity / initial_volume;

return

