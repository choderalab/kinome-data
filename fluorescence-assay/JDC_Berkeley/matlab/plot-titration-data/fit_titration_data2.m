function [Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data2(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence)
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
objective = @(protein_fluorescence, ligand_fluorescence, complex_fluorescence, ligand_quantities, Kd) sum((observed_protein_fluorescence - protein_fluorescence*protein_quantity).^2) + sum((observed_ligand_fluorescence - ligand_quantities*ligand_fluorescence).^2) + sum((observed_mixture_fluorescence - complex_fluorescence.*well_volumes.*complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - ligand_fluorescence.*well_volumes.*ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - protein_fluorescence.*well_volumes.*protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd)).^2);

% Create initial guess.
protein_quantity_estimate = initial_protein_concentration_estimate * initial_volume;
protein_fluorescence_estimate = mean(observed_protein_fluorescence) / protein_quantity_estimate;
ligand_fluorescence_estimate = observed_ligand_fluorescence(end) / ligand_quantities(end);
complex_fluorescence_estimate = observed_mixture_fluorescence(end) / max(protein_quantity_estimate, ligand_quantities(end));
Kd_estimate = 10e-9; 

% Solve for unknown parameters.
x = sqrt([protein_fluorescence_estimate ligand_fluorescence_estimate complex_fluorescence_estimate injection_concentrations(end) Kd_estimate]);
options = optimset('maxfunevals', 1e5, 'maxiter', 1e5);
for i = 1:5
  x = fminsearch(@(x) objective(x(1)^2, x(2)^2, x(3)^2, x(4)^2 / injection_concentrations(end) * ligand_quantities, x(5)^2), x, options);
x
end

% Extract solution.
protein_fluorescence = x(1)^2;
ligand_fluorescence = x(2)^2;
complex_fluorescence = x(3)^2;
ligand_quantities = x(4)^2 / injection_concentrations(end) * ligand_quantities;
Kd = x(5)^2;

% Compute model-based fluorecence to compare with observed fluorescence.
computed_protein_fluorescence = protein_fluorescence * protein_quantity * ones(1,ninjections+1);
computed_ligand_fluorescence = ligand_fluorescence * ligand_quantities;
computed_mixture_fluorescence = complex_fluorescence * well_volumes .* complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + protein_fluorescence * well_volumes .* protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + ligand_fluorescence * well_volumes .* ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd);

protein_concentration = x(4)^2;

return
