function [Kd, protein_concentration, computed_protein_fluorescence, computed_ligand_fluorescence, computed_mixture_fluorescence] = fit_titration_data_correction(initial_volume, initial_protein_concentration_estimate, injection_volumes, injection_concentrations, observed_protein_fluorescence, observed_ligand_fluorescence, observed_mixture_fluorescence)
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

% Define an objective function on observed fluorescence data.
%errorfun = @(protein_fluorescence, ligand_fluorescence, complex_fluorescence, protein_quantity, Kd) sum((observed_protein_fluorescence - protein_fluorescence*protein_quantity).^2) + sum((observed_ligand_fluorescence - ligand_fluorescence*ligand_quantities).^2) + sum((observed_mixture_fluorescence - complex_fluorescence.*well_volumes.*complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - ligand_fluorescence.*well_volumes.*ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - protein_fluorescence.*well_volumes.*protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd)).^2);

errorfun = @(protein_fluorescence, ligand_fluorescence, complex_fluorescence, protein_quantity, Kd, epsilonb) sum((observed_protein_fluorescence - protein_fluorescence*protein_quantity).^2) + sum((observed_ligand_fluorescence - ligand_fluorescence*ligand_quantities).^2) + sum((observed_mixture_fluorescence - complex_fluorescence.*well_volumes.*complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) - ligand_fluorescence*(1-exp(-epsilonb.*well_volumes.*ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd))) - protein_fluorescence.*well_volumes.*protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd)).^2);

%objective = @(x) errorfun(x(1)^2, x(2)^2, x(3)^2, x(4)^2, x(5)^2);
%objective = @(x) errorfun(exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)));
objective = @(x) errorfun(exp(x(1)), exp(x(2)), exp(x(3)), exp(x(4)), exp(x(5)), x(6));

% Create initial guess.
protein_quantity_estimate = initial_protein_concentration_estimate * initial_volume;
protein_fluorescence_estimate = mean(observed_protein_fluorescence) / protein_quantity_estimate;
ligand_fluorescence_estimate = abs(observed_ligand_fluorescence(end) / ligand_quantities(end));
complex_fluorescence_estimate = abs(observed_mixture_fluorescence(end) / max(protein_quantity_estimate, ligand_quantities(end)));
Kd_estimate = 1e-6; 
epsilonb = 1.0 * 10e-3; % epsilon (M/cm) * b

% Solve for unknown parameters.
%x = sqrt([protein_fluorescence_estimate ligand_fluorescence_estimate complex_fluorescence_estimate protein_quantity_estimate Kd_estimate]);
x = log([protein_fluorescence_estimate ligand_fluorescence_estimate complex_fluorescence_estimate protein_quantity_estimate Kd_estimate exp(epsilonb)]);
options = optimset('maxfunevals', 1e4, 'maxiter', 1e4);
for i = 1:5
  x = fminsearch(objective, x, options)
  disp(sprintf('objective = %f', objective(x)));
end


% Extract solution.
%protein_fluorescence = x(1)^2
%ligand_fluorescence = x(2)^2
%complex_fluorescence = x(3)^2
%protein_quantity = x(4)^2;
%Kd = x(5)^2;

protein_fluorescence = exp(x(1));
ligand_fluorescence = exp(x(2));
complex_fluorescence = exp(x(3));
protein_quantity = exp(x(4));
Kd = exp(x(5));
epsilonb = x(6)

% Compute model-based fluorecence to compare with observed fluorescence.
computed_protein_fluorescence = protein_fluorescence * protein_quantity * ones(1,ninjections+1);
computed_ligand_fluorescence = ligand_fluorescence * (1 - exp(-epsilonb * ligand_quantities .* well_volumes));
computed_mixture_fluorescence = complex_fluorescence * well_volumes .* complex_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + protein_fluorescence * well_volumes .* protein_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd) + ligand_fluorescence * well_volumes .* ligand_concentration(protein_quantity ./ well_volumes, ligand_quantities ./ well_volumes, Kd);

protein_concentration = protein_quantity / initial_volume;

return

