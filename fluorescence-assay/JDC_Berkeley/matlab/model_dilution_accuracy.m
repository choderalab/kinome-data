% Model accuracy of 1:2 dilution series.

desired_concentrations = 1.0 * 2.^(-(0:11));
nconcentrations = length(desired_concentrations);
single_transfer_cv = 0.017; % single-transfer coefficient of variation
dilution_factor = 0.5;
desired_volume = 100;
nreplicates = 1000;

realized_concentrations = zeros(nreplicates,nconcentrations);
realized_volumes = zeros(nreplicates,nconcentrations);
for replicate = 1:nreplicates
  realized_volumes(replicate,1) = desired_volume * (1 + single_transfer_cv * randn());
  realized_concentrations(replicate,1) = desired_concentrations(1);
  for i = 2:nconcentrations
    volume_ligand_transferred = desired_volume * dilution_factor * (1 + single_transfer_cv * randn());
    volume_buffer_transferred = desired_volume * (1-dilution_factor) * (1 + single_transfer_cv * randn());
    total_volume = volume_ligand_transferred + volume_buffer_transferred;
    realized_volumes(replicate,i) = total_volume;
    realized_concentrations(replicate,i) = volume_ligand_transferred / total_volume * realized_concentrations(replicate,i-1);
  end
end

concentration_cv = std(realized_concentrations) ./ desired_concentrations;
volume_cv = std(realized_volumes) ./ desired_volume;

clf;

subplot(2,1,1);
plot(1:12, concentration_cv, 'k.-');
hold on;
plot([0 12], single_transfer_cv * [1 1], 'r-');

subplot(2,1,2);
plot(1:12, volume_cv, 'k.-');
