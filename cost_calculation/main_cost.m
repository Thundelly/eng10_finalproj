delta_v = app.VEditField.Value;
engine_ef = app.VeEditField.Value; % engine efficiency
mass_rf = app.MoEditField.Value; % mass of rocket + max fuel
mass_r = app.MfEditField.Value; % mass of rocket

mass_f = mass_r*(exp(1)^(delta_v/engine_ef) - 1);
mass_f = num2str(mass_f);

app.ResultEditField.Value = mass_f;