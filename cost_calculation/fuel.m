%% VARIABLES
file = fopen("variables.txt");

speed = fgetl(file);
speed = str2double(speed);

g = 9.80665;

si_b = 279.3;
si_stg1 = 337.8;
si_stg2 = 450.5;

mass_b = 42630;
mass_stg1 = 284089;
mass_stg2 = 20830;

price_b = 42630*3.66;
price_stg1 = 284089*0.18;
price_stg2 = 20830*0.26;
%% APPLY VARIABLES
app.TotalchangeinvelocityEditField.Value = speed;
app.GravitationalAccelerationEditField.Value = g;

app.SpecificImpluseEditField.Value = si_b;
app.MassofthepropellantEditField.Value = mass_b;
app.PriceofthepropellantEditField.Value = price_b;

app.SpecificImpluseEditField_2.Value = si_stg1;
app.MassofthepropellantEditField_2.Value = mass_stg1;
app.PriceofthepropellantEditField_2.Value = price_stg1;

app.SpecificImpluseEditField_3.Value = si_stg2;
app.MassofthepropellantEditField_3.Value = mass_stg2;
app.PriceofthepropellantEditField_3.Value = price_stg2;