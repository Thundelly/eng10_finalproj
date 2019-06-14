global speed

disp(speed);

file = fopen("variables.txt", "w");
fprintf(file, "%s", speed);
fclose(file);