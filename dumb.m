clear; clc;

file = fopen("variables.txt");

line1 = fgetl(file)
line2 = fgetl(file)

fclose(file);