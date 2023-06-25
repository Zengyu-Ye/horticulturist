dBm_rate = xlsread('F:\Summer-Internship\raw data\p1.xlsx');
power_rate  = 10.^(dBm_rate/10-3);
