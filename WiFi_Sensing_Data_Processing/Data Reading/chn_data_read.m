clear;
clc;

csi_trace = read_bf_file('./csi/isabelle/review/supply/home/grape_08041030.dat');
csi_trace_air = read_bf_file('./csi/isabelle/review/supply/home/air_08041030.dat');
L = length(csi_trace);
L_air = length(csi_trace_air);

length=min(L,L_air);

amplitudeA = zeros(length,30);
amplitudeB = zeros(length,30);
amplitudeA_air = zeros(length,30);
amplitudeB_air = zeros(length,30);

phaseA = zeros(length,30);
phaseB = zeros(length,30);
phaseA_air = zeros(length,30);
phaseB_air = zeros(length,30);

%Extract amplitude and phase from csi
csi_test_entry=csi_trace{1};
csi_test = get_scaled_csi(csi_test_entry);
for m = 1:length
    csi_entry = csi_trace{m};
    csi = get_scaled_csi(csi_entry);
    csi_entry_air = csi_trace_air{m};
    csi_air = get_scaled_csi(csi_entry_air);
    
    csi_amp = abs(squeeze(csi));       %amplitude
    amplitudeA(m,:) = csi_amp(1,1,:);
    amplitudeB(m,:) = csi_amp(1,2,:);
    csi_amp_air = abs(squeeze(csi_air));       %air amplitude
    amplitudeA_air(m,:) = csi_amp_air(1,1,:);
    amplitudeB_air(m,:) = csi_amp_air(1,2,:);
    %
    phaseA(m,:) = squeeze(csi(1,1,:));   %phase
    phaseB(m,:) = squeeze(csi(1,2,:));
    phaseA_air(m,:) = squeeze(csi_air(1,1,:));   %air phase
    phaseB_air(m,:) = squeeze(csi_air(1,2,:));
   
end
file_name = 'original_data.mat';
save(file_name);
