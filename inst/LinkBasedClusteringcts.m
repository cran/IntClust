function LinkBasedClusteringcts(file_number) 		
M = csvread(strcat('matlabdata_',file_number,'.csv'));
S=cts(M,0.8);					
csvwrite(strcat('S_',file_number,'.csv'), S)