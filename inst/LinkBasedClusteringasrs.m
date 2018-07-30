function LinkBasedClusteringasrs(file_number) 		
M = csvread(strcat('matlabdata_',file_number,'.csv'));
S=asrs(M,0.8);					
csvwrite(strcat('S_',file_number,'.csv'), S)