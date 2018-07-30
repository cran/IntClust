function LinkBasedClusteringsrs(file_number) 		
M = csvread(strcat('matlabdata_',file_number,'.csv'));
S=srs(M,0.8,5);					
csvwrite(strcat('S_',file_number,'.csv'), S)