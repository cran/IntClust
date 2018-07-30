function CSPA(file_number)
cls = csvread(strcat('matlabdata_',file_number,'.csv'));
ClusterEnsemble=cspa(cls,file_number);
csvwrite(strcat('ClusterEnsemble_',file_number,'.csv'), ClusterEnsemble)