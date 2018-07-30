function MetisAlgorithm(Optimalk,file_number)
if (isstr(Optimalk))
  Optimalk = str2num(Optimalk);
end;
cls = csvread(strcat('matlabdata_',file_number,'.csv'));
Partition=hmetis(cls,Optimalk,file_number);
csvwrite(strcat('Partition_',file_number,'.csv'), Partition)