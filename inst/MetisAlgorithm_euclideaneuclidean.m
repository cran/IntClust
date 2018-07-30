cls = csvread('matlabdata_euclideaneuclidean.csv')
Partition_euclideaneuclidean=hmetis(cls,2,euclideaneuclidean)
csvwrite('Partition_euclideaneuclidean.csv', Partition_euclideaneuclidean)
