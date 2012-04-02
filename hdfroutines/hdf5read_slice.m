function [data] = hdf5read_slice_new(FileName, GroupName, offset, count)

file = H5F.open(FileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
datatypeID = H5T.copy('H5T_NATIVE_DOUBLE');
dataset = H5D.open(file, GroupName);

dataspace = H5D.get_space(dataset);    
H5S.select_hyperslab(dataspace, 'H5S_SELECT_SET', offset, [], count, []);

rank = length(count);
memspace = H5S.create_simple(rank, count, []);   

data = H5D.read(dataset, datatypeID, memspace, dataspace, 'H5P_DEFAULT');

H5S.close(dataspace);
H5D.close(dataset);
H5S.close(memspace);
H5F.close(file);
