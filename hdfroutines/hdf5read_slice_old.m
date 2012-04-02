function res=hdf5read_slice(filename,datasetname,dim1,dim2)

filename

  dim3=dim1+dim2;

%  [dim1 dim2 dim3]
  data=hdf5read(filename,datasetname);
%size(data)
res=data(dim1(2)+1:dim3(2),dim1(3)+1:dim3(3),dim1(1)+1:dim3(1));
clear data
