function chease=load_chease_eq(filename)

info=hdf5info(filename);

clear grid0 list1 list2
l0=length(info.GroupHierarchy.Groups(1).Groups(1).Datasets);
l1=length(info.GroupHierarchy.Groups(1).Groups(2).Datasets);
l2=length(info.GroupHierarchy.Groups(1).Groups(3).Datasets);

chease.filename=filename;
for i=1:l1
 name1=info.GroupHierarchy.Groups(1).Groups(2).Datasets(i).Name;
data=hdf5read(filename,name1);
name1 = name1(2:end);

% Remove path from names
  is = strfind(name1, '/');
  name1 = name1(is(end)+1:end);
chease=setfield(chease,name1,data);
list1{i}=name1;
end
chease=setfield(chease,'list1',list1);

for i=1:l2
  name1=info.GroupHierarchy.Groups(1).Groups(3).Datasets(i).Name;
data=hdf5read(filename,name1);
name1 = name1(2:end);

% Remove path from names
  is = strfind(name1, '/');
  name1 = name1(is(end)+1:end);
chease=setfield(chease,name1,data);
list2{i}=name1;
end
chease=setfield(chease,'list2',list2);


for i=1:l0
  name1=info.GroupHierarchy.Groups(1).Groups(1).Datasets(i).Name;
data=hdf5read(filename,name1);
name1 = name1(2:end);

% Remove path from names
  is = strfind(name1, '/');
  name1 = name1(is(end)+1:end);
chease=setfield(chease,name1,data);
grid0{i}=name1;
end

chease=setfield(chease,'grid0',grid0);

chease.psinorm=chease.PSI/max(chease.PSI);

