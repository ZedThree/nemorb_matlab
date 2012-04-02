function local_variable = nemorb_genericornot(sim, info, variable, ind)

% For backwards compatibility - previous versions had some
% parameters just in /parameters/, newer ones have them in
% /parameters/generic/'. This routine checks the former first, if
% the variable is not there, then it assumes it is in the
% latter. Once the variable is found, it is assigned to
% local_variable.

if exist('ind')==0
ind=1;
end

for i=1:length(ind)
k=ind(i);

for ii=1:length(info(k).GroupHierarchy.Groups(sim(k).ipara).Datasets)
% Could just strip off '/parameters/' from dataset name, but this way is
% more general
groupname=info(k).GroupHierarchy.Groups(sim(k).ipara).Name;
datasetname=info(k).GroupHierarchy.Groups(sim(k).ipara).Datasets(ii).Name;
datasetname=datasetname(length(groupname)+2:end);
if strcmp(datasetname,variable) == 1
local_variable = sim(k).parameters.(variable);
end
end
if exist('local_variable') == 0
local_variable = sim(k).generic.(variable);
end

end