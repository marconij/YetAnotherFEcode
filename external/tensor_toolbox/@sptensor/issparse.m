function tf = issparse(X)
tf = strcmp(X.type,'sparse');
return