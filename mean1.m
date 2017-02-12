function y = mean1(x,e,dim)
dim=1;
if nargin==1, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end

  y = sum(x);
else
  y = sum(x,dim);
end
t=size(x,dim)-e;
y=y/t;
