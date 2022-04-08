function Q = qinpfeat(inp)

load m; load A;
load Efaces;

for i = 1 : size(Efaces,2)
    temp = Efaces'*A(:,i); 
    Dfeat(:,i) = temp;  
end

[irow icol] = size(inp);

InImage = reshape(inp',irow*icol,1);

Difference = double(InImage)-m; 

Pfeat = Efaces'*Difference; 

Q = Pfeat;

return;
