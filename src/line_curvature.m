function k=line_curvature(Vertices,Lines)

if(nargin<2)
    Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
end
Na=zeros(size(Vertices,1),1); Nb=zeros(size(Vertices,1),1);
Na(Lines(:,1))=Lines(:,2); Nb(Lines(:,2))=Lines(:,1);
checkNa=Na==0; checkNb=Nb==0;
Naa=Na; Nbb=Nb;
Naa(checkNa)=find(checkNa); Nbb(checkNb)=find(checkNb);
Na(checkNa)=Nbb(Nbb(checkNa)); Nb(checkNb)=Naa(Naa(checkNb));
Ta=-sqrt(sum((Vertices-Vertices(Na,:)).^2,2));
Tb=sqrt(sum((Vertices-Vertices(Nb,:)).^2,2)); 
Ta(checkNa)=-Ta(checkNa); Tb(checkNb)=-Tb(checkNb);
x = [Vertices(Na,1) Vertices(:,1) Vertices(Nb,1)];
y = [Vertices(Na,2) Vertices(:,2) Vertices(Nb,2)];
M = [ones(size(Tb)) -Ta Ta.^2 ones(size(Tb)) zeros(size(Tb)) zeros(size(Tb)) ones(size(Tb)) -Tb Tb.^2];
invM=inverse3(M);
a(:,1)=invM(:,1,1).*x(:,1)+invM(:,2,1).*x(:,2)+invM(:,3,1).*x(:,3);
a(:,2)=invM(:,1,2).*x(:,1)+invM(:,2,2).*x(:,2)+invM(:,3,2).*x(:,3);
a(:,3)=invM(:,1,3).*x(:,1)+invM(:,2,3).*x(:,2)+invM(:,3,3).*x(:,3);
b(:,1)=invM(:,1,1).*y(:,1)+invM(:,2,1).*y(:,2)+invM(:,3,1).*y(:,3);
b(:,2)=invM(:,1,2).*y(:,1)+invM(:,2,2).*y(:,2)+invM(:,3,2).*y(:,3);
b(:,3)=invM(:,1,3).*y(:,1)+invM(:,2,3).*y(:,2)+invM(:,3,3).*y(:,3);
k = 2*(a(:,2).*b(:,3)-a(:,3).*b(:,2)) ./ ((a(:,2).^2+b(:,2).^2).^(3/2));

function  Minv  = inverse3(M)
adjM(:,1,1)=  M(:,5).*M(:,9)-M(:,8).*M(:,6);
adjM(:,1,2)=  -(M(:,4).*M(:,9)-M(:,7).*M(:,6));
adjM(:,1,3)=  M(:,4).*M(:,8)-M(:,7).*M(:,5);
adjM(:,2,1)=  -(M(:,2).*M(:,9)-M(:,8).*M(:,3));
adjM(:,2,2)=  M(:,1).*M(:,9)-M(:,7).*M(:,3);
adjM(:,2,3)=  -(M(:,1).*M(:,8)-M(:,7).*M(:,2));
adjM(:,3,1)=  M(:,2).*M(:,6)-M(:,5).*M(:,3);
adjM(:,3,2)=  -(M(:,1).*M(:,6)-M(:,4).*M(:,3));
adjM(:,3,3)=  M(:,1).*M(:,5)-M(:,4).*M(:,2);
detM=M(:,1).*M(:,5).*M(:,9)-M(:,1).*M(:,8).*M(:,6)-M(:,4).*M(:,2).*M(:,9)+M(:,4).*M(:,8).*M(:,3)+M(:,7).*M(:,2).*M(:,6)-M(:,7).*M(:,5).*M(:,3);
Minv=bsxfun(@rdivide,adjM,detM);