function X=Ch_Int(A,B,h,t)
X(1:numel(A(:,1)),1)=0;
 for i=1:numel(t)-1
    X(:,i+1)=X(:,i)+(h)*(A*X(:,i)+B');
end