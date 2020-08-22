R=[0:0.05:0.1,0.1:0.1:1]
n=size(R,2);
W=zeros(2,n);
for i=1:n
    i
    W(:,i)=NucMT(0.1,0.1,0.01,0.9,R(i),0)
end
hold;
box;

plot(R,W(1,:))
plot(R,W(2,:))

for i=1:n
    i
    W(:,i)=NucMT(0.05,0.05,0.01,0.9,R(i),0)
end
plot(R,W(1,:))
plot(R,W(2,:))

for i=1:n
    i
    W(:,i)=NucMT(0.01,0.01,0.01,0.9,R(i),0)
end
plot(R,W(1,:))
plot(R,W(2,:))
