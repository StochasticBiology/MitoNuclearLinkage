R=[0:0.05:1]
n=size(R,2);
W=zeros(2,n);
for i=1:n
    i
    W(:,i)=NucMT(0.10,0.105,0.01,0.9,0.0,R(i))
end

%hold
plot(R,W(1,:))
plot(R,W(2,:))

for i=1:n
    i
    W(:,i)=NucMT(0.105,0.1,0.01,0.9,0.1,R(i))
end

plot(R,W(1,:))
plot(R,W(2,:))

for i=1:n
    i
    W(:,i)=NucMT(0.105,0.1,0.01,0.9,0.5,R(i))
end

plot(R,W(1,:))
plot(R,W(2,:))