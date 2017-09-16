
clear all
% simulate x for given p
p = .99;
sigma_square = 1;
n = 1000;
X = zeros (n,1);
X(1) = round(rand(1));
for i = 2:n
    sample = rand(1);
    X(i) = X(i-1)*(sample < p ) + mod(X(i-1)+1, 2)*(sample > p);
end

% Simulate y conditioned on x for given ?^2.  
e = normrnd(0, sigma_square, n,1);
Y = X + e;
% dynamic programming method to estimate MAP z 
x = zeros(n,1);
S = zeros(n,2);
S(1,1) = -Y(1)^2./(2*sigma_square);
S(1,2) = -(Y(1)-1)^2./(2*sigma_square);

for k = 2:n
    for i = 0:1
        h(1) = S(k-1,1) -((Y(k)-i)^2)/(2*sigma_square) + log(p*(i==0) + (1-p)*(i~=0));
        h(2) = S(k-1,2) -((Y(k)-i)^2)/(2*sigma_square) + log(p*(i==1) + (1-p)*(i~=1));
        if  h(1) > h(2)
            x(k,i+1) = 0;
            S(k,i+1) = h(1);
        else
            x(k,i+1) = 1;
            S(k,i+1) = h(2);
        end
    end
end

z = zeros(n,1);
if S(n,1) < S(n,2)
    z(n) = 1;
else
    z(n) = 0;
end

for k = n-1:-1:1
    z(k) = x(k+1,z(k+1)+1);
end

% Plot x, y and the estimate z
plot(Y,'g')
hold on
plot(X,'LineWidth',3)
plot(z,'r','LineWidth',4)
hold off
axis([0 1000 -5 5 ])

h_legend = legend('Observed','True','Reconstructed');
set(h_legend,'FontSize',9);
title('Figure 4(a): Plot x, y & Estimate z; consider p = 0.99 & \sigma = 1','Fontsize',11)
xlabel('Time')
ylabel('value')


