clc
clear
close all

%% Given Data
K = 1000;             % Strike Price
StockPrice = 1145.45; % Current Stock Price
SIG = 62.51/100;      % Implied Volatility
r = 0.02441; 

Expiration = 255./365.25;  % Years Until Expiration
Time_Steps = 25.5/365.25;  % 1 Day Time Steps

ML = 2*StockPrice; % Maximum Stock Price
No_blocks = 10;    % Number of Elements
dS = ML/No_blocks; % Price Step

No_nodes = No_blocks + 1;    % Mesh Centered Grid
N_free   = No_nodes  - 2;    % Number of Free Nodes

S = (0:dS:ML)'; 
Tav = (0:Time_Steps:Expiration)'; 
n_time = length(Tav); 

C_1 = 0; 
C_n = ML - K*exp(-r*Tav); 
C = max(S-K, 0); 

C_total = zeros(No_nodes, n_time); 
C_total(:, 1) = C; 
C_total(end, :) = C_n; 
C = C(2:end-1); 

syms SP
dplus_T = (log(SP/K)+(r+SIG^2/2)*(Expiration))/(SIG*(Expiration)^(1/2));
dminus_T = (log(SP/K)+(r-SIG^2/2)*(Expiration))/(SIG*(Expiration)^(1/2));
Exact_Solution = SP*normcdf(dplus_T) - K*exp(-r*Expiration)*normcdf(dminus_T);

%% Initializations
A = zeros(N_free, N_free); 
B = zeros(N_free, 1); 
for t = 2:n_time
    %%  Assembly
    theta = 1/2; 
    for i = 1:N_free
        %% Options Price Constants (Theta Differencing)
        a1 = theta*S(i+1)/(2*dS)*(r-SIG^2*S(i+1)/dS);         % Constant for C(i-1, t+dt)
        a2 = 1/Time_Steps + theta*((SIG*S(i+1)/dS)^2+r);      % Constant for C(i, t+dt)
        a3 = -theta*S(i+1)/(2*dS)*(r+SIG^2*S(i+1)/dS);        % Constant for C(i+1, t+dt)
        b1 = (1-theta)*S(i+1)/(2*dS)*(-r+SIG^2*S(i+1)/dS);    % Constant for C(i-1, t)
        b2 = 1/Time_Steps - (1-theta)*((SIG*S(i+1)/dS)^2+r);  % Constant for C(i, t)
        b3 = (1-theta)*S(i+1)/(2*dS)*(r+SIG^2*S(i+1)/dS);     % Constant for C(i+1, t)

        %% Matrix Majholat
        if i == N_free
            A(i, i-1) = a1;
            A(i, i)   = a2;
        elseif i == 1
            A(i, i)   = a2;
            A(i, i+1) = a3;
        else
            A(i, i-1) = a1;
            A(i, i)   = a2;
            A(i, i+1) = a3;
        end

        %% Matrix Malomat
        if i == N_free
            B(i) = b1*C(i-1) + b2*C(i) + b3*C_n(t-1) - a3*C_n(t);
        elseif i == 1
            B(i) = (b1-a1)*C_1 + b2*C(i) + b3*C(i+1);
        else
            B(i) = b1*C(i-1) + b2*C(i) + b3*C(i+1);
        end
    end

    C = A\B;
    C_total(2:end-1, t) = C;
end
C_total

%% Rasm Sheklha
current_index = find(S == max(S(S <= StockPrice)));
Option_Price  = (C_total(current_index+1,end)-C_total(current_index,end))...
              /dS*(StockPrice-S(current_index)) + C_total(current_index,end);

%% Plotting
figure(1)
tao_days = Tav*365.25; % Change time scale to days for graphing
surf(tao_days, S, C_total)
xlabel('Time Until Expiration (Days)')
ylabel('Stock Price ($)')
zlabel('Option Price ($)')

figure(2)
fplot(SP, Exact_Solution, [0, ML], 'LineWidth', 2)
hold on
plot(S, C_total(:, end), 'r--', 'LineWidth', 2)
hold off
xlabel('Stock Price ($)')
ylabel('Option Price ($)')
legend('Black-Scholes Equation', 'Finite Difference Approximation')
ylim([0 inf])

%% Print Results
fprintf('Solution Black-Scholes Equation with Finite Difference Method is $%.4f\n', Option_Price)
fprintf('Exact Solution Black-Scholes Equation is $%.4f\n', subs(Exact_Solution, StockPrice))


%plotting error

% Define the time range
t = 0:10:1;

% Define the exact solution ce and numerical approximation cn
%ce = ExactSolution(t);
%cn = NumericalApproximation(t);

% Calculate the error
err = (SP - S) ./ S;


plot(err, t);
xlabel('Time');
ylabel('Error');
title('Error between Exact and Numerical Solutions');
