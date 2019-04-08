% Date: April 3, 2019
% Author: Stewart Pearson 
% Description: Matlab script for Assignment 9 which plots the scaling for
% the diffring and montecarlo ring and displays the respective serial
% fraction f

clear all;

% Data is formated in [#processes:time(s)]
walk_data =[32 399.219;	
            30 398.574;	
            28 406.539;
            24 431.269;	
            22 447.433;	
            21 454.279;
            18 502.090;
            16 528.647;	
            15 549.836;	
            12 642.422;	
            11 688.703;	
            10 746.157;	
            9 814.340;	
            8 908.048;	
            7 1017.082;	
            6 1153.343;	
            5 1359.299;	
            4 1641.477;	
            3 2050.384;
            2 3025.236;
            1 5860.671];

Ts = walk_data(size(walk_data,1),2);
S = Ts./walk_data(:,2);
f = (1-S(1)/32)/(S(1)-S(1)/32)
                
% Plot walkring scaling
figure;
plot(walk_data(:,1), S)
xlabel('Number of Processes')
ylabel('Speedup (S)')
title('Scaling analysis for walkring')

