clc; clear all; close all;
num = [1 0.904 0];
% k = 1
den = [1 -0.904 0.095];
sys = tf(num,den,0.1)

rlocus(sys)