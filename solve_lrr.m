function [Z,E] = solve_lrr(X,A,lambda,reg)
% This routine solves the following nuclear-norm optimization problem,
% min |Z|_*+lambda*|E|_L
% s.t., X = AZ+E
% inputs:
%        X -- D*N data matrix, D is the data dimension, and N is the number
%             of data vectors.
%        A -- D*M matrix of a dictionary, M is the size of the dictionary
%        lambda -- parameter
%        reg -- the norm chosen for characterizing E, 
%            -- reg=0 (default),                  use the l21-norm 
%            -- reg=1 (or ther values except 0),  use the l1-norm

if nargin<4 || isempty(reg)
    reg = 0;
end

Q = orth(A');
B = A*Q;

if reg==0
    [Z,E] = alm_lrr_l21(X,B,lambda);
else
    [Z,E] = alm_lrr_l1(X,B,lambda);
end
Z = Q*Z;