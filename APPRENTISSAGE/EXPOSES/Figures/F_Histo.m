function [HH, BB] = F_Histo(H, B);

B = repmat(B', 1, 2);
H = [H' zeros(length(H), 1)];
H = reshape(H', 2*length(B), 1);
B = reshape(B', 2*length(B), 1);
[HH BB] = stairs(H, B);
