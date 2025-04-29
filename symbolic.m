clear all;

syms L R C;
S_a = sym("S_a", [3, 1], 'real');
S_b = sym("S_b", [3, 1], 'real');
S_c = sym("S_c", [3, 1], 'real');
v_dc = sym("v_dc", [3, 1], 'real');
i_s = sym("i_s", [3, 1], 'real');
v_s = sym("v_s", [3, 1], 'real');

M = [ 1  0 -1  0 -1  1  0  0  0;
    0  0  0  1  0 -1  0 -1  1;
    0 -1  1  0  0  0  1  0 -1];

% Equations (4)
v_a = v_dc .* (S_a-0.5);
v_b = v_dc .* (S_b-0.5);
v_c = v_dc .* (S_c-0.5);
v_x = reshape([v_a'; v_b'; v_c'], 9, 1);

% Equation (3)
v_o = simplify(M*v_x);

% Check equations (5, 6)
collect(v_o, v_dc)

% Equation (10)
i_x = M'*i_s;

% Equation (9)
i_dc = -sum(reshape(i_x, 3, 3)' .* [S_a, S_b, S_c], 2);

% Equation (7)
i_s_dt = 1/L*v_o - R/L*i_s - 1/L*v_s;

% Equation (8)
v_dc_dt = 1/C*i_dc;

% New variables
x = [i_s; v_dc];
x_dt = [i_s_dt; v_dc_dt];

% Check equation (12)
[A, b] = equationsToMatrix(x_dt == 0, x);
b = -b;
A, b