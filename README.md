# Linear systems with the smallest infinity norm

This repo contains Matlab codes for finding the solution of the linear system ``Ax = y`` with the smallest infinity norm. Due to splitting the computation to online and offline phases, our algorithm is suited for repeated computations of this system with the same matrix ``A`` but different right-hand side ``y``. See more detailed information in our paper.

## Algorithm

The problem above can be written as a linear optimization problem. The theory of linear programming states that it is equivalent to its dual problem. Our algorithm is based on the idea that when ``A`` is fixed, the dual problem has only a finite number of possible solutions ``u``. We divide the algorithm into two phases:
- <i>Offline phase</i>: We compute all possible dual solutions and save them into a matrix ``U``.
- <i>Online phase</i>: When the right-hand side ``y`` gets known, we select the optimal solution ``u`` from the pre-computed set ``U`` and based on the complementarity conditions, we recover the primal solution ``x``.

## Application to multi-phase converters

Multi-phase converters can be written as a linear system. The goal is to compute the input voltage ``x`` for the required output voltage ``y``. Since the infinity norm amounts to the dc-link, it is natural to minimize this quantity.

We create the matrix ``A`` as a fault-tolerant system with four phases.

```
A1 = 2/3*[1, cos(2/3*pi), cos(-2/3*pi); ...
    0, sin(2/3*pi), sin(-2/3*pi); ...
    1/2, 1/2, 1/2];
A2 = [1, 0, 0, -1;...
    0, 1, 0, -1;...
    0, 0, 1, -1];
    
A = A1*A2;
[m, n] = size(A);
```

Now we perform time discretization of the interval ``[0, 0.06]`` and specify the right-hand side vectors ``y``.

```
ts = 0:100e-6:0.06;
N = length(ts);

ys = zeros(m,N);
for k = 1:N
    wt = 2*pi*50*ts(k);
    ys(:,k) = 230*sqrt(2)*[cos(wt); sin(wt); -0.4*cos(wt)];
end
```

<img src="Figures/res1.png" width="500">

The offline phase precomputes the matrix ``U``. It has only 12 elements, which means that the dual problem can be solved in 12 dot products. Our basic Matlab implementation needs only 0.05ms for one solution on a laptop.

```
U = get_u(A);
```

The online phase computes the optimal input voltage ``x`` for each realization of the output voltage ``y``.

```
xs = zeros(n,N);
for k = 1:N
    xs(:,k) = min_effort(A, ys(:,k), U);
end
```

We compare our result with the standard method minimizing the l2 norm. We depict both results next to each other with the same y axis. Our approach leads to 18% smaller dc-link, which means significant energy savings.

<img src="Figures/res2.png" width="500"><img src="Figures/res3.png" width="500">



