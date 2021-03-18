# Ax-Solver

**EEE 361**

**HW2**

**Report**

![](RackMultipart20210318-4-1qwfhwa_html_b72e489d131b72ec.png)

**Umur Gogebakan**

**21702705**

**CS**

**[Problem](#_shed7c405gw3) 3**

**[Algorithms for the Solution](#_fo8olvbpasf4) 4**

[Pseudoinverse (A+)](#_stv9pk4y07rt) 4

[Code](#_stehmp8sqfro) 4

[Comment](#_bzn2mvc8m74f) 5

[Conjugate Gradients (CG)](#_dlqaa4oxs701) 6

[Code](#_6cbn3ris1xik) 7

[Comment](#_zho8jx3k7yrj) 8

[Generalized Minimum Residual (GMRES)](#_6bfcudabpj2f) 9

[Code](#_tlyjnjuzlxuy) 11

[Comment](#_l4pzkhlpwg9r) 12

[Minimum Residual (MINRES)](#_xomo3cwz8xz7) 13

[Code](#_275bys5iyhpx) 14

[Comment](#_2n9p9b9xd7yg) 15

[Detailed Analysis of CG](#_x3b9mevz3kuu) 16

[Stopping Criterion](#_wc67uxdyg1tp) 16

[Error Metric](#_d8fe1ey05ab4) 16

[Comment](#_a9slcpcjyaju) 19

[Comment](#_r0zarb5tnb4h) 22

[Comment](#_opwa147cwzze) 25

**[Appendix](#_xvindrov989h) 26**

[Codes](#_4zdxa69rmph4) 26

[HW2\_main.m](#_oo7zs6fk1e1u) 26

[produce\_S.m](#_axfc2xjd89su) 30

[psuedo\_inverse.m](#_muk6cf90pt9x) 31

[CG.m](#_mjof169ux975) 32

[GMRES.m](#_md5v723gwvto) 33

[MINRES.m](#_nx693wm0l48r) 34

##


## Problem

We want to solve system of linear equations that has a major application area in all mathematics and engineering:

Here, many issues arise as the size, sparsity, shape etc. of A changes when we try to solve this system on a computer where we have only finite precision arithmetic.

One easy approach is simply doing a Gaussian elimination which corresponds to LU factorization. Yet, when A is large this method becomes infeasible. Moreover, if A has special properties such as being sparse, symmetric etc. these properties can be exploited for an efficient computation of the vector x.

Here, we used three methods where two of them (GMRES and CG) are iterative and the other one is not (Pseudoinverse). The other one (MINRES) is approached as a special case of GMRES when A is symmetric.

## Algorithms for the Solution

## Pseudoinverse ()

Let be the solution found for the equation,

Then can be represented as the minimum that minimized the residual sum of squares:

This yields: ; let us call as pseudoinverse of A where:

So that:

In order to find pseudoinverse of A in all cases we have to apply SVD, which yields:

### Code

1 **function** [A\_psuedo] = **psuedo\_inverse** (A)

2 [U,S,V] = svd(A);

3 A\_psuedo = V \* transpose(diag( **1**./ (S\*ones(size(S, **1** ), **1** )))) \* transpose(U);

4 **end**

##


![](RackMultipart20210318-4-1qwfhwa_html_a26b674bad880089.jpg)

**Figure 1:** ES,i values with respect to sigma for Pseudoinverse algorithm where size m ∈ {100, 500, 2500} and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

Figure 1 indicates that as the noise (σ) increases, error ES,i also increases. The trend also indicates that error is less for small tau values than for the case with tau being larger (smaller tau means matrix is more sparse and positive-definiteness). Especially for the matrix A being size of 2500 with large tau, error is significantly high. This is due to the singularity and condition of the matrix. This concludes that for relatively small matrices tau is not that effective on the result, yet for large matrices alteration becomes significant.

## Conjugate Gradients (CG)

The conjugate gradient method is using the original Krylov subspace iteration to reach the result x. This algorithm converges significantly fast and stays at the heart of computation when the matrix A is symmetric positive definite. The main different motivation behind the CG is the exploitation of the Hessenberg matrix and measuring _A-norm_ is very appropriate to measure the error after each iteration.

### ![](RackMultipart20210318-4-1qwfhwa_html_183401364826bfcc.png)
# 1

The code below is the implementation of this algorithm with added stopping criterion:

### Code

1 **function** [x\_n, i, error] =CG(A, b, maximum\_iteration\_number, tol)

2 x\_n = 0;

3 r\_n = b;

4 r\_n = b;

5 p\_n = r\_n;

6 error = zeros(1, maximum\_iteration\_number);

7 **for** i = 1:maximum\_iteration\_number

8 _%update_

9 r\_n\_prev = r\_n;

10 p\_n\_prev = p\_n;

11 x\_n\_prev = x\_n;

12 _%_

13 alpha\_n = (transpose(r\_n\_prev)\*r\_n\_prev)/(transpose(p\_n\_prev)\*A\*p\_n\_prev);

14 x\_n = x\_n\_prev + alpha\_n\*p\_n\_prev;

15 r\_n = r\_n\_prev - alpha\_n\*A\*p\_n\_prev;

16 error(i) = norm(b-A\*x\_n)^2;

17 **if** norm(r\_n) \&lt; tol

18 **return**

19 **end**

20 beta\_n = (transpose(r\_n)\*r\_n)/(transpose(r\_n\_prev)\*r\_n\_prev);

21 p\_n = r\_n + beta\_n\*p\_n\_prev;

22 **end**

23 **end**

![](RackMultipart20210318-4-1qwfhwa_html_233c3fafdeaa1e2d.jpg)

**Figure 2:** ES,i values with respect to sigma for CG algorithm where size m ∈ {100, 500, 2500} and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

Figure 2 indicates that error ES,i increases as noise increases. Similar to Figure 1 for the case where size is 2500 and tau is 0.1, error is exceptionally high and does not have a decreasing fashion. This is because of the fact that the matrix is not positive definite. Sparsity is also another issue. When m is isolated, smaller the tau smaller the error we observe. This is because of the sparsity of the matrix and especially for small matrices we see that sparsity and positive definiteness properties are highly exploited by the CG algorithm. And again similar to Figure 1, for small matrices tau seems to be less effective when it is compared to matrices with larger size. Finally, the error graph we observe here seems similar to the graph for Pseudoinverse, this algorithm computes the answer in a much faster way.

## Generalized Minimum Residual (GMRES)

The idea of GMRES is, at step n to approximate x\* by the vector xn ∈ Кn that minimizes the norm residual rn = b - Axn.

Here, we use Arnoldi iterations to construct a sequence of Krylov matrices Qn whose columns q1, q2, q3, … span the successive Krylov subspaces Кn.

So;

![](RackMultipart20210318-4-1qwfhwa_html_92c5a55f6f3fd6e4.png)
# 2

Here, finding the minimizing y can be done by using QR factorization for solving least squares:

Then it follows,

### ![](RackMultipart20210318-4-1qwfhwa_html_6657a8c0cb7d9f0.png)
# 3

The code below is the implementation of this algorithm with added stopping criterion:

Here, the built-in MATLAB function (qr) is used to find the QR factorization. One can exploit the Hessenberg structure by using Givens rotations to produce orthonormal bases. Yet, it is not implemented as it is not covered in class and not asked for.

### Code

1 **function** [x, n, error] =GMRES(A, b, maximum\_iteration\_number, tol)

2_%initial arbitrary vector q\_0_

3 q\_1 = b / norm(b);

4 n = length(A);

5 H\_tilda = zeros(n+1, n);

6 Q = zeros(length(q\_1), n+1);

7 e\_1 = zeros(n+1, 1);

8 e\_1(1,1) = 1;

9 x = 0;

10 Q(:,1) = q\_1;

11 error = zeros(1, maximum\_iteration\_number);

12 **for** n = 1:maximum\_iteration\_number

13 _%step n of arnoldi iteration_

14 _%v = A\*q\_n;_

15 v = A\*Q(:,n);

16 **for** j=1:n

17 _%h\_j\_n = q\_j&#39; \* v;_

18 H\_tilda(j,n) = Q(:,j)&#39;\* v;

19 _%v = v - h\_j\_n \* q\_j;_

20 v = v - H\_tilda(j,n) \* Q(:,j);

21 **end**

22 _%h\_n\_next\_n = norm(v);_

23 H\_tilda(n+1,n) = norm(v);

24 _%q\_n\_next = v / h\_n\_next\_n;_

25 Q(:,n+1) = v / H\_tilda(n+1,n);

26 _%end of step n of Arnoldi iteration_

27 _%find y to minimize ||r\_n|| using QR_

28 [Q\_r, R] = qr(H\_tilda(1:n+1, 1:n));

29 d = Q\_r&#39;\* (norm(b)\*e\_1(1:n+1));

30 y = R\d;

31 x = Q(:, 1:n)\*y;

32 error(n) = norm(A\*x - b);

33 **if** error(n) \&lt; tol

34 **return**

35 **end**

36 **end**

37 **end**

![](RackMultipart20210318-4-1qwfhwa_html_2b8f9acd1aeba426.jpg)

**Figure 3:** ES,i values with respect to sigma for GMRES algorithm where size m ∈ {100, 500, 2500} and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

According to Figure 3, error ES,i increases with the noise level. For the case where the size of matrix is 2500 and tau is 0.1, error is significantly high when it is compared to other cases. This is due to the less sparse structure and not being a positive definite property of the matrix. These factors result in the matrix being ill-conditioned and make it unstable close to being a singular matrix. However, unlike the CG algorithm we see that the increase in error for the largest matrix decreases after sigma = 0.01 and makes a breaking point in terms of its gradient at sigma = 0.01. Moreover, the starting error for the large matrix is less than CG. Therefore, even if the error for the large matrix with sigma = 0.01 is significantly high, this difference becomes insignificant as we go along with the sigma positive direction. This due to the fact that although matrix condition and stability of matrix is the reason for this error, GMRES does not require matrix to be positive definite and that is the reason why when sigma = 1, all algorithms gather at a single close value. Again, small tau valued matrices have less error compared to ones with larger tau, due to the sparsity and relatively well-conditioning. For small matrices, tau value is not that effective but for the large ones its effect is explicitly observable. This concludes that the plot is similar to Pseudoinverse case excluding the largest case with tau being largest, but GMRES being much faster in terms of computation time.

## Minimum Residual (MINRES)

MINRES can be taught of special case variation of GMRES as the names suggest, to be used when A is symmetric. Because of the symmetric structure of the Hessenberg matrix, the Hessenberg structure obtained is actually tridiagonal. So Arnoldi iteration can be replaced by Lanczos iteration that uses only qn-1 and qn for the computation of qn+1. That removes the nested inner second for loop and achieves better computational efficiency.

### ![](RackMultipart20210318-4-1qwfhwa_html_588ca372d92a8234.png)
# 4

The code below is the implementation of Lanczos based MINRES algorithm with added stopping criterion:

Here, the built-in MATLAB function (qr) is used to find the QR factorization. One can exploit the Hessenberg structure by using Givens rotations to produce orthonormal bases. Yet, it is not implemented as it is not covered in class and not asked for.

### Code

1 **function** [x, n, error] =MINRES(A, b, maximum\_iteration\_number, tol)

2 _%initial arbitrary vector q\_0_

3 q\_0 = 0;

4 q\_n = q\_0;

5 q\_1 = b / norm(b);

6 beta\_0 = 0;

7 beta\_n = beta\_0;

8 n = length(A);

9 Q = zeros(length(q\_1), n+1);

10 T\_tilda = zeros(n+1, n);

11 e\_1 = zeros(n+1, 1);

12 e\_1(1,1) = 1;

13 x = 0;

14 Q(:,1) = q\_1;

15 error = zeros(1, maximum\_iteration\_number);

16 **for** n = 1:maximum\_iteration\_number

17 beta\_n\_prev = beta\_n;

18 q\_n\_prev = q\_n;

19 _%Lancsoz step_

20 _%v = Aqn_

21 v = A\*Q(:,n);

22 alpha\_n = Q(:,n)&#39;\*v;

23 v = v - beta\_n\_prev\*q\_n\_prev-alpha\_n\*Q(:,n);

24 beta\_n = norm(v);

25 Q(:,n+1) = v/beta\_n;

26 T\_tilda(n,n) = alpha\_n;

27 T\_tilda(n+1,n) = beta\_n;

28 T\_tilda(n,n+1) = beta\_n;

29 _%end of step n of Lanczos iteration_

30 _%find y to minimize ||r\_n|| using QR_

31 [Q\_r, R] = qr(T\_tilda(1:n+1, 1:n));

32 d = Q\_r&#39;\* (norm(b)\*e\_1(1:n+1));

33 y = R\d;

34 x = Q(:, 1:n)\*y;

35 error(n) = norm(A\*x - b);

36 **if** error(n) \&lt; tol

37 **return**

38 **end**

39 q\_n = Q(:,n);

39 **end**

40 **end**

![](RackMultipart20210318-4-1qwfhwa_html_71622dca1bb5413e.jpg)

**Figure 4:** ES,i values with respect to sigma for MINRES algorithm where size m ∈ {100, 500, 2500} and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

Please refer to [Comment section](#_l4pzkhlpwg9r) Figure 3 as Figure 4 is indistinguishably similar.

## Detailed Analysis of CG

### Stopping Criterion

Norm of the residual is used for the stopping criterion which is:

Orthogonality of q vectors which span the Krylov subspaces can be chosen as another metric for the stopping criterion. Yet, since we do not have to store these vectors this criterion is preferable.

As the following Figures 5-13 suggests, continuing iteration after the convergence is meaningless and may start to become numerically unstable to the machine epsilon as the orthogonality () may deviate (same applies to the GMRES). The algorithm converges very fast if the matrix A is positive definite and well-conditioned. Yet, also the Figures 11-13 shows that convergence does not always have to take place due to the matrix A being not positive definite and ill-conditioned.

### Error Metric

![](RackMultipart20210318-4-1qwfhwa_html_4cd701e9c97131c5.jpg)

**Figure 5:** EO,i values with respect to n iterations for CG algorithm where size m = 100, σ = 0.0001 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

![](RackMultipart20210318-4-1qwfhwa_html_dc123a7e088fb03c.jpg)

**Figure 6:** EO,i values with respect to n iterations for CG algorithm where size m = 100, σ = 0.01 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

![](RackMultipart20210318-4-1qwfhwa_html_595d1699518b82a.jpg)

**Figure 7:** EO,i values with respect to n iterations for CG algorithm where size m = 100, σ = 1 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

As it can be seen from the Figures 5-7, noise level does not have a significant effect on convergence of CG whereas tau has. For tau = 0.01, the algorithm converges faster and reaches around machine epsilon before the 10th iteration. Yet, for tau = 0.1, the algorithm reaches machine epsilon around the 20th iteration. However, it is still an important conclusion that for both tau values and different sigma values the algorithm succeeds to minimize the residual and converge to the desired solution.

![](RackMultipart20210318-4-1qwfhwa_html_6ffcc2ce06856972.jpg)

**Figure 8:** EO,i values with respect to n iterations for CG algorithm where size m = 500, σ = 0.0001 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

![](RackMultipart20210318-4-1qwfhwa_html_e60dea26c2bbc14c.jpg)

**Figure 9:** EO,i values with respect to n iterations for CG algorithm where size m = 500, σ = 0.01 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

![](RackMultipart20210318-4-1qwfhwa_html_481f17cd71dc60cf.jpg)

**Figure 10:** EO,i values with respect to n iterations for CG algorithm where size m = 500, σ = 1 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

Figures 8-10, shows us that as the size of matrix m increases from 100 to 500, convergence occurs in later iterations. The red line showing the convergence of the case where tau = 0.01 is very similar to Fİgures 5-7 but moved a bit further in the positive x direction. So the convergence still occurs before the 10th iteration but later than the case where m was 100, being close to the 9th iteration. The blue line indicating the convergence where tau = 0.1, is more horizontal and has a smaller gradient compared to the Figures 5-7. This indicates that for m = 500, larger tau is more affected by the size and convergence slows down. Yet, it is still an important observation that both lines are decreasing and have a converging trend. Again, noise seems to be not affecting the results significantly.

![](RackMultipart20210318-4-1qwfhwa_html_61aad9c27a1c7132.jpg)

**Figure 11:** EO,i values with respect to n iterations for CG algorithm where size m = 2500, σ = 0.0001 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

![](RackMultipart20210318-4-1qwfhwa_html_3ed22365fd799a4d.jpg)

**Figure 12:** EO,i values with respect to n iterations for CG algorithm where size m = 2500, σ = 0.01 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

![](RackMultipart20210318-4-1qwfhwa_html_a4b19d94f0a7693f.jpg)

**Figure 13:** EO,i values with respect to n iterations for CG algorithm where size m = 2500, σ = 1 and τ ∈ {0.1, 0.01}; x and y axis are log scaled.

### Comment

As Figures 11-13 suggests, when matrix size m is 2500 the convergence line for tau = 0.01 is even moved to further along positive x direction. It exceeds the 10th iteration. And the most significant and important conclusion that can be drawn from Figures 11-13 is the case where tau = 0.1. As the blue line suggests in all three graphs, when the matrix size is 2500 and tau is 0.1, the CG algorithm does not converge due to the ill conditioning and not being positive definite property of the matrix. Noise level seems to not affect the result that much when tau = 0.01. Yet, it has some observable effects on the shape of the graph for the case when tau = 0.1, whereas the change does not have a meaningful interpretation as for all noise levels the algorithm does not converge when tau is large.

All these comments are adherent with the literature, stating that CG is a powerful and important method. However, it is not always the case that matrices possibly arising in practice have this kind of well-behaved spectrum as Figures 11-13 suggests.

## Appendix

## Codes

### HW2\_main.m

close all

_%generate matrix A_

arr\_A = cell(3,2);

m\_vals = [100, 500, 2500];

tau\_vals = [0.1, 0.01];

**for** m\_i = 1 : size(m\_vals, 2)

**for** tau\_j = 1 : size(tau\_vals, 2)

arr\_A{m\_i, tau\_j} = produce\_S(m\_vals(m\_i), tau\_vals(tau\_j));

**end**

**end**

sigma\_w\_values = [0.0001, 0.01, 1];

arr\_b = cell(size(sigma\_w\_values, 2),10);

arr\_x\_0 = cell(10,1);

_%%_

E\_S\_psuedo = zeros(6, 3);

E\_S\_CG = zeros(6, 3);

E\_S\_GMRES = zeros(6, 3);

E\_S\_MINRES = zeros(6, 3);

**for** A\_i = 1:6

E\_S\_sum\_psuedo = zeros(1, 3);

E\_S\_sum\_CG = zeros(1, 3);

E\_S\_sum\_GMRES = zeros(1, 3);

E\_S\_sum\_MINRES = zeros(1, 3);

A = arr\_A{idivide(int32(A\_i), int32(2),&#39;ceil&#39;),2-rem(A\_i,2)};

_%_

m = m\_vals(idivide(int32(A\_i), int32(2),&#39;ceil&#39;));

**for** i = 1:10

x\_0 = randn(m, 1);

b\_0 = A\*x\_0;

**for** sigma\_w\_j = 1:size(sigma\_w\_values, 2)

sigma\_w = sigma\_w\_values(sigma\_w\_j);

w = sigma\_w \*randn(m, 1);

arr\_b{sigma\_w\_j, i} = b\_0 + w;

**end**

arr\_x\_0{i} = x\_0;

**end**

_%_

**for** i = 1:3

**for** j = 1:10

x\_0\_j = arr\_x\_0{j};

b\_i\_j = arr\_b{i, j};

x\_hat\_i\_j\_psuedo = psuedo\_inverse(A)\*b\_i\_j;

x\_hat\_i\_j\_CG = CG(A, b\_i\_j, 25, eps);

x\_hat\_i\_j\_GMRES = GMRES(A, b\_i\_j, 25, eps);

x\_hat\_i\_j\_MINRES = MINRES(A, b\_i\_j, 25, eps);

_%error summation_

E\_S\_sum\_psuedo(i) = E\_S\_sum\_psuedo(i) + norm(-x\_hat\_i\_j\_psuedo + x\_0\_j)^2;

E\_S\_sum\_CG(i) = E\_S\_sum\_CG(i) + norm(-x\_hat\_i\_j\_CG + x\_0\_j)^2;

E\_S\_sum\_GMRES(i) = E\_S\_sum\_GMRES(i) + norm(-x\_hat\_i\_j\_GMRES + x\_0\_j)^2;

E\_S\_sum\_MINRES(i) = E\_S\_sum\_MINRES(i) + norm(-x\_hat\_i\_j\_MINRES + x\_0\_j)^2;

**end**

**end**

E\_S\_psuedo(A\_i,:) = sqrt( 0.1\* E\_S\_sum\_psuedo);

E\_S\_CG(A\_i,:) = sqrt( 0.1\* E\_S\_sum\_CG);

E\_S\_GMRES(A\_i,:) = sqrt( 0.1\* E\_S\_sum\_GMRES);

E\_S\_MINRES(A\_i,:) = sqrt( 0.1\* E\_S\_sum\_MINRES);

**end**

_%plot E\_S\_psuedo_

figure

title(&#39;E\_{S,i} with respect to \sigma\_w for Pseudoinverse&#39;);

ylabel(&#39;E\_{S,i}&#39;);

xlabel(&#39;\sigma\_w&#39;);

set(gca, &#39;YScale&#39;, &#39;log&#39;)

set(gca, &#39;XScale&#39;, &#39;log&#39;)

hold on

**for** i = 1:6

plot(sigma\_w\_values, E\_S\_psuedo(i,:));

**end**

legend(&#39;m=100, \tau=0.1&#39;,&#39;m=100, \tau=0.01&#39;, &#39;m=500, \tau=0.1&#39;, &#39;m=500, \tau=0.01&#39;, &#39;m=2500, \tau=0.1&#39;, &#39;m=2500, \tau=0.01&#39;)

hold off

_%plot E\_S\_CG_

figure

title(&#39;E\_{S,i} with respect to \sigma\_w for CG&#39;);

ylabel(&#39;E\_{S,i}&#39;);

xlabel(&#39;\sigma\_w&#39;);

set(gca, &#39;YScale&#39;, &#39;log&#39;)

set(gca, &#39;XScale&#39;, &#39;log&#39;)

hold on

**for** i = 1:6

plot(sigma\_w\_values, E\_S\_CG(i,:));

**end**

legend(&#39;m=100, \tau=0.1&#39;,&#39;m=100, \tau=0.01&#39;, &#39;m=500, \tau=0.1&#39;, &#39;m=500, \tau=0.01&#39;, &#39;m=2500, \tau=0.1&#39;, &#39;m=2500, \tau=0.01&#39;)

hold off

_%plot E\_S\_GMRES_

figure

title(&#39;E\_{S,i} with respect to \sigma\_w for GMRES&#39;);

ylabel(&#39;E\_{S,i}&#39;);

xlabel(&#39;\sigma\_w&#39;);

set(gca, &#39;YScale&#39;, &#39;log&#39;)

set(gca, &#39;XScale&#39;, &#39;log&#39;)

hold on

**for** i = 1:6

plot(sigma\_w\_values, E\_S\_GMRES(i,:));

**end**

legend(&#39;m=100, \tau=0.1&#39;,&#39;m=100, \tau=0.01&#39;, &#39;m=500, \tau=0.1&#39;, &#39;m=500, \tau=0.01&#39;, &#39;m=2500, \tau=0.1&#39;, &#39;m=2500, \tau=0.01&#39;)

hold off

_%plot E\_S\_MINRES_

figure

title(&#39;E\_{S,i} with respect to \sigma\_w for MINRES&#39;);

ylabel(&#39;E\_{S,i}&#39;);

xlabel(&#39;\sigma\_w&#39;);

set(gca, &#39;YScale&#39;, &#39;log&#39;)

set(gca, &#39;XScale&#39;, &#39;log&#39;)

hold on

**for** i = 1:6

plot(sigma\_w\_values, E\_S\_MINRES(i,:));

**end**

legend(&#39;m=100, \tau=0.1&#39;,&#39;m=100, \tau=0.01&#39;, &#39;m=500, \tau=0.1&#39;, &#39;m=500, \tau=0.01&#39;, &#39;m=2500, \tau=0.1&#39;, &#39;m=2500, \tau=0.01&#39;)

hold off

_%%_

arr\_E\_O = cell(6,1);

max\_iter\_no = 20;

**for** A\_i = 1:6

A = arr\_A{idivide(int32(A\_i), int32(2),&#39;ceil&#39;),2-rem(A\_i,2)};

_%_

m = m\_vals(idivide(int32(A\_i), int32(2),&#39;ceil&#39;));

**for** i = 1:10

x\_0 = randn(m, 1);

b\_0 = A\*x\_0;

**for** sigma\_w\_j = 1:size(sigma\_w\_values, 2)

sigma\_w = sigma\_w\_values(sigma\_w\_j);

w = sigma\_w \*randn(m, 1);

arr\_b{sigma\_w\_j, i} = b\_0 + w;

**end**

arr\_x\_0{i} = x\_0;

**end**

_%_

E\_O = zeros(3, max\_iter\_no);

**for** i = 1:3

**for** j = 1:10

b\_i\_j = arr\_b{i, j};

[x, iter, error] = CG(A, b\_i\_j, max\_iter\_no, eps);

E\_O(i,:) = E\_O(i,:) + error;

**end**

E\_O(i,:) = sqrt(0.1\* E\_O(i,:));

**end**

arr\_E\_O{A\_i} = E\_O;

**end**

**for** m\_i=1:3

**for** sigma\_i=1:3

figure

sigma = sigma\_w\_values(sigma\_i);

m = m\_vals(m\_i);

title([&#39;E\_{O,i} with respect to n for CG where m=&#39;,num2str(m),&#39; and \sigma=&#39;, num2str(sigma)]);

ylabel(&#39;E\_{O,i}&#39;);

xlabel(&#39;n&#39;);

set(gca, &#39;YScale&#39;, &#39;log&#39;)

hold on

**for** tau\_i=0:1

E\_O = arr\_E\_O{2\*m\_i-1+tau\_i};

plot((1:max\_iter\_no), E\_O(sigma\_i,:));

**end**

legend(&#39;\tau=0.1&#39;,&#39;\tau=0.01&#39;)

hold off

**end**

**end**

###


### produce\_S.m

**function** [S] =produce\_S(m,tau)

_%put 1 at each diagonal position_

d = ones(m,1); _% The diagonal values_

_%random number from a unifrom distrubtion on [-1, 1]_

a = -1;

b = 1;

t = triu(bsxfun(@min,d,d.&#39;).\*(a + (b-a).\*rand(m)),1); _% The upper trianglar random values_

S = diag(d)+t+t.&#39;; % Put them together in a symmetric matrix

_%_

S(abs(S) \&gt; tau) = 0;

S(find(eye(m))) = 1;

**end**

###


### psuedo\_inverse.m

**function** [A\_psuedo] =psuedo\_inverse(A)

[U,S,V] = svd(A);

A\_psuedo = V \* transpose(diag(1./ (S\*ones(size(S, 1),1)))) \* transpose(U);

**end**

###


### CG.m

**function** [x\_n, i, error] =CG(A, b, maximum\_iteration\_number, tol)

_% x\_0 = 0;_

_% r\_0 = b;_

_% p\_0 = r\_0;_

_%_

_% x\_n = x\_0;_

_% r\_n = r\_0;_

_% p\_n = p\_0;_

x\_n = 0;

r\_n = b;

p\_n = r\_n;

error = zeros(1, maximum\_iteration\_number);

**for** i = 1:maximum\_iteration\_number

_%update_

r\_n\_prev = r\_n;

p\_n\_prev = p\_n;

x\_n\_prev = x\_n;

_%_

alpha\_n = (transpose(r\_n\_prev)\*r\_n\_prev)/(transpose(p\_n\_prev)\*A\*p\_n\_prev);

x\_n = x\_n\_prev + alpha\_n\*p\_n\_prev;

r\_n = r\_n\_prev - alpha\_n\*A\*p\_n\_prev;

error(i) = norm(b-A\*x\_n)^2;

**if** norm(r\_n) \&lt; tol

**return**

**end**

beta\_n = (transpose(r\_n)\*r\_n)/(transpose(r\_n\_prev)\*r\_n\_prev);

p\_n = r\_n + beta\_n\*p\_n\_prev;

**end**

**end**

###


### GMRES.m

**function** [x, n, error] =GMRES(A, b, maximum\_iteration\_number, tol)

_%initial arbitrary vector q\_0_

q\_1 = b / norm(b);

n = length(A);

H\_tilda = zeros(n+1, n);

Q = zeros(length(q\_1), n+1);

e\_1 = zeros(n+1, 1);

e\_1(1,1) = 1;

x = 0;

Q(:,1) = q\_1;

error = zeros(1, maximum\_iteration\_number);

**for** n = 1:maximum\_iteration\_number

_%step n of arnoldi iteration_

_%v = A\*q\_n;_

v = A\*Q(:,n);

**for** j=1:n

_%h\_j\_n = q\_j&#39; \* v;_

H\_tilda(j,n) = Q(:,j)&#39;\* v;

_%v = v - h\_j\_n \* q\_j;_

v = v - H\_tilda(j,n) \* Q(:,j);

**end**

_%h\_n\_next\_n = norm(v);_

H\_tilda(n+1,n) = norm(v);

_%q\_n\_next = v / h\_n\_next\_n;_

Q(:,n+1) = v / H\_tilda(n+1,n);

_%end of step n of Arnoldi iteration_

_%find y to minimize ||r\_n|| using QR_

[Q\_r, R] = qr(H\_tilda(1:n+1, 1:n));

d = Q\_r&#39;\* (norm(b)\*e\_1(1:n+1));

y = R\d;

x = Q(:, 1:n)\*y;

error(n) = norm(A\*x - b);

**if** error(n) \&lt; tol

**return**

**end**

**end**

**end**

### MINRES.m

**function** [x, n, error] =MINRES(A, b, maximum\_iteration\_number, tol)

_%initial arbitrary vector q\_0_

q\_0 = 0;

q\_n = q\_0;

q\_1 = b / norm(b);

beta\_0 = 0;

beta\_n = beta\_0;

n = length(A);

Q = zeros(length(q\_1), n+1);

T\_tilda = zeros(n+1, n);

e\_1 = zeros(n+1, 1);

e\_1(1,1) = 1;

x = 0;

Q(:,1) = q\_1;

error = zeros(1, maximum\_iteration\_number);

**for** n = 1:maximum\_iteration\_number

beta\_n\_prev = beta\_n;

q\_n\_prev = q\_n;

_%Lanczos step_

_%v = Aqn_

v = A\*Q(:,n);

alpha\_n = Q(:,n)&#39;\*v;

v = v - beta\_n\_prev\*q\_n\_prev-alpha\_n\*Q(:,n);

beta\_n = norm(v);

Q(:,n+1) = v/beta\_n;

T\_tilda(n,n) = alpha\_n;

T\_tilda(n+1,n) = beta\_n;

T\_tilda(n,n+1) = beta\_n;

_%end of step n of Lanczos iteration_

_%find y to minimize ||r\_n|| using QR_

[Q\_r, R] = qr(T\_tilda(1:n+1, 1:n));

d = Q\_r&#39;\* (norm(b)\*e\_1(1:n+1));

y = R\d;

x = Q(:, 1:n)\*y;

error(n) = norm(A\*x - b);

**if** error(n) \&lt; tol

**return**

**end**

q\_n = Q(:,n);

**end**

**end**

[1](#sdfootnote1anc) Trefethen, L. N., &amp; Bau, D. I. (2000). The Conjugate Gradient Iteration. In _Numerical Linear algebra_. Philadelphia: SIAM Society for Industrial and Applied Mathematics. doi:https://mseas.mit.edu/group/References/Books/Trefethen\_Bau%20-%20Numerical%20Linear%20Algebra.pdf

[2](#sdfootnote2anc)Ibid, GMRES.

[3](#sdfootnote3anc)Ibid, Arnoldi Iteration.

[4](#sdfootnote4anc)Ibid, Lanczos Iteration.
