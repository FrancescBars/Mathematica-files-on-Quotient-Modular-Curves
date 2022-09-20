For few concrete non-hyperelliptic quotient modular curve $C:=X_0(N)/W_N$ we add a file with Mathematica code where we compute the Petri model equation for $C$ and check if $C$ has any bielliptic involution or not. The usual notation for the file corresponding to $C$ is NW_N+some details.nb
The mathematica code was running in Mathematica 11. We also attach a zip file containing all the Mathematica files in this folder.

In the file 112w16(nobielliptic).nb or .pdf (corresponding to the quotient curve $X_0(112)/\langle w_{16}\rangle$) we provide a little more explanation 
on the Mathematica code used to compute the Petri model and to decide if the quotient curve is bielliptic or not. For the other files corresponding to another quotient
modular curves, the ideas are similar only with an ad-hoc modifications in the Mathematica code.

Take a file .nb in the folder corresponding to a quotient modular curve $X_0(N)/W_N$, first we obtain its Petri model by collecting the new modular forms
that appear in the $\mathbb{Q}$-decomposition of the Jacobian of the modular curve (which is computed by Magma in another folder of this github) and lifting them to level $N$.
Secondly we apply the criteria given in Journal of Algebra paper ``Bielliptic modular curves $X_0^*(N)$" of F.Bars and J. González (Prop. 2.6) to decide if the quotient modular curve is bielliptic or not.

* First example.

Next, we explain how to decide if $X_0(90)/\langle w_5\rangle$ is bielliptic or not, explaining details of the mathematica file 90w5.nb as another example.

First we take the new modular forms that appears in the $\mathbb{Q}$-decomposition of $Jac(X_0(90)/w_5)$ and 
lift them to modular forms of level 90 by the use of the operators $B_d$ (see for example Notation and Lemma 2.1 in paper Bars-Kamel-Schweizer ``Bielliptic quotient modular curves").
In this way we obtain a $\mathbb{Q}$-basis
of the differentials for the quotient modular curve (with the usual isomorphism between weight two modular forms and differentials).

```

f1 = q - q^2 + q^3 + q^4 - q^5 - q^6 - 4*q^7 - q^8 + q^9 + q^10 + 
  q^12 + 2*q^13 + 4*q^14 - q^15 + q^16 + 6*q^17 - q^18 - 4*q^19 - 
  q^20 - 4*q^21 - q^24 + q^25 - 2*q^26 + q^27 - 4*q^28 - 6*q^29 + 
  q^30 + 8*q^31 - q^32 - 6*q^34 + 4*q^35 + q^36 + 2*q^37 + 4*q^38 + 
  2*q^39 + q^40 - 6*q^41 + 4*q^42 - 4*q^43 - q^45 + q^48 + 
  9*q^49; 
 
f2 =  q + q^2 - q^4 - q^5 - 3*q^8 - q^10 + 4*q^11 - 2*q^13 - q^16 - 
  2*q^17 + 4*q^19 + q^20 + 4*q^22 + q^25 - 2*q^26 + 2*q^29 + 5*q^32 - 
  2*q^34 - 10*q^37 + 4*q^38 + 3*q^40 - 10*q^41 + 4*q^43 - 4*q^44 - 
  8*q^47 - 7*q^49; 
  
  f3 = q + q^2 + q^4 - q^5 + 2*q^7 + q^8 - q^10 - 6*q^11 - 4*q^13 + 2*q^14 +
   q^16 + 6*q^17 - 4*q^19 - q^20 - 6*q^22 + q^25 - 4*q^26 + 2*q^28 + 
  6*q^29 - 4*q^31 + q^32 + 6*q^34 - 2*q^35 + 8*q^37 - 4*q^38 - q^40 + 
  8*q^43 - 6*q^44 - 3*q^49;
  
  
g11 = (f1 + (3*f1 /. q -> q^3)) // Expand
g12 = (f1 - (3*f1 /. q -> q^3)) // Expand
g21 = (f2 + (2*f2 /. q -> q^2)) // Expand
g22 = (f2 - (2*f2 /. q -> q^2)) // Expand
g31 = f3 // Expand

m = 44;
h1 = Series[g11, {q, 0, m}]; h1
h2 = Series[g12, {q, 0, m}]; h2
h3 = Series[g21, {q, 0, m}]; h3
h4 = Series[g22, {q, 0, m}]; h4
h5 = Series[g31, {q, 0, m}]; h5

```

Here $f1,f2,f3$ are the new modular forms associated to $\mathbb{Q}$-decomposition of the Jacobian of $X_0(90)/w_5$ $\sim E_{f_1}^2\times E_{f_2}^2\times E_{f_3}$,
and $h1,h2,h3,h4,h5$ are a basis for the modular forms of $S_2(\Gamma_0(90)\cup \langle w_5\rangle ,\mathbb{C})$, where the precision of the $q$-expansion for $h_i$'s
is until $q^{44}$. (Observe that $f_1$ has level 30, $f_2$ has level 45, and $f_3$ has level 90).

We are luckely that all automorphism of the quotient curve are defined over the rational because no quadratic twist exists between $f1,f2$ and $f3$ and no $f_i$ is a CM modular form, 
and this will help us to determine if is bielliptic or not the quotient curve because if exist should be the bielliptic involution defined over the rationals
(see an exemple bellow with quadratic twist between two modular forms). 
(We can suspect that such twist exists by comparing the coefficients of the $q$-expansion if only differs one to the other by multiplicition by $\pm 1$ or too much zeros appears)

We compute first the Petri model of the genus 5 curve from $h1,...,h5$ as follows:

```

P = {a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, 
   a15}.{x^2, y^2, z^2, t^2, s^2, x y, x z, x t, x s, y z, y t, y s, 
   z t, z s, t s}
Q = P /. {x -> h1, y -> h2, z -> h3, t -> h4, s -> h5}; 
l =  Table[ Coefficient[Q, q, i], {i, 2, 24}]; 
T =  Solve[{l == {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0}}, {a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, 
    a12, a13, a14, a15}][[1]]

QQ = P /. T // Factor // Numerator

QQ /. {x -> h1, y -> h2, z -> h3, t -> h4, s -> h5}
```

The output QQ gives us the Petri model with the free variables a1, a2 and a6, 
because Petri's model are exactly three degree two equations with variables $x,y,z,t$ and $s$ for a non-trigonal genus 5 curve.

Now in order to check if the modular curve is bielliptic for the non-repeated factors, 
we implement in Mathematica the Prop.2.6 in JA Bars-González paper ``Bielliptic modular curves X_0^*(N)".
If $E_{f_3}$ should be a bielliptic quotient (we recall that we can work all over $\mathbb{Q}$ because no quadratic twist between $f_1,f_2$ and $f_3$), by Prop.2.6
we need to check that QQ3 (bellow is zero=bielliptic or Not), because the variables $x,y,z,t,s$ follows from the $\mathbb{Q}$-decomposition of the quotient curve
and the one corresponding to $E_{f_3}$ is $s$ by construction:

```

QQ3 = (QQ - (QQ /. s -> -s))/(4 s) // Expand // Factor // Numerator
l = {Coefficient[QQ3, x, 1], Coefficient[QQ3, y, 1], 
   Coefficient[QQ3, z, 1], Coefficient[QQ3, t, 1], 
   Coefficient[QQ3, s, 1]} // Factor

```
If $E_{f_3}$ were a bielliptic quotient QQ3 should be zero, in particular Coefficient[QQ3, t, 1] is zero, but such coefficient is equal to
$-6 a6$, which imposes that $a6=0$, but this impose a condition on $a1, a2, a6$ for the general Petri  model of $X_0(90)/w_5$, and therefore $E_{f_3}$ is not a bielliptic quotient
(by Prop. 2.6 in the loc.cit. paper in JA).

In order to obtain that do not have any elliptic quotient, (because all defined over the rationals) is enought that does not have as elliptic quotient any of the elliptic curves $E_{f_i}$ with $i=1,2$.
Because $E_{f_i}$ with $i=1,2$ appears REPEATED in the Jacobian decomposition we need to consider matrices of size 2x2 (repeated 2 times in the Jacobian decomposition of the curve).

We used this ad-hoc Mathematica code to implement the result of Bars-González in JA paper ``Bielliptic modular curves $X_0^*(N)$´´ when a possible ellipic quotient appear repeated two times in the Jacobian decomposition.

```
R1 = QQ /. {x -> aa1 x + aa2 y, y -> bb1 x + bb2 y};
QQ1 = (R1 - (R1 /. x -> -x))/(4 x) // Expand // Factor // Numerator
l = {Coefficient[QQ1, x, 1], Coefficient[QQ1, y, 1], 
   Coefficient[QQ1, z, 1], Coefficient[QQ1, t, 1], 
   Coefficient[QQ1, s, 1]} // Factor
   
 l1 = l /. {aa1 -> 0, bb1 -> 1} // Factor
 l1 = l /. {aa1 -> 1} // Factor

```

where in the first we consider 2x2 matrices with aa1=0 (projective matrices, recall), and the second with =1: the result is respectively

```
{0, 3 (a6 aa2 + 2 a2 bb2), 3 a6, 0, 0} (case aa1=0)

{0, 3 (2 a1 aa2 + a6 aa2 bb1 + a6 bb2 + 2 a2 bb1 bb2), -2 a1 - 2 a2 + 
  3 a6 bb1, 0, 0} (case aa1 neq 0)
```
Thus not has $E_{f_1}$ as bielliptic quotient because with $aa1=0$ we need to impose $a6=0$, and for $aa1=1$  there is no matrix making all zeros independent of 
$a1$,$a2$ and $a6$.

For the case $E_{f_2}$ is similar, and thus conclude that $X_0(90)/w_5$ is not bielliptic.

*Another example. Bielliptic case.

We consider $X_0(90)/\langle w_9\rangle$ that has genus 5. By Magma programme in another folder of my github we obtain the Q-Jacobian decomposition of the curve, and for it we obtain
the Petri model if exists.

```mathematica

f1 = q - q^2 - q^3 - q^4 + q^5 + q^6 + 3*q^8 + q^9 - q^10 - 4*q^11 + 
  q^12 - 2*q^13 - q^15 - q^16 + 2*q^17 - q^18 + 4*q^19 - q^20 + 
  4*q^22 - 3*q^24 + q^25 + 2*q^26 - q^27 - 2*q^29; (E15a)
  
f2 = q - q^2 + q^3 + q^4 - q^5 - q^6 - 4*q^7 - q^8 + q^9 + q^10 + q^12 + 
  2*q^13 + 4*q^14 - q^15 + q^16 + 6*q^17 - q^18 - 4*q^19 - q^20 - 
  4*q^21 - q^24 + q^25 - 2*q^26 + q^27 - 4*q^28 - 6*q^29; (E30a)
  
f3 =  q - q^2 + q^4 + q^5 + 2*q^7 - q^8 - q^10 + 6*q^11 - 4*q^13 - 2*q^14 +
   q^16 - 6*q^17 - 4*q^19 + q^20 - 6*q^22 + q^25 + 4*q^26 + 2*q^28 - 
  6*q^29; (E90a)
  
f4 = q + q^2 + q^4 - q^5 + 2*q^7 + q^8 - q^10 - 6*q^11 - 4*q^13 + 2*q^14 +
   q^16 + 6*q^17 - 4*q^19 - q^20 - 6*q^22 + q^25 - 4*q^26 + 2*q^28 + 
  6*q^29; (E90b)
  
f11 = f1 + 2 (f1 /. q -> q^2);
f12 = f1 - 2 (f1 /. q -> q^2);
g1 = f11 + 3 (f11 /. q -> q^3);
g2 = f12 + 3 (f12 /. q -> q^3);
g3 = f2 - 3 (f2 /. q -> q^3);
g4 = f3;
g5 = f4;

h1 = Series[g1, {q, 0, 23}];
h2 = Series[g2, {q, 0, 23}]; 
h3 = Series[g3, {q, 0, 23}]; 
h4 = Series[g4, {q, 0, 23}]; 
h5 = Series[g5, {q, 0, 23}]; 
P0 = {a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, 
   a15}.{x^2, y^2, z^2, t^2, s^2, x y, x z, x t, x s, y z, y t, y s, 
   z t, z s, t s}
   
Q = P0 /. {x -> h1, y -> h2, z -> h3, t -> h4, s -> h5}; l = 
 Table[ Coefficient[Q, q, i], {i, 2, 25}]; T = 
 Solve[{l == {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0}}, {a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11,
     a12, a13, a14, a15}][[1]]
 QQ = (P0 /. T) // Expand // Factor // Numerator;
QQ /. {x -> h1, y -> h2, z -> h3, t -> h4, s -> h5}    
     

```

QQ gives us the Petri model equation. Observe that the Q-factorization is $(E15a)^2\times E30a\times E90a\times E90b$ where we follow Cremona's notation for elliptic curves.

In order to search if is bielliptic or not here we observe that $f3$ and $f4$ from the coefficient are candidates to become quadratic twists, and effectively:

```mathematica

l4 = Table[Coefficient[f4, q, Prime[i]], {i, 2, 10}]
l3 = Table[Coefficient[f3, q, Prime[i]], {i, 2, 10}]
Table[ l3[[i]] JacobiSymbol[-3, Prime[i + 1]] - l4[[i]], {i, 1, 9}]


```

Thus all automorphism are defined over $\mathbb{Q}(\sqrt{-3})$ and the $\mathbb{Q}(\sqrt{-3})$-factorization of the Jacobian is
$$E15a^2\times E30a\times (E90a)^2$$ because $E90a\sim_{\mathbb{Q}(\sqrt{-3})} E90b$.
The only possible bielliptic quotient are the elliptic curves listed in the decomposition, let us do the computation of all bielliptic involution in the modular quotient curve $X_0(90)/w_9$.

```mathematica
QQx = (QQ - (QQ /. z -> -z)) // Factor
lx = {Coefficient[QQx, x], Coefficient[QQx, y], Coefficient[QQx, z], 
  Coefficient[QQx, t], Coefficient[QQx, s]}
Solve[lx == {0, 0, 0, 0, 0}, {a1, a3, a12}]

QQx = (QQ - (QQ /. t -> -t)) // Factor
lx = {Coefficient[QQx, x], Coefficient[QQx, y], Coefficient[QQx, z], 
  Coefficient[QQx, t], Coefficient[QQx, s]}
Solve[lx == {0, 0, 0, 0, 0}, {a1, a3, a12}]

QQx = (QQ - (QQ /. s -> -s)) // Factor
lx = {Coefficient[QQx, x], Coefficient[QQx, y], Coefficient[QQx, z], 
  Coefficient[QQx, t], Coefficient[QQx, s]}
Solve[lx == {0, 0, 0, 0, 0}, {a1, a3, a12}]

```
We use above if $E30a$, $E90a$ or $E90b$ is a bielliptic quotient of $X_0(90)/w_9$ over the rationals. The result is:

```mathematica
{0, 0, -2 a12 t, -2 a12 z, 0}
{0, 0, -2 a12 t, -2 a12 z, 0}
{0, 2 a12 s, 0, 0, 2 a12 y}

```
Thus is not independent of the three free variables of the Petri model which ar now a1,a2 and a12.

Over the rationals remain to study the factor $(E15a)^2$ in order if any bielliptic involution or not appears.
```mathematica
R2 = QQ /. {x -> aa1 x + aa2 y, y -> bb1 x + bb2 y};
R2simx = (R2 - (R2 /. x -> -x))/(4 x) // Expand // Factor // Numerator;
l = {Coefficient[%, y, 1], Coefficient[%, z, 1], Coefficient[%, s, 1],
     Coefficient[%, t, 1]} // Factor;
l[[2]];
l1 = l /. {aa1 -> 0} // Factor

l1 = l /. {aa1 -> 1} // Factor
```
{2 a2 bb1 bb2, 0, a12 bb1, 0} (case aa1=0, thus bb1=0 and not invertible 2x2 matrix, no involution)
{2 (a1 aa2 + a2 bb1 bb2), 0, a12 bb1, 0} (case aa1=1, thus bb1=0 and aa2=0, Bielliptic involution with quotient elliptic curve E15a).

Now remains if bielliptic involutions over $\mathbb{Q}(\sqrt{-3})$ could appear, we need only to check for the variables t,s corresponding to $E90a$ and $E90b$.

```mathematica
R2 = QQ /. {t -> aa1 t + aa2 s, s -> bb1 t + bb2 s};
R2simx = (R2 - (R2 /. t -> -t))/(4 t) // Expand // Factor // Numerator;
l = {Coefficient[%, x, 1], Coefficient[%, y, 1], Coefficient[%, z, 1],
     Coefficient[%, s, 1], Coefficient[%, t, 1]} // Factor;

l1 = l /. {aa1 -> 0 } // Factor

l1 = l /. {aa1 -> 1} // Factor
```

Obtaining:
{0, a12 bb1, 0, -2 (a1 - a2) bb1 bb2, 0} (case aa1=0, thus bb1=0, no invertible matrix)
{0, a12 bb1, -a12, -2 (a1 aa2 + 5 a2 aa2 + a1 bb1 bb2 - a2 bb1 bb2),
  0} (case aa1=1, thus bb1=0 but factor -a12, thus no possibility to be zerro array, thus no new involution)
