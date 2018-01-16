# Code and data for structure of mixed weight Hilbert modular forms `(K = Q(sqrt(5)))`

## Fourier coefficients of generators
In directories `"./forms_csv/A*"`, Fourier coefficients of generators are stored as csv files. Each generator corresponds to each csv file. We fixed an order of generators. For example, `A1` has 3 generators. The first generator is of weight `(7, 9)` and some of its Fourier coefficients are stored in `"./forms_csv/A1/gen0_wt_7_9.csv"`. The second generator is of weight `(8, 10)` and the third generator is of weight `(11, 13)`. Note that this order may not be the ascending order with respect to weights of generators. In fact, the weights of generators of `A7` is ordered so that `(5, 19)`, `(4, 18)`, `(7, 21)`, `(8, 22)`. 
Each line of the csv file shows the Fourier coefficient of a generator. `a(v, u)` is corresponding to `c(v, u)` in our paper. If the corresponding value `a(v, u)` for `(v, u)` is an integer, then `a(v, u)` is corresponding to the integer itself. If `a(v, u)` is a tuple `(a, b)` then `a(v, u)` is corresponding to `(a + b sqrt(5))/2`.

## Relations
In "./relations", relations of generators are stored as text files. For example, relations of generators of `A7` is stored in "./relations/a7_rel". Each line in this file corresponds to each relation among generators. In this case, there are two relations. The first line is `[-1496880000*g5, 43200*g6 + 29*g2^3, 0, -g2]` and this shows coefficients of generators. We use the order of generators explained above. Let `F5`, `F4`, `F7` and `F8` be the generators of `A7` of weight `(5, 19)`, `(4, 18)`, `(7, 21)`, `(8, 22)` respectively. Then it satisfies `-1496880000*g5 * F5 + (43200*g6 + 29*g2^3) * F4 - g2 * F8 = 0`.