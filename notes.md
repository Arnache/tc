# v0.3 (May 2021)
Modernized the code:
- replaced my own old complex numbers library by `<complex>` standard library
- removed global scope `using namespace std;`

# v0.31
- optimized division
- added `inline` to all functions and operators
- compacted the code
- other

# v0.32
- renamed member `z` to `val` and `dz` to `der`
- added sin, cos, tan
- added complex constant I for convenience
- removed `const tc&` returning assignment operators
- tweaks in the division
- removed `typedef double real;` now the class is only for `complex<double>`
