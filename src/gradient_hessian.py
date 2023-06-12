# -*- coding: utf-8 -*-

from sympy import (symbols, diff, simplify)

a, b, c, d = symbols('a b c d')
X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11 = symbols(
    'ex[0] ex[1] ex[2] ex[3] ex[4] ex[5] ex[6] ex[7] ex[8] ex[9] ex[10] ex[11]'
)
Z0, Z1, Z2, Z3 = symbols('ez[0] ez[1] ez[2] ez[3]')

v1 = a + b * X0 + c * X1 + d * X2 - Z0
v2 = (d * d * X5
      + 2. * c * d * X4
      + (2. * b * d + c * c) * X3
      + (2. * a * d + 2. * b * c) * X2
      + (2. * a * c + b * b) * X1
      + 2. * a * b * X0
      + a * a
      - Z1)

v3 = ((d * d * d) * X8
      + (3. * c * d * d) * X7
      + (3. * b * d * d + 3. * c * c * d) * X6
      + (3. * a * d * d + 6. * b * c * d + c * c * c) * X5
      + (6. * a * c * d + 3. * b * b * d + 3. * b * c * c) * X4
      + (a * (6. * b * d + 3. * c * c) + 3. * b * b * c) * X3
      + (3. * a * a * d + 6. * a * b * c + b * b * b) * X2
      + (3. * a * a * c + 3. * a * b * b) * X1
      + 3. * a * a * b * X0
      + a * a * a
      - Z2)

v4 = ((d * d * d * d) * X11
      + (4. * c * d * d * d) * X10
      + (4. * b * d * d * d + 6. * c * c * d * d) * X9
      + (4. * a * d * d * d + 12. * b * c * d * d + 4. * c * c * c * d) * X8
      + (12. * a * c * d * d + 6. * b * b * d * d + 12. * b * c * c * d + c * c * c * c) * X7
      + (a * (12. * b * d * d + 12. * c * c * d) + 12. * b * b * c * d + 4. * b * c * c * c)
      * X6
      + (6. * a * a * d * d
         + a * (24. * b * c * d + 4. * c * c * c)
         + 4. * b * b * b * d
         + 6. * b * b * c * c)
      * X5
      + (12. * a * a * c * d + a * (12. * b * b * d + 12. * b * c * c) + 4. * b * b * b * c)
      * X4
      + (a * a * (12. * b * d + 6. * c * c) + 12. * a * b * b * c + b * b * b * b) * X3
      + (4. * a * a * a * d + 12. * a * a * b * c + 4. * a * b * b * b) * X2
      + (4. * a * a * a * c + 6. * a * a * b * b) * X1
      + (4. * a * a * a * b) * X0
      + a * a * a * a
      - Z3)
obj = v1 * v1 + v2 * v2 + v3 * v3 + v4 * v4

# gradient = [diff(obj, a), diff(obj, b), diff(obj, c), diff(obj, d)]
# for idx in range(4):
#     print(f"gradient[{idx}]={gradient[idx]}")
#
#print("=" * 50)
# hessian = [
#     [diff(obj, a, a), diff(obj, a, b), diff(obj, a, c), diff(obj, a, d)],
#     [diff(obj, b, a), diff(obj, b, b), diff(obj, b, c), diff(obj, b, d)],
#     [diff(obj, c, a), diff(obj, c, b), diff(obj, c, c), diff(obj, c, d)],
#     [diff(obj, d, a), diff(obj, d, b), diff(obj, d, c), diff(obj, d, d)]
# ]
#
# for idx in range(4):
#     for jdx in range(4):
#         print(f"hessian[{idx}][{jdx}]={hessian[idx][jdx]}")
#

jacobian = [
    [diff(v1, a), diff(v1, b), diff(v1, c), diff(v1, d)],
    [diff(v2, a), diff(v2, b), diff(v2, c), diff(v2, d)],
    [diff(v3, a), diff(v3, b), diff(v3, c), diff(v3, d)],
    [diff(v4, a), diff(v4, b), diff(v4, c), diff(v4, d)]
]

for idx in range(4):
    for jdx in range(4):
        print(f"jacobian[{idx}][{jdx}]={jacobian[idx][jdx]}")

"""
jacobian[0][0]=1
jacobian[0][1]=ex[0]
jacobian[0][2]=ex[1]
jacobian[0][3]=ex[2]
jacobian[1][0]=2*a + 2.0*b*ex[0] + 2.0*c*ex[1] + 2.0*d*ex[2]
jacobian[1][1]=2.0*a*ex[0] + 2*b*ex[1] + 2.0*c*ex[2] + 2.0*d*ex[3]
jacobian[1][2]=2.0*a*ex[1] + 2.0*b*ex[2] + 2*c*ex[3] + 2.0*d*ex[4]
jacobian[1][3]=2.0*a*ex[2] + 2.0*b*ex[3] + 2.0*c*ex[4] + 2*d*ex[5]
jacobian[2][0]=3*a**2 + 6.0*a*b*ex[0] + 6.0*c*d*ex[4] + 3.0*d**2*ex[5] + ex[1]*(6.0*a*c + 3.0*b**2) + ex[2]*(6.0*a*d + 6.0*b*c) + ex[3]*(6.0*b*d + 3.0*c**2)
jacobian[2][1]=3.0*a**2*ex[0] + 6.0*a*b*ex[1] + 6.0*c*d*ex[5] + 3.0*d**2*ex[6] + ex[2]*(6.0*a*c + 3*b**2) + ex[3]*(6.0*a*d + 6.0*b*c) + ex[4]*(6.0*b*d + 3.0*c**2)
jacobian[2][2]=3.0*a**2*ex[1] + 6.0*a*b*ex[2] + 6.0*c*d*ex[6] + 3.0*d**2*ex[7] + ex[3]*(6.0*a*c + 3.0*b**2) + ex[4]*(6.0*a*d + 6.0*b*c) + ex[5]*(6.0*b*d + 3*c**2)
jacobian[2][3]=3.0*a**2*ex[2] + 6.0*a*b*ex[3] + 6.0*c*d*ex[7] + 3*d**2*ex[8] + ex[4]*(6.0*a*c + 3.0*b**2) + ex[5]*(6.0*a*d + 6.0*b*c) + ex[6]*(6.0*b*d + 3.0*c**2)
jacobian[3][0]=4*a**3 + 12.0*a**2*b*ex[0] + 12.0*c*d**2*ex[7] + 4.0*d**3*ex[8] + ex[1]*(12.0*a**2*c + 12.0*a*b**2) + ex[2]*(12.0*a**2*d + 24.0*a*b*c + 4.0*b**3) + ex[3]*(2*a*(12.0*b*d + 6.0*c**2) + 12.0*b**2*c) + ex[4]*(24.0*a*c*d + 12.0*b**2*d + 12.0*b*c**2) + ex[5]*(12.0*a*d**2 + 24.0*b*c*d + 4.0*c**3) + ex[6]*(12.0*b*d**2 + 12.0*c**2*d)
jacobian[3][1]=4.0*a**3*ex[0] + 12.0*a**2*b*ex[1] + 12.0*c*d**2*ex[8] + 4.0*d**3*ex[9] + ex[2]*(12.0*a**2*c + 12.0*a*b**2) + ex[3]*(12.0*a**2*d + 24.0*a*b*c + 4*b**3) + ex[4]*(a*(24.0*b*d + 12.0*c**2) + 12.0*b**2*c) + ex[5]*(24.0*a*c*d + 12.0*b**2*d + 12.0*b*c**2) + ex[6]*(12.0*a*d**2 + 24.0*b*c*d + 4.0*c**3) + ex[7]*(12.0*b*d**2 + 12.0*c**2*d)
jacobian[3][2]=4.0*a**3*ex[1] + 12.0*a**2*b*ex[2] + 12.0*c*d**2*ex[9] + 4.0*d**3*ex[10] + ex[3]*(12.0*a**2*c + 12.0*a*b**2) + ex[4]*(12.0*a**2*d + 24.0*a*b*c + 4.0*b**3) + ex[5]*(a*(24.0*b*d + 12.0*c**2) + 12.0*b**2*c) + ex[6]*(24.0*a*c*d + 12.0*b**2*d + 12.0*b*c**2) + ex[7]*(12.0*a*d**2 + 24.0*b*c*d + 4*c**3) + ex[8]*(12.0*b*d**2 + 12.0*c**2*d)
jacobian[3][3]=4.0*a**3*ex[2] + 12.0*a**2*b*ex[3] + 12.0*c*d**2*ex[10] + 4*d**3*ex[11] + ex[4]*(12.0*a**2*c + 12.0*a*b**2) + ex[5]*(12.0*a**2*d + 24.0*a*b*c + 4.0*b**3) + ex[6]*(a*(24.0*b*d + 12.0*c**2) + 12.0*b**2*c) + ex[7]*(24.0*a*c*d + 12.0*b**2*d + 12.0*b*c**2) + ex[8]*(12.0*a*d**2 + 24.0*b*c*d + 4.0*c**3) + ex[9]*(12.0*b*d**2 + 12.0*c**2*d)

"""
