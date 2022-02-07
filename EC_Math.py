# -*- coding: utf-8 -*-
#
# Elliptic curve point operations (Python 2.7)
# Copyright (c) Denis Leonov <466611@gmail.com>
#

def OnCurve(x,y): # Check if the point is on the curve
    A = (y*y)%P
    B = (x*x*x)%P
    C = False
    if A == (B + 7):
        C = True
    return C

def legendre_symbol(a,p):
    ls = pow(a, (p - 1) / 2, p)
    return -1 if ls == p - 1 else ls

def modsqrt(a,p): # Square root A modulo P
    if legendre_symbol(a, p) != 1:
        return 0
    elif a == 0:
        return 0
    elif p == 2:
        return p
    elif p % 4 == 3:
        return pow(a, (p + 1) / 4, p)
    s = p - 1
    e = 0
    while s % 2 == 0:
        s /= 2
        e += 1
    n = 2
    while legendre_symbol(n, p) != -1:
        n += 1
    x = pow(a, (s + 1) / 2, p)
    b = pow(a, s, p)
    g = pow(n, s, p)
    r = e
    while True:
        t = b
        m = 0
        for m in xrange(r):
            if t == 1:
                break
            t = pow(t, 2, p)
        if m == 0:
            return x
        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m

def modinv(a,n): # Extended Euclidean Algorithm in elliptic curves 
    lm, hm = 1,0
    low, high = a%n,n
    while low > 1:
        ratio = high/low
        nm = hm - lm * ratio
        new = high - low * ratio
        hm = lm
        high = low
        lm = nm
        low = new
    return lm % n

# ECadd taken from https://github.com/wobine/blackboard101/blob/master/EllipticCurvesPart4-PrivateKeyToPublicKey.py
def ECadd(xp,yp,xq,yq): # EC point addition
    m = ((yq-yp) * modinv(xq-xp,P))%P
    xr = (m*m-xp-xq)%P
    yr = (m*(xp-xr)-yp)%P
    return (xr,yr)

# ECdouble taken from https://github.com/wobine/blackboard101/blob/master/EllipticCurvesPart4-PrivateKeyToPublicKey.py
def ECdouble(xp,yp): # EC point doubling
    LamNumer = 3*xp*xp+Acurve
    LamDenom = 2*yp
    Lam = (LamNumer * modinv(LamDenom,P)) % P
    xr = (Lam*Lam-2*xp) % P
    yr = (Lam*(xp-xr)-yp) % P
    return (xr,yr)

# ECmul taken from https://github.com/wobine/blackboard101/blob/master/EllipticCurvesPart4-PrivateKeyToPublicKey.py
def ECmul(xs,ys,Scalar): # EC point multiplication
    if Scalar == 0 or Scalar >= N: raise Exception("Invalid Scalar/Private Key")
    ScalarBin = str(bin(Scalar))[2:]
    Qx,Qy=xs,ys
    for i in range (1,len(ScalarBin)):
        Qx,Qy = ECdouble(Qx,Qy)
        if ScalarBin[i] == "1":
            Qx,Qy = ECadd(Qx,Qy,xs,ys)
    return (Qx,Qy)

def ECsub(xp,yp,xq,yq): # EC point subtraction
    X = (((yp+yq)*modinv(xq-xp,P))**2-xp-xq)%P
    A = (xp + X + xq)%P
    B = modsqrt(A,P)
    B1 = P - B
    Y = yq - (xq - X) * B
    X = X % P
    Y = Y % P
    if not OnCurve(X,Y):
        Y = yq - (xq - X) * B1
    Y = Y % P
    return X,Y

def ECdiv(Qx,Qy,Scalar): # EC point division
    A = (N-1)/Scalar
    Px,Py = ECmul(Qx,Qy,A)
    Py = P-Py
    return Px,Py

P = 2**256 - 2**32 - 2**9 - 2**8 - 2**7 - 2**6 - 2**4 - 1
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx = 55066263022277343669578718895168534326250603453777594175500187360389116729240
Gy = 32670510020758816978083085130507043184471273380659243275938904335757337482424
Acurve = 0


print (Gx,Gy)

Ax,Ay = ECadd(Gx,Gy,Gx,Gy)
print (Ax,Ay)

Bx,By = ECsub(Ax,Ay,Gx,Gy)
print (Bx,By)

Cx,Cy = ECdouble(Gx,Gy)
print (Cx,Cy)

Dx,Dy = ECdiv(Cx,Cy,2)
print (Dx,Dy)
