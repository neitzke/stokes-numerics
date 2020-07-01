import math

def latex_inf(x):
    if x > 0:
        return r'$\infty$'
    else:
        return r'$-\infty$'

def latex_real(x,prec=8):
    if math.isinf(x):
        return latex_inf(x)
    formatstring = '{:.'+str(prec)+'e}'
    return r"$\num{" + formatstring.format(x) + r"}$"

def latex_exp(x):
    if abs(x) < 0.01:
        x = 0.0
    return r'$\exp(%g)$' % x

_denoms = [ 2, 3, 4, 5, 6 ]

def tofrac(x):
    for d in _denoms:
        if abs(d*x - int(d*x)) < 1e-5 and abs(int(d*x))<15:
            return int(d*x), d
    return None
        
def pretty_float(a,frac=True):
    '''Latex a floating point number, converting to a fraction if possible'''
    if a==0:
        return []
    if a > 0:
        signum = '+'
    else:
        signum = '-'

    scalarstr = '%g' % abs(a)
    if frac and abs(a - int(a)) > 1e-5:
        r = tofrac(abs(a))
        if r:
            scalarstr = '\\frac{%d}{%d}' % r

    return [signum, scalarstr]

def pretty_monomial(k,a,frac=True):
    '''Nice string for a*z^k'''
    if a == 0:
        return []
    if k == 0:
        return pretty_float(a,frac)
    s,r = pretty_float(a,frac)
    if r == '1':
        r = ''
    if k == 1:
        return [s,r + 'z']
    else:
        return [s,r + 'z^%d'%k ]

def pretty_polynomial(coefs,frac=True):
    '''Latex a polynomial in z'''
    L = []
    for k,a in enumerate(coefs):
        L = pretty_monomial(k,a,frac) + L
    if L[0] == '+':
        L = L[1:]
    return ' '.join(L)
    
