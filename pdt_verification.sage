# PDT Core Verifications — SageMath

from sage.all import *

# ---- PDT constants ----
R = RealField(50)
t = polygen(R, 't')
rho = max(r for r in (t**3 - t - 1).roots(ring=R, multiplicities=False) if r > 0)
Q_roots = (t**4 - t - 1).roots(ring=R, multiplicities=False)
Q = max(r for r in Q_roots if r > 0)
lam3 = 1 - 1/rho
lam4 = 1 - 1/Q
chi = Q/rho
rhoQ = rho * Q

print(f"rho       = {rho}")
print(f"Q         = {Q}")
print(f"lambda_3  = {lam3}")
print(f"lambda_4  = {lam4}")
print(f"chi = Q/rho = {chi}")
print(f"rho*Q     = {rhoQ}")
print()

# ---- Elliptic curves ----
E3 = EllipticCurve('368f1')   # corresponds to x^3 - x - 1
E4 = EllipticCurve('1132b1')  # corresponds to x^4 - x - 1
print(f"E_3 = 368f1:  conductor={E3.conductor()}, rank={E3.rank()}, disc={E3.discriminant()}")
print(f"E_4 = 1132b1: conductor={E4.conductor()}, rank={E4.rank()}, disc={E4.discriminant()}")
print()

# ---- Periods ----
Om3 = E3.period_lattice().omega().real()
Om4 = E4.period_lattice().omega().real()
print(f"Omega_3 (real period E_3) = {Om3}")
print(f"Omega_4 (real period E_4) = {Om4}")
print()

# ---- Strong coupling: alpha_s = Omega_4 / 23 ----
alpha_s_pred = Om4 / 23
alpha_s_PDG = 0.1179
print(f"Omega_4 / 23       = {alpha_s_pred}")
print(f"alpha_s(M_Z) (PDG) = {alpha_s_PDG}")
print(f"Error              = {abs(float(alpha_s_pred) - alpha_s_PDG)/alpha_s_PDG * 100:.4f}%")
print()

# ---- L-values ----
L3_at_1  = E3.lseries().at1()[0]           # L(E_3, 1), rank 0
L4_prime = E4.lseries().deriv_at1()[0]     # L'(E_4, 1), rank 1
print(f"L(E_3, 1)  = {L3_at_1}")
print(f"L'(E_4, 1) = {L4_prime}")
print()

# ---- sin^2 theta_W candidates ----
sinW_PDG = 0.23122
algebraic = float(lam4 / chi**3)
Lfunc    = float(L3_at_1 * L4_prime / 23)
print("sin^2 theta_W candidates:")
print(f"  lam4 / chi^3            = {algebraic:.6f}  err {abs(algebraic - sinW_PDG)/sinW_PDG*100:.4f}%")
print(f"  L(E3,1) * L'(E4,1) / 23 = {Lfunc:.6f}  err {abs(Lfunc - sinW_PDG)/sinW_PDG*100:.4f}%")
print()

# ---- gamma_BI candidates ----
gBI_PDG = 0.23954
print("gamma_BI candidates:")
print(f"  lam4 * rho      = {float(lam4*rho):.6f}  err {abs(float(lam4*rho) - gBI_PDG)/gBI_PDG*100:.4f}%")
print(f"  Om3^2 / Om4^3   = {float(Om3**2/Om4**3):.6f}  err {abs(float(Om3**2/Om4**3) - gBI_PDG)/gBI_PDG*100:.4f}%")
print()

# ---- Number field invariants ----
x = polygen(QQ, 'x')
K3 = NumberField(x**3 - x - 1, 'a')
K4 = NumberField(x**4 - x - 1, 'b')
print(f"Q(rho): disc={K3.discriminant()}, class#={K3.class_number()}, regulator={K3.regulator()}")
print(f"Q(Q):   disc={K4.discriminant()}, class#={K4.class_number()}, regulator={K4.regulator()}")
print(f"Regulator(Q(rho)) / ln(rho) = {float(K3.regulator() / log(rho)):.10f}")
print()

# ---- Compositum and N(rho*Q) ----
result = K3.composite_fields(K4, both_maps=True)[0]
L, phi3, phi4, _ = result  # L is the compositum; phi3: K3 -> L; phi4: K4 -> L
print(f"Compositum Q(rho, Q): degree={L.degree()}, disc={L.discriminant()}")
rho_L = phi3(K3.gen())
Q_L   = phi4(K4.gen())
N_rhoQ = (rho_L * Q_L).norm()
print(f"N(rho * Q) in Q(rho, Q) = {N_rhoQ}")
if N_rhoQ == -1:
    print("*** UNIT NORM IDENTITY N(rho*Q) = -1 CONFIRMED ***")
elif N_rhoQ == 1:
    print("*** N = +1; trying sign variants for the -1 form ***")
    for s1, s2 in [(-1, 1), (1, -1), (-1, -1)]:
        N = (s1*rho_L * s2*Q_L).norm()
        print(f"  N({s1}*rho * {s2}*Q) = {N}")
        if N == -1:
            print(f"  *** FOUND: sign variant ({s1}, {s2}) gives N = -1 ***")
            break
else:
    print(f"*** N = {N_rhoQ}; not ±1 — inspect the embedding ***")
