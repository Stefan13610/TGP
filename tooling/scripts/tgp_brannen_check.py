"""
tgp_brannen_check.py - weryfikacja parametryzacji Brannena dla TGP
"""
import numpy as np
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Brannen: sqrt(m_k) = c*(1 + sqrt(2)*cos(theta + 2*pi*k/3)), k=0,1,2
# Dowod Q_K=3/2 dla DOWOLNEGO theta:
# S1 = sum(sqrt(m)) = c*sum(1+sqrt(2)*cos(phi_k)) = 3c
# S2 = sum(m) = c^2*sum((1+sqrt(2)*cos)^2)
#    = c^2*(3 + 2*sqrt(2)*0 + 2*sum(cos^2))
#    = c^2*(3 + 2*(3/2)) = 6*c^2
# Q_K = 9c^2/6c^2 = 3/2  QED

print("=== Dowod Brannen: Q_K=3/2 dla dowolnego theta ===")
print("S1=3c, S2=6c^2 => Q_K=3/2 (algebraicznie)")
print()

# Test dla kilku theta
for theta_deg in [0, 45, 90, 132.8, 180, 270, 315]:
    theta = np.radians(theta_deg)
    phi = theta + 2*np.pi*np.array([0,1,2])/3
    lam = 1.0*(1 + np.sqrt(2)*np.cos(phi))
    if np.any(lam <= 0):
        print(f"theta={theta_deg:6.1f} deg: niefizyczne (lam<=0)")
        continue
    QK = np.sum(lam)**2 / np.sum(lam**2)
    CV = np.std(lam)/np.mean(lam)
    print(f"theta={theta_deg:6.1f} deg: Q_K={QK:.8f}, CV={CV:.6f}")

print()
print("=== Wyznaczenie theta_K z wartosci TGP ===")
r21  = 105.6583755 / 0.51099895000
r31K = 3477.44
lam_obs = np.array([1.0, np.sqrt(r21), np.sqrt(r31K)])
c_K = np.mean(lam_obs)
print(f"c_K = {c_K:.6f}")

cosvals = (lam_obs/c_K - 1) / np.sqrt(2)
print(f"cos-wartosci: [{cosvals[0]:.6f}, {cosvals[1]:.6f}, {cosvals[2]:.6f}]")
print(f"suma cos: {np.sum(cosvals):.2e}  (oczekiwane 0)")

theta_K = np.arccos(np.clip(cosvals[0], -1, 1))
theta_K_deg = np.degrees(theta_K)
print(f"theta_K = {theta_K_deg:.6f} deg = {theta_K:.8f} rad")
print(f"theta_K/pi = {theta_K/np.pi:.8f}")

# weryfikacja k=1,2
for k in [1,2]:
    pred = np.cos(theta_K + 2*np.pi*k/3)
    print(f"  k={k}: pred={pred:.8f}, obs={cosvals[k]:.8f}, blad={abs(pred-cosvals[k]):.2e}")

# Proste ulamki
print()
print("Szukam prostej postaci theta_K:")
fracs = [(1,6),(1,5),(1,4),(2,7),(1,3),(3,8),(2,5),(5,12),(1,2),
         (7,12),(2,3),(3,4),(5,6),(11,12),(1,1)]
for p,q in fracs:
    cand = p*np.pi/q
    cand_deg = np.degrees(cand)
    diff = abs(theta_K_deg - cand_deg)
    if diff < 3.0:
        print(f"  {p}/{q}*pi = {cand_deg:.4f} deg  (diff={diff:.4f} deg)")

# Sprawdz cos(theta_K)
print(f"cos(theta_K) = {cosvals[0]:.8f}")
print(f"Czy cos(theta_K) = -1/sqrt(2)+eps? cos(135)={np.cos(np.radians(135)):.6f}")
print(f"theta_K vs 3*pi/4 = 135 deg: diff={abs(theta_K_deg-135):.4f} deg")
