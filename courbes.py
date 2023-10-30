from functions import *
# φ =  (porosity), σ =  (resistivity), αh = (tortuosity).

dictionnaire = {
    "foam": [0.99, 14000, 1.02],
    "pine": [0.54, 35540, 1.22],
    "oak": [0.67, 197525, 1.36]
}

W=np.linspace(120,30000,10)
phi,sigma,alpha_h=dictionnaire["foam"]
Re1,Im1=minimise(g2,W)
phi,sigma,alpha_h=dictionnaire["pine"]
Re2,Im2=minimise(g2,W)
phi,sigma,alpha_h=dictionnaire["oak"]
Re3,Im3=minimise(g2,W)
print(phi)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

# Tracer les parties réelles
ax1.plot(W, Re1, label='Re(g1)', color='blue')
ax1.plot(W , Re2, label='Re(g2)', color='green')
ax1.plot(W, Re3, label='Re(g3)', color='red')
ax1.set_xlabel('w')
ax1.set_ylabel('Re(gi)')
ax1.legend()

# Tracer les parties imaginaires
ax2.plot(W, Im1, label='Im(g1)', color='blue')
ax2.plot(W, Im2, label='Im(g2)', color='green')
ax2.plot(W, Im3, label='Im(g3)', color='red')
ax2.set_xlabel('w')
ax2.set_ylabel('Im(gi)')
ax2.legend()

plt.tight_layout()
plt.show()




plt.plot(W,Re)
plt.show()
