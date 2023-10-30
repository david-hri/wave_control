from functions import *



W=np.linspace(120,30000,100) #la liste des w 
Re1,Im1=minimise(g3,W) #on obtient les listes des parties réelles et imaginaires grace à la fonction de minimisation 

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

# Tracé des parties réelles exemple ici pour la fonction g2 avec le matériau oak
ax1.plot(W, Re1, label='Re(g3)', color='blue')

ax1.set_xlabel('w')
ax1.set_ylabel('Re(gi)')
ax1.legend()
ax1.set_title('Partie Réelle de g3 en fonction de w pour le oak')

# Tracé des parties imaginaires exemple ici pour la fonction g2 avec le matériau oak
ax2.plot(W, Im1, label='Im(g3)', color='blue')
ax2.set_xlabel('w')
ax2.set_ylabel('Im(gi)')
ax2.legend()
ax2.set_title('Partie Imaginaire de g3 en fonction de w pour le oak')

plt.tight_layout()
plt.show()



