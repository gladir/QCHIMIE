# Exemple de script PML pour PyMOL Clone - QChimie
# Charge un tripeptide, configure la vue et applique des couleurs.

load ALANINE.PDB

# Mode de rendu
show ball-and-stick

# Coloration
color cyan all
color yellow hetatm
color red water

# Selection nommee
select bb backbone

# Vue
bg_color black
zoom

# Active la rotation automatique
spin on

echo Script charge avec succes.
