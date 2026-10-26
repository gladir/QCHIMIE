# QChimie
Liste de commandes écrite en Pascal (Turbo Pascal ou Free Pascal) sur la chimie.

![image](https://user-images.githubusercontent.com/11842176/173241903-be1d3f65-5db1-4675-b16a-ba13d25c9800.png)

Aperçu de la commande PyMOL :

<img width="638" height="477" alt="image" src="https://github.com/user-attachments/assets/1b5824e1-6cc2-4585-bce9-04c317f3e912" />


Aperçu de la commande ChemDraw :

<img width="640" height="480" alt="ChemDraw" src="SAMPLES/CHEMDRAW/APERCU.PNG" />

Aperçu de la commande VCCTL :

<img width="640" height="480" alt="VCCTL" src="SAMPLES/VCCTL/APERCU.PNG" />

Aperçu de la commande ORCA :

<img width="640" height="480" alt="ORCA" src="SAMPLES/ORCA/APERCU.PNG" />

<h2>Liste des fichiers</h3>

Voici la liste des différents fichiers proposés dans QChimie :

<table>
  <tr>
    <th>Nom</th>
    <th>Description</th>
  </tr>
  <tr>
      <td><b>ELEMENTS.PAS</b></td>
      <td>Cette commande permet d'afficher le tableau des éléments chimiques ou un élément chimique.</td>
  </tr>
  <tr>
      <td><b>FORMULA.PAS</b></td>
      <td>Cette commande permet d'afficher la formule d'un élément chimique ou d'un composé chimique.</td>
  </tr>	
  <tr>
      <td><b>FINDELEM.PAS</b></td>
      <td>Cette commande permet de détecter les éléments chimiques mentionné dans un texte.</td>
  </tr>
  <tr>
	  <td><b>PYMOL.PAS</b></td>
	  <td>Cette commande permet de visualiser, manipuler et analyser des structures moléculaires tridimensionnelles, not amment des protéines, des acides nucléiques et des composés chimiques, à l'aide du logiciel PyMOL.</td>
  </tr>
  <tr>
	  <td><b>CHEMDRAW.PAS</b></td>
	  <td>Cette commande permet de dessiner des structures chimiques en deux dimensions (molécules, réactions, flèches, cycles, groupes d'atomes), de les analyser (formule, poids moléculaire, masse exacte, spectre de masse, prédiction RMN ¹H et ¹³C, nom du composé, SMILES, InChI, stéréochimie R/S et E/Z) et de les importer ou exporter aux formats CDT, MOL, SMILES, InChI, SVG et BMP, à la manière du logiciel ChemDraw.</td>
  </tr>
  <tr>
	  <td><b>VCCTL.PAS</b></td>
	  <td>Cette commande permet de simuler, à la manière du logiciel VCCTL (Virtual Cement and Concrete Testing Laboratory), la microstructure 3D et l'hydratation d'une pâte de ciment (ciment Portland, cendres volantes, laitier, calcaire, fumée de silice) par automate cellulaire, puis d'estimer ses propriétés : degré d'hydratation, chaleur dégagée, retrait chimique, porosité, percolation, distribution des pores, modules élastiques, résistance à la compression, conductivité et diffusivité, zone de transition interfaciale. Elle comprend aussi une bibliothèque de matériaux, le calculateur de Bogue, l'empilement 3D de granulats, le calibrage de la cinétique et un mode lot (<code>VCCTL /B</code>).</td>
  </tr>
  <tr>
	  <td><b>ORCA.PAS</b></td>
	  <td>Cette commande permet d'effectuer, à la manière du logiciel ORCA, des calculs de chimie quantique sur des molécules à partir d'un fichier d'entrée au format ORCA (lignes <code>! mots-clés</code>, blocs <code>%scf</code>, <code>%geom</code>, <code>%tddft</code>, <code>%freq</code> et coordonnées <code>* xyz</code>, <code>* int</code> ou <code>* xyzfile</code>) : énergie Hartree-Fock (RHF/UHF), MP2, états excités CIS et DFT (LDA, BLYP, B3LYP, PBE, PBE0) avec les bases STO-3G, 3-21G et 6-31G, optimisation de géométrie (avec contraintes), balayage relaxé de surface, fréquences harmoniques, spectre infrarouge et thermochimie. Elle produit un fichier de sortie <code>.OUT</code> (énergies, orbitales, charges de Mulliken et de Löwdin, ordres de liaison de Mayer, moment dipolaire) et propose un éditeur d'entrée, la visualisation 3D de la molécule, le diagramme et les coupes des orbitales moléculaires ainsi que les graphiques SCF, optimisation, IR, UV et surface d'énergie (<code>ORCA /B fichier.INP</code> en mode lot).</td>
  </tr>
</table>

<h2>Compilation</h2>
	
Les fichiers Pascal n'ont aucune dépendances, il suffit de télécharger le fichier désiré et de le compiler avec Free Pascal avec la syntaxe de commande  :

<pre><b>fpc</b> <i>LEFICHIER.PAS</i></pre>
	
Sinon, vous pouvez également le compiler avec le Turbo Pascal à l'aide de la syntaxe de commande suivante :	

<pre><b>tpc</b> <i>LEFICHIER.PAS</i></pre>
	
Par exemple, si vous voulez compiler ELEMENTS.PAS, vous devrez tapez la commande suivante :

<pre><b>fpc</b> ELEMENTS.PAS</pre>

<h2>Remarque</h2>

<ul>
	<li>La commande CHEMDRAW accepte en ligne de commande un fichier à ouvrir (CDT, MOL, SDF ou SMI) ou une structure SMILES avec le paramètre -s (par exemple : <i>CHEMDRAW -s CC(=O)Oc1ccccc1C(O)=O</i>). Des exemples sont fournis dans le dossier SAMPLES/CHEMDRAW. Le format CDT est le format texte natif de la commande.</li>
	<li>L'InChI produit par la commande CHEMDRAW (menu Analyse) couvre la formule, la connectivité, les hydrogènes (y compris mobiles), la charge et la stéréochimie (couches /b, /t, /m et /s, calculées à partir des coins pleins, des hachures et de la géométrie des doubles liaisons) d'une molécule connexe, sans isotopes. Sur un échantillon de 4 900 molécules sans stéréochimie, il est identique à celui de la bibliothèque InChI officielle dans 98,8 % des cas ; sur 3 000 stéréoisomères la couche stéréochimique est identique dans tous les cas où la structure est correctement lue, hormis les doubles liaisons de grands cycles et quelques cations conjugués. Les fragments multiples et les organométalliques ne sont pas pris en charge.</li>
	<li>La stéréochimie de la commande CHEMDRAW (menu Affichage &gt; Stéréochimie, ou menu Analyse &gt; Stéréochimie) attribue les descripteurs CIP R/S aux centres dessinés avec un coin plein ou hachuré (extrémité fine sur le centre) et les descripteurs E/Z aux doubles liaisons substituées. Les priorités sont calculées avec un digraphe hiérarchique (atomes dupliqués, hydrogènes implicites, doublets libres) selon les règles 1 à 5 (numéro atomique, nombre de masse, Z avant E, unités chirales avant pseudo-asymétriques, like avant unlike, r avant s, R avant S) : les centres pseudo-asymétriques (r/s), les alditols symétriques et les cycles symétriques 1,3- ou 1,4-disubstitués sont donc étiquetés. Les symboles @ / @@ et / \ du SMILES sont lus (ils deviennent des coins pleins ou hachurés) et écrits, y compris pour les doubles liaisons des grands cycles (jusqu'à 40 atomes, dessinés avec des coudes pour les liaisons trans) ; la stéréochimie est conservée par la commande Nettoyer. Sur 3 001 stéréoisomères tirés d'une collection de molécules, 2 982 (99,4 %) reçoivent exactement les mêmes descripteurs que la bibliothèque RDKit ; les écarts restants concernent des systèmes pontés ou fusionnés (bicycles, azote en tête de pont), des complexes métalliques et quelques cycles symétriques couplés. L'InChI exporté contient la couche stéréochimique (/b, /t, /m, /s).</li>
	<li>La commande ORCA accepte en ligne de commande un fichier d'entrée à ouvrir (INP) ou une géométrie XYZ, ou le paramètre <i>/B</i> pour exécuter un calcul sans interface (par exemple : <i>ORCA /B EAU.INP EAU.OUT</i>). Des exemples d'entrées (avec deux sorties) sont fournis dans le dossier SAMPLES/ORCA. Le moteur de calcul est entièrement écrit en Pascal : intégrales gaussiennes (schéma de McMurchie-Davidson), SCF avec DIIS, MP2 à cœur gelé, CIS, DFT sur grille de Becke et dérivées numériques pour les gradients, les fréquences et les intensités IR. Sur les molécules testées (eau, OH, NH3, éthane, benzène, O2…), les énergies HF, MP2, CIS, B3LYP, PBE et PBE0, les charges, les dipôles, les fréquences harmoniques et les valeurs de &lt;S²&gt; concordent avec celles de la bibliothèque PySCF à environ 10<sup>-6</sup> Eh près (quelques 10<sup>-6</sup> Eh pour la DFT, à cause de la grille d'intégration). Les bases se limitent à STO-3G (H à Ne), 3-21G et 6-31G (H, He, C à F) ; les autres bases d'ORCA (def2-SVP, etc.) sont remplacées par 6-31G avec un avertissement. Les calculs DFT et les optimisations sont plus lents qu'avec ORCA (pas de gradient analytique, pas d'approximation RI) : prévoyez quelques dizaines de secondes pour optimiser une petite molécule en B3LYP.</li>
	<li>La commande FORMULA a été bonifié du paramètre /FORMAT:HTML à partir de sa version 1.1.</li>
</ul>

<h2>Licence</h2>
<ul>
 <li>Le code source est publié sous la licence <a href="https://github.com/gladir/QCHIMIE/blob/main/LICENSE">MIT</a>.</li>
 <li>Le paquet original est publié sous la licence <a href="https://github.com/gladir/QCHIMIE/blob/main/LICENSE">MIT</a>.</li>
</ul>
