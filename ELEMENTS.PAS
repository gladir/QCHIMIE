{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2022
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal, Free Pascal)
}

Program Elements;

Uses Crt;

Type
 Element=Record
  Groupe,Periode:Byte;
  ElementName:String[20];
  NumeroAtomique:Byte;
  SymboleChimique:String[2];
  MasseAtomique:Real;
  Categorie:(NonClasses,Alcalins,Alcalinoterreux,Lanthanides,
  Actinides,MetauxDeTransition,MetauxPauvres,Metalloides,
  AutresNonMetaux,Halogenes,GazNobles);
  Etat:(Gaz,Liquide,Solide);
 End;

Const
 ElementsChimiques:Array[0..117]of Element=(
  (Groupe:1;
   Periode:1;
   ElementName:'Hydrogäne';
   NumeroAtomique:1;
   SymboleChimique:'H';
   MasseAtomique:1.007975;
   Categorie:AutresNonMetaux;
   Etat:Gaz),
  (Groupe:18;
   Periode:1;
   ElementName:'HÇlium';
   NumeroAtomique:2;
   SymboleChimique:'He';
   MasseAtomique:4.002602;
   Categorie:GazNobles;
   Etat:Gaz),
  (Groupe:1;
   Periode:2;
   ElementName:'Lithium';
   NumeroAtomique:3;
   SymboleChimique:'Li';
   MasseAtomique:6.9395;
   Categorie:Alcalins;
   Etat:Solide),
  (Groupe:2;
   Periode:2;
   ElementName:'BÇryllium';
   NumeroAtomique:4;
   SymboleChimique:'Be';
   MasseAtomique:9.0121831;
   Categorie:Alcalinoterreux;
   Etat:Solide),
  (Groupe:13;
   Periode:2;
   ElementName:'Bore';
   NumeroAtomique:5;
   SymboleChimique:'B';
   MasseAtomique:10.8135;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:14;
   Periode:2;
   ElementName:'Carbone';
   NumeroAtomique:6;
   SymboleChimique:'C';
   MasseAtomique:12.0106;
   Categorie:AutresNonMetaux;
   Etat:Solide),
  (Groupe:15;
   Periode:2;
   ElementName:'Azote';
   NumeroAtomique:7;
   SymboleChimique:'N';
   MasseAtomique:14.006855;
   Categorie:AutresNonMetaux;
   Etat:Gaz),
  (Groupe:16;
   Periode:2;
   ElementName:'Oxygäne';
   NumeroAtomique:8;
   SymboleChimique:'O';
   MasseAtomique:15.99940;
   Categorie:AutresNonMetaux;
   Etat:Gaz),
  (Groupe:17;
   Periode:2;
   ElementName:'Fluor';
   NumeroAtomique:9;
   SymboleChimique:'F';
   MasseAtomique:18.99840316;
   Categorie:Halogenes;
   Etat:Gaz),
  (Groupe:18;
   Periode:2;
   ElementName:'NÇon';
   NumeroAtomique:10;
   SymboleChimique:'Ne';
   MasseAtomique:20.1797;
   Categorie:GazNobles;
   Etat:Gaz),
  (Groupe:1;
   Periode:3;
   ElementName:'Sodium';
   NumeroAtomique:11;
   SymboleChimique:'Na';
   MasseAtomique:22.98976928;
   Categorie:Alcalins;
   Etat:Solide),
  (Groupe:2;
   Periode:3;
   ElementName:'MagnÇsium';
   NumeroAtomique:12;
   SymboleChimique:'Mg';
   MasseAtomique:24.3055;
   Categorie:Alcalinoterreux;
   Etat:Solide),
  (Groupe:13;
   Periode:3;
   ElementName:'Aluminium';
   NumeroAtomique:13;
   SymboleChimique:'Al';
   MasseAtomique:24.3055;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:14;
   Periode:3;
   ElementName:'Silicium';
   NumeroAtomique:14;
   SymboleChimique:'Si';
   MasseAtomique:28.085;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:15;
   Periode:3;
   ElementName:'Phosphore';
   NumeroAtomique:15;
   SymboleChimique:'P';
   MasseAtomique:30.9737620;
   Categorie:AutresNonMetaux;
   Etat:Solide),
  (Groupe:16;
   Periode:3;
   ElementName:'Soufre';
   NumeroAtomique:16;
   SymboleChimique:'S';
   MasseAtomique:32.0675;
   Categorie:AutresNonMetaux;
   Etat:Solide),
  (Groupe:17;
   Periode:3;
   ElementName:'Chlore';
   NumeroAtomique:17;
   SymboleChimique:'Cl';
   MasseAtomique:35.4515;
   Categorie:Halogenes;
   Etat:Gaz),
  (Groupe:18;
   Periode:3;
   ElementName:'Argon';
   NumeroAtomique:18;
   SymboleChimique:'Ar';
   MasseAtomique:39.948;
   Categorie:GazNobles;
   Etat:Gaz),
  (Groupe:1;
   Periode:4;
   ElementName:'Potassium';
   NumeroAtomique:19;
   SymboleChimique:'K';
   MasseAtomique:39.0983;
   Categorie:Alcalins;
   Etat:Solide),
  (Groupe:2;
   Periode:4;
   ElementName:'Calcium';
   NumeroAtomique:20;
   SymboleChimique:'Ca';
   MasseAtomique:40.078;
   Categorie:Alcalinoterreux;
   Etat:Solide),
  (Groupe:3;
   Periode:4;
   ElementName:'Scandium';
   NumeroAtomique:21;
   SymboleChimique:'Sc';
   MasseAtomique:44.955908;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:4;
   Periode:4;
   ElementName:'Titane';
   NumeroAtomique:22;
   SymboleChimique:'Ti';
   MasseAtomique:47.867;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:5;
   Periode:4;
   ElementName:'Vanadium';
   NumeroAtomique:23;
   SymboleChimique:'V';
   MasseAtomique:50.9415;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:6;
   Periode:4;
   ElementName:'Chrome';
   NumeroAtomique:24;
   SymboleChimique:'Cr';
   MasseAtomique:51.9961;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:7;
   Periode:4;
   ElementName:'Manganäse';
   NumeroAtomique:25;
   SymboleChimique:'Mn';
   MasseAtomique:54.938044;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:8;
   Periode:4;
   ElementName:'Fer';
   NumeroAtomique:26;
   SymboleChimique:'Fe';
   MasseAtomique:55.845;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:9;
   Periode:4;
   ElementName:'Cobalt';
   NumeroAtomique:27;
   SymboleChimique:'Co';
   MasseAtomique:58.933194;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:10;
   Periode:4;
   ElementName:'Nickel';
   NumeroAtomique:28;
   SymboleChimique:'Ni';
   MasseAtomique:58.6934;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:11;
   Periode:4;
   ElementName:'Cuivre';
   NumeroAtomique:29;
   SymboleChimique:'Cu';
   MasseAtomique:63.546;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:12;
   Periode:4;
   ElementName:'Zinc';
   NumeroAtomique:30;
   SymboleChimique:'Zn';
   MasseAtomique:63.38;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:13;
   Periode:4;
   ElementName:'Gallium';
   NumeroAtomique:31;
   SymboleChimique:'Ga';
   MasseAtomique:69.723;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:14;
   Periode:4;
   ElementName:'Germanium';
   NumeroAtomique:32;
   SymboleChimique:'Ge';
   MasseAtomique:72.630;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:15;
   Periode:4;
   ElementName:'Arsenic';
   NumeroAtomique:33;
   SymboleChimique:'As';
   MasseAtomique:74.921595;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:16;
   Periode:4;
   ElementName:'SÇlÇmium';
   NumeroAtomique:34;
   SymboleChimique:'Se';
   MasseAtomique:78.971;
   Categorie:AutresNonMetaux;
   Etat:Solide),
  (Groupe:17;
   Periode:4;
   ElementName:'Brome';
   NumeroAtomique:35;
   SymboleChimique:'Br';
   MasseAtomique:79.904;
   Categorie:Halogenes;
   Etat:Liquide),
  (Groupe:18;
   Periode:4;
   ElementName:'Krypton';
   NumeroAtomique:36;
   SymboleChimique:'Kr';
   MasseAtomique:83.798;
   Categorie:GazNobles;
   Etat:Gaz),
  (Groupe:1;
   Periode:5;
   ElementName:'Rubidium';
   NumeroAtomique:37;
   SymboleChimique:'Rb';
   MasseAtomique:85.4678;
   Categorie:Alcalins;
   Etat:Solide),
  (Groupe:2;
   Periode:5;
   ElementName:'Strontium';
   NumeroAtomique:38;
   SymboleChimique:'Sr';
   MasseAtomique:87.62;
   Categorie:Alcalinoterreux;
   Etat:Solide),
  (Groupe:3;
   Periode:5;
   ElementName:'Yttrium';
   NumeroAtomique:39;
   SymboleChimique:'Y';
   MasseAtomique:88.90584;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:4;
   Periode:5;
   ElementName:'Zirconium';
   NumeroAtomique:40;
   SymboleChimique:'Zr';
   MasseAtomique:91.224;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:5;
   Periode:5;
   ElementName:'Niobium';
   NumeroAtomique:41;
   SymboleChimique:'Nb';
   MasseAtomique:92.90637;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:6;
   Periode:5;
   ElementName:'Molybdäne';
   NumeroAtomique:42;
   SymboleChimique:'Mo';
   MasseAtomique:95.95;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:7;
   Periode:5;
   ElementName:'TechnÇtium';
   NumeroAtomique:43;
   SymboleChimique:'Tc';
   MasseAtomique:98.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:8;
   Periode:5;
   ElementName:'RuthÇnium';
   NumeroAtomique:44;
   SymboleChimique:'Ru';
   MasseAtomique:101.07;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:9;
   Periode:5;
   ElementName:'Rhodium';
   NumeroAtomique:45;
   SymboleChimique:'Rh';
   MasseAtomique:102.90550;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:10;
   Periode:5;
   ElementName:'Palladium';
   NumeroAtomique:46;
   SymboleChimique:'Pd';
   MasseAtomique:106.42;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:11;
   Periode:5;
   ElementName:'Argent';
   NumeroAtomique:47;
   SymboleChimique:'Ag';
   MasseAtomique:107.8682;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:12;
   Periode:5;
   ElementName:'Cadmium';
   NumeroAtomique:48;
   SymboleChimique:'Cd';
   MasseAtomique:112.414;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:13;
   Periode:5;
   ElementName:'Indium';
   NumeroAtomique:49;
   SymboleChimique:'In';
   MasseAtomique:114.818;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:14;
   Periode:5;
   ElementName:'êtain';
   NumeroAtomique:50;
   SymboleChimique:'Sn';
   MasseAtomique:118.710;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:15;
   Periode:5;
   ElementName:'Antimoine';
   NumeroAtomique:51;
   SymboleChimique:'Sb';
   MasseAtomique:121.760;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:16;
   Periode:5;
   ElementName:'Tellure';
   NumeroAtomique:52;
   SymboleChimique:'Te';
   MasseAtomique:127.60;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:17;
   Periode:5;
   ElementName:'Iode';
   NumeroAtomique:53;
   SymboleChimique:'I';
   MasseAtomique:126.90447;
   Categorie:Halogenes;
   Etat:Solide),
  (Groupe:18;
   Periode:5;
   ElementName:'XÇnon';
   NumeroAtomique:54;
   SymboleChimique:'Xe';
   MasseAtomique:131.293;
   Categorie:GazNobles;
   Etat:Gaz),
  (Groupe:1;
   Periode:6;
   ElementName:'CÇsium';
   NumeroAtomique:55;
   SymboleChimique:'Cs';
   MasseAtomique:132.905452;
   Categorie:Alcalins;
   Etat:Solide),
  (Groupe:2;
   Periode:6;
   ElementName:'Baryum';
   NumeroAtomique:56;
   SymboleChimique:'Ba';
   MasseAtomique:137.327;
   Categorie:Alcalinoterreux;
   Etat:Solide),

    { Lanthanides }
  (Groupe:3;
   Periode:6;
   ElementName:'Lanthane';
   NumeroAtomique:57;
   SymboleChimique:'La';
   MasseAtomique:138.90547;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'CÇrium';
   NumeroAtomique:58;
   SymboleChimique:'Ce';
   MasseAtomique:140.116;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'PrasÇodyme';
   NumeroAtomique:59;
   SymboleChimique:'Pr';
   MasseAtomique:140.90766;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'NÇodyme';
   NumeroAtomique:60;
   SymboleChimique:'Nd';
   MasseAtomique:144.242;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'PromÇthium';
   NumeroAtomique:61;
   SymboleChimique:'Pm';
   MasseAtomique:145.0;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Samarium';
   NumeroAtomique:62;
   SymboleChimique:'Sm';
   MasseAtomique:150.36;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Europium';
   NumeroAtomique:63;
   SymboleChimique:'Eu';
   MasseAtomique:151.964;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Gadolinium';
   NumeroAtomique:64;
   SymboleChimique:'Gd';
   MasseAtomique:157.25;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Terbium';
   NumeroAtomique:65;
   SymboleChimique:'Tb';
   MasseAtomique:158.92535;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Dysprosium';
   NumeroAtomique:66;
   SymboleChimique:'Dy';
   MasseAtomique:162.500;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Holmium';
   NumeroAtomique:67;
   SymboleChimique:'Ho';
   MasseAtomique:164.93033;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Erbium';
   NumeroAtomique:68;
   SymboleChimique:'Er';
   MasseAtomique:167.259;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Thulium';
   NumeroAtomique:69;
   SymboleChimique:'Tm';
   MasseAtomique:168.93422;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'Ytterbium';
   NumeroAtomique:70;
   SymboleChimique:'Yb';
   MasseAtomique:173.045;
   Categorie:Lanthanides;
   Etat:Solide),
  (Groupe:3;
   Periode:6;
   ElementName:'LutÇcium';
   NumeroAtomique:71;
   SymboleChimique:'Lu';
   MasseAtomique:174.9668;
   Categorie:Lanthanides;
   Etat:Solide),
    { Fin des Lanthanides }

  (Groupe:4;
   Periode:6;
   ElementName:'Hafnium';
   NumeroAtomique:72;
   SymboleChimique:'Hf';
   MasseAtomique:178.49;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:5;
   Periode:6;
   ElementName:'Tantale';
   NumeroAtomique:73;
   SymboleChimique:'Ta';
   MasseAtomique:180.94788;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:6;
   Periode:6;
   ElementName:'Tungstäne';
   NumeroAtomique:74;
   SymboleChimique:'W';
   MasseAtomique:183.84;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:7;
   Periode:6;
   ElementName:'RhÇnium';
   NumeroAtomique:75;
   SymboleChimique:'Re';
   MasseAtomique:186.207;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:8;
   Periode:6;
   ElementName:'Osmium';
   NumeroAtomique:76;
   SymboleChimique:'Os';
   MasseAtomique:190.23;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:9;
   Periode:6;
   ElementName:'Iridium';
   NumeroAtomique:77;
   SymboleChimique:'Ir';
   MasseAtomique:192.217;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:10;
   Periode:6;
   ElementName:'Platine';
   NumeroAtomique:78;
   SymboleChimique:'Pt';
   MasseAtomique:195.084;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:11;
   Periode:6;
   ElementName:'Or';
   NumeroAtomique:79;
   SymboleChimique:'Au';
   MasseAtomique:196.966569;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:12;
   Periode:6;
   ElementName:'Mercure';
   NumeroAtomique:80;
   SymboleChimique:'Hg';
   MasseAtomique:200.592;
   Categorie:MetauxPauvres;
   Etat:Liquide),
  (Groupe:13;
   Periode:6;
   ElementName:'Thallium';
   NumeroAtomique:81;
   SymboleChimique:'Tl';
   MasseAtomique:204.3835;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:14;
   Periode:6;
   ElementName:'Plomb';
   NumeroAtomique:82;
   SymboleChimique:'Pb';
   MasseAtomique:207.2;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:15;
   Periode:6;
   ElementName:'Bismuth';
   NumeroAtomique:83;
   SymboleChimique:'Bi';
   MasseAtomique:208.98040;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:16;
   Periode:6;
   ElementName:'Polonium';
   NumeroAtomique:84;
   SymboleChimique:'Po';
   MasseAtomique:209.0;
   Categorie:MetauxPauvres;
   Etat:Solide),
  (Groupe:17;
   Periode:6;
   ElementName:'Astate';
   NumeroAtomique:85;
   SymboleChimique:'At';
   MasseAtomique:210.0;
   Categorie:Metalloides;
   Etat:Solide),
  (Groupe:18;
   Periode:6;
   ElementName:'Radon';
   NumeroAtomique:86;
   SymboleChimique:'Rn';
   MasseAtomique:222.0;
   Categorie:GazNobles;
   Etat:Gaz),
  (Groupe:1;
   Periode:7;
   ElementName:'Francium';
   NumeroAtomique:87;
   SymboleChimique:'Fr';
   MasseAtomique:223.0;
   Categorie:Alcalins;
   Etat:Solide),
  (Groupe:2;
   Periode:7;
   ElementName:'Radium';
   NumeroAtomique:88;
   SymboleChimique:'Ra';
   MasseAtomique:226.0;
   Categorie:Alcalinoterreux;
   Etat:Solide),

    { Actinides }
  (Groupe:3;
   Periode:7;
   ElementName:'Actinium';
   NumeroAtomique:89;
   SymboleChimique:'Ac';
   MasseAtomique:227.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Thorium';
   NumeroAtomique:90;
   SymboleChimique:'Th';
   MasseAtomique:232.0377;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Protactinium';
   NumeroAtomique:91;
   SymboleChimique:'Pa';
   MasseAtomique:231.03588;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Uranium';
   NumeroAtomique:92;
   SymboleChimique:'U';
   MasseAtomique:238.02891;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Neptunium';
   NumeroAtomique:93;
   SymboleChimique:'Np';
   MasseAtomique:237.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Plutonium';
   NumeroAtomique:94;
   SymboleChimique:'Pu';
   MasseAtomique:244.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'AmÇricium';
   NumeroAtomique:95;
   SymboleChimique:'Am';
   MasseAtomique:243.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Carium';
   NumeroAtomique:96;
   SymboleChimique:'Cm';
   MasseAtomique:247.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'BerkÇlium';
   NumeroAtomique:97;
   SymboleChimique:'Bk';
   MasseAtomique:247.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Californium';
   NumeroAtomique:98;
   SymboleChimique:'Cf';
   MasseAtomique:251.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Einsteinium';
   NumeroAtomique:99;
   SymboleChimique:'Es';
   MasseAtomique:252.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Fermium';
   NumeroAtomique:100;
   SymboleChimique:'Fm';
   MasseAtomique:257.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'MendÇlÇvium';
   NumeroAtomique:101;
   SymboleChimique:'Md';
   MasseAtomique:258.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'NobÇlium';
   NumeroAtomique:102;
   SymboleChimique:'No';
   MasseAtomique:259.0;
   Categorie:Actinides;
   Etat:Solide),
  (Groupe:3;
   Periode:7;
   ElementName:'Lawrencium';
   NumeroAtomique:103;
   SymboleChimique:'Lr';
   MasseAtomique:266.0;
   Categorie:Actinides;
   Etat:Solide),
    { Fin du Actinides }

  (Groupe:4;
   Periode:7;
   ElementName:'Rutherfordium';
   NumeroAtomique:104;
   SymboleChimique:'Rf';
   MasseAtomique:267.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:5;
   Periode:7;
   ElementName:'Dubnium';
   NumeroAtomique:105;
   SymboleChimique:'Db';
   MasseAtomique:268.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:6;
   Periode:7;
   ElementName:'Seaborgium';
   NumeroAtomique:106;
   SymboleChimique:'Sg';
   MasseAtomique:269.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:7;
   Periode:7;
   ElementName:'Bohrium';
   NumeroAtomique:107;
   SymboleChimique:'Bh';
   MasseAtomique:270.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:8;
   Periode:7;
   ElementName:'Hassium';
   NumeroAtomique:108;
   SymboleChimique:'Hs';
   MasseAtomique:277.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:9;
   Periode:7;
   ElementName:'MeitnÇrium';
   NumeroAtomique:109;
   SymboleChimique:'Mt';
   MasseAtomique:278.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:10;
   Periode:7;
   ElementName:'Darmstaditium';
   NumeroAtomique:110;
   SymboleChimique:'Ds';
   MasseAtomique:281.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:11;
   Periode:7;
   ElementName:'Roentgenium';
   NumeroAtomique:111;
   SymboleChimique:'Rg';
   MasseAtomique:282.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:12;
   Periode:7;
   ElementName:'Copernicium';
   NumeroAtomique:112;
   SymboleChimique:'Cn';
   MasseAtomique:285.0;
   Categorie:MetauxDeTransition;
   Etat:Solide),
  (Groupe:13;
   Periode:7;
   ElementName:'Nihonium';
   NumeroAtomique:113;
   SymboleChimique:'Nh';
   MasseAtomique:286.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:14;
   Periode:7;
   ElementName:'FlÇrovium';
   NumeroAtomique:114;
   SymboleChimique:'Fl';
   MasseAtomique:289.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:15;
   Periode:7;
   ElementName:'Moscovium';
   NumeroAtomique:115;
   SymboleChimique:'Mc';
   MasseAtomique:289.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:16;
   Periode:7;
   ElementName:'Livermorium';
   NumeroAtomique:116;
   SymboleChimique:'Lv';
   MasseAtomique:293.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:17;
   Periode:7;
   ElementName:'Tennesse';
   NumeroAtomique:117;
   SymboleChimique:'Ts';
   MasseAtomique:294.0;
   Categorie:NonClasses;
   Etat:Solide),
  (Groupe:18;
   Periode:7;
   ElementName:'Oganesson';
   NumeroAtomique:118;
   SymboleChimique:'Og';
   MasseAtomique:294.0;
   Categorie:NonClasses;
   Etat:Solide)

 );

Var
 Groupe,Periode,I,J:Integer;
 Found:Boolean;
 S:String;

Function StrToUpper(S:String):String;
Var
 I:Byte;
Begin
 For I:=1 to Length(S)do Begin
  If S[I] in['a'..'z']Then S[I]:=Chr(Ord(S[I])-32);
 End;
 StrToUpper:=S;
End;

Function PadRight(S:String;Space:Byte):String;
Var
 I:Byte;
Begin
 If Length(S)<Space Then For I:=Length(S)+1 to Space do S:=S+' ';
 PadRight:=S;
End;

Procedure SetColorElement(NumeroAtomique:Byte);
Var
 J:Integer;
Begin
 TextColor(7);
 TextBackground(0);
 For J:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
  If(ElementsChimiques[J].NumeroAtomique=NumeroAtomique)Then Begin
   Case ElementsChimiques[J].Etat of
    Solide:TextColor(Black);
    Gaz:TextColor(LightRed);
    Liquide:TextColor(LightBlue);
   End;
   Case ElementsChimiques[J].Categorie of
    NonClasses:TextBackground(White);
    Alcalins:TextBackground(Red);
    Alcalinoterreux:TextBackground(Brown);
    Lanthanides:TextBackground(Magenta);
    Actinides:TextBackground(Magenta);
    MetauxDeTransition:TextBackground(Magenta);
    MetauxPauvres:TextBackground(LightGray);
    Metalloides:TextBackground(Brown);
    AutresNonMetaux:TextBackground(Green);
    Halogenes:TextBackground(Yellow);
    GazNobles:TextBackground(Cyan);
   End;
   Exit;
  End;
 End;
End;

BEGIN
 If(ParamStr(1)='/?')or(ParamStr(1)='--help')or(ParamStr(1)='-h')or
   (ParamStr(1)='/h')or(ParamStr(1)='/H')Then Begin
  WriteLn('ELEMENTS : Cette commande permet d''afficher ',
                      'le tableau periodique des elements chimiques.');
  WriteLn;
  WriteLn('Syntaxe : ELEMENTS');
  WriteLn('Syntaxe : ELEMENTS nomelement');
 End
  Else
 If ParamCount>0Then Begin
  For I:=1 to ParamCount do Begin
   Found:=False;
   For J:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If StrToUpper(ElementsChimiques[J].ElementName)=StrToUpper(ParamStr(I))Then Begin
     WriteLn('Nom : ',ElementsChimiques[J].ElementName);
     WriteLn('Groupe : ',ElementsChimiques[J].Groupe);
     WriteLn('Periode : ',ElementsChimiques[J].Periode);
     WriteLn('Numero atomique : ',ElementsChimiques[J].NumeroAtomique);
     WriteLn('Symbole chimique : ',ElementsChimiques[J].SymboleChimique);
     WriteLn('Masse atomique : ',ElementsChimiques[J].MasseAtomique:1:4);
     Write('Categorie : ');
     Case ElementsChimiques[J].Categorie of
      NonClasses:WriteLn('Non classes');
      Alcalins:WriteLn('Alcalins');
      Alcalinoterreux:WriteLn('Alcalinoterreux');
      Lanthanides:WriteLn('Lanthanides');
      Actinides:WriteLn('Actinides');
      MetauxDeTransition:WriteLn('Metaux de transition');
      MetauxPauvres:WriteLn('Metaux pauvres');
      Metalloides:WriteLn('Metalloides');
      AutresNonMetaux:WriteLn('Autres non metaux');
      Halogenes:WriteLn('Halogenes');
      GazNobles:WriteLn('Gaz nobles');
     End;
     Write('Etat a 0 Celcius et 101,3 kPa : ');
     Case ElementsChimiques[J].Etat of
      Gaz:WriteLn('Gaz');
      Liquide:WriteLn('Liquide');
      Solide:WriteLn('Solide');
     End;
     Found:=True;
     Break;
    End;
   End;
   If Not(Found)Then Begin
    WriteLn('Element non reconnu');
    Halt;
   End;
  End;
 End
  Else
 Begin
  For Periode:=1 to 7 do Begin
    { Nom }
   For Groupe:=1 to 18 do Begin
    Found:=False;
    SetColorElement(0);
    For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
     If(ElementsChimiques[I].Periode=Periode)and
       (ElementsChimiques[I].Groupe=Groupe)Then Begin
      SetColorElement(ElementsChimiques[I].NumeroAtomique);
      If(ElementsChimiques[I].Periode=6)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(Copy('Lanthanides',1,4));
      End
       Else
      If(ElementsChimiques[I].Periode=7)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(Copy('Actinides',1,4));
      End
       Else
      Begin
       Write(PadRight(Copy(ElementsChimiques[I].ElementName,1,4),4));
      End;
      Found:=True;
      Break;
     End;
    End;
    If Not(Found)Then Write(' ':4);
   End;
   SetColorElement(0);
   WriteLn;
    { NumÇro atomique }
   For Groupe:=1 to 18 do Begin
    Found:=False;
    SetColorElement(0);
    For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
     If(ElementsChimiques[I].Periode=Periode)and
       (ElementsChimiques[I].Groupe=Groupe)Then Begin
      SetColorElement(ElementsChimiques[I].NumeroAtomique);
      If(ElementsChimiques[I].Periode=6)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write('57':4);
      End
       Else
      If(ElementsChimiques[I].Periode=7)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(Copy('89 Ö ',1,4));
      End
       Else
      Begin
       Write(' ',ElementsChimiques[I].NumeroAtomique:3);
      End;
      Found:=True;
      Break;
     End;
    End;
    If Not(Found)Then Write(' ':4);
   End;
   SetColorElement(0);
   WriteLn;
    { Symbole chimique }
   For Groupe:=1 to 18 do Begin
    Found:=False;
    SetColorElement(0);
    For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
     If(ElementsChimiques[I].Periode=Periode)and
       (ElementsChimiques[I].Groupe=Groupe)Then Begin
      SetColorElement(ElementsChimiques[I].NumeroAtomique);
      If(ElementsChimiques[I].Periode=6)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(Copy('Ö 71',1,4));
      End
       Else
     If(ElementsChimiques[I].Periode=7)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(PadRight(Copy('103',1,4),4));
      End
       Else
      Begin
       Write(PadRight(' '+ElementsChimiques[I].SymboleChimique,4));
      End;
      Found:=True;
      Break;
     End;
    End;
    If Not(Found)Then Write(' ':4);
   End;
   SetColorElement(0);
   WriteLn;
    { Masse atomique }
   For Groupe:=1 to 18 do Begin
    Found:=False;
    SetColorElement(0);
    For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
     If(ElementsChimiques[I].Periode=Periode)and
       (ElementsChimiques[I].Groupe=Groupe)Then Begin
      SetColorElement(ElementsChimiques[I].NumeroAtomique);
       If(ElementsChimiques[I].Periode=6)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(' ':4);
      End
       Else
      If(ElementsChimiques[I].Periode=7)and
         (ElementsChimiques[I].Groupe=3)Then Begin
       Write(' ':4);
      End
       Else
      Begin
       Str(ElementsChimiques[I].MasseAtomique:1:2,S);
       Write(Copy(S,1,4));
      End;
      Found:=True;
      Break;
     End;
    End;
    If Not(Found)Then Write(' ':4);
   End;
   SetColorElement(0);
   WriteLn;
  End;
  WriteLn;
   { Nom }
  For J:=57 to 71 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Write(PadRight(Copy(ElementsChimiques[I].ElementName,1,4),4));
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
   { NumÇro atomique }
  For J:=57 to 71 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Write(' ',ElementsChimiques[I].NumeroAtomique:3);
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
   { Symbole chimique }
  For J:=57 to 71 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Write(PadRight(' '+ElementsChimiques[I].SymboleChimique,4));
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
   { Masse atomique }
  For J:=57 to 71 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Str(ElementsChimiques[I].MasseAtomique:1:2,S);
     Write(Copy(S,1,4));
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
  WriteLn;
   { Nom }
  For J:=89 to 103 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Write(PadRight(Copy(ElementsChimiques[I].ElementName,1,4),4));
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
   { NumÇro atomique }
  For J:=89 to 103 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Write(' ',ElementsChimiques[I].NumeroAtomique:3);
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
   { Symbole chimique }
  For J:=89 to 103 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Write(PadRight(' '+ElementsChimiques[I].SymboleChimique,4));
     Break;
    End;
   End;
  End;
  SetColorElement(0);
  WriteLn;
   { Masse atomique }
  For J:=89 to 103 do Begin
   For I:=Low(ElementsChimiques)to High(ElementsChimiques)do Begin
    If ElementsChimiques[I].NumeroAtomique=J Then Begin
     SetColorElement(J);
     Str(ElementsChimiques[I].MasseAtomique:1:2,S);
     Write(Copy(S,1,4));
     Break;
    End;
   End;
  End;
 End;
END.