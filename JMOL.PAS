{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2026
  @version: 1.00
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal 7, Free Pascal)
  @description: Clone de l'application JMOL - visualiseur moleculaire 3D
                (PDB, XYZ, MOL, SDF, MOL2, CIF) avec menus en francais.

  Ligne de commande :
    JMOL [fichier.PDB|.XYZ|.MOL|.SDF|.MOL2|.CIF|.SPT]

  Fonctionnalites :
    - Barre de titre, menus deroulants, barre d'outils, barre d'etat
    - Menu contextuel (bouton droit de la souris) avec sous-menus
    - Styles : remplissage CPK, boules et batonnets, batonnets, fil de fer,
      points de surface, squelette, trace, rubans, cartoon, fusees, brins
    - Couleurs : element, chaine, residu, structure, arc-en-ciel, temperature
    - Selection (expressions a la Jmol), etiquettes, axes, boite englobante,
      cellule unitaire, liaisons hydrogene, halos de selection
    - Mesures de distances, d'angles et de diedres
    - Modeles multiples (animation) et modes de vibration (panneau Vibration)
    - Console de commandes (langage de script de type Jmol)
    - Enregistrement PDB / XYZ / MOL / image BMP

  Inspirations :
    - Menus / dialogues : CLONEDRAW\PAINT.PAS
    - Chargeurs de fichiers moleculaires et rendu : PYMOL.PAS
}

Program JMOL;

{$IFDEF FPC}
 {$MODE TP}
 {$Q-}
 {$R-}
{$ENDIF}

Uses {$IFDEF FPC}
      DOS,Crt,PtcGraph,PtcCrt,PtcMouse
     {$ELSE}
      DOS,Crt,Graph
     {$ENDIF};

{ --------------------------------------------------------------------- }
{ Constantes generales                                                  }
{ --------------------------------------------------------------------- }

Const
 PROG_NAME     = 'Jmol';
 PROG_TITLE    = 'Jmol Clone 1.00 - QChimie';
 PROG_VERSION  = '1.00';

 {$IFDEF FPC}
 MAX_ATOMS     = 8192;
 MAX_BONDS     = 16384;
 MAX_RES       = 4096;
 MAX_PRIM      = 32000;
 MAX_QUAD      = 30000;
 MAX_HB        = 4096;
 MAX_MODELS    = 400;
 MAX_RIB       = 12000;
 MAX_PAIRS     = 30000;
 {$ELSE}
 MAX_ATOMS     = 700;
 MAX_BONDS     = 1400;
 MAX_RES       = 400;
 MAX_PRIM      = 4000;
 MAX_QUAD      = 1800;
 MAX_HB        = 400;
 MAX_MODELS    = 80;
 MAX_RIB       = 1200;
 MAX_PAIRS     = 2500;
 {$ENDIF}
 MAX_MEAS      = 40;
 {$IFDEF FPC}
 MAX_TEXT      = 200;      { lignes d'une fenetre de texte }
 MAX_CONLN     = 100;      { lignes du tampon de la console }
 FD_MAX        = 300;      { entrees de la boite d'ouverture }
 {$ELSE}
 MAX_TEXT      = 100;
 MAX_CONLN     = 50;
 FD_MAX        = 150;
 {$ENDIF}
 MAX_LBLFMT    = 12;
 MAX_PM        = 56;       { nombre maximal de menus }
 MAX_PMI       = 26;       { nombre maximal d'elements par menu }
 MAX_DOTS      = 112;      { points par atome pour la surface }
 DOT_BYTES     = 14;       { MAX_DOTS / 8 }

 { Disposition de la fenetre }
 TITLE_H       = 18;
 MENU_H        = 18;
 TOOL_H        = 28;
 STAT_H        = 18;
 CONS_H        = 112;

 { Types d'atome }
 ATOMK_HETERO  = 0;
 ATOMK_PROT    = 1;
 ATOMK_NUCL    = 2;
 ATOMK_WATER   = 3;
 ATOMK_LIGAND  = 4;

 { Styles de structure secondaire }
 RIB_NONE      = 0;
 RIB_BACKBONE  = 1;
 RIB_TRACE     = 2;
 RIB_RIBBON    = 3;
 RIB_CARTOON   = 4;
 RIB_ROCKET    = 5;
 RIB_STRAND    = 6;

 { Types de liaison }
 BOND_HBOND    = 9;

 { Types de primitives }
 PRIM_ATOM     = 1;
 PRIM_BOND     = 2;
 PRIM_QUAD     = 3;
 PRIM_TUBE     = 4;
 PRIM_DOTS     = 5;
 PRIM_HBOND    = 6;
 PRIM_VECTOR   = 7;

 { Modes de selection a la souris (Choisir) }
 PICK_IDENT    = 0;
 PICK_DIST     = 1;
 PICK_ANGLE    = 2;
 PICK_TORSION  = 3;
 PICK_SELATOM  = 4;
 PICK_SELRES   = 5;
 PICK_SELCHAIN = 6;
 PICK_CENTER   = 7;
 PICK_LABEL    = 8;

 PI_VAL        = 3.14159265358979;
 DEG2RAD       = 0.01745329251994;

{ Codes de commandes (menus, barre d'outils, clavier) }
Const
 CM_NONE        = 0;
 CM_OPEN        = 100;
 CM_RELOAD      = 101;
 CM_SAVEPDB     = 102;
 CM_SAVEXYZ     = 103;
 CM_SAVEMOL     = 104;
 CM_SAVEBMP     = 105;
 CM_RUNSCRIPT   = 106;
 CM_FILEINFO    = 107;
 CM_QUIT        = 108;

 CM_SELALL      = 120;
 CM_SELNONE     = 121;
 CM_SELINV      = 122;
 CM_HIDESEL     = 123;
 CM_SHOWALL     = 124;
 CM_CENTERSEL   = 125;
 CM_CLEARMEAS   = 126;
 CM_CLEARLABELS = 127;
 CM_SELEXPR     = 128;
 CM_DELETESEL   = 129;

 CM_SCH_CPK     = 140;
 CM_SCH_BALL    = 141;
 CM_SCH_STICK   = 142;
 CM_SCH_WIRE    = 143;

 CM_TOG_HYDRO   = 150;
 CM_TOG_HBOND   = 151;
 CM_TOG_HALO    = 152;
 CM_TOG_PERSP   = 153;
 CM_TOG_AXES    = 154;
 CM_TOG_BBOX    = 155;
 CM_TOG_CELL    = 156;
 CM_TOG_VECTOR  = 157;
 CM_TOG_WATER   = 158;
 CM_TOG_CONSOLE = 159;
 CM_TOG_SPIN    = 160;
 CM_TOG_VIB     = 161;
 CM_TOG_ANIM    = 162;

 CM_VIEW_FRONT  = 170;
 CM_VIEW_BACK   = 171;
 CM_VIEW_LEFT   = 172;
 CM_VIEW_RIGHT  = 173;
 CM_VIEW_TOP    = 174;
 CM_VIEW_BOTTOM = 175;
 CM_VIEW_RESET  = 176;
 CM_ZOOM_IN     = 177;
 CM_ZOOM_OUT    = 178;
 CM_ZOOM_FIT    = 179;

 CM_SPIN_X      = 180;
 CM_SPIN_Y      = 181;
 CM_SPIN_Z      = 182;
 CM_SPIN_SLOW   = 183;
 CM_SPIN_MED    = 184;
 CM_SPIN_FAST   = 185;

 CM_FR_FIRST    = 190;
 CM_FR_PREV     = 191;
 CM_FR_NEXT     = 192;
 CM_FR_LAST     = 193;
 CM_FR_ALL      = 194;
 CM_AN_ONCE     = 195;
 CM_AN_LOOP     = 196;
 CM_AN_PALIN    = 197;
 CM_AN_SLOWER   = 198;
 CM_AN_FASTER   = 199;

 CM_VIB_PANEL   = 200;
 CM_VIB_BIGGER  = 201;
 CM_VIB_SMALLER = 202;
 CM_VIB_NEXT    = 203;
 CM_VIB_PREV    = 204;

 CM_MEAS_DIST   = 210;
 CM_MEAS_ANGLE  = 211;
 CM_MEAS_TORS   = 212;
 CM_DEL_DIST    = 213;
 CM_DEL_ANGLE   = 214;
 CM_DEL_TORS    = 215;
 CM_MEAS_LIST   = 216;
 CM_UNIT_ANG    = 217;
 CM_UNIT_NM     = 218;
 CM_UNIT_PM     = 219;

 CM_SURF_SEL    = 230;
 CM_SURF_ALL    = 231;
 CM_SURF_OFF    = 232;
 CM_SURF_DENSE  = 233;
 CM_SURF_LIGHT  = 234;

 CM_SHOW_FILE   = 240;
 CM_SHOW_MODEL  = 241;
 CM_SHOW_ATOMS  = 242;
 CM_SHOW_SEQ    = 243;
 CM_SHOW_CELL   = 244;
 CM_SHOW_ORIENT = 245;
 CM_SHOW_SEL    = 246;

 CM_CALC_FORM   = 250;
 CM_CALC_MASS   = 251;
 CM_CALC_COM    = 252;
 CM_CALC_GYR    = 253;
 CM_CALC_DIM    = 254;
 CM_CALC_MIN    = 255;
 CM_CALC_HBOND  = 256;

 CM_HELP        = 270;
 CM_HELP_MOUSE  = 271;
 CM_HELP_CMDS   = 272;
 CM_ABOUT       = 273;
 CM_LANG_FR     = 274;

 CM_STRUCT_BASE = 300;   { + RIB_xxx }
 CM_PICK_BASE   = 320;   { + PICK_xxx }
 CM_SEL_BASE    = 1000;  { + indice de l'expression de selection }
 CM_COLSCH_BASE = 1050;  { + indice du schema de couleurs }
 CM_COLATOM_BASE= 1100;  { + indice de couleur }
 CM_COLRIB_BASE = 1150;
 CM_COLBG_BASE  = 1200;
 CM_ATOMPCT_BASE= 1250;
 CM_BONDRAD_BASE= 1300;
 CM_ZOOM_BASE   = 1350;
 CM_LABEL_BASE  = 1400;
 CM_MACRO_BASE  = 1450;
 CM_ELEM_BASE   = 1500;  { + numero atomique }
 CM_CHAIN_BASE  = 1700;  { + indice de chaine }
 CM_COLRIBSCH_BASE = 1750;

{ --------------------------------------------------------------------- }
{ Types                                                                 }
{ --------------------------------------------------------------------- }

Type
 TAtom = Record
  X, Y, Z   : Real;          { coordonnees (A) }
  Temp      : Real;          { facteur de temperature }
  Serial    : LongInt;
  Name      : String[4];     { nom de l'atome (CA, N, O...) }
  Res       : String[4];     { nom du residu }
  Chain     : Char;
  ResSeq    : Integer;
  ICode     : Char;
  Model     : Word;          { numero du modele (1..NModels) }
  ElemZ     : Byte;          { numero atomique }
  Kind      : Byte;          { ATOMK_xxx }
  Charge    : ShortInt;
  Sf        : Integer;       { >=0 : % vdW ; <0 : rayon en 1/100 A }
  Col       : Byte;          { 255 = couleur de l'element }
  Lbl       : Byte;          { indice du format d'etiquette (0 = aucune) }
  ResIdx    : Integer;       { indice dans Resi }
  Sel       : Boolean;
  Hid       : Boolean;
  Dots      : Boolean;
 End;
 TAtomArr = Array[0..MAX_ATOMS-1] Of TAtom;
 PAtomArr = ^TAtomArr;

 TBond = Record
  A, B  : Integer;
  Ord   : Byte;            { 1, 2, 3, 4 (aromatique) }
  BR    : Byte;            { 0 = cachee, 1 = fil, n = rayon en 1/100 A }
 End;
 TBondArr = Array[0..MAX_BONDS-1] Of TBond;
 PBondArr = ^TBondArr;

 TResidue = Record
  Chain  : Char;
  Seq    : Integer;
  ICode  : Char;
  Name   : String[4];
  First, Last : Integer;
  CA, CC, OO  : Integer;   { atomes guides : CA (ou P), C, O ; -1 si absent }
  Model  : Word;
  SS     : Char;           { H=helice, E=brin, T=coude, L=boucle }
  Rib    : Byte;           { style de structure (RIB_xxx) }
  RCol   : Byte;           { couleur du ruban (255 = selon la structure) }
  IsAA, IsNA : Boolean;
 End;
 TResArr = Array[0..MAX_RES-1] Of TResidue;
 PResArr = ^TResArr;

 TModelInfo = Record
  Name  : String[40];
  Freq  : Real;
  Vib   : Boolean;
  First : Integer;         { premier et dernier atome du modele }
  Last  : Integer;
 End;
 TModelArr = Array[1..MAX_MODELS] Of TModelInfo;
 PModelArr = ^TModelArr;

 THBond = Record A, B : Integer; End;
 THBArr = Array[0..MAX_HB-1] Of THBond;
 PHBArr = ^THBArr;

 TMeas = Record
  N : Byte;
  A : Array[0..3] Of Integer;
 End;

 TVec3 = Record X, Y, Z : Real; End;
 TBoxPts = Array[0..7] Of TVec3;
 TVecArr = Array[0..MAX_ATOMS-1] Of TVec3;
 PVecArr = ^TVecArr;

 TRealArr = Array[0..MAX_ATOMS-1] Of Real;
 PRealArr = ^TRealArr;
 TIntArr  = Array[0..MAX_ATOMS-1] Of Integer;
 PIntArr  = ^TIntArr;
 TByteArr = Array[0..MAX_ATOMS-1] Of Byte;
 PByteArr = ^TByteArr;
 TRibReal = Array[0..MAX_RIB-1] Of Real;
 PRibReal = ^TRibReal;
 TRibByte = Array[0..MAX_RIB-1] Of Byte;
 PRibByte = ^TRibByte;
 TRibInt  = Array[0..MAX_RIB-1] Of Integer;
 PRibInt  = ^TRibInt;
 TDotMask = Array[0..MAX_ATOMS-1] Of Array[0..DOT_BYTES-1] Of Byte;
 PDotMask = ^TDotMask;

 TPrim = Record
  Key  : LongInt;
  Kind : Byte;
  A, B : Integer;
 End;
 TPrimArr = Array[0..MAX_PRIM-1] Of TPrim;
 PPrimArr = ^TPrimArr;

 TQuad = Record
  X, Y : Array[0..3] Of Integer;
  Col, Dk : Byte;
 End;
 TQuadArr = Array[0..MAX_QUAD-1] Of TQuad;
 PQuadArr = ^TQuadArr;

 TPMItem = Record
  Cap  : String[34];       { 'texte|raccourci' ; '-' = separateur }
  Cmd  : Integer;
  Sub  : ShortInt;         { indice du sous-menu, -1 si aucun }
  Kind : Byte;             { 0 normal, 1 case a cocher, 2 bouton radio }
 End;
 TPMenu = Record
  N     : Integer;
  Title : String[14];
  It    : Array[0..MAX_PMI-1] Of TPMItem;
 End;
 TPMArr = Array[0..MAX_PM-1] Of TPMenu;
 PPMArr = ^TPMArr;

 TElementInfo = Record
  Z       : Byte;
  Sym     : String[2];
  Name    : String[16];
  Color   : Byte;
  Radius  : Real;          { rayon covalent (A) }
  VdW     : Real;          { rayon de van der Waals (A) }
  Mass    : Real;          { masse atomique }
 End;
 TElemDB = Array[0..103] Of TElementInfo;

 TTokens = Array[0..23] Of String[60];

 TLvl = Record
  Menu   : Integer;
  X, Y, W, H : Integer;
  Sel    : Integer;
  Buf    : Pointer;
  BufSz  : LongInt;
 End;

{ --------------------------------------------------------------------- }
{ Table d'elements chimiques                                            }
{ --------------------------------------------------------------------- }

Const
 ElemDB : TElemDB = (
 (Z:0; Sym:'?'; Name:'Inconnu'; Color:LightGray; Radius:0.70; VdW:1.50; Mass:0.0),
 (Z:1; Sym:'H'; Name:'Hydrogene'; Color:White; Radius:0.31; VdW:1.20; Mass:1.008),
 (Z:2; Sym:'He'; Name:'Helium'; Color:LightCyan; Radius:0.28; VdW:1.40; Mass:4.0026),
 (Z:3; Sym:'Li'; Name:'Lithium'; Color:LightMagenta; Radius:1.28; VdW:1.82; Mass:6.94),
 (Z:4; Sym:'Be'; Name:'Beryllium'; Color:LightGreen; Radius:0.96; VdW:1.53; Mass:9.0122),
 (Z:5; Sym:'B'; Name:'Bore'; Color:Brown; Radius:0.84; VdW:1.92; Mass:10.81),
 (Z:6; Sym:'C'; Name:'Carbone'; Color:LightGray; Radius:0.76; VdW:1.70; Mass:12.011),
 (Z:7; Sym:'N'; Name:'Azote'; Color:LightBlue; Radius:0.71; VdW:1.55; Mass:14.007),
 (Z:8; Sym:'O'; Name:'Oxygene'; Color:LightRed; Radius:0.66; VdW:1.52; Mass:15.999),
 (Z:9; Sym:'F'; Name:'Fluor'; Color:LightGreen; Radius:0.57; VdW:1.47; Mass:18.998),
 (Z:10; Sym:'Ne'; Name:'Neon'; Color:LightCyan; Radius:0.58; VdW:1.54; Mass:20.180),
 (Z:11; Sym:'Na'; Name:'Sodium'; Color:LightMagenta; Radius:1.66; VdW:2.27; Mass:22.990),
 (Z:12; Sym:'Mg'; Name:'Magnesium'; Color:Green; Radius:1.41; VdW:1.73; Mass:24.305),
 (Z:13; Sym:'Al'; Name:'Aluminium'; Color:LightGray; Radius:1.21; VdW:1.84; Mass:26.982),
 (Z:14; Sym:'Si'; Name:'Silicium'; Color:Brown; Radius:1.11; VdW:2.10; Mass:28.085),
 (Z:15; Sym:'P'; Name:'Phosphore'; Color:Brown; Radius:1.07; VdW:1.80; Mass:30.974),
 (Z:16; Sym:'S'; Name:'Soufre'; Color:Yellow; Radius:1.05; VdW:1.80; Mass:32.06),
 (Z:17; Sym:'Cl'; Name:'Chlore'; Color:Green; Radius:1.02; VdW:1.75; Mass:35.45),
 (Z:18; Sym:'Ar'; Name:'Argon'; Color:LightCyan; Radius:1.06; VdW:1.88; Mass:39.948),
 (Z:19; Sym:'K'; Name:'Potassium'; Color:LightMagenta; Radius:2.03; VdW:2.75; Mass:39.098),
 (Z:20; Sym:'Ca'; Name:'Calcium'; Color:LightGreen; Radius:1.76; VdW:2.31; Mass:40.078),
 (Z:21; Sym:'Sc'; Name:'Scandium'; Color:LightGray; Radius:1.70; VdW:2.11; Mass:44.956),
 (Z:22; Sym:'Ti'; Name:'Titane'; Color:LightGray; Radius:1.60; VdW:2.00; Mass:47.867),
 (Z:23; Sym:'V'; Name:'Vanadium'; Color:LightGray; Radius:1.53; VdW:2.00; Mass:50.942),
 (Z:24; Sym:'Cr'; Name:'Chrome'; Color:LightGray; Radius:1.39; VdW:2.00; Mass:51.996),
 (Z:25; Sym:'Mn'; Name:'Manganese'; Color:Magenta; Radius:1.39; VdW:2.00; Mass:54.938),
 (Z:26; Sym:'Fe'; Name:'Fer'; Color:Brown; Radius:1.32; VdW:2.00; Mass:55.845),
 (Z:27; Sym:'Co'; Name:'Cobalt'; Color:LightMagenta; Radius:1.26; VdW:2.00; Mass:58.933),
 (Z:28; Sym:'Ni'; Name:'Nickel'; Color:Green; Radius:1.24; VdW:1.63; Mass:58.693),
 (Z:29; Sym:'Cu'; Name:'Cuivre'; Color:Brown; Radius:1.32; VdW:1.40; Mass:63.546),
 (Z:30; Sym:'Zn'; Name:'Zinc'; Color:DarkGray; Radius:1.22; VdW:1.39; Mass:65.38),
 (Z:31; Sym:'Ga'; Name:'Gallium'; Color:Brown; Radius:1.22; VdW:1.87; Mass:69.723),
 (Z:32; Sym:'Ge'; Name:'Germanium'; Color:LightGray; Radius:1.20; VdW:2.11; Mass:72.630),
 (Z:33; Sym:'As'; Name:'Arsenic'; Color:Magenta; Radius:1.19; VdW:1.85; Mass:74.922),
 (Z:34; Sym:'Se'; Name:'Selenium'; Color:Brown; Radius:1.20; VdW:1.90; Mass:78.971),
 (Z:35; Sym:'Br'; Name:'Brome'; Color:Brown; Radius:1.20; VdW:1.85; Mass:79.904),
 (Z:36; Sym:'Kr'; Name:'Krypton'; Color:LightCyan; Radius:1.16; VdW:2.02; Mass:83.798),
 (Z:37; Sym:'Rb'; Name:'Rubidium'; Color:Magenta; Radius:2.20; VdW:3.03; Mass:85.468),
 (Z:38; Sym:'Sr'; Name:'Strontium'; Color:LightGreen; Radius:1.95; VdW:2.49; Mass:87.62),
 (Z:39; Sym:'Y'; Name:'Yttrium'; Color:LightCyan; Radius:1.90; VdW:2.00; Mass:88.906),
 (Z:40; Sym:'Zr'; Name:'Zirconium'; Color:LightCyan; Radius:1.75; VdW:2.00; Mass:91.224),
 (Z:41; Sym:'Nb'; Name:'Niobium'; Color:LightCyan; Radius:1.64; VdW:2.00; Mass:92.906),
 (Z:42; Sym:'Mo'; Name:'Molybdene'; Color:LightCyan; Radius:1.54; VdW:2.00; Mass:95.95),
 (Z:43; Sym:'Tc'; Name:'Technetium'; Color:LightCyan; Radius:1.47; VdW:2.00; Mass:98),
 (Z:44; Sym:'Ru'; Name:'Ruthenium'; Color:LightCyan; Radius:1.46; VdW:2.00; Mass:101.07),
 (Z:45; Sym:'Rh'; Name:'Rhodium'; Color:LightCyan; Radius:1.42; VdW:2.00; Mass:102.91),
 (Z:46; Sym:'Pd'; Name:'Palladium'; Color:LightGray; Radius:1.39; VdW:1.63; Mass:106.42),
 (Z:47; Sym:'Ag'; Name:'Argent'; Color:LightGray; Radius:1.45; VdW:1.72; Mass:107.87),
 (Z:48; Sym:'Cd'; Name:'Cadmium'; Color:LightGray; Radius:1.44; VdW:1.58; Mass:112.41),
 (Z:49; Sym:'In'; Name:'Indium'; Color:LightGray; Radius:1.42; VdW:1.93; Mass:114.82),
 (Z:50; Sym:'Sn'; Name:'Etain'; Color:LightGray; Radius:1.39; VdW:2.17; Mass:118.71),
 (Z:51; Sym:'Sb'; Name:'Antimoine'; Color:LightGray; Radius:1.39; VdW:2.06; Mass:121.76),
 (Z:52; Sym:'Te'; Name:'Tellure'; Color:Brown; Radius:1.38; VdW:2.06; Mass:127.60),
 (Z:53; Sym:'I'; Name:'Iode'; Color:Magenta; Radius:1.39; VdW:1.98; Mass:126.90),
 (Z:54; Sym:'Xe'; Name:'Xenon'; Color:LightCyan; Radius:1.40; VdW:2.16; Mass:131.29),
 (Z:55; Sym:'Cs'; Name:'Cesium'; Color:Magenta; Radius:2.44; VdW:3.43; Mass:132.91),
 (Z:56; Sym:'Ba'; Name:'Baryum'; Color:LightGreen; Radius:2.15; VdW:2.68; Mass:137.33),
 (Z:57; Sym:'La'; Name:'Lanthane'; Color:LightCyan; Radius:2.07; VdW:2.40; Mass:138.91),
 (Z:58; Sym:'Ce'; Name:'Cerium'; Color:LightCyan; Radius:2.04; VdW:2.35; Mass:140.12),
 (Z:59; Sym:'Pr'; Name:'Praseodyme'; Color:LightCyan; Radius:2.03; VdW:2.39; Mass:140.91),
 (Z:60; Sym:'Nd'; Name:'Neodyme'; Color:LightCyan; Radius:2.01; VdW:2.29; Mass:144.24),
 (Z:61; Sym:'Pm'; Name:'Promethium'; Color:LightCyan; Radius:1.99; VdW:2.36; Mass:145),
 (Z:62; Sym:'Sm'; Name:'Samarium'; Color:LightCyan; Radius:1.98; VdW:2.29; Mass:150.36),
 (Z:63; Sym:'Eu'; Name:'Europium'; Color:LightCyan; Radius:1.98; VdW:2.33; Mass:151.96),
 (Z:64; Sym:'Gd'; Name:'Gadolinium'; Color:LightCyan; Radius:1.96; VdW:2.37; Mass:157.25),
 (Z:65; Sym:'Tb'; Name:'Terbium'; Color:LightCyan; Radius:1.94; VdW:2.21; Mass:158.93),
 (Z:66; Sym:'Dy'; Name:'Dysprosium'; Color:LightCyan; Radius:1.92; VdW:2.29; Mass:162.50),
 (Z:67; Sym:'Ho'; Name:'Holmium'; Color:LightCyan; Radius:1.92; VdW:2.16; Mass:164.93),
 (Z:68; Sym:'Er'; Name:'Erbium'; Color:LightCyan; Radius:1.89; VdW:2.35; Mass:167.26),
 (Z:69; Sym:'Tm'; Name:'Thulium'; Color:LightCyan; Radius:1.90; VdW:2.27; Mass:168.93),
 (Z:70; Sym:'Yb'; Name:'Ytterbium'; Color:LightCyan; Radius:1.87; VdW:2.42; Mass:173.05),
 (Z:71; Sym:'Lu'; Name:'Lutetium'; Color:LightCyan; Radius:1.87; VdW:2.21; Mass:174.97),
 (Z:72; Sym:'Hf'; Name:'Hafnium'; Color:LightCyan; Radius:1.75; VdW:2.12; Mass:178.49),
 (Z:73; Sym:'Ta'; Name:'Tantale'; Color:LightCyan; Radius:1.70; VdW:2.17; Mass:180.95),
 (Z:74; Sym:'W'; Name:'Tungstene'; Color:LightCyan; Radius:1.62; VdW:2.10; Mass:183.84),
 (Z:75; Sym:'Re'; Name:'Rhenium'; Color:LightCyan; Radius:1.51; VdW:2.17; Mass:186.21),
 (Z:76; Sym:'Os'; Name:'Osmium'; Color:LightCyan; Radius:1.44; VdW:2.16; Mass:190.23),
 (Z:77; Sym:'Ir'; Name:'Iridium'; Color:LightCyan; Radius:1.41; VdW:2.02; Mass:192.22),
 (Z:78; Sym:'Pt'; Name:'Platine'; Color:LightGray; Radius:1.36; VdW:1.75; Mass:195.08),
 (Z:79; Sym:'Au'; Name:'Or'; Color:Yellow; Radius:1.36; VdW:1.66; Mass:196.97),
 (Z:80; Sym:'Hg'; Name:'Mercure'; Color:LightGray; Radius:1.32; VdW:1.55; Mass:200.59),
 (Z:81; Sym:'Tl'; Name:'Thallium'; Color:LightGray; Radius:1.45; VdW:1.96; Mass:204.38),
 (Z:82; Sym:'Pb'; Name:'Plomb'; Color:DarkGray; Radius:1.46; VdW:2.02; Mass:207.2),
 (Z:83; Sym:'Bi'; Name:'Bismuth'; Color:Magenta; Radius:1.48; VdW:2.07; Mass:208.98),
 (Z:84; Sym:'Po'; Name:'Polonium'; Color:Brown; Radius:1.40; VdW:1.97; Mass:209),
 (Z:85; Sym:'At'; Name:'Astate'; Color:Brown; Radius:1.50; VdW:2.02; Mass:210),
 (Z:86; Sym:'Rn'; Name:'Radon'; Color:LightCyan; Radius:1.50; VdW:2.20; Mass:222),
 (Z:87; Sym:'Fr'; Name:'Francium'; Color:Magenta; Radius:2.60; VdW:3.48; Mass:223),
 (Z:88; Sym:'Ra'; Name:'Radium'; Color:LightGreen; Radius:2.21; VdW:2.83; Mass:226),
 (Z:89; Sym:'Ac'; Name:'Actinium'; Color:LightCyan; Radius:2.15; VdW:2.00; Mass:227),
 (Z:90; Sym:'Th'; Name:'Thorium'; Color:LightCyan; Radius:2.06; VdW:2.43; Mass:232.04),
 (Z:91; Sym:'Pa'; Name:'Protactinium'; Color:LightCyan; Radius:2.00; VdW:2.43; Mass:231.04),
 (Z:92; Sym:'U'; Name:'Uranium'; Color:LightCyan; Radius:1.96; VdW:2.41; Mass:238.03),
 (Z:93; Sym:'Np'; Name:'Neptunium'; Color:LightCyan; Radius:1.90; VdW:2.39; Mass:237),
 (Z:94; Sym:'Pu'; Name:'Plutonium'; Color:LightCyan; Radius:1.87; VdW:2.43; Mass:244),
 (Z:95; Sym:'Am'; Name:'Americium'; Color:LightCyan; Radius:1.80; VdW:2.44; Mass:243),
 (Z:96; Sym:'Cm'; Name:'Curium'; Color:LightCyan; Radius:1.69; VdW:2.45; Mass:247),
 (Z:97; Sym:'Bk'; Name:'Berkelium'; Color:LightCyan; Radius:1.68; VdW:2.44; Mass:247),
 (Z:98; Sym:'Cf'; Name:'Californium'; Color:LightCyan; Radius:1.68; VdW:2.45; Mass:251),
 (Z:99; Sym:'Es'; Name:'Einsteinium'; Color:LightCyan; Radius:1.65; VdW:2.45; Mass:252),
 (Z:100; Sym:'Fm'; Name:'Fermium'; Color:LightCyan; Radius:1.67; VdW:2.45; Mass:257),
 (Z:101; Sym:'Md'; Name:'Mendelevium'; Color:LightCyan; Radius:1.73; VdW:2.46; Mass:258),
 (Z:102; Sym:'No'; Name:'Nobelium'; Color:LightCyan; Radius:1.76; VdW:2.46; Mass:259),
 (Z:103; Sym:'Lr'; Name:'Lawrencium'; Color:LightCyan; Radius:1.61; VdW:2.46; Mass:262)
 );

 { Couleur sombre et claire associees a chaque couleur de base pour
   l'ombrage des spheres et des rubans }
 ShadeDark : Array[0..15] Of Byte = (
  Black, Black, Black, Blue, Black, Black, Black, DarkGray,
  Black, Blue, Green, Cyan, Red, Magenta, Brown, LightGray);
 ShadeLight : Array[0..15] Of Byte = (
  DarkGray, LightBlue, LightGreen, LightCyan, LightRed, LightMagenta, Yellow,
  White, LightGray, LightCyan, White, White, White, White, White, White);

 { Couleurs nommees (commande color, menus) }
 NUM_COLORS = 12;
 ColIdx : Array[0..NUM_COLORS-1] Of Byte = (
  Black, White, LightGray, LightRed, Brown, Yellow, LightGreen, LightCyan,
  LightBlue, Magenta, LightMagenta, DarkGray);
 ColNameFr : Array[0..NUM_COLORS-1] Of String[12] = (
  'Noir', 'Blanc', 'Gris', 'Rouge', 'Orange', 'Jaune', 'Vert', 'Cyan',
  'Bleu', 'Violet', 'Rose', 'Gris fonce');
 ColNameEn : Array[0..NUM_COLORS-1] Of String[12] = (
  'black', 'white', 'gray', 'red', 'orange', 'yellow', 'green', 'cyan',
  'blue', 'purple', 'pink', 'darkgray');

 { Couleurs arc-en-ciel (N-terminal -> C-terminal) }
 NUM_RAINBOW = 9;
 Rainbow : Array[0..NUM_RAINBOW-1] Of Byte = (
  Blue, LightBlue, Cyan, Green, LightGreen, Yellow, Brown, LightRed, Red);

 { Couleurs des chaines }
 ChainCols : Array[0..7] Of Byte = (
  LightGreen, LightCyan, LightMagenta, Yellow, LightRed, LightBlue, Brown,
  Cyan);

 { Expressions de selection du menu Selection }
 NUM_SELEXPR = 14;
 SelExpr : Array[0..NUM_SELEXPR-1] Of String[16] = (
  'all', 'none', 'protein', 'nucleic', 'water', 'hetero', 'hydrogen',
  'backbone', 'sidechain', 'helix', 'sheet', 'turn', 'ligand', 'carbon');

{ --------------------------------------------------------------------- }
{ Variables globales                                                    }
{ --------------------------------------------------------------------- }

Var
 { Ecran }
 SW, SH               : Integer;
 CvX1, CvY1, CvX2, CvY2 : Integer;     { zone de dessin 3D (absolue) }
 CvW, CvH             : Integer;
 ConsoleOn            : Boolean;

 { Donnees moleculaires }
 Atoms   : PAtomArr;
 Bonds   : PBondArr;
 Resi    : PResArr;
 Models  : PModelArr;
 HBs     : PHBArr;
 Vib     : PVecArr;
 NAtoms, NBonds, NRes, NModels, NHB : Integer;
 CurModel   : Integer;                 { 0 = tous les modeles }
 FileNameCur, MolTitle, MolHeader, MolFormat : String;
 Loaded     : Boolean;
 HasVibData : Boolean;
 AtomLimit  : Boolean;
 HasCell    : Boolean;
 CellP      : Array[0..5] Of Real;     { a, b, c, alpha, beta, gamma }
 CellPts    : TBoxPts;

 { Tableaux de travail (une image) }
 WX, WY, WZ : PRealArr;                { positions courantes }
 PX, PY, PR : PIntArr;                 { ecran + rayon en pixels }
 PZ         : PRealArr;                { profondeur }
 PVis       : PByteArr;                { 1 = atome dessine }
 Act        : PIntArr;                 { atomes du modele courant }
 NAct       : Integer;
 DotMask    : PDotMask;
 DotsDirty  : Boolean;
 Prims      : PPrimArr;
 NPrim      : Integer;
 Quads      : PQuadArr;
 NQuad      : Integer;
 RibX, RibY, RibZ : PRibReal;          { points de la spline }
 RibWX, RibWY, RibWZ : PRibReal;       { vecteurs de largeur }
 RibCol     : PRibByte;
 RibKind    : PRibByte;
 RibHW      : PRibByte;
 RibRes     : PRibInt;
 RibSX, RibSY : PRibInt;
 RibV, RibF : PRibReal;
 NRib       : Integer;
 UX, UY, UZ : Array[0..MAX_DOTS-1] Of Real;

 { Vue }
 RM         : Array[0..2, 0..2] Of Real;
 CenX, CenY, CenZ, MolRad : Real;
 ZoomF      : Real;
 PanX, PanY : Integer;
 Persp      : Boolean;
 BgCol      : Byte;
 BaseScale, CamD, ViewScale : Real;
 CxMid, CyMid : Integer;               { centre de la zone 3D (relatif) }

 { Options d'affichage }
 ShowHydro, ShowWater, ShowHB, ShowHalo : Boolean;
 ShowAxes, ShowBBox, ShowCell, ShowVec  : Boolean;
 SpinOn     : Boolean;
 SpinAxis   : Integer;                 { 0=X 1=Y 2=Z }
 SpinSpeed  : Real;
 AnimOn     : Boolean;
 AnimMode   : Integer;                 { 0=une fois 1=boucle 2=aller-retour }
 AnimDir    : Integer;
 AnimDelay  : Integer;                 { ms entre deux images }
 VibOn      : Boolean;
 VibPhase   : Real;
 VibScale   : Real;
 VibDelay   : Real;                    { 0..1 s }
 VibPanel   : Boolean;
 VibPaused  : Boolean;
 VecScale   : Real;
 PickMode   : Integer;
 PickQ      : Array[0..3] Of Integer;
 PickN      : Integer;
 MeasUnit   : Integer;                 { 0=A 1=nm 2=pm }
 Meas       : Array[1..MAX_MEAS] Of TMeas;
 NMeas      : Integer;
 LblFmt     : Array[1..MAX_LBLFMT] Of String[24];
 NLblFmt    : Integer;

 { Menus }
 PM         : PPMArr;

 { Souris }
 MouseX, MouseY, MouseB : LongInt;
 OldMB      : LongInt;
 MouseOk    : Boolean;
 DragMode   : Byte;                    { 0 rien, 1 rotation, 2 translation, 3 zoom }
 DragX, DragY : Integer;
 DragMoved  : Boolean;
 PressX, PressY : Integer;
 HoverAtom  : Integer;
 ToolHover  : Integer;

 { Console }
 ConLines   : Array[0..MAX_CONLN-1] Of String[80];
 ConCount   : Integer;
 ConInput   : String;
 ConHist    : Array[0..15] Of String[80];
 ConHistN, ConHistPos : Integer;

 { Etat de l'application }
 Running    : Boolean;
 NeedScene  : Boolean;
 NeedChrome : Boolean;
 NeedStatus : Boolean;
 NeedLayout : Boolean;
 StatMsg    : String;
 RenderMs   : LongInt;
 ScriptDepth: Integer;
 ScriptAbort: Boolean;
 PendingScript: String;
 TextBuf    : Array[0..MAX_TEXT-1] Of String[78];
 TextN      : Integer;

{ --------------------------------------------------------------------- }
{ Utilitaires : chaines, conversions, math                              }
{ --------------------------------------------------------------------- }

Function IntToStr(N: LongInt): String;
Var S: String;
Begin
 Str(N, S);
 IntToStr := S;
End;

Function RealToStr(R: Real; Decim: Integer): String;
Var S: String;
Begin
 Str(R:0:Decim, S);
 RealToStr := S;
End;

Function UpCh(C: Char): Char;
Begin
 If (C >= 'a') And (C <= 'z') Then UpCh := Chr(Ord(C) - 32) Else UpCh := C;
End;

Function LoCh(C: Char): Char;
Begin
 If (C >= 'A') And (C <= 'Z') Then LoCh := Chr(Ord(C) + 32) Else LoCh := C;
End;

Function UpStr(Const S: String): String;
Var I: Integer; R: String;
Begin
 R := S;
 For I := 1 To Length(R) Do R[I] := UpCh(R[I]);
 UpStr := R;
End;

Function LoStr(Const S: String): String;
Var I: Integer; R: String;
Begin
 R := S;
 For I := 1 To Length(R) Do R[I] := LoCh(R[I]);
 LoStr := R;
End;

Function MinI(A, B: Integer): Integer;
Begin If A < B Then MinI := A Else MinI := B; End;

Function MaxI(A, B: Integer): Integer;
Begin If A > B Then MaxI := A Else MaxI := B; End;

Function MinR(A, B: Real): Real;
Begin If A < B Then MinR := A Else MinR := B; End;

Function MaxR(A, B: Real): Real;
Begin If A > B Then MaxR := A Else MaxR := B; End;

Function ClampI(V, Lo, Hi: Integer): Integer;
Begin
 If V < Lo Then V := Lo;
 If V > Hi Then V := Hi;
 ClampI := V;
End;

Function ClampR(V, Lo, Hi: Real): Real;
Begin
 If V < Lo Then V := Lo;
 If V > Hi Then V := Hi;
 ClampR := V;
End;

Function AbsR(R: Real): Real;
Begin If R < 0 Then AbsR := -R Else AbsR := R; End;

Function RoundI(R: Real): Integer;
Begin
 If R > 30000 Then R := 30000;
 If R < -30000 Then R := -30000;
 RoundI := Round(R);
End;

Function TrimStr(Const S: String): String;
Var I, J: Integer;
Begin
 I := 1; J := Length(S);
 While (I <= J) And ((S[I] = ' ') Or (S[I] = #9) Or (S[I] = #13)) Do Inc(I);
 While (J >= I) And ((S[J] = ' ') Or (S[J] = #9) Or (S[J] = #13)) Do Dec(J);
 TrimStr := Copy(S, I, J - I + 1);
End;

Function StartsWith(Const S, Prefix: String): Boolean;
Begin
 StartsWith := Copy(S, 1, Length(Prefix)) = Prefix;
End;

Function StartsWithI(Const S, Prefix: String): Boolean;
Begin
 StartsWithI := UpStr(Copy(S, 1, Length(Prefix))) = UpStr(Prefix);
End;

Function StrToIntDef(Const S: String; Def: LongInt): LongInt;
Var V: LongInt; Code: Integer; T: String;
Begin
 T := TrimStr(S);
 If T = '' Then Begin StrToIntDef := Def; Exit; End;
 Val(T, V, Code);
 If Code = 0 Then StrToIntDef := V Else StrToIntDef := Def;
End;

Function StrToReal(Const S: String; Var V: Real): Boolean;
Var Code: Integer; T: String; R: Real;
Begin
 T := TrimStr(S);
 If T = '' Then Begin StrToReal := False; Exit; End;
 If T[1] = '+' Then Delete(T, 1, 1);
 Val(T, R, Code);
 If Code = 0 Then Begin V := R; StrToReal := True; End
 Else StrToReal := False;
End;

{ Sous-chaine 1-based [A..B] (tolere les chaines trop courtes) }
Function FieldS(Const S: String; A, B: Integer): String;
Begin
 If A > Length(S) Then Begin FieldS := ''; Exit; End;
 If B > Length(S) Then B := Length(S);
 FieldS := Copy(S, A, B - A + 1);
End;

Function FieldT(Const S: String; A, B: Integer): String;
Begin
 FieldT := TrimStr(FieldS(S, A, B));
End;

{ Decoupe une chaine en mots (separateurs : espace, tabulation) }
Procedure Tokenize(Const S: String; Var T: TTokens; Var N: Integer);
Var I, L: Integer; Cur: String;
Begin
 N := 0;
 Cur := '';
 L := Length(S);
 For I := 1 To L + 1 Do Begin
  If (I > L) Or (S[I] = ' ') Or (S[I] = #9) Or (S[I] = #13) Then Begin
   If Cur <> '' Then Begin
    If N <= High(T) Then Begin T[N] := Cur; Inc(N); End;
    Cur := '';
   End;
  End
  Else If Length(Cur) < 60 Then Cur := Cur + S[I];
 End;
End;

Function FixedStr(Const S: String; W: Integer): String;
Var R: String;
Begin
 R := S;
 While Length(R) < W Do R := R + ' ';
 FixedStr := Copy(R, 1, W);
End;

Function PadL(Const S: String; W: Integer): String;
Var R: String;
Begin
 R := S;
 While Length(R) < W Do R := ' ' + R;
 PadL := R;
End;

Function RealFmt(R: Real; W, D: Integer): String;
Var S: String;
Begin
 Str(R:W:D, S);
 RealFmt := S;
End;

{ Millisecondes depuis minuit (precision 10 ms) }
Function TickMs: LongInt;
Var H, M, S, S100: Word;
Begin
 GetTime(H, M, S, S100);
 TickMs := ((LongInt(H) * 60 + M) * 60 + S) * 1000 + LongInt(S100) * 10;
End;

Function ElapsedMs(T0: LongInt): LongInt;
Var D: LongInt;
Begin
 D := TickMs - T0;
 If D < 0 Then D := D + 86400000;
 ElapsedMs := D;
End;

{ --------------------------------------------------------------------- }
{ Gestion de la souris                                                  }
{ --------------------------------------------------------------------- }

{$IFNDEF FPC}
 Function MouseDriverFound: Boolean;
 Var
  Regs: Registers;
 Begin
  Regs.AX := 0;
  Intr($33, Regs);
  MouseDriverFound := Regs.AX = $FFFF;
 End;

 Procedure ShowMouse;
 Var
  Regs: Registers;
 Begin
  Regs.AX := $0001;
  Intr($33, Regs);
 End;

 Procedure HideMouse;
 Var
  Regs: Registers;
 Begin
  Regs.AX := $0002;
  Intr($33, Regs);
 End;

 Procedure GetMouseState(Var X, Y, Button: LongInt);
 Var
  Regs: Registers;
 Begin
  Regs.AX := $0003;
  Intr($33, Regs);
  Button := Regs.BX;
  X := Regs.CX;
  Y := Regs.DX;
 End;
{$ELSE}
 Function MouseDriverFound: Boolean;
 Begin
  PtcMouse.InitMouse;
  MouseDriverFound := True;
 End;

 Procedure ShowMouse;
 Begin
  PtcMouse.ShowMouse;
 End;

 Procedure HideMouse;
 Begin
  PtcMouse.HideMouse;
 End;

 Procedure GetMouseState(Var X, Y, Button: LongInt);
 Begin
  PtcMouse.GetMouseState(X, Y, Button);
 End;
{$ENDIF}

Procedure PollMouse;
Begin
 If MouseOk Then GetMouseState(MouseX, MouseY, MouseB)
 Else Begin MouseX := 0; MouseY := 0; MouseB := 0; End;
End;

Function PtIn(X, Y, X1, Y1, X2, Y2: Integer): Boolean;
Begin
 PtIn := (X >= X1) And (X <= X2) And (Y >= Y1) And (Y <= Y2);
End;

{ --------------------------------------------------------------------- }
{ Allocation de la memoire (tas) - tous les tableaux < 64 Ko            }
{ --------------------------------------------------------------------- }

Procedure AllocAll;
Begin
 New(Atoms);  New(Bonds);  New(Resi);  New(Models); New(HBs);  New(Vib);
 New(WX); New(WY); New(WZ);
 New(PX); New(PY); New(PR); New(PZ); New(PVis); New(Act); New(DotMask);
 New(Prims); New(Quads);
 New(RibX); New(RibY); New(RibZ);
 New(RibWX); New(RibWY); New(RibWZ);
 New(RibCol); New(RibKind); New(RibHW); New(RibRes);
 New(RibSX); New(RibSY); New(RibV); New(RibF);
 New(PM);
End;

{ --------------------------------------------------------------------- }
{ Initialisation graphique et primitives de dessin                      }
{ --------------------------------------------------------------------- }

Procedure InitializeGraphics;
Var
 Driver, Mode: Integer;
 ErrCode: Integer;
Begin
 {$IFDEF FPC}
  Driver:=VGA;
  Mode:=VGAHi;
 {$ELSE}
  Driver:=Detect;
  Mode:=VGAHi;
 {$ENDIF}
 InitGraph(Driver,Mode,'');
 ErrCode:=GraphResult;
 If ErrCode=grOk Then Begin
  SetColor(White);
  SetLineStyle(0,0,1);
  SW := GetMaxX + 1;
  SH := GetMaxY + 1;
 End
  Else
 Begin
  WriteLn('Erreur graphique : ',GraphErrorMsg(ErrCode));
  Halt;
 End;
End;

Procedure Txt(X, Y: Integer; Const S: String; C: Byte);
Begin
 SetColor(C);
 OutTextXY(X, Y, S);
End;

{ Cadre 3D : Raised=True releve, False enfonce }
Procedure Frame3D(X1, Y1, X2, Y2: Integer; Raised: Boolean);
Var C1, C2: Byte;
Begin
 If Raised Then Begin C1 := White; C2 := DarkGray; End
 Else Begin C1 := DarkGray; C2 := White; End;
 SetColor(C1);
 Line(X1, Y1, X2 - 1, Y1);
 Line(X1, Y1, X1, Y2 - 1);
 SetColor(C2);
 Line(X2, Y1, X2, Y2);
 Line(X1, Y2, X2, Y2);
End;

Procedure PanelBar(X1, Y1, X2, Y2: Integer; Col: Byte);
Begin
 SetFillStyle(SolidFill, Col);
 Bar(X1, Y1, X2, Y2);
End;

Procedure Button3D(X1, Y1, X2, Y2: Integer; Down: Boolean);
Begin
 PanelBar(X1, Y1, X2, Y2, LightGray);
 Frame3D(X1, Y1, X2, Y2, Not Down);
 If Not Down Then Begin
  SetColor(Black);
  Line(X2, Y1, X2, Y2);
  Line(X1, Y2, X2, Y2);
 End;
End;

{ Sphere ombree : disque sombre, corps, zone claire et reflet }
Procedure DrawBall(Cx, Cy, R: Integer; Col: Byte);
Var
 Dk, Lt, R1, R2, R3, O1, O2, O3: Integer;
Begin
 If R <= 1 Then Begin
  SetColor(Col);
  If R <= 0 Then PutPixel(Cx, Cy, Col)
  Else Begin
   SetFillStyle(SolidFill, Col);
   FillEllipse(Cx, Cy, R, R);
  End;
  Exit;
 End;
 Dk := ShadeDark[Col And 15];
 Lt := ShadeLight[Col And 15];
 If R <= 3 Then Begin
  SetColor(Dk); SetFillStyle(SolidFill, Dk);
  FillEllipse(Cx, Cy, R, R);
  SetColor(Col); SetFillStyle(SolidFill, Col);
  FillEllipse(Cx - 1, Cy - 1, R - 1, R - 1);
  PutPixel(Cx - R Div 2, Cy - R Div 2, White);
  Exit;
 End;
 R1 := (R * 7) Div 8; O1 := R Div 10;
 R2 := (R * 12) Div 25; O2 := (R * 3) Div 10;
 R3 := R Div 6; If R3 < 1 Then R3 := 1;
 O3 := (R * 4) Div 10;
 SetColor(Dk); SetFillStyle(SolidFill, Dk);
 FillEllipse(Cx, Cy, R, R);
 SetColor(Col); SetFillStyle(SolidFill, Col);
 FillEllipse(Cx - O1, Cy - O1, R1, R1);
 If Lt <> Col Then Begin
  SetColor(Lt); SetFillStyle(SolidFill, Lt);
  FillEllipse(Cx - O2, Cy - O2, R2, R2);
 End;
 SetColor(White); SetFillStyle(SolidFill, White);
 FillEllipse(Cx - O3, Cy - O3, R3, R3);
End;

{ Bande parallele a un segment : decalage Off (pixels), largeur Wid }
Procedure BandPoly(X1, Y1, X2, Y2: Integer; Nx, Ny, Off, Wid: Real; Col: Byte);
Var
 P: Array[0..3] Of PointType;
 Hw: Real;
Begin
 Hw := Wid * 0.5;
 P[0].X := RoundI(X1 + Nx * (Off + Hw)); P[0].Y := RoundI(Y1 + Ny * (Off + Hw));
 P[1].X := RoundI(X2 + Nx * (Off + Hw)); P[1].Y := RoundI(Y2 + Ny * (Off + Hw));
 P[2].X := RoundI(X2 + Nx * (Off - Hw)); P[2].Y := RoundI(Y2 + Ny * (Off - Hw));
 P[3].X := RoundI(X1 + Nx * (Off - Hw)); P[3].Y := RoundI(Y1 + Ny * (Off - Hw));
 SetColor(Col);
 SetFillStyle(SolidFill, Col);
 FillPoly(4, P);
End;

{ Cylindre (batonnet) de largeur totale W pixels ombre }
Procedure DrawStick(X1, Y1, X2, Y2, W: Integer; Col: Byte);
Var
 Dx, Dy, Len, Nx, Ny, Hw, Sg: Real;
 Dk, Lt: Byte;
Begin
 If W <= 1 Then Begin
  SetColor(Col);
  SetLineStyle(SolidLn, 0, NormWidth);
  Line(X1, Y1, X2, Y2);
  Exit;
 End;
 If W <= 3 Then Begin
  SetColor(Col);
  SetLineStyle(SolidLn, 0, ThickWidth);
  Line(X1, Y1, X2, Y2);
  SetLineStyle(SolidLn, 0, NormWidth);
  Exit;
 End;
 Dx := X2 - X1; Dy := Y2 - Y1;
 Len := Sqrt(Dx * Dx + Dy * Dy);
 If Len < 0.5 Then Exit;
 Nx := -Dy / Len; Ny := Dx / Len;
 { orienter la normale vers la source lumineuse (haut gauche) }
 Sg := 1.0;
 If (Nx * (-0.7) + Ny * (-0.7)) < 0 Then Sg := -1.0;
 Nx := Nx * Sg; Ny := Ny * Sg;
 Hw := W * 0.5;
 Dk := ShadeDark[Col And 15];
 Lt := ShadeLight[Col And 15];
 BandPoly(X1, Y1, X2, Y2, Nx, Ny, 0, W, Dk);
 BandPoly(X1, Y1, X2, Y2, Nx, Ny, Hw * 0.15, W * 0.72, Col);
 If Lt <> Col Then
  BandPoly(X1, Y1, X2, Y2, Nx, Ny, Hw * 0.38, W * 0.26, Lt);
End;

{ Ligne en pointilles }
Procedure DashLine(X1, Y1, X2, Y2: Integer; Col: Byte; Gap: Integer);
Var
 Dx, Dy, Len, T, Stp: Real;
 N, I: Integer;
Begin
 SetColor(Col);
 Dx := X2 - X1; Dy := Y2 - Y1;
 Len := Sqrt(Dx * Dx + Dy * Dy);
 If Len < 1 Then Exit;
 If Gap < 2 Then Gap := 2;
 N := Trunc(Len / Gap);
 If N < 1 Then Begin Line(X1, Y1, X2, Y2); Exit; End;
 Stp := 1.0 / N;
 For I := 0 To N - 1 Do Begin
  T := I * Stp;
  If (I And 1) = 0 Then
   Line(X1 + RoundI(Dx * T), Y1 + RoundI(Dy * T),
        X1 + RoundI(Dx * (T + Stp)), Y1 + RoundI(Dy * (T + Stp)));
 End;
End;

{ Fleche : segment + tete }
Procedure DrawArrow(X1, Y1, X2, Y2: Integer; Col: Byte);
Var
 Dx, Dy, Len, Ux, Uy: Real;
 P: Array[0..2] Of PointType;
Begin
 SetColor(Col);
 SetLineStyle(SolidLn, 0, NormWidth);
 Line(X1, Y1, X2, Y2);
 Dx := X2 - X1; Dy := Y2 - Y1;
 Len := Sqrt(Dx * Dx + Dy * Dy);
 If Len < 6 Then Exit;
 Ux := Dx / Len; Uy := Dy / Len;
 P[0].X := X2; P[0].Y := Y2;
 P[1].X := RoundI(X2 - Ux * 7 - Uy * 3); P[1].Y := RoundI(Y2 - Uy * 7 + Ux * 3);
 P[2].X := RoundI(X2 - Ux * 7 + Uy * 3); P[2].Y := RoundI(Y2 - Uy * 7 - Ux * 3);
 SetFillStyle(SolidFill, Col);
 FillPoly(3, P);
End;

{ Quadrilatere plein avec contour }
Procedure DrawQuad(Const Q: TQuad);
Var P: Array[0..3] Of PointType; I: Integer;
Begin
 For I := 0 To 3 Do Begin P[I].X := Q.X[I]; P[I].Y := Q.Y[I]; End;
 SetColor(Q.Col);
 SetFillStyle(SolidFill, Q.Col);
 FillPoly(4, P);
 SetColor(Q.Dk);
 SetLineStyle(SolidLn, 0, NormWidth);
 Line(Q.X[0], Q.Y[0], Q.X[1], Q.Y[1]);
 Line(Q.X[2], Q.Y[2], Q.X[3], Q.Y[3]);
End;

{ --------------------------------------------------------------------- }
{ Reconnaissance des elements et des residus                            }
{ --------------------------------------------------------------------- }

Function ElemBySymbol(Const Sym: String): Byte;
Var
 Cleaned : String;
 I : Integer;
 C : Char;
 R : Byte;
Begin
 Cleaned := '';
 For I := 1 To Length(Sym) Do Begin
  C := Sym[I];
  If ((C >= 'A') And (C <= 'Z')) Or ((C >= 'a') And (C <= 'z')) Then
   Cleaned := Cleaned + C;
 End;
 If Cleaned = '' Then Begin ElemBySymbol := 0; Exit; End;
 If Length(Cleaned) >= 2 Then Begin
  Cleaned[1] := UpCh(Cleaned[1]);
  Cleaned[2] := LoCh(Cleaned[2]);
  R := 0;
  For I := 1 To 103 Do
   If ElemDB[I].Sym = Copy(Cleaned, 1, 2) Then Begin R := I; Break; End;
  If R <> 0 Then Begin ElemBySymbol := R; Exit; End;
 End;
 Cleaned[1] := UpCh(Cleaned[1]);
 R := 0;
 For I := 1 To 103 Do
  If ElemDB[I].Sym = Copy(Cleaned, 1, 1) Then Begin R := I; Break; End;
 ElemBySymbol := R;
End;

{ Deduit l'element a partir du nom d'atome PDB }
Function ElemFromAtomName(Const Nm: String): Byte;
Var
 S, Two : String;
 I : Integer;
 R : Byte;
Begin
 S := TrimStr(Nm);
 ElemFromAtomName := 0;
 If S = '' Then Exit;
 I := 1;
 While (I <= Length(S)) And (S[I] >= '0') And (S[I] <= '9') Do Inc(I);
 If I > Length(S) Then Exit;
 If I + 1 <= Length(S) Then Begin
  Two := UpCh(S[I]) + LoCh(S[I + 1]);
  If (Two = 'He') Or (Two = 'Li') Or (Two = 'Be') Or (Two = 'Ne') Or
     (Two = 'Na') Or (Two = 'Mg') Or (Two = 'Al') Or (Two = 'Si') Or
     (Two = 'Cl') Or (Two = 'Ar') Or (Two = 'Ca') Or (Two = 'Sc') Or
     (Two = 'Ti') Or (Two = 'Cr') Or (Two = 'Mn') Or (Two = 'Fe') Or
     (Two = 'Co') Or (Two = 'Ni') Or (Two = 'Cu') Or (Two = 'Zn') Or
     (Two = 'Br') Or (Two = 'Se') Or (Two = 'Mo') Or (Two = 'Ag') Or
     (Two = 'Cd') Or (Two = 'Sn') Or (Two = 'Au') Or (Two = 'Hg') Or
     (Two = 'Pb') Or (Two = 'Bi') Or (Two = 'Ra') Or (Two = 'Pt') Then Begin
   R := ElemBySymbol(Two);
   If R <> 0 Then Begin ElemFromAtomName := R; Exit; End;
  End;
 End;
 ElemFromAtomName := ElemBySymbol(UpCh(S[I]));
End;

Function IsAminoAcidName(Const R: String): Boolean;
Var U: String;
Begin
 U := UpStr(TrimStr(R));
 IsAminoAcidName :=
  (U = 'ALA') Or (U = 'ARG') Or (U = 'ASN') Or (U = 'ASP') Or
  (U = 'CYS') Or (U = 'GLN') Or (U = 'GLU') Or (U = 'GLY') Or
  (U = 'HIS') Or (U = 'ILE') Or (U = 'LEU') Or (U = 'LYS') Or
  (U = 'MET') Or (U = 'PHE') Or (U = 'PRO') Or (U = 'SER') Or
  (U = 'THR') Or (U = 'TRP') Or (U = 'TYR') Or (U = 'VAL') Or
  (U = 'SEC') Or (U = 'PYL') Or (U = 'MSE') Or (U = 'HYP') Or
  (U = 'ASX') Or (U = 'GLX') Or (U = 'CSE') Or (U = 'CSO') Or
  (U = 'HID') Or (U = 'HIE') Or (U = 'HIP') Or (U = 'CYX');
End;

Function IsNucleicName(Const R: String): Boolean;
Var U: String;
Begin
 U := UpStr(TrimStr(R));
 IsNucleicName :=
  (U = 'A')  Or (U = 'G')  Or (U = 'C')  Or (U = 'U')  Or (U = 'T') Or
  (U = 'DA') Or (U = 'DG') Or (U = 'DC') Or (U = 'DT') Or (U = 'DU') Or
  (U = 'ADE') Or (U = 'GUA') Or (U = 'CYT') Or (U = 'THY') Or
  (U = 'URA') Or (U = '+A') Or (U = '+G') Or (U = '+C') Or (U = '+T');
End;

Function IsWaterName(Const R: String): Boolean;
Var U: String;
Begin
 U := UpStr(TrimStr(R));
 IsWaterName := (U = 'HOH') Or (U = 'WAT') Or (U = 'TIP') Or
                (U = 'TIP3') Or (U = 'SOL') Or (U = 'H2O') Or (U = 'DOD');
End;

Function AtomKindFor(Const ResName: String; HetAtm: Boolean): Byte;
Begin
 If IsWaterName(ResName) Then AtomKindFor := ATOMK_WATER
 Else If IsAminoAcidName(ResName) Then AtomKindFor := ATOMK_PROT
 Else If IsNucleicName(ResName) Then AtomKindFor := ATOMK_NUCL
 Else If HetAtm Then AtomKindFor := ATOMK_LIGAND
 Else AtomKindFor := ATOMK_HETERO;
End;

{ --------------------------------------------------------------------- }
{ Manipulation de la molecule                                           }
{ --------------------------------------------------------------------- }

Procedure InitAtom(Var A: TAtom);
Begin
 A.X := 0; A.Y := 0; A.Z := 0; A.Temp := 0;
 A.Serial := 0;
 A.Name := ''; A.Res := 'UNL'; A.Chain := ' ';
 A.ResSeq := 1; A.ICode := ' ';
 A.Model := 1; A.ElemZ := 0; A.Kind := ATOMK_LIGAND;
 A.Charge := 0; A.Sf := 23; A.Col := 255; A.Lbl := 0;
 A.ResIdx := -1;
 A.Sel := False; A.Hid := False; A.Dots := False;
End;

Procedure ClearMol;
Begin
 NAtoms := 0; NBonds := 0; NRes := 0; NModels := 1; NHB := 0;
 CurModel := 1;
 FileNameCur := ''; MolTitle := ''; MolHeader := ''; MolFormat := '';
 Loaded := False; HasVibData := False; HasCell := False; AtomLimit := False;
 NMeas := 0; PickN := 0; NAct := 0;
 NLblFmt := 0;
 DotsDirty := True;
 VibOn := False;
 HoverAtom := -1;
End;

Function AddAtom(Const A: TAtom): Integer;
Begin
 If NAtoms >= MAX_ATOMS Then Begin AddAtom := -1; AtomLimit := True; Exit; End;
 Atoms^[NAtoms] := A;
 AddAtom := NAtoms;
 Inc(NAtoms);
End;

Function FindBond(A, B: Integer): Integer;
Var I: Integer;
Begin
 FindBond := -1;
 For I := 0 To NBonds - 1 Do
  If ((Bonds^[I].A = A) And (Bonds^[I].B = B)) Or
     ((Bonds^[I].A = B) And (Bonds^[I].B = A)) Then Begin
   FindBond := I; Exit;
  End;
End;

{ Ajoute une liaison sans verifier les doublons }
Function AddBondFast(A, B: Integer; Ord: Byte): Integer;
Begin
 If (NBonds >= MAX_BONDS) Or (A = B) Then Begin AddBondFast := -1; Exit; End;
 Bonds^[NBonds].A := A;
 Bonds^[NBonds].B := B;
 Bonds^[NBonds].Ord := Ord;
 Bonds^[NBonds].BR := 15;
 AddBondFast := NBonds;
 Inc(NBonds);
End;

Function AddBond(A, B: Integer; Ord: Byte): Integer;
Var K: Integer;
Begin
 AddBond := -1;
 If (A < 0) Or (B < 0) Or (A >= NAtoms) Or (B >= NAtoms) Or (A = B) Then Exit;
 K := FindBond(A, B);
 If K >= 0 Then Begin AddBond := K; Exit; End;
 AddBond := AddBondFast(A, B, Ord);
End;

{ Tri rapide d'indices d'atomes selon X }
Procedure QSortX(Var Idx: TIntArr; Lo, Hi: Integer);
Var
 I, J, T: Integer;
 Pv: Real;
Begin
 While Lo < Hi Do Begin
  I := Lo; J := Hi;
  Pv := Atoms^[Idx[(Lo + Hi) Div 2]].X;
  Repeat
   While Atoms^[Idx[I]].X < Pv Do Inc(I);
   While Atoms^[Idx[J]].X > Pv Do Dec(J);
   If I <= J Then Begin
    T := Idx[I]; Idx[I] := Idx[J]; Idx[J] := T;
    Inc(I); Dec(J);
   End;
  Until I > J;
  If (J - Lo) < (Hi - I) Then Begin
   QSortX(Idx, Lo, J);
   Lo := I;
  End
  Else Begin
   QSortX(Idx, I, Hi);
   Hi := J;
  End;
 End;
End;

{ Detection des liaisons covalentes entre les atomes First..Last }
Procedure DetectBonds(First, Last: Integer);
Var
 N, I, J, K, Ia, Ib : Integer;
 Dx, Dy, Dz, D2, Lim, R : Real;
 Za, Zb : Byte;
Begin
 N := 0;
 For I := First To Last Do Begin Act^[N] := I; Inc(N); End;
 If N < 2 Then Exit;
 QSortX(Act^, 0, N - 1);
 For I := 0 To N - 2 Do Begin
  Ia := Act^[I];
  Za := Atoms^[Ia].ElemZ; If Za > 103 Then Za := 0;
  For K := I + 1 To N - 1 Do Begin
   Ib := Act^[K];
   Dx := Atoms^[Ib].X - Atoms^[Ia].X;
   If Dx > 3.6 Then Break;
   Dy := Atoms^[Ib].Y - Atoms^[Ia].Y;
   If (Dy > 3.6) Or (Dy < -3.6) Then Continue;
   Dz := Atoms^[Ib].Z - Atoms^[Ia].Z;
   If (Dz > 3.6) Or (Dz < -3.6) Then Continue;
   Zb := Atoms^[Ib].ElemZ; If Zb > 103 Then Zb := 0;
   If (Za = 1) And (Zb = 1) Then Continue;
   D2 := Dx * Dx + Dy * Dy + Dz * Dz;
   If D2 < 0.16 Then Continue;
   R := ElemDB[Za].Radius + ElemDB[Zb].Radius + 0.45;
   Lim := R * R;
   If D2 <= Lim Then Begin
    If Ia < Ib Then J := AddBondFast(Ia, Ib, 1) Else J := AddBondFast(Ib, Ia, 1);
    If J < 0 Then Exit;
   End;
  End;
 End;
End;

{ Ordre de liaison estime d'apres la longueur (fichiers sans ordres) }
Procedure GuessBondOrders;
Var
 I : Integer;
 Za, Zb, T : Byte;
 D : Real;
Begin
 For I := 0 To NBonds - 1 Do Begin
  Za := Atoms^[Bonds^[I].A].ElemZ;
  Zb := Atoms^[Bonds^[I].B].ElemZ;
  If Za > Zb Then Begin T := Za; Za := Zb; Zb := T; End;
  D := Sqrt(Sqr(Atoms^[Bonds^[I].A].X - Atoms^[Bonds^[I].B].X) +
            Sqr(Atoms^[Bonds^[I].A].Y - Atoms^[Bonds^[I].B].Y) +
            Sqr(Atoms^[Bonds^[I].A].Z - Atoms^[Bonds^[I].B].Z));
  Bonds^[I].Ord := 1;
  If (Za = 6) And (Zb = 6) Then Begin
   If D < 1.25 Then Bonds^[I].Ord := 3
   Else If D < 1.36 Then Bonds^[I].Ord := 2;
  End
  Else If (Za = 6) And (Zb = 7) Then Begin
   If D < 1.20 Then Bonds^[I].Ord := 3
   Else If D < 1.31 Then Bonds^[I].Ord := 2;
  End
  Else If (Za = 6) And (Zb = 8) Then Begin
   If D < 1.30 Then Bonds^[I].Ord := 2;
  End
  Else If (Za = 7) And (Zb = 8) Then Begin
   If D < 1.28 Then Bonds^[I].Ord := 2;
  End
  Else If ((Za = 8) And (Zb = 15)) Or ((Za = 8) And (Zb = 16)) Then Begin
   If D < 1.52 Then Bonds^[I].Ord := 2;
  End;
 End;
End;

Procedure ComputeCellPoints;
Var
 A, B, C, Al, Be, Ga, Ca, Cb, Cg, Sg, Cx, Cy, Cz: Real;
 I: Integer;
 Va, Vb, Vc: TVec3;
Begin
 If Not HasCell Then Exit;
 A := CellP[0]; B := CellP[1]; C := CellP[2];
 Al := CellP[3] * DEG2RAD; Be := CellP[4] * DEG2RAD; Ga := CellP[5] * DEG2RAD;
 Ca := Cos(Al); Cb := Cos(Be); Cg := Cos(Ga); Sg := Sin(Ga);
 If Sg < 0.01 Then Sg := 1.0;
 Va.X := A; Va.Y := 0; Va.Z := 0;
 Vb.X := B * Cg; Vb.Y := B * Sg; Vb.Z := 0;
 Cx := C * Cb;
 Cy := C * (Ca - Cb * Cg) / Sg;
 Cz := 1.0 - Cb * Cb - Sqr((Ca - Cb * Cg) / Sg);
 If Cz < 0 Then Cz := 0;
 Cz := C * Sqrt(Cz);
 Vc.X := Cx; Vc.Y := Cy; Vc.Z := Cz;
 For I := 0 To 7 Do Begin
  CellPts[I].X := (I And 1) * Va.X + ((I Shr 1) And 1) * Vb.X + ((I Shr 2) And 1) * Vc.X;
  CellPts[I].Y := (I And 1) * Va.Y + ((I Shr 1) And 1) * Vb.Y + ((I Shr 2) And 1) * Vc.Y;
  CellPts[I].Z := (I And 1) * Va.Z + ((I Shr 1) And 1) * Vb.Z + ((I Shr 2) And 1) * Vc.Z;
 End;
End;

{ Convertit des coordonnees fractionnaires en cartesiennes }
Procedure FracToCart(Fx, Fy, Fz: Real; Var X, Y, Z: Real);
Begin
 X := Fx * (CellPts[1].X - CellPts[0].X) + Fy * (CellPts[2].X - CellPts[0].X) +
      Fz * (CellPts[4].X - CellPts[0].X);
 Y := Fx * (CellPts[1].Y - CellPts[0].Y) + Fy * (CellPts[2].Y - CellPts[0].Y) +
      Fz * (CellPts[4].Y - CellPts[0].Y);
 Z := Fx * (CellPts[1].Z - CellPts[0].Z) + Fy * (CellPts[2].Z - CellPts[0].Z) +
      Fz * (CellPts[4].Z - CellPts[0].Z);
End;

{ --------------------------------------------------------------------- }
{ Residus, structure secondaire                                         }
{ --------------------------------------------------------------------- }

Type
 TSSERng = Record
  Kind     : Char;
  Chain    : Char;
  StartSeq : Integer;
  EndSeq   : Integer;
 End;

Var
 SSEs  : Array[0..199] Of TSSERng;
 NSSE  : Integer;

Function AtomDist(I, J: Integer): Real;
Begin
 AtomDist := Sqrt(Sqr(Atoms^[I].X - Atoms^[J].X) +
                  Sqr(Atoms^[I].Y - Atoms^[J].Y) +
                  Sqr(Atoms^[I].Z - Atoms^[J].Z));
End;

Procedure BuildResidues;
Var
 I, R : Integer;
 LastCh, LastIC : Char;
 LastSeq : Integer;
 LastName : String[4];
 LastModel : Word;
 Nm : String;
Begin
 NRes := 0;
 If NAtoms = 0 Then Exit;
 LastCh := #0; LastSeq := -32768; LastIC := #0; LastName := ''; LastModel := 0;
 R := -1;
 For I := 0 To NAtoms - 1 Do Begin
  If (Atoms^[I].Chain <> LastCh) Or (Atoms^[I].ResSeq <> LastSeq) Or
     (Atoms^[I].ICode <> LastIC) Or (Atoms^[I].Res <> LastName) Or
     (Atoms^[I].Model <> LastModel) Then Begin
   If R >= 0 Then Resi^[R].Last := I - 1;
   If NRes < MAX_RES Then Begin
    R := NRes;
    Resi^[R].Chain := Atoms^[I].Chain;
    Resi^[R].Seq := Atoms^[I].ResSeq;
    Resi^[R].ICode := Atoms^[I].ICode;
    Resi^[R].Name := Atoms^[I].Res;
    Resi^[R].First := I;
    Resi^[R].Last := I;
    Resi^[R].CA := -1; Resi^[R].CC := -1; Resi^[R].OO := -1;
    Resi^[R].Model := Atoms^[I].Model;
    Resi^[R].SS := 'L';
    Resi^[R].Rib := RIB_NONE;
    Resi^[R].RCol := 255;
    Resi^[R].IsAA := IsAminoAcidName(Atoms^[I].Res);
    Resi^[R].IsNA := IsNucleicName(Atoms^[I].Res);
    Inc(NRes);
   End
   Else R := -1;
   LastCh := Atoms^[I].Chain; LastSeq := Atoms^[I].ResSeq;
   LastIC := Atoms^[I].ICode; LastName := Atoms^[I].Res;
   LastModel := Atoms^[I].Model;
  End;
  If R >= 0 Then Begin
   Atoms^[I].ResIdx := R;
   Resi^[R].Last := I;
   Nm := UpStr(TrimStr(Atoms^[I].Name));
   If Resi^[R].IsAA Then Begin
    If (Nm = 'CA') And (Atoms^[I].ElemZ = 6) Then Resi^[R].CA := I
    Else If (Nm = 'C') And (Atoms^[I].ElemZ = 6) Then Resi^[R].CC := I
    Else If (Nm = 'O') Then Resi^[R].OO := I;
   End
   Else If Resi^[R].IsNA Then Begin
    If (Nm = 'P') Or (Nm = 'C4*') Or (Nm = 'C4''') Then Begin
     If (Resi^[R].CA < 0) Or (Nm = 'P') Then Resi^[R].CA := I;
    End
    Else If (Nm = 'O5*') Or (Nm = 'O5''') Or (Nm = 'OP1') Then
     If Resi^[R].OO < 0 Then Resi^[R].OO := I;
   End;
  End
  Else Atoms^[I].ResIdx := -1;
 End;
End;

{ Residus consecutifs d'une meme chaine ? }
Function ResLinked(R1, R2: Integer): Boolean;
Var D: Real;
Begin
 ResLinked := False;
 If (R1 < 0) Or (R2 >= NRes) Or (R1 >= NRes) Or (R2 < 0) Then Exit;
 If (Resi^[R1].Model <> Resi^[R2].Model) Or (Resi^[R1].Chain <> Resi^[R2].Chain) Then Exit;
 If (Resi^[R1].CA < 0) Or (Resi^[R2].CA < 0) Then Exit;
 If Resi^[R1].IsAA <> Resi^[R2].IsAA Then Exit;
 D := AtomDist(Resi^[R1].CA, Resi^[R2].CA);
 If Resi^[R1].IsAA Then ResLinked := D < 4.4 Else ResLinked := D < 8.0;
End;

{ Structure secondaire : enregistrements HELIX/SHEET ou heuristique CA }
Procedure AssignStructure;
Var
 I, J, K, RunS, RunL : Integer;
 D3, D4, D2 : Real;
 Ok : Boolean;
Begin
 If NRes = 0 Then Exit;
 For I := 0 To NRes - 1 Do Resi^[I].SS := 'L';
 If NSSE > 0 Then Begin
  For J := 0 To NSSE - 1 Do
   For I := 0 To NRes - 1 Do
    If ((SSEs[J].Chain = ' ') Or (Resi^[I].Chain = SSEs[J].Chain)) And
       (Resi^[I].Seq >= SSEs[J].StartSeq) And (Resi^[I].Seq <= SSEs[J].EndSeq) Then
     Resi^[I].SS := SSEs[J].Kind;
  Exit;
 End;
 { Heuristique basee sur les distances entre carbones alpha }
 For I := 0 To NRes - 1 Do Begin
  If Not Resi^[I].IsAA Then Continue;
  If Not (ResLinked(I, I + 1) And ResLinked(I + 1, I + 2) And
          ResLinked(I + 2, I + 3)) Then Continue;
  D3 := AtomDist(Resi^[I].CA, Resi^[I + 3].CA);
  Ok := (D3 > 4.4) And (D3 < 5.9);
  If Ok And ResLinked(I + 3, I + 4) Then Begin
   D4 := AtomDist(Resi^[I].CA, Resi^[I + 4].CA);
   Ok := (D4 > 5.5) And (D4 < 6.9);
  End;
  If Ok Then
   For K := I To I + 3 Do Resi^[K].SS := 'H';
 End;
 { brins : residus etendus }
 For I := 0 To NRes - 1 Do Begin
  If (Resi^[I].SS <> 'L') Or Not Resi^[I].IsAA Then Continue;
  If ResLinked(I - 1, I) And ResLinked(I, I + 1) Then Begin
   D2 := AtomDist(Resi^[I - 1].CA, Resi^[I + 1].CA);
   If D2 > 6.3 Then Resi^[I].SS := 'e';
  End;
 End;
 { une suite d'au moins 3 'e' devient un brin ; sinon boucle }
 I := 0;
 While I < NRes Do Begin
  If Resi^[I].SS = 'e' Then Begin
   RunS := I; RunL := 0;
   While (I < NRes) And (Resi^[I].SS = 'e') Do Begin Inc(RunL); Inc(I); End;
   For K := RunS To RunS + RunL - 1 Do
    If RunL >= 3 Then Resi^[K].SS := 'E' Else Resi^[K].SS := 'L';
  End
  Else Inc(I);
 End;
 { coudes : boucles dont les i,i+3 sont proches }
 For I := 0 To NRes - 1 Do Begin
  If (Resi^[I].SS <> 'L') Or Not Resi^[I].IsAA Then Continue;
  If ResLinked(I, I + 1) And ResLinked(I + 1, I + 2) And ResLinked(I + 2, I + 3) Then Begin
   D3 := AtomDist(Resi^[I].CA, Resi^[I + 3].CA);
   If D3 < 7.0 Then
    For K := I To I + 1 Do
     If Resi^[K].SS = 'L' Then Resi^[K].SS := 'T';
  End;
 End;
End;

{ Liaisons hydrogene (N/O ... N/O) entre residus eloignes }
Procedure ComputeHBonds;
Var
 I, J, K, N, Ia, Ib : Integer;
 D2, Dx, Dy, Dz : Real;
 Ok : Boolean;
Begin
 NHB := 0;
 N := 0;
 For I := 0 To NAtoms - 1 Do
  If ((Atoms^[I].ElemZ = 7) Or (Atoms^[I].ElemZ = 8)) And
     (Atoms^[I].Kind <> ATOMK_WATER) Then Begin Act^[N] := I; Inc(N); End;
 If N < 2 Then Exit;
 QSortX(Act^, 0, N - 1);
 For I := 0 To N - 2 Do Begin
  Ia := Act^[I];
  For K := I + 1 To N - 1 Do Begin
   Ib := Act^[K];
   Dx := Atoms^[Ib].X - Atoms^[Ia].X;
   If Dx > 3.3 Then Break;
   If Atoms^[Ia].Model <> Atoms^[Ib].Model Then Continue;
   Dy := Atoms^[Ib].Y - Atoms^[Ia].Y;
   If (Dy > 3.3) Or (Dy < -3.3) Then Continue;
   Dz := Atoms^[Ib].Z - Atoms^[Ia].Z;
   If (Dz > 3.3) Or (Dz < -3.3) Then Continue;
   D2 := Dx * Dx + Dy * Dy + Dz * Dz;
   If (D2 < 6.25) Or (D2 > 10.89) Then Continue;      { 2.5 .. 3.3 A }
   If (Atoms^[Ia].ElemZ = 8) And (Atoms^[Ib].ElemZ = 8) Then Continue;
   { pas de liaison H entre atomes proches dans la sequence }
   Ok := True;
   If (Atoms^[Ia].ResIdx >= 0) And (Atoms^[Ib].ResIdx >= 0) Then Begin
    If (Atoms^[Ia].Chain = Atoms^[Ib].Chain) And
       (Abs(Atoms^[Ia].ResIdx - Atoms^[Ib].ResIdx) < 3) Then Ok := False;
   End;
   If Ok Then Begin
    J := FindBond(Ia, Ib);
    If J >= 0 Then Ok := False;
   End;
   If Ok And (NHB < MAX_HB) Then Begin
    HBs^[NHB].A := Ia; HBs^[NHB].B := Ib;
    Inc(NHB);
   End;
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Finalisation d'un chargement                                          }
{ --------------------------------------------------------------------- }

Procedure ComputeBounds;
Var
 I, N : Integer;
 Mnx, Mxx, Mny, Mxy, Mnz, Mxz, D : Real;
 M1 : Integer;
Begin
 If NAtoms = 0 Then Begin
  CenX := 0; CenY := 0; CenZ := 0; MolRad := 5.0;
  Exit;
 End;
 M1 := Models^[1].Last;
 If (M1 < 0) Or (M1 >= NAtoms) Then M1 := NAtoms - 1;
 N := M1 + 1;
 If HasVibData And (NModels > 1) Then N := NAtoms;
 Mnx := Atoms^[0].X; Mxx := Mnx; Mny := Atoms^[0].Y; Mxy := Mny;
 Mnz := Atoms^[0].Z; Mxz := Mnz;
 For I := 1 To N - 1 Do Begin
  If Atoms^[I].X < Mnx Then Mnx := Atoms^[I].X;
  If Atoms^[I].X > Mxx Then Mxx := Atoms^[I].X;
  If Atoms^[I].Y < Mny Then Mny := Atoms^[I].Y;
  If Atoms^[I].Y > Mxy Then Mxy := Atoms^[I].Y;
  If Atoms^[I].Z < Mnz Then Mnz := Atoms^[I].Z;
  If Atoms^[I].Z > Mxz Then Mxz := Atoms^[I].Z;
 End;
 CenX := (Mnx + Mxx) * 0.5; CenY := (Mny + Mxy) * 0.5; CenZ := (Mnz + Mxz) * 0.5;
 MolRad := 1.0;
 For I := 0 To N - 1 Do Begin
  D := Sqrt(Sqr(Atoms^[I].X - CenX) + Sqr(Atoms^[I].Y - CenY) +
            Sqr(Atoms^[I].Z - CenZ));
  If D > MolRad Then MolRad := D;
 End;
 MolRad := MolRad + 1.6;
End;

Procedure ApplyDefaultStyle;
Var I: Integer;
Begin
 For I := 0 To NAtoms - 1 Do Begin
  Atoms^[I].Col := 255;
  Atoms^[I].Lbl := 0;
  Atoms^[I].Dots := False;
  If NAtoms > 3500 Then Atoms^[I].Sf := 0 Else Atoms^[I].Sf := 23;
 End;
 For I := 0 To NBonds - 1 Do
  If NAtoms > 3500 Then Bonds^[I].BR := 1 Else Bonds^[I].BR := 15;
 For I := 0 To NRes - 1 Do Begin
  Resi^[I].Rib := RIB_NONE;
  Resi^[I].RCol := 255;
 End;
End;

Procedure FinishLoad(Detect, Guess: Boolean; Const Fmt: String);
Var
 M, I, K, Off, Cnt1, Cnt, B0 : Integer;
Begin
 If NAtoms = 0 Then Exit;
 { informations de modele }
 If NModels < 1 Then NModels := 1;
 For M := 1 To NModels Do Begin
  Models^[M].First := -1; Models^[M].Last := -1;
 End;
 For I := 0 To NAtoms - 1 Do Begin
  M := Atoms^[I].Model;
  If (M >= 1) And (M <= NModels) Then Begin
   If Models^[M].First < 0 Then Models^[M].First := I;
   Models^[M].Last := I;
  End;
 End;
 { liaisons }
 If Detect Then Begin
  If Models^[1].First >= 0 Then DetectBonds(Models^[1].First, Models^[1].Last);
 End;
 B0 := NBonds;
 Cnt1 := 0;
 If Models^[1].First >= 0 Then Cnt1 := Models^[1].Last - Models^[1].First + 1;
 For M := 2 To NModels Do Begin
  If Models^[M].First < 0 Then Continue;
  Cnt := Models^[M].Last - Models^[M].First + 1;
  If (Cnt = Cnt1) And (B0 > 0) Then Begin
   Off := Models^[M].First - Models^[1].First;
   For K := 0 To B0 - 1 Do
    If Bonds^[K].B < Models^[1].Last + 1 Then
     AddBondFast(Bonds^[K].A + Off, Bonds^[K].B + Off, Bonds^[K].Ord);
  End
  Else If Detect Then DetectBonds(Models^[M].First, Models^[M].Last);
 End;
 If Guess Then GuessBondOrders;
 BuildResidues;
 AssignStructure;
 ComputeBounds;
 ApplyDefaultStyle;
 Loaded := True;
 DotsDirty := True;
 HasVibData := False;
 For M := 1 To NModels Do If Models^[M].Vib Then HasVibData := True;
 CurModel := 1;
 MolFormat := Fmt;
End;

Procedure BeginLoad(Const FN: String);
Var M: Integer;
Begin
 ClearMol;
 FileNameCur := FN;
 NSSE := 0;
 NModels := 1;
 For M := 1 To MAX_MODELS Do Begin
  Models^[M].Name := ''; Models^[M].Freq := 0; Models^[M].Vib := False;
  Models^[M].First := -1; Models^[M].Last := -1;
 End;
End;

Function OpenText(Const FN: String; Var F: Text): Boolean;
Begin
 {$I-}
 Assign(F, FN);
 Reset(F);
 {$I+}
 OpenText := IOResult = 0;
End;

Procedure ConsoleSay(Const S: String); Forward;
Procedure ConsoleErr(Const S: String); Forward;

{ --------------------------------------------------------------------- }
{ Chargeur PDB                                                          }
{ --------------------------------------------------------------------- }

Function ParseAtomLine(Const Line: String; HetAtm: Boolean; Var A: TAtom): Boolean;
Var
 X, Y, Z, Tf : Real;
 Sym : String;
 Z2 : Byte;
Begin
 ParseAtomLine := False;
 If Length(Line) < 54 Then Exit;
 InitAtom(A);
 A.Serial := StrToIntDef(FieldT(Line, 7, 11), 0);
 A.Name := FieldT(Line, 13, 16);
 A.Res := FieldT(Line, 18, 20);
 A.Chain := ' ';
 If Length(Line) >= 22 Then A.Chain := Line[22];
 A.ResSeq := StrToIntDef(FieldT(Line, 23, 26), 0);
 A.ICode := ' ';
 If Length(Line) >= 27 Then A.ICode := Line[27];
 If Not StrToReal(FieldT(Line, 31, 38), X) Then Exit;
 If Not StrToReal(FieldT(Line, 39, 46), Y) Then Exit;
 If Not StrToReal(FieldT(Line, 47, 54), Z) Then Exit;
 A.X := X; A.Y := Y; A.Z := Z;
 Tf := 0.0;
 If Length(Line) >= 66 Then StrToReal(FieldT(Line, 61, 66), Tf);
 A.Temp := Tf;
 Sym := '';
 If Length(Line) >= 78 Then Sym := FieldT(Line, 77, 78);
 Z2 := 0;
 If Sym <> '' Then Z2 := ElemBySymbol(Sym);
 If Z2 = 0 Then Z2 := ElemFromAtomName(A.Name);
 If (Z2 = 20) And (UpStr(A.Name) = 'CA') And IsAminoAcidName(A.Res) Then Z2 := 6;
 A.ElemZ := Z2;
 If Length(Line) >= 80 Then Begin
  Sym := FieldT(Line, 79, 80);
  If Sym <> '' Then A.Charge := StrToIntDef(Sym, 0);
 End;
 A.Kind := AtomKindFor(A.Res, HetAtm);
 ParseAtomLine := True;
End;

Function LoadPDB(Const FN: String): Boolean;
Var
 F : Text;
 Line, Rec : String;
 A : TAtom;
 I, K, IxA, IxB, CurM, Added : Integer;
 Sn : LongInt;
 Alt : Char;
 HaveCell, DetDone : Boolean;
 V : Real;
 HdrClass : String;
Begin
 LoadPDB := False;
 HdrClass := '';
 DetDone := False;
 If Not OpenText(FN, F) Then Exit;
 BeginLoad(FN);
 CurM := 1;
 HaveCell := False;
 Added := 0;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  If Line = '' Then Continue;
  Rec := UpStr(Copy(Line, 1, 6));
  If Rec = 'HEADER' Then Begin
   MolHeader := TrimStr(FieldS(Line, 63, 66));
   HdrClass := TrimStr(FieldS(Line, 11, 50));
  End
  Else If Rec = 'TITLE ' Then Begin
   If MolTitle <> '' Then MolTitle := MolTitle + ' ';
   MolTitle := MolTitle + TrimStr(FieldS(Line, 11, 80));
  End
  Else If Rec = 'CRYST1' Then Begin
   HaveCell := True;
   For K := 0 To 5 Do Begin
    V := 0;
    Case K Of
     0: StrToReal(FieldT(Line, 7, 15), V);
     1: StrToReal(FieldT(Line, 16, 24), V);
     2: StrToReal(FieldT(Line, 25, 33), V);
     3: StrToReal(FieldT(Line, 34, 40), V);
     4: StrToReal(FieldT(Line, 41, 47), V);
     5: StrToReal(FieldT(Line, 48, 54), V);
    End;
    CellP[K] := V;
   End;
   If (CellP[0] < 1.5) Or (CellP[1] < 1.5) Or (CellP[2] < 1.5) Then HaveCell := False;
   If CellP[3] < 1 Then CellP[3] := 90;
   If CellP[4] < 1 Then CellP[4] := 90;
   If CellP[5] < 1 Then CellP[5] := 90;
  End
  Else If Rec = 'HELIX ' Then Begin
   If NSSE < 200 Then Begin
    SSEs[NSSE].Kind := 'H';
    If Length(Line) >= 20 Then SSEs[NSSE].Chain := Line[20] Else SSEs[NSSE].Chain := ' ';
    SSEs[NSSE].StartSeq := StrToIntDef(FieldT(Line, 22, 25), 0);
    SSEs[NSSE].EndSeq   := StrToIntDef(FieldT(Line, 34, 37), 0);
    Inc(NSSE);
   End;
  End
  Else If Rec = 'SHEET ' Then Begin
   If NSSE < 200 Then Begin
    SSEs[NSSE].Kind := 'E';
    If Length(Line) >= 22 Then SSEs[NSSE].Chain := Line[22] Else SSEs[NSSE].Chain := ' ';
    SSEs[NSSE].StartSeq := StrToIntDef(FieldT(Line, 23, 26), 0);
    SSEs[NSSE].EndSeq   := StrToIntDef(FieldT(Line, 34, 37), 0);
    Inc(NSSE);
   End;
  End
  Else If Rec = 'MODEL ' Then Begin
   CurM := NModels + 1;
   If Added = 0 Then CurM := 1;
   If CurM > MAX_MODELS Then CurM := MAX_MODELS;
   If CurM > NModels Then NModels := CurM;
   Models^[CurM].Name := 'Modele ' + IntToStr(CurM);
  End
  Else If (Rec = 'ATOM  ') Or (Rec = 'HETATM') Then Begin
   If ParseAtomLine(Line, Rec = 'HETATM', A) Then Begin
    Alt := ' ';
    If Length(Line) >= 17 Then Alt := Line[17];
    If (Alt = ' ') Or (Alt = 'A') Or (Alt = '1') Then Begin
     A.Model := CurM;
     If AddAtom(A) >= 0 Then Inc(Added);
    End;
   End;
  End
  Else If Rec = 'CONECT' Then Begin
   If Not DetDone Then Begin
    DetDone := True;
    K := 0;
    For I := 0 To NAtoms - 1 Do If Atoms^[I].Model = 1 Then K := I;
    DetectBonds(0, K);
   End;
   Sn := StrToIntDef(FieldT(Line, 7, 11), 0);
   IxA := -1;
   For I := 0 To NAtoms - 1 Do
    If (Atoms^[I].Serial = Sn) And (Atoms^[I].Model = 1) Then Begin IxA := I; Break; End;
   If IxA >= 0 Then
    For K := 0 To 3 Do Begin
     Sn := StrToIntDef(FieldT(Line, 12 + K * 5, 16 + K * 5), 0);
     If Sn = 0 Then Continue;
     IxB := -1;
     For I := 0 To NAtoms - 1 Do
      If (Atoms^[I].Serial = Sn) And (Atoms^[I].Model = 1) Then Begin IxB := I; Break; End;
     If IxB >= 0 Then AddBond(IxA, IxB, 1);
    End;
  End
  Else If (Rec = 'END   ') Then Break;
 End;
 Close(F);
 If MolTitle = '' Then MolTitle := HdrClass;
 If MolTitle = '' Then MolTitle := FN;
 HasCell := HaveCell;
 If HasCell Then ComputeCellPoints;
 { Les CONECT du fichier sont conserves ; les autres liaisons sont detectees }
 FinishLoad(Not DetDone, True, 'PDB');
 LoadPDB := NAtoms > 0;
End;

{ --------------------------------------------------------------------- }
{ Chargeur XYZ (images multiples, vecteurs de vibration facultatifs)    }
{ --------------------------------------------------------------------- }

Function LoadXYZ(Const FN: String): Boolean;
Var
 F : Text;
 Line, Cmt : String;
 N, I, NT, Frame : Integer;
 T : TTokens;
 A : TAtom;
 X, Y, Z, Vx, Vy, Vz, Fq : Real;
 Z2 : Byte;
 Ix, J : Integer;
 HasV : Boolean;
Begin
 LoadXYZ := False;
 If Not OpenText(FN, F) Then Exit;
 BeginLoad(FN);
 NModels := 0;
 Frame := 0;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  If Line = '' Then Continue;
  N := StrToIntDef(Line, -1);
  If N <= 0 Then Break;
  If Eof(F) Then Break;
  ReadLn(F, Cmt);
  Cmt := TrimStr(Cmt);
  Inc(Frame);
  If Frame > MAX_MODELS Then Break;
  NModels := Frame;
  Models^[Frame].Name := Copy(Cmt, 1, 40);
  If (Frame = 1) And (Cmt <> '') Then MolTitle := Cmt;
  HasV := False;
  Fq := 0;
  Tokenize(Cmt, T, NT);
  For J := 0 To NT - 1 Do
   If StrToReal(T[J], Fq) And (Fq <> 0) Then Break
   Else Fq := 0;
  For I := 1 To N Do Begin
   If Eof(F) Then Break;
   ReadLn(F, Line);
   Tokenize(Line, T, NT);
   If NT < 4 Then Continue;
   If Not StrToReal(T[1], X) Then Continue;
   If Not StrToReal(T[2], Y) Then Continue;
   If Not StrToReal(T[3], Z) Then Continue;
   InitAtom(A);
   Z2 := ElemBySymbol(T[0]);
   If (Z2 = 0) And (StrToIntDef(T[0], 0) > 0) Then Z2 := StrToIntDef(T[0], 0);
   If Z2 > 103 Then Z2 := 0;
   A.ElemZ := Z2;
   A.Name := ElemDB[Z2].Sym;
   A.Serial := I;
   A.X := X; A.Y := Y; A.Z := Z;
   A.Model := Frame;
   A.ResSeq := 1;
   Ix := AddAtom(A);
   If Ix >= 0 Then Begin Vib^[Ix].X := 0; Vib^[Ix].Y := 0; Vib^[Ix].Z := 0; End;
   If (Ix >= 0) And (NT >= 7) Then Begin
    If StrToReal(T[4], Vx) And StrToReal(T[5], Vy) And StrToReal(T[6], Vz) Then Begin
     Vib^[Ix].X := Vx; Vib^[Ix].Y := Vy; Vib^[Ix].Z := Vz;
     HasV := True;
    End;
   End;
  End;
  If HasV Then Begin
   Models^[Frame].Vib := True;
   Models^[Frame].Freq := Fq;
  End;
 End;
 Close(F);
 If NModels < 1 Then NModels := 1;
 If MolTitle = '' Then MolTitle := FN;
 FinishLoad(True, True, 'XYZ');
 LoadXYZ := NAtoms > 0;
End;

{ --------------------------------------------------------------------- }
{ Chargeur MOL / SDF (V2000)                                            }
{ --------------------------------------------------------------------- }

Function LoadMOLFile(Const FN: String; Multi: Boolean): Boolean;
Var
 F : Text;
 Line, Counts, Sym : String;
 NA, NB, I, A1, A2, Tp, Base, Frame : Integer;
 X, Y, Z : Real;
 A : TAtom;
Begin
 LoadMOLFile := False;
 If Not OpenText(FN, F) Then Exit;
 BeginLoad(FN);
 NModels := 0;
 Frame := 0;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  Inc(Frame);
  If Frame > MAX_MODELS Then Break;
  NModels := Frame;
  Models^[Frame].Name := Copy(Line, 1, 40);
  If Frame = 1 Then MolTitle := Line;
  If Eof(F) Then Break;
  ReadLn(F, Line);
  ReadLn(F, Line);
  If Eof(F) Then Break;
  ReadLn(F, Counts);
  NA := StrToIntDef(TrimStr(FieldS(Counts, 1, 3)), 0);
  NB := StrToIntDef(TrimStr(FieldS(Counts, 4, 6)), 0);
  If NA <= 0 Then Begin
   Dec(Frame);
   NModels := Frame;
   Break;
  End;
  Base := NAtoms;
  For I := 1 To NA Do Begin
   If Eof(F) Then Break;
   ReadLn(F, Line);
   Sym := TrimStr(FieldS(Line, 32, 34));
   If Sym = '' Then Sym := TrimStr(FieldS(Line, 31, 34));
   X := 0; Y := 0; Z := 0;
   StrToReal(FieldT(Line, 1, 10), X);
   StrToReal(FieldT(Line, 11, 20), Y);
   StrToReal(FieldT(Line, 21, 30), Z);
   InitAtom(A);
   A.ElemZ := ElemBySymbol(Sym);
   A.Name := Copy(Sym, 1, 4);
   A.Serial := I;
   A.X := X; A.Y := Y; A.Z := Z;
   A.Model := Frame;
   AddAtom(A);
  End;
  For I := 1 To NB Do Begin
   If Eof(F) Then Break;
   ReadLn(F, Line);
   A1 := StrToIntDef(FieldT(Line, 1, 3), 0);
   A2 := StrToIntDef(FieldT(Line, 4, 6), 0);
   Tp := StrToIntDef(FieldT(Line, 7, 9), 1);
   If Tp > 4 Then Tp := 1;
   If (Frame = 1) And (A1 >= 1) And (A2 >= 1) And (A1 <= NA) And (A2 <= NA) Then
    AddBondFast(Base + A1 - 1, Base + A2 - 1, Tp);
  End;
  { sauter jusqu'a la fin de l'enregistrement }
  While Not Eof(F) Do Begin
   ReadLn(F, Line);
   If StartsWith(Line, '$$$$') Then Break;
   If Not Multi Then Break;
  End;
  If Not Multi Then Break;
 End;
 Close(F);
 If NModels < 1 Then NModels := 1;
 If MolTitle = '' Then MolTitle := FN;
 FinishLoad(NBonds = 0, NBonds = 0, 'MDL MOL');
 LoadMOLFile := NAtoms > 0;
End;

{ --------------------------------------------------------------------- }
{ Chargeur MOL2 (Tripos)                                                }
{ --------------------------------------------------------------------- }

Function LoadMOL2(Const FN: String): Boolean;
Var
 F : Text;
 Line, Sym, Bt : String;
 Section, MolLn, Idx, A1, A2 : Integer;
 T : TTokens;
 NT : Integer;
 X, Y, Z : Real;
 A : TAtom;
 O : Byte;
Begin
 LoadMOL2 := False;
 If Not OpenText(FN, F) Then Exit;
 BeginLoad(FN);
 Section := 0; MolLn := 0;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  If Line = '' Then Continue;
  If StartsWith(Line, '@<TRIPOS>') Then Begin
   If UpStr(Line) = '@<TRIPOS>MOLECULE' Then Begin Section := 3; MolLn := 0; End
   Else If UpStr(Line) = '@<TRIPOS>ATOM' Then Section := 1
   Else If UpStr(Line) = '@<TRIPOS>BOND' Then Section := 2
   Else Section := 0;
   Continue;
  End;
  If Section = 3 Then Begin
   Inc(MolLn);
   If MolLn = 1 Then MolTitle := Line;
  End
  Else If Section = 1 Then Begin
   Tokenize(Line, T, NT);
   If NT < 6 Then Continue;
   If Not StrToReal(T[2], X) Then Continue;
   If Not StrToReal(T[3], Y) Then Continue;
   If Not StrToReal(T[4], Z) Then Continue;
   InitAtom(A);
   A.Serial := StrToIntDef(T[0], NAtoms + 1);
   A.Name := Copy(T[1], 1, 4);
   Sym := T[5];
   Idx := Pos('.', Sym);
   If Idx > 0 Then Sym := Copy(Sym, 1, Idx - 1);
   A.ElemZ := ElemBySymbol(Sym);
   If A.ElemZ = 0 Then A.ElemZ := ElemFromAtomName(A.Name);
   A.X := X; A.Y := Y; A.Z := Z;
   If NT >= 7 Then A.ResSeq := StrToIntDef(T[6], 1);
   If NT >= 8 Then A.Res := Copy(T[7], 1, 3);
   A.Kind := AtomKindFor(A.Res, True);
   AddAtom(A);
  End
  Else If Section = 2 Then Begin
   Tokenize(Line, T, NT);
   If NT < 4 Then Continue;
   A1 := StrToIntDef(T[1], 0) - 1;
   A2 := StrToIntDef(T[2], 0) - 1;
   Bt := UpStr(T[3]);
   O := 1;
   If Bt = '2' Then O := 2
   Else If Bt = '3' Then O := 3
   Else If (Bt = 'AR') Or (Bt = 'AM') Then O := 4;
   If (A1 >= 0) And (A2 >= 0) And (A1 < NAtoms) And (A2 < NAtoms) Then
    AddBondFast(A1, A2, O);
  End;
 End;
 Close(F);
 If MolTitle = '' Then MolTitle := FN;
 FinishLoad(NBonds = 0, NBonds = 0, 'Tripos MOL2');
 LoadMOL2 := NAtoms > 0;
End;

{ --------------------------------------------------------------------- }
{ Chargeur CIF (mmCIF/PDBx et petites molecules)                        }
{ --------------------------------------------------------------------- }

Function LoadCIF(Const FN: String): Boolean;
Var
 F : Text;
 Line, Key : String;
 InLoop, IsAtomLoop, GotCols : Boolean;
 NCols, NT, I, Cur : Integer;
 T : TTokens;
 CGrp, CLabel, CComp, CSeq, CAsym, CAlt, CX, CY, CZ, CB, CEl, CIns, CMod,
 CAuthAsym, CAuthSeq, CFx, CFy, CFz : Integer;
 A : TAtom;
 X, Y, Z, V : Real;
 FirstModel, Mv : Integer;
 Frac : Boolean;
 SeenCell : Integer;

 Procedure ResetCols;
 Begin
  CGrp := -1; CLabel := -1; CComp := -1; CSeq := -1; CAsym := -1; CAlt := -1;
  CX := -1; CY := -1; CZ := -1; CB := -1; CEl := -1; CIns := -1; CMod := -1;
  CAuthAsym := -1; CAuthSeq := -1; CFx := -1; CFy := -1; CFz := -1;
 End;

Begin
 LoadCIF := False;
 If Not OpenText(FN, F) Then Exit;
 BeginLoad(FN);
 InLoop := False; IsAtomLoop := False; GotCols := False;
 NCols := 0; ResetCols; FirstModel := -1; SeenCell := 0;
 For I := 0 To 5 Do CellP[I] := 0;
 CellP[3] := 90; CellP[4] := 90; CellP[5] := 90;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  If Line = '' Then Continue;
  If Line[1] = '#' Then Continue;
  If StartsWithI(Line, 'data_') Then Begin
   If MolTitle = '' Then MolTitle := Copy(Line, 6, Length(Line) - 5);
   Continue;
  End;
  Key := LoStr(Line);
  If StartsWith(Key, '_struct.title') Then Begin
   If Pos(' ', Line) > 0 Then MolTitle := TrimStr(Copy(Line, Pos(' ', Line), 200));
   If (Length(MolTitle) > 0) And (MolTitle[1] = '''') Then MolTitle := Copy(MolTitle, 2, Length(MolTitle) - 2);
   Continue;
  End;
  If (Key[1] = '_') And (Not InLoop Or GotCols) Then Begin
   Tokenize(Line, T, NT);
   If NT >= 2 Then Begin
    If (StartsWith(Key, '_cell_length_a')) Or (StartsWith(Key, '_cell.length_a')) Then Begin
     If StrToReal(T[1], V) Then CellP[0] := V; Inc(SeenCell);
    End
    Else If (StartsWith(Key, '_cell_length_b')) Or (StartsWith(Key, '_cell.length_b')) Then Begin
     If StrToReal(T[1], V) Then CellP[1] := V; Inc(SeenCell);
    End
    Else If (StartsWith(Key, '_cell_length_c')) Or (StartsWith(Key, '_cell.length_c')) Then Begin
     If StrToReal(T[1], V) Then CellP[2] := V; Inc(SeenCell);
    End
    Else If (StartsWith(Key, '_cell_angle_alpha')) Or (StartsWith(Key, '_cell.angle_alpha')) Then Begin
     If StrToReal(T[1], V) Then CellP[3] := V;
    End
    Else If (StartsWith(Key, '_cell_angle_beta')) Or (StartsWith(Key, '_cell.angle_beta')) Then Begin
     If StrToReal(T[1], V) Then CellP[4] := V;
    End
    Else If (StartsWith(Key, '_cell_angle_gamma')) Or (StartsWith(Key, '_cell.angle_gamma')) Then Begin
     If StrToReal(T[1], V) Then CellP[5] := V;
    End;
   End;
  End;
  If Key = 'loop_' Then Begin
   InLoop := True; IsAtomLoop := False; GotCols := False;
   NCols := 0; ResetCols;
   Continue;
  End;
  If InLoop And (Line[1] = '_') And Not GotCols Then Begin
   If NCols < 40 Then Begin
    If StartsWith(Key, '_atom_site') Then Begin
     IsAtomLoop := True;
     If Key = '_atom_site.group_pdb' Then CGrp := NCols
     Else If (Key = '_atom_site.label_atom_id') Or (Key = '_atom_site_label') Then CLabel := NCols
     Else If Key = '_atom_site.label_comp_id' Then CComp := NCols
     Else If Key = '_atom_site.label_seq_id' Then CSeq := NCols
     Else If Key = '_atom_site.label_asym_id' Then CAsym := NCols
     Else If Key = '_atom_site.label_alt_id' Then CAlt := NCols
     Else If Key = '_atom_site.cartn_x' Then CX := NCols
     Else If Key = '_atom_site.cartn_y' Then CY := NCols
     Else If Key = '_atom_site.cartn_z' Then CZ := NCols
     Else If Key = '_atom_site.b_iso_or_equiv' Then CB := NCols
     Else If (Key = '_atom_site.type_symbol') Or (Key = '_atom_site_type_symbol') Then CEl := NCols
     Else If Key = '_atom_site.pdbx_pdb_ins_code' Then CIns := NCols
     Else If Key = '_atom_site.pdbx_pdb_model_num' Then CMod := NCols
     Else If Key = '_atom_site.auth_asym_id' Then CAuthAsym := NCols
     Else If Key = '_atom_site.auth_seq_id' Then CAuthSeq := NCols
     Else If Key = '_atom_site_fract_x' Then CFx := NCols
     Else If Key = '_atom_site_fract_y' Then CFy := NCols
     Else If Key = '_atom_site_fract_z' Then CFz := NCols;
    End;
    Inc(NCols);
   End;
   Continue;
  End;
  If InLoop Then Begin
   If Line[1] = '_' Then Begin InLoop := False; Continue; End;
   GotCols := True;
   If Not IsAtomLoop Then Continue;
   Frac := (CFx >= 0) And (CX < 0);
   If (CX < 0) And Not Frac Then Continue;
   Tokenize(Line, T, NT);
   If NT < NCols Then Continue;
   If CMod >= 0 Then Begin
    Mv := StrToIntDef(T[CMod], 1);
    If FirstModel < 0 Then FirstModel := Mv;
    Cur := Mv - FirstModel + 1;
    If (Cur < 1) Or (Cur > MAX_MODELS) Then Continue;
    If Cur > NModels Then Begin
     NModels := Cur;
     Models^[Cur].Name := 'Modele ' + IntToStr(Cur);
    End;
   End
   Else Cur := 1;
   InitAtom(A);
   A.Serial := NAtoms + 1;
   A.Model := Cur;
   If CLabel >= 0 Then A.Name := Copy(T[CLabel], 1, 4) Else A.Name := 'X';
   If (Length(A.Name) > 1) And (A.Name[1] = '"') Then A.Name := Copy(A.Name, 2, Length(A.Name) - 2);
   If (CAlt >= 0) And (T[CAlt] <> '.') And (T[CAlt] <> '?') And
      (T[CAlt] <> 'A') Then Continue;
   If CComp >= 0 Then A.Res := Copy(T[CComp], 1, 3);
   If (CAuthAsym >= 0) Then A.Chain := T[CAuthAsym][1]
   Else If CAsym >= 0 Then A.Chain := T[CAsym][1];
   If CAuthSeq >= 0 Then A.ResSeq := StrToIntDef(T[CAuthSeq], 0)
   Else If CSeq >= 0 Then A.ResSeq := StrToIntDef(T[CSeq], 0);
   If (CIns >= 0) And (T[CIns] <> '?') And (T[CIns] <> '.') Then A.ICode := T[CIns][1];
   If Frac Then Begin
    If Not StrToReal(T[CFx], X) Then Continue;
    If Not StrToReal(T[CFy], Y) Then Continue;
    If Not StrToReal(T[CFz], Z) Then Continue;
    If SeenCell >= 3 Then Begin
     HasCell := True; ComputeCellPoints;
     FracToCart(X, Y, Z, X, Y, Z);
    End;
   End
   Else Begin
    If Not StrToReal(T[CX], X) Then Continue;
    If Not StrToReal(T[CY], Y) Then Continue;
    If Not StrToReal(T[CZ], Z) Then Continue;
   End;
   A.X := X; A.Y := Y; A.Z := Z;
   If (CB >= 0) And StrToReal(T[CB], V) Then A.Temp := V;
   If CEl >= 0 Then A.ElemZ := ElemBySymbol(T[CEl])
   Else A.ElemZ := ElemFromAtomName(A.Name);
   If A.ElemZ = 0 Then A.ElemZ := ElemFromAtomName(A.Name);
   If (A.ElemZ = 20) And (UpStr(A.Name) = 'CA') And IsAminoAcidName(A.Res) Then A.ElemZ := 6;
   A.Kind := AtomKindFor(A.Res, (CGrp >= 0) And (UpStr(T[CGrp]) = 'HETATM'));
   AddAtom(A);
  End;
 End;
 Close(F);
 If (SeenCell >= 3) And (CellP[0] > 1.5) Then Begin
  HasCell := True;
  ComputeCellPoints;
 End;
 If MolTitle = '' Then MolTitle := FN;
 FinishLoad(True, True, 'CIF');
 LoadCIF := NAtoms > 0;
End;

{ --------------------------------------------------------------------- }
{ Detection du format / chargement                                      }
{ --------------------------------------------------------------------- }

Function ExtractExt(Const FN: String): String;
Var I: Integer;
Begin
 ExtractExt := '';
 For I := Length(FN) DownTo 1 Do Begin
  If FN[I] = '.' Then Begin
   ExtractExt := UpStr(Copy(FN, I + 1, Length(FN) - I));
   Exit;
  End;
  If (FN[I] = '\') Or (FN[I] = '/') Then Exit;
 End;
End;

Function ExtractName(Const FN: String): String;
Var I: Integer;
Begin
 ExtractName := FN;
 For I := Length(FN) DownTo 1 Do
  If (FN[I] = '\') Or (FN[I] = '/') Then Begin
   ExtractName := Copy(FN, I + 1, Length(FN) - I);
   Exit;
  End;
End;

Function LoadMolecule(Const FN: String): Boolean;
Var E: String; R: Boolean;
Begin
 E := ExtractExt(FN);
 R := False;
 If (E = 'PDB') Or (E = 'ENT') Then R := LoadPDB(FN)
 Else If E = 'XYZ' Then R := LoadXYZ(FN)
 Else If E = 'MOL' Then R := LoadMOLFile(FN, False)
 Else If (E = 'SDF') Or (E = 'SD') Then R := LoadMOLFile(FN, True)
 Else If E = 'MOL2' Then R := LoadMOL2(FN)
 Else If (E = 'CIF') Or (E = 'MMCIF') Then R := LoadCIF(FN)
 Else Begin
  R := LoadPDB(FN);
  If Not R Then R := LoadXYZ(FN);
  If Not R Then R := LoadMOLFile(FN, False);
 End;
 LoadMolecule := R;
End;

{ --------------------------------------------------------------------- }
{ Vue : rotation, zoom, projection                                      }
{ --------------------------------------------------------------------- }

Procedure IdentityRM;
Var I, J: Integer;
Begin
 For I := 0 To 2 Do
  For J := 0 To 2 Do
   If I = J Then RM[I, J] := 1.0 Else RM[I, J] := 0.0;
End;

{ Re-orthonormalise la matrice de rotation }
Procedure OrthoRM;
Var
 L : Real;
 I : Integer;
 D : Real;
Begin
 L := Sqrt(Sqr(RM[0, 0]) + Sqr(RM[0, 1]) + Sqr(RM[0, 2]));
 If L < 1E-6 Then Begin IdentityRM; Exit; End;
 For I := 0 To 2 Do RM[0, I] := RM[0, I] / L;
 D := RM[0, 0] * RM[1, 0] + RM[0, 1] * RM[1, 1] + RM[0, 2] * RM[1, 2];
 For I := 0 To 2 Do RM[1, I] := RM[1, I] - D * RM[0, I];
 L := Sqrt(Sqr(RM[1, 0]) + Sqr(RM[1, 1]) + Sqr(RM[1, 2]));
 If L < 1E-6 Then Begin IdentityRM; Exit; End;
 For I := 0 To 2 Do RM[1, I] := RM[1, I] / L;
 RM[2, 0] := RM[0, 1] * RM[1, 2] - RM[0, 2] * RM[1, 1];
 RM[2, 1] := RM[0, 2] * RM[1, 0] - RM[0, 0] * RM[1, 2];
 RM[2, 2] := RM[0, 0] * RM[1, 1] - RM[0, 1] * RM[1, 0];
End;

{ Rotation incrementale autour d'un axe de l'ecran (0=X 1=Y 2=Z) }
Procedure RotateView(Ax: Integer; Deg: Real);
Var
 R, N : Array[0..2, 0..2] Of Real;
 C, S : Real;
 I, J, K : Integer;
Begin
 C := Cos(Deg * DEG2RAD); S := Sin(Deg * DEG2RAD);
 For I := 0 To 2 Do For J := 0 To 2 Do
  If I = J Then R[I, J] := 1.0 Else R[I, J] := 0.0;
 Case Ax Of
  0: Begin R[1, 1] := C; R[1, 2] := -S; R[2, 1] := S; R[2, 2] := C; End;
  1: Begin R[0, 0] := C; R[0, 2] := S; R[2, 0] := -S; R[2, 2] := C; End;
  2: Begin R[0, 0] := C; R[0, 1] := -S; R[1, 0] := S; R[1, 1] := C; End;
 End;
 For I := 0 To 2 Do
  For J := 0 To 2 Do Begin
   N[I, J] := 0;
   For K := 0 To 2 Do N[I, J] := N[I, J] + R[I, K] * RM[K, J];
  End;
 For I := 0 To 2 Do For J := 0 To 2 Do RM[I, J] := N[I, J];
End;

Procedure SetupView;
Begin
 CvW := CvX2 - CvX1 + 1;
 CvH := CvY2 - CvY1 + 1;
 CxMid := CvW Div 2;
 CyMid := CvH Div 2;
 If MolRad < 1.0 Then MolRad := 1.0;
 BaseScale := MinR(CvW, CvH) / (2.0 * MolRad) * 0.96;
 CamD := 5.0 * MolRad;
 ViewScale := BaseScale * ZoomF;
End;

Procedure ResetView;
Begin
 IdentityRM;
 ZoomF := 1.0;
 PanX := 0; PanY := 0;
 SetupView;
 NeedScene := True;
End;

Procedure SetOrientation(Ax, Deg: Integer);
Begin
 IdentityRM;
 If Deg <> 0 Then RotateView(Ax, Deg);
 NeedScene := True;
End;

Function FgCol: Byte;
Begin
 If (BgCol = White) Or (BgCol = LightGray) Or (BgCol = Yellow) Or
    (BgCol = LightCyan) Or (BgCol = LightGreen) Then FgCol := Black
 Else FgCol := White;
End;

Function AtomColor(I: Integer): Byte;
Var Z: Byte;
Begin
 If Atoms^[I].Col <> 255 Then Begin AtomColor := Atoms^[I].Col; Exit; End;
 Z := Atoms^[I].ElemZ; If Z > 103 Then Z := 0;
 AtomColor := ElemDB[Z].Color;
End;

Function AtomVdW(I: Integer): Real;
Var Z: Byte;
Begin
 Z := Atoms^[I].ElemZ; If Z > 103 Then Z := 0;
 AtomVdW := ElemDB[Z].VdW;
End;

Function FOfZ(Z: Real): Real;
Var F: Real;
Begin
 If Persp Then Begin
  F := CamD - Z;
  If F < 1.0 Then F := 1.0;
  FOfZ := CamD / F;
 End
 Else FOfZ := 1.0;
End;

Procedure ProjPoint(X, Y, Z: Real; Var Sx, Sy: Integer; Var Zv, F: Real);
Var
 Dx, Dy, Dz, Vx, Vy : Real;
 Tx, Ty : LongInt;
Begin
 Dx := X - CenX; Dy := Y - CenY; Dz := Z - CenZ;
 Vx := RM[0, 0] * Dx + RM[0, 1] * Dy + RM[0, 2] * Dz;
 Vy := RM[1, 0] * Dx + RM[1, 1] * Dy + RM[1, 2] * Dz;
 Zv := RM[2, 0] * Dx + RM[2, 1] * Dy + RM[2, 2] * Dz;
 F := FOfZ(Zv);
 Tx := CxMid + PanX + RoundI(Vx * ViewScale * F);
 Ty := CyMid + PanY - RoundI(Vy * ViewScale * F);
 If Tx > 30000 Then Tx := 30000;
 If Tx < -30000 Then Tx := -30000;
 If Ty > 30000 Then Ty := 30000;
 If Ty < -30000 Then Ty := -30000;
 Sx := Tx; Sy := Ty;
End;

{ --------------------------------------------------------------------- }
{ Preparation d'une image : atomes actifs, positions, projection        }
{ --------------------------------------------------------------------- }

Var
 VibNorm : Real;
 FrFirst, FrLast : Integer;

Function AtomShown(I: Integer): Boolean;
Begin
 AtomShown := Not Atoms^[I].Hid And
              (ShowHydro Or (Atoms^[I].ElemZ <> 1)) And
              (ShowWater Or (Atoms^[I].Kind <> ATOMK_WATER));
End;

Function InFrame(I: Integer): Boolean;
Begin
 InFrame := (I >= FrFirst) And (I <= FrLast);
End;

Procedure PrepareFrame;
Var
 I : Integer;
 Mg, Mx, Sn : Real;
Begin
 FrFirst := 0; FrLast := NAtoms - 1;
 If (CurModel >= 1) And (CurModel <= NModels) And (Models^[CurModel].First >= 0) Then Begin
  FrFirst := Models^[CurModel].First;
  FrLast := Models^[CurModel].Last;
 End;
 NAct := 0;
 For I := FrFirst To FrLast Do Begin
  Act^[NAct] := I; Inc(NAct);
  WX^[I] := Atoms^[I].X; WY^[I] := Atoms^[I].Y; WZ^[I] := Atoms^[I].Z;
 End;
 If VibOn And HasVibData Then Begin
  Mx := 0;
  For I := FrFirst To FrLast Do Begin
   Mg := Sqrt(Sqr(Vib^[I].X) + Sqr(Vib^[I].Y) + Sqr(Vib^[I].Z));
   If Mg > Mx Then Mx := Mg;
  End;
  If Mx > 1E-6 Then VibNorm := 1.0 / Mx Else VibNorm := 0;
  Sn := Sin(VibPhase) * VibScale * VibNorm;
  For I := FrFirst To FrLast Do Begin
   WX^[I] := WX^[I] + Vib^[I].X * Sn;
   WY^[I] := WY^[I] + Vib^[I].Y * Sn;
   WZ^[I] := WZ^[I] + Vib^[I].Z * Sn;
  End;
 End;
End;

Procedure ProjectAtoms;
Var
 K, I, Sx, Sy : Integer;
 Zv, F, Rad : Real;
Begin
 FillChar(PVis^, SizeOf(PVis^), 0);
 For K := 0 To NAct - 1 Do Begin
  I := Act^[K];
  ProjPoint(WX^[I], WY^[I], WZ^[I], Sx, Sy, Zv, F);
  PX^[I] := Sx; PY^[I] := Sy; PZ^[I] := Zv;
  If AtomShown(I) Then Begin
   PVis^[I] := 1;
   If Atoms^[I].Sf >= 0 Then Rad := AtomVdW(I) * Atoms^[I].Sf / 100.0
   Else Rad := -Atoms^[I].Sf / 100.0;
   If Rad > 0 Then Begin
    PR^[I] := RoundI(Rad * ViewScale * F);
    If PR^[I] < 1 Then PR^[I] := 1;
    If PR^[I] > 2000 Then PR^[I] := 2000;
   End
   Else PR^[I] := 0;
  End
  Else PR^[I] := 0;
 End;
End;

{ --------------------------------------------------------------------- }
{ Surface de points (van der Waals)                                     }
{ --------------------------------------------------------------------- }

Procedure InitUnitSphere;
Var I: Integer; Y, R, Phi: Real;
Begin
 For I := 0 To MAX_DOTS - 1 Do Begin
  Y := 1.0 - 2.0 * (I + 0.5) / MAX_DOTS;
  R := Sqrt(1.0 - Y * Y);
  Phi := I * 2.399963229728653;
  UX[I] := R * Cos(Phi); UY[I] := Y; UZ[I] := R * Sin(Phi);
 End;
End;

Procedure ComputeDotMask;
Var
 K, I, J, D, NN : Integer;
 Ri, Rj, Dx, Dy, Dz, Px, Py, Pz, Lim : Real;
 Nb : Array[0..39] Of Integer;
 Hid : Boolean;
 Mask : Array[0..DOT_BYTES-1] Of Byte;
Begin
 For K := 0 To NAct - 1 Do Begin
  I := Act^[K];
  FillChar(DotMask^[I], DOT_BYTES, 0);
  If Not Atoms^[I].Dots Then Continue;
  Ri := AtomVdW(I);
  NN := 0;
  For J := 0 To NAct - 1 Do Begin
   If Act^[J] = I Then Continue;
   Dx := Atoms^[Act^[J]].X - Atoms^[I].X;
   If AbsR(Dx) > 5 Then Continue;
   Dy := Atoms^[Act^[J]].Y - Atoms^[I].Y;
   If AbsR(Dy) > 5 Then Continue;
   Dz := Atoms^[Act^[J]].Z - Atoms^[I].Z;
   If AbsR(Dz) > 5 Then Continue;
   Rj := AtomVdW(Act^[J]);
   If Dx * Dx + Dy * Dy + Dz * Dz < Sqr(Ri + Rj) Then
    If NN < 40 Then Begin Nb[NN] := Act^[J]; Inc(NN); End;
  End;
  FillChar(Mask, DOT_BYTES, 0);
  For D := 0 To MAX_DOTS - 1 Do Begin
   Px := Atoms^[I].X + UX[D] * Ri;
   Py := Atoms^[I].Y + UY[D] * Ri;
   Pz := Atoms^[I].Z + UZ[D] * Ri;
   Hid := False;
   For J := 0 To NN - 1 Do Begin
    Lim := AtomVdW(Nb[J]);
    If Sqr(Px - Atoms^[Nb[J]].X) + Sqr(Py - Atoms^[Nb[J]].Y) +
       Sqr(Pz - Atoms^[Nb[J]].Z) < Lim * Lim Then Begin Hid := True; Break; End;
   End;
   If Not Hid Then Mask[D Shr 3] := Mask[D Shr 3] Or (1 Shl (D And 7));
  End;
  Move(Mask, DotMask^[I], DOT_BYTES);
 End;
 DotsDirty := False;
End;

{ --------------------------------------------------------------------- }
{ Primitives triees en profondeur                                       }
{ --------------------------------------------------------------------- }

Procedure AddPrim(Z: Real; Kind: Byte; A, B: Integer);
Begin
 If NPrim >= MAX_PRIM Then Exit;
 Prims^[NPrim].Key := Round(Z * 64.0);
 Prims^[NPrim].Kind := Kind;
 Prims^[NPrim].A := A;
 Prims^[NPrim].B := B;
 Inc(NPrim);
End;

Procedure SortPrims(Lo, Hi: Integer);
Var
 I, J : Integer;
 Pv : LongInt;
 T : TPrim;
Begin
 While Lo < Hi Do Begin
  I := Lo; J := Hi;
  Pv := Prims^[(LongInt(Lo) + Hi) Div 2].Key;
  Repeat
   While Prims^[I].Key < Pv Do Inc(I);
   While Prims^[J].Key > Pv Do Dec(J);
   If I <= J Then Begin
    T := Prims^[I]; Prims^[I] := Prims^[J]; Prims^[J] := T;
    Inc(I); Dec(J);
   End;
  Until I > J;
  If (J - Lo) < (Hi - I) Then Begin
   SortPrims(Lo, J);
   Lo := I;
  End
  Else Begin
   SortPrims(I, Hi);
   Hi := J;
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Rubans (squelette, trace, cartoon...)                                 }
{ --------------------------------------------------------------------- }

Function RibColorOf(R: Integer): Byte;
Begin
 If Resi^[R].RCol <> 255 Then Begin RibColorOf := Resi^[R].RCol; Exit; End;
 Case Resi^[R].SS Of
  'H': RibColorOf := Magenta;
  'E': RibColorOf := Yellow;
  'T': RibColorOf := LightBlue;
 Else
  RibColorOf := White;
 End;
End;

Function RibShown(R: Integer): Boolean;
Begin
 RibShown := (Resi^[R].Rib <> RIB_NONE) And (Resi^[R].CA >= 0) And
             InFrame(Resi^[R].CA);
End;

Procedure GuidePt(S0, N, K: Integer; Var X, Y, Z: Real);
Var A, B: Integer;
Begin
 If K < 0 Then Begin
  A := Resi^[S0].CA; B := Resi^[S0 + 1].CA;
  X := 2 * WX^[A] - WX^[B]; Y := 2 * WY^[A] - WY^[B]; Z := 2 * WZ^[A] - WZ^[B];
 End
 Else If K >= N Then Begin
  A := Resi^[S0 + N - 1].CA; B := Resi^[S0 + N - 2].CA;
  X := 2 * WX^[A] - WX^[B]; Y := 2 * WY^[A] - WY^[B]; Z := 2 * WZ^[A] - WZ^[B];
 End
 Else Begin
  A := Resi^[S0 + K].CA;
  X := WX^[A]; Y := WY^[A]; Z := WZ^[A];
 End;
End;

{ Vecteur de largeur brut du ruban au residu K du segment }
Procedure RawWidth(S0, N, K: Integer; Var Dx, Dy, Dz: Real);
Var
 R, C, O : Integer;
 X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, L : Real;
 Tx, Ty, Tz, Cx, Cy, Cz : Real;
Begin
 R := S0 + K;
 C := Resi^[R].CC; O := Resi^[R].OO;
 If Resi^[R].IsAA And (C >= 0) And (O >= 0) Then Begin
  Dx := WX^[O] - WX^[C]; Dy := WY^[O] - WY^[C]; Dz := WZ^[O] - WZ^[C];
 End
 Else Begin
  { sans atomes C et O : tangente x courbure (axe d'une helice, normale au
    plan en zigzag d'un brin) }
  GuidePt(S0, N, K - 1, X0, Y0, Z0);
  GuidePt(S0, N, K, X1, Y1, Z1);
  GuidePt(S0, N, K + 1, X2, Y2, Z2);
  Tx := X2 - X0; Ty := Y2 - Y0; Tz := Z2 - Z0;
  Cx := X0 - 2 * X1 + X2; Cy := Y0 - 2 * Y1 + Y2; Cz := Z0 - 2 * Z1 + Z2;
  Dx := Ty * Cz - Tz * Cy; Dy := Tz * Cx - Tx * Cz; Dz := Tx * Cy - Ty * Cx;
 End;
 L := Sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
 If L < 1E-6 Then Begin Dx := 0; Dy := 0; Dz := 1; End
 Else Begin Dx := Dx / L; Dy := Dy / L; Dz := Dz / L; End;
End;

{ Vecteur de largeur lisse : les plans peptidiques des brins alternent,
  on moyenne donc avec les voisins dans les feuillets }
Procedure WidthVec(S0, N, K: Integer; Var Dx, Dy, Dz: Real);
Var
 Bx, By, Bz, L : Real;
 J : Integer;
Begin
 RawWidth(S0, N, K, Dx, Dy, Dz);
 If Resi^[S0 + K].SS <> 'E' Then Exit;
 For J := -1 To 1 Do Begin
  If (J = 0) Or (K + J < 0) Or (K + J >= N) Then Continue;
  RawWidth(S0, N, K + J, Bx, By, Bz);
  If Dx * Bx + Dy * By + Dz * Bz < 0 Then Begin Bx := -Bx; By := -By; Bz := -Bz; End;
  Dx := Dx + Bx; Dy := Dy + By; Dz := Dz + Bz;
 End;
 L := Sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
 If L < 1E-6 Then Begin Dx := 0; Dy := 0; Dz := 1; End
 Else Begin Dx := Dx / L; Dy := Dy / L; Dz := Dz / L; End;
End;

Procedure PushRibPoint(X, Y, Z, Wx, Wy, Wz: Real; Col, Kind, HW: Byte;
                       R: Integer);
Begin
 If NRib >= MAX_RIB Then Exit;
 RibX^[NRib] := X; RibY^[NRib] := Y; RibZ^[NRib] := Z;
 RibWX^[NRib] := Wx; RibWY^[NRib] := Wy; RibWZ^[NRib] := Wz;
 RibCol^[NRib] := Col; RibKind^[NRib] := Kind; RibHW^[NRib] := HW;
 RibRes^[NRib] := R;
 Inc(NRib);
End;

{ Style de profil d'un residu : Kind (0 fil, 1 tube, 2 ruban plat, 3 brins),
  demi-largeur en 0.02 A }
Procedure ProfileOf(R: Integer; Var Kind, HW: Byte);
Begin
 Case Resi^[R].Rib Of
  RIB_BACKBONE: Begin Kind := 0; HW := 0; End;
  RIB_TRACE   : Begin Kind := 1; HW := 9; End;
  RIB_RIBBON  : Begin Kind := 2; HW := 40; End;
  RIB_STRAND  : Begin Kind := 3; HW := 45; End;
 Else
  { cartoon / fusees }
  If Resi^[R].SS = 'H' Then Begin Kind := 2; HW := 72; End
  Else If Resi^[R].SS = 'E' Then Begin Kind := 2; HW := 55; End
  Else Begin Kind := 1; HW := 12; End;
 End;
 If (Resi^[R].Rib = RIB_STRAND) And (Resi^[R].SS <> 'E') Then Begin Kind := 1; HW := 9; End;
End;

{ Spline de Catmull-Rom entre B et C (A et D sont les points voisins) }
Function CatRom(A, B, C, D, T: Real): Real;
Var T2, T3: Real;
Begin
 T2 := T * T; T3 := T2 * T;
 CatRom := 0.5 * (2 * B + (C - A) * T + (2 * A - 5 * B + 4 * C - D) * T2 +
                  (3 * B - A - 3 * C + D) * T3);
End;

Procedure BuildSegment(S0, N: Integer);
Var
 K, S, NS : Integer;
 T, Pxx, Pyy, Pzz, Tx, Ty, Tz, L, Dt : Real;
 P0, P1, P2, P3 : TVec3;
 Dk, Dk1, Wv : TVec3;
 R : Integer;
 Kind, HW : Byte;
 Dp : Real;
 Arrow : Boolean;
Begin
 If N < 2 Then Exit;
 NS := 8;
 If NRes > 300 Then NS := 6;
 If NRes > 600 Then NS := 4;
 If NRes > 1000 Then NS := 3;
 If NRes > 2500 Then NS := 2;
 WidthVec(S0, N, 0, Dk.X, Dk.Y, Dk.Z);
 For K := 0 To N - 2 Do Begin
  GuidePt(S0, N, K - 1, P0.X, P0.Y, P0.Z);
  GuidePt(S0, N, K, P1.X, P1.Y, P1.Z);
  GuidePt(S0, N, K + 1, P2.X, P2.Y, P2.Z);
  GuidePt(S0, N, K + 2, P3.X, P3.Y, P3.Z);
  WidthVec(S0, N, K + 1, Dk1.X, Dk1.Y, Dk1.Z);
  If Dk.X * Dk1.X + Dk.Y * Dk1.Y + Dk.Z * Dk1.Z < 0 Then Begin
   Dk1.X := -Dk1.X; Dk1.Y := -Dk1.Y; Dk1.Z := -Dk1.Z;
  End;
  Arrow := (Resi^[S0 + K].Rib = RIB_CARTOON) And (Resi^[S0 + K].SS = 'E') And
           (Resi^[S0 + K + 1].SS <> 'E');
  If Resi^[S0 + K].Rib = RIB_RIBBON Then Arrow := False;
  For S := 0 To NS - 1 Do Begin
   T := S / NS;
   Pxx := CatRom(P0.X, P1.X, P2.X, P3.X, T);
   Pyy := CatRom(P0.Y, P1.Y, P2.Y, P3.Y, T);
   Pzz := CatRom(P0.Z, P1.Z, P2.Z, P3.Z, T);
   Tx := P2.X - P0.X; Ty := P2.Y - P0.Y; Tz := P2.Z - P0.Z;
   L := Sqrt(Tx * Tx + Ty * Ty + Tz * Tz);
   If L < 1E-6 Then L := 1;
   Tx := Tx / L; Ty := Ty / L; Tz := Tz / L;
   Wv.X := Dk.X + (Dk1.X - Dk.X) * T;
   Wv.Y := Dk.Y + (Dk1.Y - Dk.Y) * T;
   Wv.Z := Dk.Z + (Dk1.Z - Dk.Z) * T;
   Dp := Wv.X * Tx + Wv.Y * Ty + Wv.Z * Tz;
   Wv.X := Wv.X - Dp * Tx; Wv.Y := Wv.Y - Dp * Ty; Wv.Z := Wv.Z - Dp * Tz;
   L := Sqrt(Wv.X * Wv.X + Wv.Y * Wv.Y + Wv.Z * Wv.Z);
   If L < 1E-6 Then Begin Wv.X := 0; Wv.Y := 1; Wv.Z := 0; End
   Else Begin Wv.X := Wv.X / L; Wv.Y := Wv.Y / L; Wv.Z := Wv.Z / L; End;
   If T < 0.5 Then R := S0 + K Else R := S0 + K + 1;
   ProfileOf(R, Kind, HW);
   If Arrow Then Begin
    Dt := 110 * (1.0 - T);
    HW := RoundI(Dt);
    Kind := 2;
   End;
   PushRibPoint(Pxx, Pyy, Pzz, Wv.X, Wv.Y, Wv.Z, RibColorOf(R), Kind, HW, R);
   If (S = 0) And (K = 0) Then RibKind^[NRib - 1] := RibKind^[NRib - 1] Or 128;
  End;
  If Arrow Then Begin
   { pointe de la fleche }
   PushRibPoint(P2.X, P2.Y, P2.Z, Dk1.X, Dk1.Y, Dk1.Z, RibColorOf(S0 + K), 2, 0, S0 + K);
  End;
  Dk := Dk1;
 End;
 { dernier point }
 R := S0 + N - 1;
 ProfileOf(R, Kind, HW);
 GuidePt(S0, N, N - 1, P1.X, P1.Y, P1.Z);
 PushRibPoint(P1.X, P1.Y, P1.Z, Dk.X, Dk.Y, Dk.Z, RibColorOf(R), Kind, HW, R);
End;

Procedure BuildRibbons;
Var
 R, S0, N : Integer;
 Any : Boolean;
Begin
 NRib := 0;
 Any := False;
 For R := 0 To NRes - 1 Do If Resi^[R].Rib <> RIB_NONE Then Begin Any := True; Break; End;
 If Not Any Then Exit;
 R := 0;
 While R < NRes Do Begin
  If Not RibShown(R) Then Begin Inc(R); Continue; End;
  S0 := R; N := 1;
  While (R + 1 < NRes) And RibShown(R + 1) And ResLinked(R, R + 1) Do Begin
   Inc(R); Inc(N);
  End;
  Inc(R);
  If N >= 2 Then BuildSegment(S0, N);
 End;
End;

{ Transforme les points de ruban en primitives }
Procedure RibbonPrims;
Var
 I, A, B, K : Integer;
 Kind : Byte;
 Hwa, Hwb, Sh, Nx, Ny, Nz, Tx, Ty, Tz, Mx, My, Mz, L, Vx, Vy, Vz : Real;
 Zv, F, Za : Real;
 Qi, Sx, Sy : Integer;
 Col : Byte;
Begin
 If NRib < 2 Then Exit;
 For I := 0 To NRib - 1 Do Begin
  ProjPoint(RibX^[I], RibY^[I], RibZ^[I], Sx, Sy, Zv, F);
  RibSX^[I] := Sx; RibSY^[I] := Sy; RibV^[I] := Zv; RibF^[I] := F;
 End;
 For B := 1 To NRib - 1 Do Begin
  If (RibKind^[B] And 128) <> 0 Then Continue;
  A := B - 1;
  Kind := RibKind^[B] And 127;
  Col := RibCol^[B];
  If NQuad >= MAX_QUAD - 3 Then Exit;
  Za := (RibV^[A] + RibV^[B]) * 0.5;
  F := (RibF^[A] + RibF^[B]) * 0.5;
  Case Kind Of
   0: Begin
    Qi := NQuad; Inc(NQuad);
    Quads^[Qi].X[0] := RibSX^[A]; Quads^[Qi].Y[0] := RibSY^[A];
    Quads^[Qi].X[1] := RibSX^[B]; Quads^[Qi].Y[1] := RibSY^[B];
    Quads^[Qi].X[2] := 1; Quads^[Qi].Col := Col;
    AddPrim(Za, PRIM_TUBE, Qi, 0);
   End;
   1: Begin
    Qi := NQuad; Inc(NQuad);
    Quads^[Qi].X[0] := RibSX^[A]; Quads^[Qi].Y[0] := RibSY^[A];
    Quads^[Qi].X[1] := RibSX^[B]; Quads^[Qi].Y[1] := RibSY^[B];
    Quads^[Qi].X[2] := MaxI(2, RoundI(2.0 * RibHW^[B] * 0.02 * ViewScale * F));
    Quads^[Qi].Col := Col;
    AddPrim(Za, PRIM_TUBE, Qi, 0);
   End;
   3: Begin
    For K := -1 To 1 Do Begin
     Qi := NQuad; Inc(NQuad);
     Hwa := K * RibHW^[B] * 0.02;
     ProjPoint(RibX^[A] + RibWX^[A] * Hwa, RibY^[A] + RibWY^[A] * Hwa,
               RibZ^[A] + RibWZ^[A] * Hwa, Sx, Sy, Zv, F);
     Quads^[Qi].X[0] := Sx; Quads^[Qi].Y[0] := Sy;
     ProjPoint(RibX^[B] + RibWX^[B] * Hwa, RibY^[B] + RibWY^[B] * Hwa,
               RibZ^[B] + RibWZ^[B] * Hwa, Sx, Sy, Zv, F);
     Quads^[Qi].X[1] := Sx; Quads^[Qi].Y[1] := Sy;
     Quads^[Qi].X[2] := 2; Quads^[Qi].Col := Col;
     AddPrim(Za, PRIM_TUBE, Qi, 0);
    End;
   End;
  Else
   Begin
    { ruban plat }
    Hwa := RibHW^[A] * 0.02; Hwb := RibHW^[B] * 0.02;
    Qi := NQuad; Inc(NQuad);
    ProjPoint(RibX^[A] + RibWX^[A] * Hwa, RibY^[A] + RibWY^[A] * Hwa,
              RibZ^[A] + RibWZ^[A] * Hwa, Sx, Sy, Zv, F);
    Quads^[Qi].X[0] := Sx; Quads^[Qi].Y[0] := Sy;
    ProjPoint(RibX^[B] + RibWX^[B] * Hwb, RibY^[B] + RibWY^[B] * Hwb,
              RibZ^[B] + RibWZ^[B] * Hwb, Sx, Sy, Zv, F);
    Quads^[Qi].X[1] := Sx; Quads^[Qi].Y[1] := Sy;
    ProjPoint(RibX^[B] - RibWX^[B] * Hwb, RibY^[B] - RibWY^[B] * Hwb,
              RibZ^[B] - RibWZ^[B] * Hwb, Sx, Sy, Zv, F);
    Quads^[Qi].X[2] := Sx; Quads^[Qi].Y[2] := Sy;
    ProjPoint(RibX^[A] - RibWX^[A] * Hwa, RibY^[A] - RibWY^[A] * Hwa,
              RibZ^[A] - RibWZ^[A] * Hwa, Sx, Sy, Zv, F);
    Quads^[Qi].X[3] := Sx; Quads^[Qi].Y[3] := Sy;
    { normale = T x W (en coordonnees de la molecule) puis passage en vue }
    Tx := RibX^[B] - RibX^[A]; Ty := RibY^[B] - RibY^[A]; Tz := RibZ^[B] - RibZ^[A];
    Mx := RibWX^[A] + RibWX^[B]; My := RibWY^[A] + RibWY^[B]; Mz := RibWZ^[A] + RibWZ^[B];
    Nx := Ty * Mz - Tz * My; Ny := Tz * Mx - Tx * Mz; Nz := Tx * My - Ty * Mx;
    L := Sqrt(Nx * Nx + Ny * Ny + Nz * Nz);
    If L < 1E-9 Then Sh := 0.5
    Else Begin
     Nx := Nx / L; Ny := Ny / L; Nz := Nz / L;
     Vx := RM[0, 0] * Nx + RM[0, 1] * Ny + RM[0, 2] * Nz;
     Vy := RM[1, 0] * Nx + RM[1, 1] * Ny + RM[1, 2] * Nz;
     Vz := RM[2, 0] * Nx + RM[2, 1] * Ny + RM[2, 2] * Nz;
     Sh := AbsR(-0.41 * Vx + 0.51 * Vy + 0.76 * Vz);
    End;
    Quads^[Qi].Dk := ShadeDark[Col And 15];
    If Sh < 0.25 Then Quads^[Qi].Col := ShadeDark[Col And 15]
    Else If Sh < 0.62 Then Quads^[Qi].Col := Col
    Else Quads^[Qi].Col := ShadeLight[Col And 15];
    If Quads^[Qi].Col = Quads^[Qi].Dk Then Quads^[Qi].Dk := Col;
    AddPrim(Za, PRIM_QUAD, Qi, 0);
   End;
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Construction de la liste de primitives                                }
{ --------------------------------------------------------------------- }

Procedure BuildPrims;
Var
 K, I, B, A1, A2 : Integer;
Begin
 NPrim := 0; NQuad := 0;
 For K := 0 To NAct - 1 Do Begin
  I := Act^[K];
  If PVis^[I] = 0 Then Continue;
  If (PR^[I] > 0) Or (ShowHalo And Atoms^[I].Sel) Then
   AddPrim(PZ^[I], PRIM_ATOM, I, 0);
  If Atoms^[I].Dots Then AddPrim(PZ^[I], PRIM_DOTS, I, 0);
 End;
 For B := 0 To NBonds - 1 Do Begin
  If Bonds^[B].BR = 0 Then Continue;
  A1 := Bonds^[B].A; A2 := Bonds^[B].B;
  If Not (InFrame(A1) And InFrame(A2)) Then Continue;
  If (PVis^[A1] = 0) Or (PVis^[A2] = 0) Then Continue;
  AddPrim(PZ^[A1] + (PZ^[A2] - PZ^[A1]) * 0.25, PRIM_BOND, B, 0);
  AddPrim(PZ^[A2] + (PZ^[A1] - PZ^[A2]) * 0.25, PRIM_BOND, B, 1);
 End;
 If ShowHB Then
  For B := 0 To NHB - 1 Do Begin
   A1 := HBs^[B].A; A2 := HBs^[B].B;
   If Not (InFrame(A1) And InFrame(A2)) Then Continue;
   If (PVis^[A1] = 0) Or (PVis^[A2] = 0) Then Continue;
   AddPrim((PZ^[A1] + PZ^[A2]) * 0.5, PRIM_HBOND, B, 0);
  End;
 If ShowVec And HasVibData Then
  For K := 0 To NAct - 1 Do Begin
   I := Act^[K];
   If PVis^[I] = 0 Then Continue;
   If (Vib^[I].X <> 0) Or (Vib^[I].Y <> 0) Or (Vib^[I].Z <> 0) Then
    AddPrim(PZ^[I] + 0.5, PRIM_VECTOR, I, 0);
  End;
 BuildRibbons;
 RibbonPrims;
 If NPrim > 1 Then SortPrims(0, NPrim - 1);
End;

{ --------------------------------------------------------------------- }
{ Dessin des primitives                                                 }
{ --------------------------------------------------------------------- }

Procedure DrawAtomPrim(I: Integer);
Var R: Integer;
Begin
 R := PR^[I];
 If ShowHalo And Atoms^[I].Sel Then Begin
  SetColor(Yellow); SetFillStyle(SolidFill, Yellow);
  If R < 3 Then FillEllipse(PX^[I], PY^[I], 5, 5)
  Else FillEllipse(PX^[I], PY^[I], R + 3, R + 3);
 End;
 If R > 0 Then DrawBall(PX^[I], PY^[I], R, AtomColor(I));
End;

Procedure DrawBondPrim(B, H: Integer);
Var
 A1, A2, X1, Y1, X2, Y2, W, N, I, Ord : Integer;
 Dx, Dy, Len, Nx, Ny, Off, Fz, Sp : Real;
 Col : Byte;
 Ox, Oy : Integer;
Begin
 If H = 0 Then Begin A1 := Bonds^[B].A; A2 := Bonds^[B].B; End
 Else Begin A1 := Bonds^[B].B; A2 := Bonds^[B].A; End;
 X1 := PX^[A1]; Y1 := PY^[A1];
 X2 := (PX^[A1] + PX^[A2]) Div 2; Y2 := (PY^[A1] + PY^[A2]) Div 2;
 Col := AtomColor(A1);
 Ord := Bonds^[B].Ord;
 If Ord > 3 Then Ord := 1;
 If Ord < 1 Then Ord := 1;
 Dx := X2 - X1; Dy := Y2 - Y1;
 Len := Sqrt(Dx * Dx + Dy * Dy);
 If Len < 0.5 Then Exit;
 Nx := -Dy / Len; Ny := Dx / Len;
 Fz := FOfZ((PZ^[A1] + PZ^[A2]) * 0.5);
 If Bonds^[B].BR <= 1 Then Begin
  { fil de fer }
  Sp := 3.0;
  For I := 0 To Ord - 1 Do Begin
   Off := (I - (Ord - 1) * 0.5) * Sp;
   Ox := RoundI(Nx * Off); Oy := RoundI(Ny * Off);
   DrawStick(X1 + Ox, Y1 + Oy, X2 + Ox, Y2 + Oy, 1, Col);
  End;
  Exit;
 End;
 W := RoundI(Bonds^[B].BR * 0.02 * ViewScale * Fz);
 If W < 2 Then W := 2;
 If Ord = 1 Then DrawStick(X1, Y1, X2, Y2, W, Col)
 Else Begin
  N := Ord;
  Sp := W * 0.85;
  W := RoundI(W * (1.5 / N));
  If W < 2 Then W := 2;
  For I := 0 To N - 1 Do Begin
   Off := (I - (N - 1) * 0.5) * Sp;
   Ox := RoundI(Nx * Off); Oy := RoundI(Ny * Off);
   DrawStick(X1 + Ox, Y1 + Oy, X2 + Ox, Y2 + Oy, W, Col);
  End;
 End;
End;

Procedure DrawDotsPrim(I: Integer);
Var
 D, Sx, Sy : Integer;
 Vx, Vy, Vz, R, F : Real;
 Col, Dk : Byte;
Begin
 Col := AtomColor(I);
 Dk := ShadeDark[Col And 15];
 If Dk = BgCol Then Dk := Col;
 R := AtomVdW(I);
 F := FOfZ(PZ^[I]);
 For D := 0 To MAX_DOTS - 1 Do Begin
  If (DotMask^[I][D Shr 3] And (1 Shl (D And 7))) = 0 Then Continue;
  Vx := RM[0, 0] * UX[D] + RM[0, 1] * UY[D] + RM[0, 2] * UZ[D];
  Vy := RM[1, 0] * UX[D] + RM[1, 1] * UY[D] + RM[1, 2] * UZ[D];
  Vz := RM[2, 0] * UX[D] + RM[2, 1] * UY[D] + RM[2, 2] * UZ[D];
  Sx := PX^[I] + RoundI(Vx * R * ViewScale * F);
  Sy := PY^[I] - RoundI(Vy * R * ViewScale * F);
  If Vz >= 0 Then PutPixel(Sx, Sy, Col) Else PutPixel(Sx, Sy, Dk);
 End;
End;

Procedure DrawPrim(Var P: TPrim);
Var
 A, B, W : Integer;
 Vx, Vy, Sc : Real;
Begin
 Case P.Kind Of
  PRIM_ATOM : DrawAtomPrim(P.A);
  PRIM_BOND : DrawBondPrim(P.A, P.B);
  PRIM_QUAD : DrawQuad(Quads^[P.A]);
  PRIM_TUBE : Begin
   W := Quads^[P.A].X[2];
   DrawStick(Quads^[P.A].X[0], Quads^[P.A].Y[0], Quads^[P.A].X[1],
             Quads^[P.A].Y[1], W, Quads^[P.A].Col);
   If W >= 5 Then Begin
    SetColor(Quads^[P.A].Col); SetFillStyle(SolidFill, Quads^[P.A].Col);
    FillEllipse(Quads^[P.A].X[1], Quads^[P.A].Y[1], W Div 2 - 1, W Div 2 - 1);
   End;
  End;
  PRIM_DOTS : DrawDotsPrim(P.A);
  PRIM_HBOND: Begin
   A := HBs^[P.A].A; B := HBs^[P.A].B;
   DashLine(PX^[A], PY^[A], PX^[B], PY^[B], LightCyan, 4);
  End;
  PRIM_VECTOR: Begin
   A := P.A;
   Sc := VibNorm * VecScale * ViewScale;
   Vx := RM[0, 0] * Vib^[A].X + RM[0, 1] * Vib^[A].Y + RM[0, 2] * Vib^[A].Z;
   Vy := RM[1, 0] * Vib^[A].X + RM[1, 1] * Vib^[A].Y + RM[1, 2] * Vib^[A].Z;
   DrawArrow(PX^[A], PY^[A], PX^[A] + RoundI(Vx * Sc), PY^[A] - RoundI(Vy * Sc), FgCol);
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Mesures                                                               }
{ --------------------------------------------------------------------- }

Function DistW(I, J: Integer): Real;
Begin
 DistW := Sqrt(Sqr(WX^[I] - WX^[J]) + Sqr(WY^[I] - WY^[J]) + Sqr(WZ^[I] - WZ^[J]));
End;

Function AngleW(I, J, K: Integer): Real;
Var
 Ax, Ay, Az, Bx, By, Bz, La, Lb, C: Real;
Begin
 Ax := WX^[I] - WX^[J]; Ay := WY^[I] - WY^[J]; Az := WZ^[I] - WZ^[J];
 Bx := WX^[K] - WX^[J]; By := WY^[K] - WY^[J]; Bz := WZ^[K] - WZ^[J];
 La := Sqrt(Ax * Ax + Ay * Ay + Az * Az);
 Lb := Sqrt(Bx * Bx + By * By + Bz * Bz);
 If (La < 1E-9) Or (Lb < 1E-9) Then Begin AngleW := 0; Exit; End;
 C := (Ax * Bx + Ay * By + Az * Bz) / (La * Lb);
 If C > 1 Then C := 1;
 If C < -1 Then C := -1;
 If C < -1E-12 Then AngleW := 180.0 - ArcTan(Sqrt(1 - C * C) / (-C)) / DEG2RAD
 Else If C > 1E-12 Then AngleW := ArcTan(Sqrt(1 - C * C) / C) / DEG2RAD
 Else AngleW := 90.0;
End;

Function TorsionW(I, J, K, L: Integer): Real;
Var
 B1x, B1y, B1z, B2x, B2y, B2z, B3x, B3y, B3z : Real;
 N1x, N1y, N1z, N2x, N2y, N2z, M1x, M1y, M1z, X, Y, Lb : Real;
Begin
 B1x := WX^[J] - WX^[I]; B1y := WY^[J] - WY^[I]; B1z := WZ^[J] - WZ^[I];
 B2x := WX^[K] - WX^[J]; B2y := WY^[K] - WY^[J]; B2z := WZ^[K] - WZ^[J];
 B3x := WX^[L] - WX^[K]; B3y := WY^[L] - WY^[K]; B3z := WZ^[L] - WZ^[K];
 N1x := B1y * B2z - B1z * B2y; N1y := B1z * B2x - B1x * B2z; N1z := B1x * B2y - B1y * B2x;
 N2x := B2y * B3z - B2z * B3y; N2y := B2z * B3x - B2x * B3z; N2z := B2x * B3y - B2y * B3x;
 Lb := Sqrt(B2x * B2x + B2y * B2y + B2z * B2z);
 If Lb < 1E-9 Then Begin TorsionW := 0; Exit; End;
 M1x := N1y * (B2z / Lb) - N1z * (B2y / Lb);
 M1y := N1z * (B2x / Lb) - N1x * (B2z / Lb);
 M1z := N1x * (B2y / Lb) - N1y * (B2x / Lb);
 X := N1x * N2x + N1y * N2y + N1z * N2z;
 Y := M1x * N2x + M1y * N2y + M1z * N2z;
 If (X = 0) And (Y = 0) Then Begin TorsionW := 0; Exit; End;
 { atan2(-Y, X) }
 If X > 0 Then TorsionW := ArcTan(-Y / X) / DEG2RAD
 Else If X < 0 Then Begin
  If -Y >= 0 Then TorsionW := ArcTan(-Y / X) / DEG2RAD + 180.0
  Else TorsionW := ArcTan(-Y / X) / DEG2RAD - 180.0;
 End
 Else If -Y > 0 Then TorsionW := 90.0 Else TorsionW := -90.0;
End;

Function UnitStr: String;
Begin
 Case MeasUnit Of
  1: UnitStr := ' nm';
  2: UnitStr := ' pm';
 Else
  UnitStr := ' ' + #143;
 End;
End;

Function DistStr(D: Real): String;
Begin
 Case MeasUnit Of
  1: DistStr := RealToStr(D / 10.0, 4) + ' nm';
  2: DistStr := RealToStr(D * 100.0, 1) + ' pm';
 Else
  DistStr := RealToStr(D, 3) + ' ' + #143;
 End;
End;

Function MeasValueStr(Const M: TMeas): String;
Begin
 Case M.N Of
  2: MeasValueStr := DistStr(DistW(M.A[0], M.A[1]));
  3: MeasValueStr := RealToStr(AngleW(M.A[0], M.A[1], M.A[2]), 1) + #248;
 Else
  MeasValueStr := RealToStr(TorsionW(M.A[0], M.A[1], M.A[2], M.A[3]), 1) + #248;
 End;
End;

Function MeasInFrame(Const M: TMeas): Boolean;
Var I: Integer;
Begin
 MeasInFrame := False;
 For I := 0 To M.N - 1 Do
  If (M.A[I] < 0) Or (M.A[I] >= NAtoms) Or Not InFrame(M.A[I]) Then Exit;
 MeasInFrame := True;
End;

Procedure DrawMeasures;
Var
 I, J, X, Y : Integer;
 S : String;
 C : Byte;
Begin
 C := Yellow;
 If (BgCol = Yellow) Or (BgCol = White) Then C := Blue;
 For I := 1 To NMeas Do Begin
  If Not MeasInFrame(Meas[I]) Then Continue;
  For J := 0 To Meas[I].N - 2 Do
   DashLine(PX^[Meas[I].A[J]], PY^[Meas[I].A[J]],
            PX^[Meas[I].A[J + 1]], PY^[Meas[I].A[J + 1]], C, 4);
  S := MeasValueStr(Meas[I]);
  If Meas[I].N = 2 Then Begin
   X := (PX^[Meas[I].A[0]] + PX^[Meas[I].A[1]]) Div 2;
   Y := (PY^[Meas[I].A[0]] + PY^[Meas[I].A[1]]) Div 2;
  End
  Else If Meas[I].N = 3 Then Begin
   X := PX^[Meas[I].A[1]]; Y := PY^[Meas[I].A[1]];
  End
  Else Begin
   X := (PX^[Meas[I].A[1]] + PX^[Meas[I].A[2]]) Div 2;
   Y := (PY^[Meas[I].A[1]] + PY^[Meas[I].A[2]]) Div 2;
  End;
  SetColor(C);
  SetTextStyle(DefaultFont, HorizDir, 1);
  OutTextXY(X + 3, Y - 10, S);
 End;
End;

{ --------------------------------------------------------------------- }
{ Etiquettes                                                            }
{ --------------------------------------------------------------------- }

Function LabelText(I: Integer; Const Fmt: String): String;
Var
 S : String;
 K : Integer;
Begin
 S := '';
 K := 1;
 While K <= Length(Fmt) Do Begin
  If (Fmt[K] = '%') And (K < Length(Fmt)) Then Begin
   Inc(K);
   Case Fmt[K] Of
    'a': S := S + TrimStr(Atoms^[I].Name);
    'e': S := S + ElemDB[Atoms^[I].ElemZ].Sym;
    'i': S := S + IntToStr(Atoms^[I].Serial);
    'r': S := S + TrimStr(Atoms^[I].Res);
    'n': S := S + IntToStr(Atoms^[I].ResSeq);
    'c': S := S + Atoms^[I].Chain;
    'x': S := S + RealToStr(Atoms^[I].X, 2);
    'y': S := S + RealToStr(Atoms^[I].Y, 2);
    'z': S := S + RealToStr(Atoms^[I].Z, 2);
    '%': S := S + '%';
   End;
  End
  Else S := S + Fmt[K];
  Inc(K);
 End;
 LabelText := S;
End;

Procedure DrawLabels;
Var
 K, I : Integer;
 S : String;
 C : Byte;
Begin
 C := FgCol;
 SetTextStyle(DefaultFont, HorizDir, 1);
 For K := 0 To NAct - 1 Do Begin
  I := Act^[K];
  If (PVis^[I] = 0) Or (Atoms^[I].Lbl = 0) Or (Atoms^[I].Lbl > NLblFmt) Then Continue;
  S := LabelText(I, LblFmt[Atoms^[I].Lbl]);
  SetColor(C);
  OutTextXY(PX^[I] + PR^[I] + 2, PY^[I] - 4, S);
 End;
End;

{ --------------------------------------------------------------------- }
{ Axes, boite englobante, cellule unitaire                              }
{ --------------------------------------------------------------------- }

Procedure Line3D(X1, Y1, Z1, X2, Y2, Z2: Real; C: Byte; Dashed: Boolean);
Var
 Sx1, Sy1, Sx2, Sy2 : Integer;
 Zv, F : Real;
Begin
 ProjPoint(X1, Y1, Z1, Sx1, Sy1, Zv, F);
 ProjPoint(X2, Y2, Z2, Sx2, Sy2, Zv, F);
 SetColor(C);
 SetLineStyle(SolidLn, 0, NormWidth);
 If Dashed Then DashLine(Sx1, Sy1, Sx2, Sy2, C, 4)
 Else Line(Sx1, Sy1, Sx2, Sy2);
End;

Procedure DrawBoxFrame(Var P: TBoxPts; C: Byte);
Var I, J: Integer;
Begin
 For I := 0 To 7 Do
  For J := 0 To 2 Do
   If (I And (1 Shl J)) = 0 Then
    Line3D(P[I].X, P[I].Y, P[I].Z, P[I Or (1 Shl J)].X, P[I Or (1 Shl J)].Y,
           P[I Or (1 Shl J)].Z, C, False);
End;

Procedure DrawOverlays3D;
Var
 I, K, Sx, Sy : Integer;
 Mn, Mx : TVec3;
 B : TBoxPts;
 L, Zv, F : Real;
Begin
 If ShowBBox And (NAct > 0) Then Begin
  Mn.X := 1E9; Mn.Y := 1E9; Mn.Z := 1E9; Mx.X := -1E9; Mx.Y := -1E9; Mx.Z := -1E9;
  For K := 0 To NAct - 1 Do Begin
   I := Act^[K];
   If WX^[I] < Mn.X Then Mn.X := WX^[I]; If WX^[I] > Mx.X Then Mx.X := WX^[I];
   If WY^[I] < Mn.Y Then Mn.Y := WY^[I]; If WY^[I] > Mx.Y Then Mx.Y := WY^[I];
   If WZ^[I] < Mn.Z Then Mn.Z := WZ^[I]; If WZ^[I] > Mx.Z Then Mx.Z := WZ^[I];
  End;
  For I := 0 To 7 Do Begin
   If (I And 1) = 0 Then B[I].X := Mn.X Else B[I].X := Mx.X;
   If (I And 2) = 0 Then B[I].Y := Mn.Y Else B[I].Y := Mx.Y;
   If (I And 4) = 0 Then B[I].Z := Mn.Z Else B[I].Z := Mx.Z;
  End;
  DrawBoxFrame(B, Brown);
 End;
 If ShowCell And HasCell Then Begin
  DrawBoxFrame(CellPts, LightBlue);
 End;
 If ShowAxes Then Begin
  L := MolRad * 0.9;
  Line3D(CenX, CenY, CenZ, CenX + L, CenY, CenZ, LightRed, False);
  Line3D(CenX, CenY, CenZ, CenX, CenY + L, CenZ, LightGreen, False);
  Line3D(CenX, CenY, CenZ, CenX, CenY, CenZ + L, LightBlue, False);
  SetTextStyle(DefaultFont, HorizDir, 1);
  ProjPoint(CenX + L, CenY, CenZ, Sx, Sy, Zv, F);
  Txt(Sx + 3, Sy - 4, 'X', LightRed);
  ProjPoint(CenX, CenY + L, CenZ, Sx, Sy, Zv, F);
  Txt(Sx + 3, Sy - 4, 'Y', LightGreen);
  ProjPoint(CenX, CenY, CenZ + L, Sx, Sy, Zv, F);
  Txt(Sx + 3, Sy - 4, 'Z', LightBlue);
 End;
End;

{ --------------------------------------------------------------------- }
{ Rendu complet de la zone 3D                                           }
{ --------------------------------------------------------------------- }

Procedure DrawPickMarks;
Var I, K: Integer;
Begin
 For K := 0 To PickN - 1 Do Begin
  I := PickQ[K];
  If (I >= 0) And (I < NAtoms) And InFrame(I) And (PVis^[I] <> 0) Then Begin
   SetColor(LightGreen);
   SetLineStyle(SolidLn, 0, NormWidth);
   Circle(PX^[I], PY^[I], PR^[I] + 5);
   Circle(PX^[I], PY^[I], PR^[I] + 6);
  End;
 End;
End;

Procedure RenderCanvas;
Var
 K : Integer;
 T0 : LongInt;
Begin
 T0 := TickMs;
 HideMouse;
 SetViewPort(CvX1, CvY1, CvX2, CvY2, True);
 SetFillStyle(SolidFill, BgCol);
 Bar(0, 0, CvW - 1, CvH - 1);
 SetLineStyle(SolidLn, 0, NormWidth);
 If Loaded And (NAtoms > 0) Then Begin
  SetupView;
  PrepareFrame;
  ProjectAtoms;
  If DotsDirty Then ComputeDotMask;
  DrawOverlays3D;
  BuildPrims;
  For K := 0 To NPrim - 1 Do DrawPrim(Prims^[K]);
  DrawLabels;
  DrawMeasures;
  DrawPickMarks;
 End
 Else Begin
  SetColor(DarkGray);
  SetTextStyle(DefaultFont, HorizDir, 1);
  OutTextXY((CvW - 8 * 40) Div 2, CvH Div 2 - 10, 'Aucun fichier charge.');
  OutTextXY((CvW - 8 * 40) Div 2, CvH Div 2 + 4, 'Fichier > Ouvrir... (Ctrl+O) pour charger une molecule.');
 End;
 SetColor(DarkGray);
 SetTextStyle(DefaultFont, HorizDir, 2);
 OutTextXY(CvW - 66, CvH - 22, PROG_NAME);
 SetTextStyle(DefaultFont, HorizDir, 1);
 If VibOn And HasVibData And (CurModel >= 1) And (Models^[CurModel].Freq <> 0) Then Begin
  SetColor(FgCol);
  OutTextXY(6, 6, 'Vibration : ' + RealToStr(Models^[CurModel].Freq, 2) + ' cm^-1');
 End;
 SetViewPort(0, 0, SW - 1, SH - 1, True);
 ShowMouse;
 RenderMs := ElapsedMs(T0);
End;

{ --------------------------------------------------------------------- }
{ Console (tampon de lignes)                                            }
{ --------------------------------------------------------------------- }

Procedure ConsolePush(Const S: String);
Var I: Integer;
Begin
 If ConCount >= MAX_CONLN Then Begin
  For I := 0 To MAX_CONLN - 2 Do ConLines[I] := ConLines[I + 1];
  ConCount := MAX_CONLN - 1;
 End;
 ConLines[ConCount] := Copy(S, 1, 80);
 Inc(ConCount);
 NeedChrome := True;
End;

Procedure ConsoleSay(Const S: String);
Var T: String;
Begin
 T := S;
 While Length(T) > 76 Do Begin
  ConsolePush(' ' + Copy(T, 1, 76));
  Delete(T, 1, 76);
 End;
 ConsolePush(' ' + T);
 StatMsg := S;
 NeedStatus := True;
End;

Procedure ConsoleErr(Const S: String);
Begin
 ConsolePush('!' + Copy('Erreur : ' + S, 1, 79));
 StatMsg := 'Erreur : ' + S;
 NeedStatus := True;
End;

Procedure ConsoleEcho(Const S: String);
Begin
 ConsolePush('>' + Copy(S, 1, 79));
End;

{ --------------------------------------------------------------------- }
{ Entrees clavier / souris pour les boites de dialogue                  }
{ --------------------------------------------------------------------- }

Procedure FlushKeys;
Var K: Char;
Begin
 While KeyPressed Do Begin
  K := ReadKey;
  If K = #0 Then If KeyPressed Then K := ReadKey;
 End;
End;

{ Lit une touche : Ext=True pour les touches etendues (fleches, F1...) }
Function GetKey(Var Ext: Boolean): Char;
Var K: Char;
Begin
 Ext := False;
 K := ReadKey;
 If K = #0 Then Begin
  Ext := True;
  If KeyPressed Then K := ReadKey Else K := #0;
 End;
 GetKey := K;
End;

Procedure WaitRelease;
Begin
 Repeat
  PollMouse;
  Delay(10);
 Until MouseB = 0;
 OldMB := 0;
End;

{ --------------------------------------------------------------------- }
{ Dialogues : cadre, boutons, messages                                  }
{ --------------------------------------------------------------------- }

Procedure DrawDialogBox(X1, Y1, X2, Y2: Integer; Const Title: String);
Begin
 PanelBar(X1, Y1, X2, Y2, LightGray);
 Frame3D(X1, Y1, X2, Y2, True);
 Frame3D(X1 + 1, Y1 + 1, X2 - 1, Y2 - 1, True);
 PanelBar(X1 + 3, Y1 + 3, X2 - 3, Y1 + 17, Blue);
 Txt(X1 + 8, Y1 + 6, Title, White);
 SetColor(Black);
End;

Procedure DrawButton(X1, Y1, X2, Y2: Integer; Const Caption: String; Focused: Boolean);
Var TX, TY: Integer;
Begin
 Button3D(X1, Y1, X2, Y2, False);
 TX := (X1 + X2 - TextWidth(Caption)) Div 2;
 TY := (Y1 + Y2 - TextHeight(Caption)) Div 2;
 Txt(TX, TY, Caption, Black);
 If Focused Then Begin
  SetColor(Black);
  SetLineStyle(DottedLn, 0, NormWidth);
  Rectangle(X1 + 3, Y1 + 3, X2 - 3, Y2 - 3);
  SetLineStyle(SolidLn, 0, NormWidth);
 End;
End;

{ Message multi-lignes : les lignes sont separees par '|' }
Procedure ShowMessage(Const Title, Msg: String);
Const BW = 74; BH = 22;
Var
 X1, Y1, X2, Y2, BX1, BY1, W, H, I, N, Y : Integer;
 K : Char;
 Ext, Done : Boolean;
 Ln : Array[0..9] Of String[78];
 Cur : String;
Begin
 N := 0; Cur := '';
 For I := 1 To Length(Msg) + 1 Do Begin
  If (I > Length(Msg)) Or (Msg[I] = '|') Then Begin
   If N < 10 Then Begin Ln[N] := Copy(Cur, 1, 78); Inc(N); End;
   Cur := '';
  End
  Else Cur := Cur + Msg[I];
 End;
 W := TextWidth(Title) + 40;
 For I := 0 To N - 1 Do If TextWidth(Ln[I]) + 40 > W Then W := TextWidth(Ln[I]) + 40;
 If W < 260 Then W := 260;
 If W > SW - 20 Then W := SW - 20;
 H := 70 + N * 12;
 X1 := (SW - W) Div 2; Y1 := (SH - H) Div 2;
 X2 := X1 + W; Y2 := Y1 + H;
 BX1 := (X1 + X2 - BW) Div 2; BY1 := Y2 - BH - 8;
 HideMouse;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 Y := Y1 + 28;
 For I := 0 To N - 1 Do Begin
  Txt(X1 + 20, Y, Ln[I], Black);
  Inc(Y, 12);
 End;
 DrawButton(BX1, BY1, BX1 + BW, BY1 + BH, 'OK', True);
 ShowMouse;
 FlushKeys;
 WaitRelease;
 Done := False;
 While Not Done Do Begin
  If KeyPressed Then Begin
   K := GetKey(Ext);
   If Not Ext And ((K = #13) Or (K = #27) Or (K = ' ')) Then Done := True;
  End;
  PollMouse;
  If ((MouseB And 1) <> 0) And PtIn(MouseX, MouseY, BX1, BY1, BX1 + BW, BY1 + BH) Then Done := True;
  Delay(15);
 End;
 WaitRelease;
 NeedScene := True; NeedChrome := True;
End;

Function InputString(Const Title, Prompt: String; Var Value: String): Boolean;
Const BW = 80; BH = 22;
Var
 X1, Y1, X2, Y2, BOk1, BCa1, BBY : Integer;
 K : Char;
 Ext, Done, OK : Boolean;
 Cur : String;
Begin
 Cur := Value;
 X1 := (SW - 440) Div 2; Y1 := (SH - 130) Div 2;
 X2 := X1 + 440; Y2 := Y1 + 130;
 BBY := Y2 - BH - 8;
 BOk1 := X1 + 220 - BW - 10; BCa1 := X1 + 220 + 10;
 HideMouse;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 Txt(X1 + 15, Y1 + 28, Prompt, Black);
 DrawButton(BOk1, BBY, BOk1 + BW, BBY + BH, 'OK', False);
 DrawButton(BCa1, BBY, BCa1 + BW, BBY + BH, 'Annuler', False);
 ShowMouse;
 FlushKeys;
 WaitRelease;
 Done := False; OK := False;
 While Not Done Do Begin
  HideMouse;
  PanelBar(X1 + 15, Y1 + 46, X2 - 15, Y1 + 66, White);
  Frame3D(X1 + 14, Y1 + 45, X2 - 14, Y1 + 67, False);
  If Length(Cur) > 50 Then Txt(X1 + 18, Y1 + 52, Copy(Cur, Length(Cur) - 49, 50) + '_', Black)
  Else Txt(X1 + 18, Y1 + 52, Cur + '_', Black);
  ShowMouse;
  Repeat
   PollMouse;
   If (MouseB And 1) <> 0 Then Begin
    If PtIn(MouseX, MouseY, BOk1, BBY, BOk1 + BW, BBY + BH) Then Begin Done := True; OK := True; End
    Else If PtIn(MouseX, MouseY, BCa1, BBY, BCa1 + BW, BBY + BH) Then Begin Done := True; OK := False; End;
   End;
   If Done Or KeyPressed Then Break;
   Delay(15);
  Until False;
  If Done Then Break;
  K := GetKey(Ext);
  If Ext Then Begin
   { ignorer }
  End
  Else If K = #13 Then Begin Done := True; OK := True; End
  Else If K = #27 Then Begin Done := True; OK := False; End
  Else If K = #8 Then Begin
   If Length(Cur) > 0 Then Delete(Cur, Length(Cur), 1);
  End
  Else If (Ord(K) >= 32) And (Length(Cur) < 120) Then Cur := Cur + K;
 End;
 WaitRelease;
 If OK Then Value := Cur;
 InputString := OK;
 NeedScene := True; NeedChrome := True;
End;

{ Affiche TextBuf[0..TextN-1] dans une fenetre defilante }
Procedure ShowTextDlg(Const Title: String);
Const BW = 80; BH = 22;
Var
 X1, Y1, X2, Y2, BX1, BY1, LX1, LY1, LX2, LY2, MaxLn, Top, I, Y : Integer;
 K : Char;
 Ext, Done, Redraw : Boolean;
Begin
 X1 := 30; Y1 := 24; X2 := SW - 30; Y2 := SH - 24;
 LX1 := X1 + 10; LY1 := Y1 + 24; LX2 := X2 - 10; LY2 := Y2 - 40;
 MaxLn := (LY2 - LY1 - 4) Div 10;
 BX1 := (X1 + X2 - BW) Div 2; BY1 := Y2 - BH - 8;
 Top := 0;
 FlushKeys;
 WaitRelease;
 Done := False; Redraw := True;
 While Not Done Do Begin
  If Redraw Then Begin
   HideMouse;
   DrawDialogBox(X1, Y1, X2, Y2, Title);
   PanelBar(LX1, LY1, LX2, LY2, White);
   Frame3D(LX1, LY1, LX2, LY2, False);
   For I := 0 To MaxLn - 1 Do Begin
    If Top + I >= TextN Then Break;
    Y := LY1 + 3 + I * 10;
    Txt(LX1 + 4, Y, TextBuf[Top + I], Black);
   End;
   If TextN > MaxLn Then Begin
    Txt(LX2 - 90, LY2 + 4, IntToStr(Top + 1) + '-' +
        IntToStr(MinI(Top + MaxLn, TextN)) + '/' + IntToStr(TextN), Black);
   End;
   DrawButton(BX1, BY1, BX1 + BW, BY1 + BH, 'Fermer', True);
   ShowMouse;
   Redraw := False;
  End;
  If KeyPressed Then Begin
   K := GetKey(Ext);
   If Ext Then Begin
    Case K Of
     #72: If Top > 0 Then Begin Dec(Top); Redraw := True; End;
     #80: If Top + MaxLn < TextN Then Begin Inc(Top); Redraw := True; End;
     #73: Begin Top := MaxI(0, Top - MaxLn); Redraw := True; End;
     #81: Begin
      If Top + MaxLn < TextN Then Top := MinI(TextN - MaxLn, Top + MaxLn);
      If Top < 0 Then Top := 0;
      Redraw := True;
     End;
     #71: Begin Top := 0; Redraw := True; End;
     #79: Begin Top := MaxI(0, TextN - MaxLn); Redraw := True; End;
    End;
   End
   Else If (K = #13) Or (K = #27) Then Done := True;
  End;
  PollMouse;
  If (MouseB And 1) <> 0 Then Begin
   If PtIn(MouseX, MouseY, BX1, BY1, BX1 + BW, BY1 + BH) Then Done := True
   Else If PtIn(MouseX, MouseY, LX1, LY1, LX2, LY1 + (LY2 - LY1) Div 2) And (Top > 0) Then Begin
    Top := MaxI(0, Top - 1); Redraw := True; Delay(60);
   End
   Else If PtIn(MouseX, MouseY, LX1, LY1 + (LY2 - LY1) Div 2, LX2, LY2) And (Top + MaxLn < TextN) Then Begin
    Inc(Top); Redraw := True; Delay(60);
   End;
  End;
  Delay(15);
 End;
 WaitRelease;
 NeedScene := True; NeedChrome := True;
End;

{ --------------------------------------------------------------------- }
{ Boite d'ouverture de fichier                                          }
{ --------------------------------------------------------------------- }

Type
 TFDEntry = Record
  Name  : String[40];
  IsDir : Boolean;
 End;

Var
 FDFiles : Array[0..FD_MAX-1] Of TFDEntry;
 FDCount, FDTop, FDSel : Integer;
 FDMask  : String;

Procedure FD_Scan(Const Mask: String);
Var
 Sr : SearchRec;
 I, J : Integer;
 Tmp : TFDEntry;
Begin
 FDCount := 0;
 FindFirst('*', AnyFile, Sr);
 While (DosError = 0) And (FDCount < FD_MAX) Do Begin
  If ((Sr.Attr And Directory) <> 0) And (Sr.Name <> '.') And
     ((Sr.Name = '..') Or (Sr.Name[1] <> '.')) Then Begin
   FDFiles[FDCount].Name := Sr.Name;
   FDFiles[FDCount].IsDir := True;
   Inc(FDCount);
  End;
  FindNext(Sr);
 End;
 {$IFDEF FPC}FindClose(Sr);{$ENDIF}
 FindFirst(Mask, Archive Or ReadOnly, Sr);
 While (DosError = 0) And (FDCount < FD_MAX) Do Begin
  If (Sr.Attr And Directory) = 0 Then Begin
   FDFiles[FDCount].Name := Sr.Name;
   FDFiles[FDCount].IsDir := False;
   Inc(FDCount);
  End;
  FindNext(Sr);
 End;
 {$IFDEF FPC}FindClose(Sr);{$ENDIF}
 For I := 1 To FDCount - 1 Do Begin
  Tmp := FDFiles[I];
  J := I - 1;
  While (J >= 0) And ((FDFiles[J].IsDir < Tmp.IsDir) Or
        ((FDFiles[J].IsDir = Tmp.IsDir) And (UpStr(FDFiles[J].Name) > UpStr(Tmp.Name)))) Do Begin
   FDFiles[J + 1] := FDFiles[J];
   Dec(J);
  End;
  FDFiles[J + 1] := Tmp;
 End;
End;

Procedure FD_DrawList(X1, Y1, X2, Y2: Integer);
Var
 I, Y, MaxLn : Integer;
 S : String;
Begin
 PanelBar(X1, Y1, X2, Y2, White);
 Frame3D(X1, Y1, X2, Y2, False);
 MaxLn := (Y2 - Y1 - 4) Div 10;
 For I := 0 To MaxLn - 1 Do Begin
  If FDTop + I >= FDCount Then Break;
  Y := Y1 + 3 + I * 10;
  If FDTop + I = FDSel Then Begin
   PanelBar(X1 + 2, Y, X2 - 2, Y + 9, Blue);
   SetColor(White);
  End
  Else SetColor(Black);
  If FDFiles[FDTop + I].IsDir Then S := '[' + FDFiles[FDTop + I].Name + ']'
  Else S := FDFiles[FDTop + I].Name;
  OutTextXY(X1 + 6, Y + 1, S);
 End;
End;

Function FileOpenDialog(Const Title: String; Var FileName: String;
                        Const InitialMask: String): Boolean;
Const BW = 84; BH = 22;
Var
 X1, Y1, X2, Y2, LX1, LY1, LX2, LY2, BOu1, BAn1, BMa1, BBY, MaxLn, I : Integer;
 K : Char;
 Ext, Done, Res, Redraw, Full : Boolean;
 OldDir, CurDir, NewMask, Typed : String;
 T0, T1 : LongInt;
 LastSel : Integer;
 Ac : Integer;
Begin
 FileOpenDialog := False;
 GetDir(0, OldDir);
 FDMask := InitialMask;
 FD_Scan(FDMask);
 FDTop := 0; FDSel := 0;
 X1 := 50; Y1 := 40; X2 := SW - 50; Y2 := SH - 40;
 LX1 := X1 + 10; LY1 := Y1 + 52; LX2 := X2 - 10; LY2 := Y2 - 76;
 MaxLn := (LY2 - LY1 - 4) Div 10;
 BBY := Y2 - BH - 8;
 BOu1 := X2 - 3 * (BW + 8) - 4; BAn1 := X2 - 2 * (BW + 8) - 4; BMa1 := X2 - (BW + 8) - 4;
 Typed := '';
 Done := False; Res := False; Full := True; Redraw := True;
 LastSel := -1; T0 := 0;
 FlushKeys;
 WaitRelease;
 While Not Done Do Begin
  If Full Then Begin
   HideMouse;
   DrawDialogBox(X1, Y1, X2, Y2, Title);
   GetDir(0, CurDir);
   Txt(X1 + 10, Y1 + 24, 'Dossier : ' + Copy(CurDir, 1, 60), Black);
   Txt(X1 + 10, Y1 + 38, 'Masque  : ' + FDMask, Black);
   Txt(X1 + 10, Y2 - 64, 'Fichier :', Black);
   DrawButton(BOu1, BBY, BOu1 + BW, BBY + BH, 'Ouvrir', True);
   DrawButton(BAn1, BBY, BAn1 + BW, BBY + BH, 'Annuler', False);
   DrawButton(BMa1, BBY, BMa1 + BW, BBY + BH, 'Masque...', False);
   Full := False; Redraw := True;
   ShowMouse;
  End;
  If Redraw Then Begin
   HideMouse;
   FD_DrawList(LX1, LY1, LX2, LY2);
   PanelBar(X1 + 70, Y2 - 67, X2 - 10, Y2 - 52, White);
   Frame3D(X1 + 69, Y2 - 68, X2 - 9, Y2 - 51, False);
   Txt(X1 + 74, Y2 - 64, Typed + '_', Black);
   ShowMouse;
   Redraw := False;
  End;
  Ac := 0;
  PollMouse;
  If (MouseB And 1) <> 0 Then Begin
   If PtIn(MouseX, MouseY, BOu1, BBY, BOu1 + BW, BBY + BH) Then Ac := 1
   Else If PtIn(MouseX, MouseY, BAn1, BBY, BAn1 + BW, BBY + BH) Then Ac := 2
   Else If PtIn(MouseX, MouseY, BMa1, BBY, BMa1 + BW, BBY + BH) Then Ac := 3
   Else If PtIn(MouseX, MouseY, LX1, LY1, LX2, LY2) Then Begin
    I := FDTop + (MouseY - LY1 - 3) Div 10;
    If (I >= 0) And (I < FDCount) Then Begin
     T1 := TickMs;
     If (I = LastSel) And (ElapsedMs(T0) < 500) Then Ac := 1;
     FDSel := I; LastSel := I; T0 := T1;
     Typed := '';
     Redraw := True;
    End;
    WaitRelease;
   End;
   If Ac <> 0 Then WaitRelease;
  End;
  If KeyPressed Then Begin
   K := GetKey(Ext);
   If Ext Then Begin
    Case K Of
     #72: Begin If FDSel > 0 Then Dec(FDSel); If FDSel < FDTop Then FDTop := FDSel; Redraw := True; End;
     #80: Begin
      If FDSel < FDCount - 1 Then Inc(FDSel);
      If FDSel >= FDTop + MaxLn Then FDTop := FDSel - MaxLn + 1;
      Redraw := True;
     End;
     #73: Begin Dec(FDSel, MaxLn); If FDSel < 0 Then FDSel := 0; FDTop := FDSel; Redraw := True; End;
     #81: Begin
      Inc(FDSel, MaxLn);
      If FDSel >= FDCount Then FDSel := FDCount - 1;
      FDTop := FDSel - MaxLn + 1;
      If FDTop < 0 Then FDTop := 0;
      Redraw := True;
     End;
     #71: Begin FDSel := 0; FDTop := 0; Redraw := True; End;
     #79: Begin FDSel := FDCount - 1; FDTop := FDSel - MaxLn + 1; If FDTop < 0 Then FDTop := 0; Redraw := True; End;
    End;
   End
   Else If K = #13 Then Ac := 1
   Else If K = #27 Then Ac := 2
   Else If K = #9 Then Ac := 3
   Else If K = #8 Then Begin
    If Length(Typed) > 0 Then Delete(Typed, Length(Typed), 1);
    Redraw := True;
   End
   Else If (Ord(K) >= 32) And (Length(Typed) < 70) Then Begin
    Typed := Typed + K;
    Redraw := True;
   End;
  End;
  Case Ac Of
   1: Begin
    If Typed <> '' Then Begin
     FileName := Typed;
     Res := True; Done := True;
    End
    Else If (FDSel >= 0) And (FDSel < FDCount) Then Begin
     If FDFiles[FDSel].IsDir Then Begin
      {$I-} ChDir(FDFiles[FDSel].Name); {$I+}
      If IOResult <> 0 Then ;
      FD_Scan(FDMask);
      FDTop := 0; FDSel := 0; LastSel := -1;
      Full := True;
     End
     Else Begin
      FileName := FDFiles[FDSel].Name;
      Res := True; Done := True;
     End;
    End;
   End;
   2: Begin Res := False; Done := True; End;
   3: Begin
    NewMask := FDMask;
    If InputString('Masque', 'Entrer le masque (ex: *.PDB) :', NewMask) Then Begin
     FDMask := NewMask;
     FD_Scan(FDMask);
     FDTop := 0; FDSel := 0;
    End;
    Full := True;
   End;
  End;
  Delay(15);
 End;
 If Not Res Then Begin {$I-} ChDir(OldDir); {$I+} If IOResult <> 0 Then ; End;
 FileOpenDialog := Res;
 WaitRelease;
 NeedScene := True; NeedChrome := True;
End;

Function FileSaveDialog(Const Title: String; Var FileName: String): Boolean;
Var Cur: String;
Begin
 Cur := FileName;
 FileSaveDialog := False;
 If InputString(Title, 'Nom du fichier a enregistrer :', Cur) And (Cur <> '') Then Begin
  FileName := Cur;
  FileSaveDialog := True;
 End;
End;

{ --------------------------------------------------------------------- }
{ Menus : definition                                                    }
{ --------------------------------------------------------------------- }

Const
 PM_ROOT    = 0;
 PM_SEL     = 1;
 PM_ELEM    = 2;
 PM_CHAIN   = 3;
 PM_VIEW    = 4;
 PM_STYLE   = 5;
 PM_SCHEME  = 6;
 PM_ATOMS   = 7;
 PM_BONDS   = 8;
 PM_STRUCT  = 9;
 PM_LABEL   = 10;
 PM_COLOR   = 11;
 PM_COLSCH  = 12;
 PM_COLATOM = 13;
 PM_COLRIB  = 14;
 PM_COLBG   = 15;
 PM_SURF    = 16;
 PM_SYM     = 17;
 PM_ZOOM    = 18;
 PM_SPIN    = 19;
 PM_VIB     = 20;
 PM_ANIM    = 21;
 PM_MEAS    = 22;
 PM_PICK    = 23;
 PM_SHOW    = 24;
 PM_FILE    = 25;
 PM_CALC    = 26;
 PM_LANG    = 27;
 PM_ABOUT   = 28;
 PM_MODEL   = 29;
 PM_UNITS   = 30;
 PM_B_FILE  = 31;
 PM_B_EDIT  = 32;
 PM_B_DISP  = 33;
 PM_B_VIEW  = 34;
 PM_B_MEAS  = 35;
 PM_B_TOOLS = 36;
 PM_B_MACRO = 37;
 PM_B_HELP  = 38;
 NUM_BAR    = 8;

 AtomPct  : Array[0..9] Of Integer = (0, 10, 20, 30, 40, 50, 60, 80, 100, 120);
 BondRad  : Array[0..6] Of Integer = (0, 1, 5, 10, 15, 20, 30);
 ZoomPct  : Array[0..6] Of Integer = (25, 50, 100, 150, 200, 400, 800);
 NUM_LBLPRE = 7;
 LblPre   : Array[0..NUM_LBLPRE-1] Of String[8] = (
  '', '%a', '%e', '%i', '%r%n', '%a %i', '%e%i');

Var
 BarX, BarW : Array[0..NUM_BAR-1] Of Integer;
 BarCur     : Integer;
 ChainList  : Array[0..25] Of Char;
 NChainList : Integer;
 Lv         : Array[0..3] Of TLvl;
 NLv        : Integer;
 NSelCount  : Integer;

Procedure MNew(Id: Integer; Const Title: String);
Begin
 PM^[Id].N := 0;
 PM^[Id].Title := Copy(Title, 1, 14);
End;

Procedure MAdd(Id: Integer; Const Cap: String; Cmd, Sub: Integer; Kind: Byte);
Begin
 If PM^[Id].N >= MAX_PMI Then Exit;
 PM^[Id].It[PM^[Id].N].Cap := Cap;
 PM^[Id].It[PM^[Id].N].Cmd := Cmd;
 PM^[Id].It[PM^[Id].N].Sub := Sub;
 PM^[Id].It[PM^[Id].N].Kind := Kind;
 Inc(PM^[Id].N);
End;

Procedure MI(Id: Integer; Const Cap: String; Cmd: Integer);
Begin MAdd(Id, Cap, Cmd, -1, 0); End;

Procedure MC(Id: Integer; Const Cap: String; Cmd: Integer);
Begin MAdd(Id, Cap, Cmd, -1, 1); End;

Procedure MR(Id: Integer; Const Cap: String; Cmd: Integer);
Begin MAdd(Id, Cap, Cmd, -1, 2); End;

Procedure MSub(Id: Integer; Const Cap: String; Sub: Integer);
Begin MAdd(Id, Cap, 0, Sub, 0); End;

Procedure MSep(Id: Integer);
Begin MAdd(Id, '-', 0, -1, 0); End;

Procedure ColorItems(Id: Integer; Base: Integer);
Var I: Integer;
Begin
 For I := 0 To NUM_COLORS - 1 Do MI(Id, ColNameFr[I], Base + I);
End;

Procedure InitMenus;
Var I: Integer;
Begin
 For I := 0 To MAX_PM - 1 Do Begin PM^[I].N := 0; PM^[I].Title := ''; End;

 { --- Menu contextuel principal --- }
 MNew(PM_ROOT, 'Jmol');
 MSub(PM_ROOT, 'modele 1/1', PM_MODEL);
 MSub(PM_ROOT, 'Configurations', -1);
 MSub(PM_ROOT, 'Selection (0)', PM_SEL);
 MSub(PM_ROOT, 'Vue', PM_VIEW);
 MSub(PM_ROOT, 'Style', PM_STYLE);
 MSub(PM_ROOT, 'Couleur', PM_COLOR);
 MSep(PM_ROOT);
 MSub(PM_ROOT, 'Surfaces', PM_SURF);
 MSub(PM_ROOT, 'Symetrie', PM_SYM);
 MSep(PM_ROOT);
 MSub(PM_ROOT, 'Zoom', PM_ZOOM);
 MSub(PM_ROOT, 'Rotation', PM_SPIN);
 MSub(PM_ROOT, 'Vibration', PM_VIB);
 MSub(PM_ROOT, 'Animation', PM_ANIM);
 MSep(PM_ROOT);
 MSub(PM_ROOT, 'Mesures', PM_MEAS);
 MSub(PM_ROOT, 'Choisir', PM_PICK);
 MSep(PM_ROOT);
 MI(PM_ROOT, 'Console', CM_TOG_CONSOLE);
 MSub(PM_ROOT, 'Afficher', PM_SHOW);
 MSub(PM_ROOT, 'Fichier', PM_FILE);
 MSub(PM_ROOT, 'Calcul', PM_CALC);
 MSep(PM_ROOT);
 MSub(PM_ROOT, 'Langue', PM_LANG);
 MSub(PM_ROOT, 'A propos', PM_ABOUT);

 { --- Modele --- }
 MNew(PM_MODEL, 'Modele');
 MI(PM_MODEL, 'Premier modele', CM_FR_FIRST);
 MI(PM_MODEL, 'Modele precedent', CM_FR_PREV);
 MI(PM_MODEL, 'Modele suivant', CM_FR_NEXT);
 MI(PM_MODEL, 'Dernier modele', CM_FR_LAST);
 MSep(PM_MODEL);
 MC(PM_MODEL, 'Tous les modeles', CM_FR_ALL);

 { --- Selection --- }
 MNew(PM_SEL, 'Selection');
 MI(PM_SEL, 'Tout', CM_SEL_BASE + 0);
 MI(PM_SEL, 'Aucun', CM_SEL_BASE + 1);
 MI(PM_SEL, 'Inverser la selection', CM_SELINV);
 MSep(PM_SEL);
 MI(PM_SEL, 'Proteine', CM_SEL_BASE + 2);
 MI(PM_SEL, 'Acides nucleiques', CM_SEL_BASE + 3);
 MI(PM_SEL, 'Eau', CM_SEL_BASE + 4);
 MI(PM_SEL, 'Hetero-atomes', CM_SEL_BASE + 5);
 MI(PM_SEL, 'Hydrogenes', CM_SEL_BASE + 6);
 MI(PM_SEL, 'Squelette', CM_SEL_BASE + 7);
 MI(PM_SEL, 'Chaines laterales', CM_SEL_BASE + 8);
 MSep(PM_SEL);
 MI(PM_SEL, 'Helices', CM_SEL_BASE + 9);
 MI(PM_SEL, 'Feuillets', CM_SEL_BASE + 10);
 MI(PM_SEL, 'Coudes', CM_SEL_BASE + 11);
 MSep(PM_SEL);
 MSub(PM_SEL, 'Elements', PM_ELEM);
 MSub(PM_SEL, 'Chaines', PM_CHAIN);
 MNew(PM_ELEM, 'Elements');
 MNew(PM_CHAIN, 'Chaines');

 { --- Vue --- }
 MNew(PM_VIEW, 'Vue');
 MI(PM_VIEW, 'Avant', CM_VIEW_FRONT);
 MI(PM_VIEW, 'Arriere', CM_VIEW_BACK);
 MI(PM_VIEW, 'Gauche', CM_VIEW_LEFT);
 MI(PM_VIEW, 'Droite', CM_VIEW_RIGHT);
 MI(PM_VIEW, 'Dessus', CM_VIEW_TOP);
 MI(PM_VIEW, 'Dessous', CM_VIEW_BOTTOM);
 MSep(PM_VIEW);
 MI(PM_VIEW, 'Reinitialiser la vue|Home', CM_VIEW_RESET);
 MI(PM_VIEW, 'Centrer sur la selection', CM_CENTERSEL);

 { --- Style --- }
 MNew(PM_STYLE, 'Style');
 MSub(PM_STYLE, 'Schema', PM_SCHEME);
 MSub(PM_STYLE, 'Atomes', PM_ATOMS);
 MSub(PM_STYLE, 'Liaisons', PM_BONDS);
 MSub(PM_STYLE, 'Structure', PM_STRUCT);
 MSub(PM_STYLE, 'Etiquettes', PM_LABEL);
 MSep(PM_STYLE);
 MC(PM_STYLE, 'Hydrogenes|H', CM_TOG_HYDRO);
 MC(PM_STYLE, 'Molecules d''eau', CM_TOG_WATER);
 MC(PM_STYLE, 'Liaisons hydrogene', CM_TOG_HBOND);
 MC(PM_STYLE, 'Halos de selection', CM_TOG_HALO);
 MC(PM_STYLE, 'Perspective', CM_TOG_PERSP);
 MSep(PM_STYLE);
 MC(PM_STYLE, 'Axes|A', CM_TOG_AXES);
 MC(PM_STYLE, 'Boite englobante|B', CM_TOG_BBOX);
 MC(PM_STYLE, 'Cellule unitaire', CM_TOG_CELL);
 MC(PM_STYLE, 'Vecteurs de vibration', CM_TOG_VECTOR);

 MNew(PM_SCHEME, 'Schema');
 MI(PM_SCHEME, 'Remplissage CPK', CM_SCH_CPK);
 MI(PM_SCHEME, 'Boules et batonnets', CM_SCH_BALL);
 MI(PM_SCHEME, 'Batonnets', CM_SCH_STICK);
 MI(PM_SCHEME, 'Fil de fer', CM_SCH_WIRE);

 MNew(PM_ATOMS, 'Atomes');
 For I := 0 To 9 Do
  If I = 0 Then MI(PM_ATOMS, 'Aucun', CM_ATOMPCT_BASE)
  Else MI(PM_ATOMS, IntToStr(AtomPct[I]) + '% van der Waals', CM_ATOMPCT_BASE + I);

 MNew(PM_BONDS, 'Liaisons');
 MI(PM_BONDS, 'Aucune', CM_BONDRAD_BASE + 0);
 MI(PM_BONDS, 'Fil de fer', CM_BONDRAD_BASE + 1);
 For I := 2 To 6 Do
  MI(PM_BONDS, RealToStr(BondRad[I] / 100.0, 2) + ' ' + #143, CM_BONDRAD_BASE + I);

 MNew(PM_STRUCT, 'Structure');
 MI(PM_STRUCT, 'Aucune', CM_STRUCT_BASE + RIB_NONE);
 MI(PM_STRUCT, 'Squelette', CM_STRUCT_BASE + RIB_BACKBONE);
 MI(PM_STRUCT, 'Trace', CM_STRUCT_BASE + RIB_TRACE);
 MI(PM_STRUCT, 'Rubans', CM_STRUCT_BASE + RIB_RIBBON);
 MI(PM_STRUCT, 'Cartoon', CM_STRUCT_BASE + RIB_CARTOON);
 MI(PM_STRUCT, 'Brins', CM_STRUCT_BASE + RIB_STRAND);

 MNew(PM_LABEL, 'Etiquettes');
 MI(PM_LABEL, 'Aucune', CM_LABEL_BASE + 0);
 MI(PM_LABEL, 'Nom de l''atome', CM_LABEL_BASE + 1);
 MI(PM_LABEL, 'Symbole', CM_LABEL_BASE + 2);
 MI(PM_LABEL, 'Numero', CM_LABEL_BASE + 3);
 MI(PM_LABEL, 'Residu', CM_LABEL_BASE + 4);
 MI(PM_LABEL, 'Nom et numero', CM_LABEL_BASE + 5);
 MI(PM_LABEL, 'Symbole et numero', CM_LABEL_BASE + 6);

 { --- Couleur --- }
 MNew(PM_COLOR, 'Couleur');
 MSub(PM_COLOR, 'Schema', PM_COLSCH);
 MSub(PM_COLOR, 'Atomes', PM_COLATOM);
 MSub(PM_COLOR, 'Structure', PM_COLRIB);
 MSub(PM_COLOR, 'Fond', PM_COLBG);

 MNew(PM_COLSCH, 'Schema');
 MI(PM_COLSCH, 'Selon l''element (CPK)', CM_COLSCH_BASE + 0);
 MI(PM_COLSCH, 'Selon la chaine', CM_COLSCH_BASE + 1);
 MI(PM_COLSCH, 'Selon le residu', CM_COLSCH_BASE + 2);
 MI(PM_COLSCH, 'Selon la structure', CM_COLSCH_BASE + 3);
 MI(PM_COLSCH, 'Arc-en-ciel (groupe)', CM_COLSCH_BASE + 4);
 MI(PM_COLSCH, 'Selon la temperature', CM_COLSCH_BASE + 5);

 MNew(PM_COLATOM, 'Atomes');
 ColorItems(PM_COLATOM, CM_COLATOM_BASE);

 MNew(PM_COLRIB, 'Structure');
 MI(PM_COLRIB, 'Selon la structure', CM_COLRIBSCH_BASE + 0);
 MI(PM_COLRIB, 'Arc-en-ciel', CM_COLRIBSCH_BASE + 1);
 MI(PM_COLRIB, 'Selon la chaine', CM_COLRIBSCH_BASE + 2);
 MSep(PM_COLRIB);
 ColorItems(PM_COLRIB, CM_COLRIB_BASE);

 MNew(PM_COLBG, 'Fond');
 ColorItems(PM_COLBG, CM_COLBG_BASE);

 { --- Surfaces, symetrie --- }
 MNew(PM_SURF, 'Surfaces');
 MI(PM_SURF, 'Points - atomes choisis', CM_SURF_SEL);
 MI(PM_SURF, 'Points - tous les atomes', CM_SURF_ALL);
 MI(PM_SURF, 'Aucun point', CM_SURF_OFF);

 MNew(PM_SYM, 'Symetrie');
 MC(PM_SYM, 'Boite englobante', CM_TOG_BBOX);
 MC(PM_SYM, 'Cellule unitaire', CM_TOG_CELL);
 MC(PM_SYM, 'Axes', CM_TOG_AXES);

 { --- Zoom, rotation, vibration, animation --- }
 MNew(PM_ZOOM, 'Zoom');
 MI(PM_ZOOM, 'Zoom avant|+', CM_ZOOM_IN);
 MI(PM_ZOOM, 'Zoom arriere|-', CM_ZOOM_OUT);
 MSep(PM_ZOOM);
 For I := 0 To 6 Do MI(PM_ZOOM, IntToStr(ZoomPct[I]) + '%', CM_ZOOM_BASE + I);
 MSep(PM_ZOOM);
 MI(PM_ZOOM, 'Ajuster a la fenetre', CM_ZOOM_FIT);

 MNew(PM_SPIN, 'Rotation');
 MC(PM_SPIN, 'Rotation automatique|Espace', CM_TOG_SPIN);
 MSep(PM_SPIN);
 MR(PM_SPIN, 'Axe X', CM_SPIN_X);
 MR(PM_SPIN, 'Axe Y', CM_SPIN_Y);
 MR(PM_SPIN, 'Axe Z', CM_SPIN_Z);
 MSep(PM_SPIN);
 MR(PM_SPIN, 'Vitesse lente', CM_SPIN_SLOW);
 MR(PM_SPIN, 'Vitesse moyenne', CM_SPIN_MED);
 MR(PM_SPIN, 'Vitesse rapide', CM_SPIN_FAST);

 MNew(PM_VIB, 'Vibration');
 MC(PM_VIB, 'Activer la vibration', CM_TOG_VIB);
 MI(PM_VIB, 'Panneau de vibration...', CM_VIB_PANEL);
 MC(PM_VIB, 'Vecteurs', CM_TOG_VECTOR);
 MSep(PM_VIB);
 MI(PM_VIB, 'Amplitude +', CM_VIB_BIGGER);
 MI(PM_VIB, 'Amplitude -', CM_VIB_SMALLER);
 MI(PM_VIB, 'Mode suivant', CM_VIB_NEXT);
 MI(PM_VIB, 'Mode precedent', CM_VIB_PREV);

 MNew(PM_ANIM, 'Animation');
 MC(PM_ANIM, 'Lecture', CM_TOG_ANIM);
 MSep(PM_ANIM);
 MI(PM_ANIM, 'Premiere image', CM_FR_FIRST);
 MI(PM_ANIM, 'Image precedente', CM_FR_PREV);
 MI(PM_ANIM, 'Image suivante', CM_FR_NEXT);
 MI(PM_ANIM, 'Derniere image', CM_FR_LAST);
 MSep(PM_ANIM);
 MR(PM_ANIM, 'Une fois', CM_AN_ONCE);
 MR(PM_ANIM, 'Boucle', CM_AN_LOOP);
 MR(PM_ANIM, 'Aller-retour', CM_AN_PALIN);
 MSep(PM_ANIM);
 MI(PM_ANIM, 'Ralentir', CM_AN_SLOWER);
 MI(PM_ANIM, 'Accelerer', CM_AN_FASTER);

 { --- Mesures --- }
 MNew(PM_MEAS, 'Mesures');
 MI(PM_MEAS, 'Distance', CM_MEAS_DIST);
 MI(PM_MEAS, 'Angle', CM_MEAS_ANGLE);
 MI(PM_MEAS, 'Diedre', CM_MEAS_TORS);
 MSep(PM_MEAS);
 MI(PM_MEAS, 'Supprimer les mesures', CM_CLEARMEAS);
 MI(PM_MEAS, 'Liste des mesures...', CM_MEAS_LIST);
 MSep(PM_MEAS);
 MSub(PM_MEAS, 'Unites', PM_UNITS);

 MNew(PM_UNITS, 'Unites');
 MR(PM_UNITS, 'Angstrom', CM_UNIT_ANG);
 MR(PM_UNITS, 'Nanometre', CM_UNIT_NM);
 MR(PM_UNITS, 'Picometre', CM_UNIT_PM);

 MNew(PM_PICK, 'Choisir');
 MR(PM_PICK, 'Identifier l''atome', CM_PICK_BASE + PICK_IDENT);
 MR(PM_PICK, 'Mesurer une distance', CM_PICK_BASE + PICK_DIST);
 MR(PM_PICK, 'Mesurer un angle', CM_PICK_BASE + PICK_ANGLE);
 MR(PM_PICK, 'Mesurer un diedre', CM_PICK_BASE + PICK_TORSION);
 MR(PM_PICK, 'Selectionner un atome', CM_PICK_BASE + PICK_SELATOM);
 MR(PM_PICK, 'Selectionner un residu', CM_PICK_BASE + PICK_SELRES);
 MR(PM_PICK, 'Selectionner une chaine', CM_PICK_BASE + PICK_SELCHAIN);
 MR(PM_PICK, 'Centrer sur l''atome', CM_PICK_BASE + PICK_CENTER);
 MR(PM_PICK, 'Etiqueter l''atome', CM_PICK_BASE + PICK_LABEL);

 { --- Afficher, fichier, calcul, langue, a propos --- }
 MNew(PM_SHOW, 'Afficher');
 MI(PM_SHOW, 'Informations du fichier', CM_SHOW_FILE);
 MI(PM_SHOW, 'Informations du modele', CM_SHOW_MODEL);
 MI(PM_SHOW, 'Liste des atomes', CM_SHOW_ATOMS);
 MI(PM_SHOW, 'Sequence', CM_SHOW_SEQ);
 MI(PM_SHOW, 'Cellule unitaire', CM_SHOW_CELL);
 MI(PM_SHOW, 'Orientation', CM_SHOW_ORIENT);
 MI(PM_SHOW, 'Selection', CM_SHOW_SEL);

 MNew(PM_FILE, 'Fichier');
 MI(PM_FILE, 'Ouvrir...', CM_OPEN);
 MI(PM_FILE, 'Recharger', CM_RELOAD);
 MSep(PM_FILE);
 MI(PM_FILE, 'Enregistrer en PDB...', CM_SAVEPDB);
 MI(PM_FILE, 'Enregistrer en XYZ...', CM_SAVEXYZ);
 MI(PM_FILE, 'Enregistrer en MOL...', CM_SAVEMOL);
 MI(PM_FILE, 'Exporter l''image (BMP)...', CM_SAVEBMP);

 MNew(PM_CALC, 'Calcul');
 MI(PM_CALC, 'Formule moleculaire', CM_CALC_FORM);
 MI(PM_CALC, 'Masse molaire', CM_CALC_MASS);
 MI(PM_CALC, 'Centre de masse', CM_CALC_COM);
 MI(PM_CALC, 'Rayon de giration', CM_CALC_GYR);
 MI(PM_CALC, 'Dimensions', CM_CALC_DIM);
 MSep(PM_CALC);
 MI(PM_CALC, 'Liaisons hydrogene', CM_CALC_HBOND);
 MI(PM_CALC, 'Minimiser l''energie', CM_CALC_MIN);

 MNew(PM_LANG, 'Langue');
 MR(PM_LANG, 'Francais', CM_LANG_FR);

 MNew(PM_ABOUT, 'A propos');
 MI(PM_ABOUT, 'A propos de Jmol...', CM_ABOUT);
 MI(PM_ABOUT, 'Aide...', CM_HELP);

 { --- Barre de menus --- }
 MNew(PM_B_FILE, 'Fichier');
 MI(PM_B_FILE, 'Ouvrir...|Ctrl+O', CM_OPEN);
 MI(PM_B_FILE, 'Recharger|F3', CM_RELOAD);
 MSep(PM_B_FILE);
 MI(PM_B_FILE, 'Enregistrer en PDB...|Ctrl+S', CM_SAVEPDB);
 MI(PM_B_FILE, 'Enregistrer en XYZ...', CM_SAVEXYZ);
 MI(PM_B_FILE, 'Enregistrer en MOL...', CM_SAVEMOL);
 MI(PM_B_FILE, 'Exporter l''image (BMP)...', CM_SAVEBMP);
 MSep(PM_B_FILE);
 MI(PM_B_FILE, 'Executer un script...', CM_RUNSCRIPT);
 MI(PM_B_FILE, 'Informations du fichier', CM_FILEINFO);
 MSep(PM_B_FILE);
 MI(PM_B_FILE, 'Quitter|Ctrl+Q', CM_QUIT);

 MNew(PM_B_EDIT, 'Edition');
 MI(PM_B_EDIT, 'Tout selectionner|Ctrl+A', CM_SELALL);
 MI(PM_B_EDIT, 'Ne rien selectionner', CM_SELNONE);
 MI(PM_B_EDIT, 'Inverser la selection', CM_SELINV);
 MI(PM_B_EDIT, 'Selection par expression...', CM_SELEXPR);
 MSep(PM_B_EDIT);
 MI(PM_B_EDIT, 'Masquer la selection', CM_HIDESEL);
 MI(PM_B_EDIT, 'Afficher tout', CM_SHOWALL);
 MSep(PM_B_EDIT);
 MI(PM_B_EDIT, 'Centrer sur la selection', CM_CENTERSEL);
 MSep(PM_B_EDIT);
 MI(PM_B_EDIT, 'Supprimer les mesures', CM_CLEARMEAS);
 MI(PM_B_EDIT, 'Supprimer les etiquettes', CM_CLEARLABELS);

 MNew(PM_B_DISP, 'Affichage');
 MSub(PM_B_DISP, 'Schema', PM_SCHEME);
 MSub(PM_B_DISP, 'Atomes', PM_ATOMS);
 MSub(PM_B_DISP, 'Liaisons', PM_BONDS);
 MSub(PM_B_DISP, 'Structure', PM_STRUCT);
 MSub(PM_B_DISP, 'Etiquettes', PM_LABEL);
 MSep(PM_B_DISP);
 MC(PM_B_DISP, 'Hydrogenes|H', CM_TOG_HYDRO);
 MC(PM_B_DISP, 'Molecules d''eau', CM_TOG_WATER);
 MC(PM_B_DISP, 'Liaisons hydrogene', CM_TOG_HBOND);
 MC(PM_B_DISP, 'Halos de selection', CM_TOG_HALO);
 MSep(PM_B_DISP);
 MC(PM_B_DISP, 'Perspective', CM_TOG_PERSP);
 MC(PM_B_DISP, 'Axes|A', CM_TOG_AXES);
 MC(PM_B_DISP, 'Boite englobante|B', CM_TOG_BBOX);
 MC(PM_B_DISP, 'Cellule unitaire', CM_TOG_CELL);
 MC(PM_B_DISP, 'Vecteurs de vibration', CM_TOG_VECTOR);

 MNew(PM_B_VIEW, 'Vue');
 MSub(PM_B_VIEW, 'Orientation', PM_VIEW);
 MSub(PM_B_VIEW, 'Zoom', PM_ZOOM);
 MSub(PM_B_VIEW, 'Rotation', PM_SPIN);
 MSep(PM_B_VIEW);
 MI(PM_B_VIEW, 'Reinitialiser la vue|Home', CM_VIEW_RESET);
 MI(PM_B_VIEW, 'Centrer sur la selection', CM_CENTERSEL);
 MSep(PM_B_VIEW);
 MSub(PM_B_VIEW, 'Modele', PM_MODEL);
 MSub(PM_B_VIEW, 'Animation', PM_ANIM);
 MSub(PM_B_VIEW, 'Vibration', PM_VIB);

 MNew(PM_B_MEAS, 'Mesure');
 MI(PM_B_MEAS, 'Distance...', CM_MEAS_DIST);
 MI(PM_B_MEAS, 'Angle...', CM_MEAS_ANGLE);
 MI(PM_B_MEAS, 'Diedre...', CM_MEAS_TORS);
 MSep(PM_B_MEAS);
 MI(PM_B_MEAS, 'Supprimer les distances', CM_DEL_DIST);
 MI(PM_B_MEAS, 'Supprimer les angles', CM_DEL_ANGLE);
 MI(PM_B_MEAS, 'Supprimer les diedres', CM_DEL_TORS);
 MSep(PM_B_MEAS);
 MI(PM_B_MEAS, 'Supprimer toutes les mesures', CM_CLEARMEAS);
 MSep(PM_B_MEAS);
 MI(PM_B_MEAS, 'Voir la liste des mesures...', CM_MEAS_LIST);
 MSep(PM_B_MEAS);
 MSub(PM_B_MEAS, 'Unites', PM_UNITS);

 MNew(PM_B_TOOLS, 'Outils');
 MI(PM_B_TOOLS, 'Console|F4', CM_TOG_CONSOLE);
 MI(PM_B_TOOLS, 'Panneau de vibration...', CM_VIB_PANEL);
 MSep(PM_B_TOOLS);
 MSub(PM_B_TOOLS, 'Couleur', PM_COLOR);
 MSub(PM_B_TOOLS, 'Surfaces', PM_SURF);
 MSub(PM_B_TOOLS, 'Choisir', PM_PICK);
 MSub(PM_B_TOOLS, 'Afficher', PM_SHOW);
 MSub(PM_B_TOOLS, 'Calcul', PM_CALC);

 MNew(PM_B_MACRO, 'Macros');
 MI(PM_B_MACRO, 'Boules et batonnets', CM_MACRO_BASE + 0);
 MI(PM_B_MACRO, 'Remplissage CPK', CM_MACRO_BASE + 1);
 MI(PM_B_MACRO, 'Fil de fer', CM_MACRO_BASE + 2);
 MI(PM_B_MACRO, 'Cartoon par structure', CM_MACRO_BASE + 3);
 MI(PM_B_MACRO, 'Cartoon arc-en-ciel', CM_MACRO_BASE + 4);
 MI(PM_B_MACRO, 'Cartoon et chaines laterales', CM_MACRO_BASE + 5);
 MI(PM_B_MACRO, 'Surface de points (vdW)', CM_MACRO_BASE + 6);
 MI(PM_B_MACRO, 'Cartoon par chaine', CM_MACRO_BASE + 7);
 MI(PM_B_MACRO, 'Liaisons hydrogene', CM_MACRO_BASE + 8);
 MI(PM_B_MACRO, 'Visite en rotation', CM_MACRO_BASE + 9);
 MSep(PM_B_MACRO);
 MI(PM_B_MACRO, 'Tout reinitialiser', CM_MACRO_BASE + 10);

 MNew(PM_B_HELP, 'Aide');
 MI(PM_B_HELP, 'Aide Jmol...|F1', CM_HELP);
 MI(PM_B_HELP, 'Souris et clavier...', CM_HELP_MOUSE);
 MI(PM_B_HELP, 'Commandes de la console...', CM_HELP_CMDS);
 MSep(PM_B_HELP);
 MI(PM_B_HELP, 'A propos de Jmol...', CM_ABOUT);
End;

{ --------------------------------------------------------------------- }
{ Etat des elements de menu                                             }
{ --------------------------------------------------------------------- }

Function CountSelected: Integer;
Var I, N: Integer;
Begin
 N := 0;
 For I := 0 To NAtoms - 1 Do If Atoms^[I].Sel Then Inc(N);
 CountSelected := N;
End;

Function ItemChecked(Cmd: Integer): Boolean;
Begin
 ItemChecked := False;
 Case Cmd Of
  CM_TOG_HYDRO : ItemChecked := ShowHydro;
  CM_TOG_WATER : ItemChecked := ShowWater;
  CM_TOG_HBOND : ItemChecked := ShowHB;
  CM_TOG_HALO  : ItemChecked := ShowHalo;
  CM_TOG_PERSP : ItemChecked := Persp;
  CM_TOG_AXES  : ItemChecked := ShowAxes;
  CM_TOG_BBOX  : ItemChecked := ShowBBox;
  CM_TOG_CELL  : ItemChecked := ShowCell;
  CM_TOG_VECTOR: ItemChecked := ShowVec;
  CM_TOG_SPIN  : ItemChecked := SpinOn;
  CM_TOG_VIB   : ItemChecked := VibOn;
  CM_TOG_ANIM  : ItemChecked := AnimOn;
  CM_SPIN_X    : ItemChecked := SpinAxis = 0;
  CM_SPIN_Y    : ItemChecked := SpinAxis = 1;
  CM_SPIN_Z    : ItemChecked := SpinAxis = 2;
  CM_SPIN_SLOW : ItemChecked := SpinSpeed < 2.0;
  CM_SPIN_MED  : ItemChecked := (SpinSpeed >= 2.0) And (SpinSpeed < 5.0);
  CM_SPIN_FAST : ItemChecked := SpinSpeed >= 5.0;
  CM_AN_ONCE   : ItemChecked := AnimMode = 0;
  CM_AN_LOOP   : ItemChecked := AnimMode = 1;
  CM_AN_PALIN  : ItemChecked := AnimMode = 2;
  CM_FR_ALL    : ItemChecked := CurModel = 0;
  CM_UNIT_ANG  : ItemChecked := MeasUnit = 0;
  CM_UNIT_NM   : ItemChecked := MeasUnit = 1;
  CM_UNIT_PM   : ItemChecked := MeasUnit = 2;
  CM_LANG_FR   : ItemChecked := True;
  CM_VIB_PANEL : ItemChecked := VibPanel;
 Else
  If (Cmd >= CM_PICK_BASE) And (Cmd <= CM_PICK_BASE + 9) Then
   ItemChecked := PickMode = Cmd - CM_PICK_BASE;
 End;
End;

Function ItemEnabled(M, I: Integer): Boolean;
Var Cmd, Sub: Integer;
Begin
 Cmd := PM^[M].It[I].Cmd;
 Sub := PM^[M].It[I].Sub;
 ItemEnabled := True;
 If (Cmd = 0) And (Sub < 0) Then Begin ItemEnabled := False; Exit; End;
 If Sub = PM_VIB Then ItemEnabled := HasVibData
 Else If Sub = PM_ANIM Then ItemEnabled := NModels > 1
 Else If Sub = PM_MODEL Then ItemEnabled := NModels > 1;
 Case Cmd Of
  CM_TOG_CELL, CM_SHOW_CELL: ItemEnabled := HasCell;
  CM_TOG_VIB, CM_VIB_PANEL, CM_TOG_VECTOR, CM_VIB_BIGGER, CM_VIB_SMALLER,
  CM_VIB_NEXT, CM_VIB_PREV: ItemEnabled := HasVibData;
  CM_TOG_ANIM, CM_FR_FIRST, CM_FR_PREV, CM_FR_NEXT, CM_FR_LAST, CM_FR_ALL,
  CM_AN_ONCE, CM_AN_LOOP, CM_AN_PALIN, CM_AN_SLOWER, CM_AN_FASTER:
   ItemEnabled := NModels > 1;
  CM_RELOAD, CM_SAVEPDB, CM_SAVEXYZ, CM_SAVEMOL, CM_SAVEBMP, CM_FILEINFO:
   ItemEnabled := Loaded;
  CM_TOG_HBOND: ItemEnabled := Loaded;
 End;
 If (Sub >= 0) And Not Loaded And (Sub <> PM_FILE) And (Sub <> PM_LANG) And
    (Sub <> PM_ABOUT) And (Sub <> PM_SHOW) Then ItemEnabled := False;
End;

{ Met a jour les titres et les listes dynamiques avant l'affichage }
Procedure RefreshMenus;
Var
 I, Z, K, C : Integer;
 Seen : Array[0..103] Of Integer;
 S : String;
 Found : Boolean;
Begin
 NSelCount := CountSelected;
 PM^[PM_ROOT].It[0].Cap := 'modele ' + IntToStr(MaxI(CurModel, 1)) + '/' + IntToStr(NModels);
 If CurModel = 0 Then PM^[PM_ROOT].It[0].Cap := 'modele tous/' + IntToStr(NModels);
 PM^[PM_ROOT].It[2].Cap := 'Selection (' + IntToStr(NSelCount) + ')';
 PM^[PM_ROOT].It[1].Sub := -1;
 { elements presents }
 PM^[PM_ELEM].N := 0;
 FillChar(Seen, SizeOf(Seen), 0);
 For I := 0 To NAtoms - 1 Do Begin
  Z := Atoms^[I].ElemZ; If Z > 103 Then Z := 0;
  If (Atoms^[I].Model = MaxI(CurModel, 1)) Or (CurModel = 0) Then Inc(Seen[Z]);
 End;
 For Z := 0 To 103 Do
  If (Seen[Z] > 0) And (PM^[PM_ELEM].N < MAX_PMI) Then
   MI(PM_ELEM, ElemDB[Z].Name + ' (' + ElemDB[Z].Sym + ') : ' + IntToStr(Seen[Z]), CM_ELEM_BASE + Z);
 { chaines presentes }
 PM^[PM_CHAIN].N := 0;
 NChainList := 0;
 For I := 0 To NRes - 1 Do Begin
  Found := False;
  For K := 0 To NChainList - 1 Do If ChainList[K] = Resi^[I].Chain Then Found := True;
  If (Not Found) And (NChainList < 26) And (Resi^[I].IsAA Or Resi^[I].IsNA) Then Begin
   ChainList[NChainList] := Resi^[I].Chain;
   Inc(NChainList);
  End;
 End;
 For K := 0 To NChainList - 1 Do Begin
  C := 0;
  For I := 0 To NRes - 1 Do If Resi^[I].Chain = ChainList[K] Then Inc(C);
  S := ChainList[K];
  If S = ' ' Then S := '(sans nom)';
  MI(PM_CHAIN, 'Chaine ' + S + ' : ' + IntToStr(C) + ' residus', CM_CHAIN_BASE + K);
 End;
End;

{ --------------------------------------------------------------------- }
{ Menus : affichage et moteur modal                                     }
{ --------------------------------------------------------------------- }

Function ItemH(M, I: Integer): Integer;
Begin
 If PM^[M].It[I].Cap = '-' Then ItemH := 8 Else ItemH := 17;
End;

Function ItemTop(M, I: Integer): Integer;
Var K, Y: Integer;
Begin
 Y := 0;
 For K := 0 To I - 1 Do Y := Y + ItemH(M, K);
 ItemTop := Y;
End;

Function CapLeft(Const Cap: String): String;
Var P: Integer;
Begin
 P := Pos('|', Cap);
 If P > 0 Then CapLeft := Copy(Cap, 1, P - 1) Else CapLeft := Cap;
End;

Function CapRight(Const Cap: String): String;
Var P: Integer;
Begin
 P := Pos('|', Cap);
 If P > 0 Then CapRight := Copy(Cap, P + 1, 30) Else CapRight := '';
End;

Procedure SaveBack(Var L: TLvl);
Var Sz: LongInt;
Begin
 L.Buf := Nil; L.BufSz := 0;
 Sz := ImageSize(L.X, L.Y, L.X + L.W + 2, L.Y + L.H + 2);
 {$IFNDEF FPC}
 If (Sz <= 0) Or (Sz > 65000) Then Exit;
 GetMem(L.Buf, Word(Sz));
 {$ELSE}
 If Sz <= 0 Then Exit;
 GetMem(L.Buf, Sz);
 {$ENDIF}
 L.BufSz := Sz;
 GetImage(L.X, L.Y, L.X + L.W + 2, L.Y + L.H + 2, L.Buf^);
End;

Procedure RestoreBack(Var L: TLvl);
Begin
 If L.Buf <> Nil Then Begin
  PutImage(L.X, L.Y, L.Buf^, NormalPut);
  {$IFNDEF FPC}
  FreeMem(L.Buf, Word(L.BufSz));
  {$ELSE}
  FreeMem(L.Buf, L.BufSz);
  {$ENDIF}
  L.Buf := Nil; L.BufSz := 0;
 End;
End;

Procedure DrawLevelItem(Li, I: Integer);
Var
 M, Y, X, W, H, Cx, Cy : Integer;
 Sel, En : Boolean;
 Cap, Sc : String;
 P : Array[0..2] Of PointType;
 Fg : Byte;
Begin
 M := Lv[Li].Menu;
 X := Lv[Li].X; W := Lv[Li].W;
 Y := Lv[Li].Y + 2 + ItemTop(M, I);
 H := ItemH(M, I);
 If PM^[M].It[I].Cap = '-' Then Begin
  PanelBar(X + 2, Y, X + W - 2, Y + H - 1, LightGray);
  SetColor(DarkGray); Line(X + 4, Y + 3, X + W - 4, Y + 3);
  SetColor(White); Line(X + 4, Y + 4, X + W - 4, Y + 4);
  Exit;
 End;
 Sel := Lv[Li].Sel = I;
 En := ItemEnabled(M, I);
 Cap := CapLeft(PM^[M].It[I].Cap);
 Sc := CapRight(PM^[M].It[I].Cap);
 If Sel And En Then Begin
  PanelBar(X + 2, Y, X + W - 2, Y + H - 1, Blue);
  Fg := White;
 End
 Else Begin
  PanelBar(X + 2, Y, X + W - 2, Y + H - 1, LightGray);
  If En Then Fg := Black Else Fg := DarkGray;
 End;
 If Not En And Not Sel Then Txt(X + 23, Y + 6, Cap, White);
 Txt(X + 22, Y + 5, Cap, Fg);
 If Sc <> '' Then Txt(X + W - 8 - TextWidth(Sc), Y + 5, Sc, Fg);
 If PM^[M].It[I].Kind <> 0 Then
  If ItemChecked(PM^[M].It[I].Cmd) Then Begin
   If PM^[M].It[I].Kind = 1 Then Txt(X + 7, Y + 5, #251, Fg)
   Else Txt(X + 7, Y + 5, #7, Fg);
  End;
 If PM^[M].It[I].Sub >= 0 Then Begin
  Cx := X + W - 9; Cy := Y + H Div 2;
  P[0].X := Cx - 3; P[0].Y := Cy - 4;
  P[1].X := Cx - 3; P[1].Y := Cy + 4;
  P[2].X := Cx + 2; P[2].Y := Cy;
  SetColor(Fg); SetFillStyle(SolidFill, Fg);
  FillPoly(3, P);
 End;
End;

Procedure DrawLevel(Li: Integer);
Var I, M: Integer;
Begin
 M := Lv[Li].Menu;
 HideMouse;
 SetColor(Black);
 Line(Lv[Li].X + Lv[Li].W + 1, Lv[Li].Y + 1, Lv[Li].X + Lv[Li].W + 1, Lv[Li].Y + Lv[Li].H + 1);
 Line(Lv[Li].X + 1, Lv[Li].Y + Lv[Li].H + 1, Lv[Li].X + Lv[Li].W + 1, Lv[Li].Y + Lv[Li].H + 1);
 PanelBar(Lv[Li].X, Lv[Li].Y, Lv[Li].X + Lv[Li].W, Lv[Li].Y + Lv[Li].H, LightGray);
 Frame3D(Lv[Li].X, Lv[Li].Y, Lv[Li].X + Lv[Li].W, Lv[Li].Y + Lv[Li].H, True);
 For I := 0 To PM^[M].N - 1 Do DrawLevelItem(Li, I);
 ShowMouse;
End;

Function MenuWidth(M: Integer): Integer;
Var I, Tw, Wmax: Integer;
Begin
 Wmax := 100;
 For I := 0 To PM^[M].N - 1 Do Begin
  Tw := TextWidth(CapLeft(PM^[M].It[I].Cap)) + 42;
  If CapRight(PM^[M].It[I].Cap) <> '' Then Tw := Tw + TextWidth(CapRight(PM^[M].It[I].Cap)) + 16;
  If Tw > Wmax Then Wmax := Tw;
 End;
 MenuWidth := Wmax;
End;

Procedure OpenLevel(M, X, Y: Integer);
Var
 I, W, H : Integer;
Begin
 H := 4;
 For I := 0 To PM^[M].N - 1 Do H := H + ItemH(M, I);
 W := MenuWidth(M);
 If X + W + 3 > SW Then X := SW - W - 3;
 If X < 0 Then X := 0;
 If Y + H + 3 > SH Then Y := SH - H - 3;
 If Y < 0 Then Y := 0;
 If NLv >= 4 Then Exit;
 Lv[NLv].Menu := M; Lv[NLv].X := X; Lv[NLv].Y := Y;
 Lv[NLv].W := W; Lv[NLv].H := H; Lv[NLv].Sel := -1;
 HideMouse;
 SaveBack(Lv[NLv]);
 ShowMouse;
 Inc(NLv);
 DrawLevel(NLv - 1);
End;

Procedure CloseLevelsFrom(L: Integer);
Begin
 HideMouse;
 While NLv > L Do Begin
  Dec(NLv);
  RestoreBack(Lv[NLv]);
 End;
 ShowMouse;
End;

Function HitItem(Li, Mx, My: Integer): Integer;
Var I, Y, M: Integer;
Begin
 HitItem := -1;
 M := Lv[Li].Menu;
 If (Mx < Lv[Li].X + 2) Or (Mx > Lv[Li].X + Lv[Li].W - 2) Then Exit;
 If (My < Lv[Li].Y + 2) Or (My > Lv[Li].Y + Lv[Li].H - 2) Then Exit;
 Y := Lv[Li].Y + 2;
 For I := 0 To PM^[M].N - 1 Do Begin
  If (My >= Y) And (My < Y + ItemH(M, I)) Then Begin HitItem := I; Exit; End;
  Y := Y + ItemH(M, I);
 End;
End;

Procedure OpenSubOf(Li: Integer);
Var M, I, Sub, Y, X : Integer;
Begin
 M := Lv[Li].Menu; I := Lv[Li].Sel;
 If I < 0 Then Exit;
 Sub := PM^[M].It[I].Sub;
 If (Sub < 0) Or Not ItemEnabled(M, I) Then Exit;
 If PM^[Sub].N = 0 Then Exit;
 Y := Lv[Li].Y + 2 + ItemTop(M, I) - 3;
 X := Lv[Li].X + Lv[Li].W - 3;
 If X + MenuWidth(Sub) + 3 > SW Then X := Lv[Li].X - MenuWidth(Sub) + 3;
 OpenLevel(Sub, X, Y);
End;

Procedure BarLayout;
Var I, X: Integer;
Begin
 X := 6;
 For I := 0 To NUM_BAR - 1 Do Begin
  BarX[I] := X;
  BarW[I] := TextWidth(PM^[PM_B_FILE + I].Title) + 14;
  X := X + BarW[I];
 End;
End;

{ Boucle modale. Renvoie le code de commande choisi (0 si annule).
  BarIdx >= 0 : ouverture depuis la barre de menus. }
Function RunMenu(Root, X, Y, BarIdx: Integer): Integer;
Var
 Result_, I, HL, HI, Moved0X, Moved0Y, K, M, Cmd, NewBar : Integer;
 Done, Moved, Ext : Boolean;
 Ch : Char;
 HoverT : LongInt;
 OldB : LongInt;
 Top : Integer;
Begin
 Result_ := 0;
 RefreshMenus;
 NLv := 0;
 BarCur := BarIdx;
 OpenLevel(Root, X, Y);
 PollMouse;
 Moved0X := MouseX; Moved0Y := MouseY;
 OldB := 0;
 If BarIdx < 0 Then OldB := MouseB;
 If BarIdx >= 0 Then OldB := MouseB;
 Moved := False;
 Done := False;
 HoverT := TickMs;
 FlushKeys;
 While Not Done Do Begin
  PollMouse;
  If (Abs(MouseX - Moved0X) > 3) Or (Abs(MouseY - Moved0Y) > 3) Then Moved := True;
  { survol de la barre de menus : changer de menu }
  If (BarIdx >= 0) And (MouseY >= TITLE_H) And (MouseY < TITLE_H + MENU_H) Then Begin
   NewBar := -1;
   For I := 0 To NUM_BAR - 1 Do
    If (MouseX >= BarX[I]) And (MouseX < BarX[I] + BarW[I]) Then NewBar := I;
   If (NewBar >= 0) And (NewBar <> BarCur) And Moved Then Begin
    CloseLevelsFrom(0);
    BarCur := NewBar;
    BarIdx := NewBar;
    OpenLevel(PM_B_FILE + NewBar, BarX[NewBar], TITLE_H + MENU_H);
    HoverT := TickMs;
   End;
  End;
  { niveau survole }
  HL := -1; HI := -1;
  For K := NLv - 1 DownTo 0 Do Begin
   HI := HitItem(K, MouseX, MouseY);
   If (HI >= 0) Or PtIn(MouseX, MouseY, Lv[K].X, Lv[K].Y, Lv[K].X + Lv[K].W, Lv[K].Y + Lv[K].H) Then Begin
    HL := K; Break;
   End;
  End;
  If HL >= 0 Then Begin
   M := Lv[HL].Menu;
   If (HI >= 0) And (PM^[M].It[HI].Cap = '-') Then HI := -1;
   If Lv[HL].Sel <> HI Then Begin
    If NLv > HL + 1 Then CloseLevelsFrom(HL + 1);
    I := Lv[HL].Sel;
    Lv[HL].Sel := HI;
    If I >= 0 Then DrawLevelItem(HL, I);
    If HI >= 0 Then DrawLevelItem(HL, HI);
    HoverT := TickMs;
   End;
   { ouverture differee du sous-menu }
   If (HI >= 0) And (NLv = HL + 1) And (PM^[M].It[HI].Sub >= 0) And
      (ElapsedMs(HoverT) > 120) Then Begin
    HideMouse; OpenSubOf(HL); ShowMouse;
   End;
  End;
  { clic }
  If (MouseB <> 0) And (OldB = 0) Then Begin
   If (HL >= 0) And (HI >= 0) Then Begin
    M := Lv[HL].Menu;
    If ItemEnabled(M, HI) Then Begin
     If PM^[M].It[HI].Sub >= 0 Then Begin
      If NLv = HL + 1 Then Begin HideMouse; OpenSubOf(HL); ShowMouse; End;
     End
     Else Begin Result_ := PM^[M].It[HI].Cmd; Done := True; End;
    End;
   End
   Else If HL < 0 Then Begin
    Done := True;
   End;
  End;
  { relachement apres un deplacement : choix }
  If (MouseB = 0) And (OldB <> 0) And Moved Then Begin
   If (HL >= 0) And (HI >= 0) Then Begin
    M := Lv[HL].Menu;
    If ItemEnabled(M, HI) And (PM^[M].It[HI].Sub < 0) Then Begin
     Result_ := PM^[M].It[HI].Cmd; Done := True;
    End;
   End;
  End;
  OldB := MouseB;
  { clavier }
  If KeyPressed And Not Done Then Begin
   Ch := GetKey(Ext);
   Top := NLv - 1;
   If Ext Then Begin
    M := Lv[Top].Menu;
    Case Ch Of
     #80, #72: Begin
      K := Lv[Top].Sel;
      I := 0;
      Repeat
       If Ch = #80 Then Begin Inc(K); If K >= PM^[M].N Then K := 0; End
       Else Begin Dec(K); If K < 0 Then K := PM^[M].N - 1; End;
       Inc(I);
      Until ((PM^[M].It[K].Cap <> '-') And ItemEnabled(M, K)) Or (I > PM^[M].N);
      I := Lv[Top].Sel;
      Lv[Top].Sel := K;
      If I >= 0 Then DrawLevelItem(Top, I);
      DrawLevelItem(Top, K);
     End;
     #77: Begin
      K := Lv[Top].Sel;
      If (K >= 0) And (PM^[M].It[K].Sub >= 0) Then Begin
       HideMouse; OpenSubOf(Top); ShowMouse;
       If NLv > Top + 1 Then Begin
        Lv[NLv - 1].Sel := 0;
        While (Lv[NLv - 1].Sel < PM^[Lv[NLv - 1].Menu].N - 1) And
              ((PM^[Lv[NLv - 1].Menu].It[Lv[NLv - 1].Sel].Cap = '-') Or
               Not ItemEnabled(Lv[NLv - 1].Menu, Lv[NLv - 1].Sel)) Do
         Inc(Lv[NLv - 1].Sel);
        DrawLevel(NLv - 1);
       End;
      End
      Else If BarIdx >= 0 Then Begin
       CloseLevelsFrom(0);
       BarCur := (BarCur + 1) Mod NUM_BAR; BarIdx := BarCur;
       OpenLevel(PM_B_FILE + BarCur, BarX[BarCur], TITLE_H + MENU_H);
      End;
     End;
     #75: Begin
      If NLv > 1 Then CloseLevelsFrom(NLv - 1)
      Else If BarIdx >= 0 Then Begin
       CloseLevelsFrom(0);
       BarCur := (BarCur + NUM_BAR - 1) Mod NUM_BAR; BarIdx := BarCur;
       OpenLevel(PM_B_FILE + BarCur, BarX[BarCur], TITLE_H + MENU_H);
      End;
     End;
    End;
   End
   Else If Ch = #27 Then Begin
    If NLv > 1 Then CloseLevelsFrom(NLv - 1) Else Done := True;
   End
   Else If Ch = #13 Then Begin
    M := Lv[Top].Menu; K := Lv[Top].Sel;
    If K >= 0 Then Begin
     If PM^[M].It[K].Sub >= 0 Then Begin
      HideMouse; OpenSubOf(Top); ShowMouse;
     End
     Else If ItemEnabled(M, K) Then Begin
      Result_ := PM^[M].It[K].Cmd; Done := True;
     End;
    End;
   End;
  End;
  If Not Done Then Delay(10);
 End;
 CloseLevelsFrom(0);
 WaitRelease;
 Cmd := Result_;
 RunMenu := Cmd;
End;

{ --------------------------------------------------------------------- }
{ Selection d'atomes : expressions de type Jmol                         }
{ --------------------------------------------------------------------- }

Const
 NUM_ENAMES = 25;
 EnNames : Array[1..NUM_ENAMES] Of String[12] = (
  'hydrogen', 'helium', 'lithium', 'beryllium', 'boron', 'carbon',
  'nitrogen', 'oxygen', 'fluorine', 'neon', 'sodium', 'magnesium',
  'aluminum', 'silicon', 'phosphorus', 'sulfur', 'chlorine', 'argon',
  'potassium', 'calcium', 'iron', 'copper', 'zinc', 'bromine', 'iodine');
 EnZ : Array[1..NUM_ENAMES] Of Byte = (
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
  26, 29, 30, 35, 53);

Var
 SStk   : Array[0..7] Of PByteArr;
 SExpr  : String;
 SPos   : Integer;
 STok   : String;
 STokK  : Integer;       { 0 fin, 1 mot, 2 nombre, 3 symbole }
 SErr   : String;

Function InCur(I: Integer): Boolean;
Begin
 InCur := (CurModel = 0) Or (Atoms^[I].Model = CurModel);
End;

Procedure SNext;
Var C: Char;
Begin
 While (SPos <= Length(SExpr)) And (SExpr[SPos] = ' ') Do Inc(SPos);
 STok := ''; STokK := 0;
 If SPos > Length(SExpr) Then Exit;
 C := SExpr[SPos];
 If ((C >= 'a') And (C <= 'z')) Or (C = '_') Or ((C >= 'A') And (C <= 'Z')) Then Begin
  STokK := 1;
  While (SPos <= Length(SExpr)) And
        (((SExpr[SPos] >= 'a') And (SExpr[SPos] <= 'z')) Or
         ((SExpr[SPos] >= 'A') And (SExpr[SPos] <= 'Z')) Or
         ((SExpr[SPos] >= '0') And (SExpr[SPos] <= '9')) Or (SExpr[SPos] = '_')) Do Begin
   STok := STok + LoCh(SExpr[SPos]); Inc(SPos);
  End;
 End
 Else If (C >= '0') And (C <= '9') Then Begin
  STokK := 2;
  While (SPos <= Length(SExpr)) And (((SExpr[SPos] >= '0') And (SExpr[SPos] <= '9')) Or
        (SExpr[SPos] = '.')) Do Begin
   STok := STok + SExpr[SPos]; Inc(SPos);
  End;
 End
 Else Begin
  STokK := 3; STok := C; Inc(SPos);
  If (SPos <= Length(SExpr)) And (SExpr[SPos] = '=') And
     ((C = '<') Or (C = '>') Or (C = '!')) Then Begin
   STok := STok + '='; Inc(SPos);
  End;
 End;
End;

Procedure ParseOr(D: Integer); Forward;

Function MatchResRange(Const Ch: String; N1, N2: Integer; D: Integer): Boolean;
Var I: Integer;
Begin
 For I := 0 To NAtoms - 1 Do
  SStk[D]^[I] := Ord((Atoms^[I].ResSeq >= N1) And (Atoms^[I].ResSeq <= N2) And
                     ((Ch = '') Or (UpCh(Atoms^[I].Chain) = UpCh(Ch[1]))));
 MatchResRange := True;
End;

Function SSOf(I: Integer): Char;
Begin
 If Atoms^[I].ResIdx >= 0 Then SSOf := Resi^[Atoms^[I].ResIdx].SS Else SSOf := ' ';
End;

Function IsBackboneAtom(I: Integer): Boolean;
Var N: String;
Begin
 N := UpStr(TrimStr(Atoms^[I].Name));
 IsBackboneAtom := (Atoms^[I].Kind = ATOMK_PROT) And
                   ((N = 'N') Or (N = 'CA') Or (N = 'C') Or (N = 'O') Or (N = 'OXT'));
End;

Function ParseCompare(Const Prop: String; D: Integer): Boolean;
Var
 Op : String;
 V : Real;
 I : Integer;
 Sv : String;
 Ok, Res : Boolean;
 X : Real;
Begin
 ParseCompare := True;
 SNext;
 Op := STok;
 If (Op = '=') Then Begin
  { '==' tolere }
  If (SPos <= Length(SExpr)) And (SExpr[SPos] = '=') Then Inc(SPos);
 End;
 SNext;
 Sv := STok;
 If (Prop = 'chain') Or (Prop = 'atomname') Or (Prop = 'resname') Or (Prop = 'element') Then Begin
  For I := 0 To NAtoms - 1 Do Begin
   Res := False;
   If Prop = 'chain' Then Res := UpStr(Atoms^[I].Chain) = UpStr(Sv)
   Else If Prop = 'atomname' Then Res := UpStr(TrimStr(Atoms^[I].Name)) = UpStr(Sv)
   Else If Prop = 'resname' Then Res := UpStr(TrimStr(Atoms^[I].Res)) = UpStr(Sv)
   Else Res := UpStr(ElemDB[Atoms^[I].ElemZ].Sym) = UpStr(Sv);
   If Op = '!=' Then Res := Not Res;
   SStk[D]^[I] := Ord(Res);
  End;
  Exit;
 End;
 { propriete numerique }
 Ok := True;
 Sv := Sv;
 If Not StrToReal(STok, V) Then Begin Ok := False; V := 0; End;
 If Not Ok Then Begin SErr := 'nombre attendu apres ' + Prop; ParseCompare := False; Exit; End;
 For I := 0 To NAtoms - 1 Do Begin
  If Prop = 'atomno' Then X := Atoms^[I].Serial
  Else If Prop = 'resno' Then X := Atoms^[I].ResSeq
  Else If Prop = 'elemno' Then X := Atoms^[I].ElemZ
  Else If Prop = 'model' Then X := Atoms^[I].Model
  Else If Prop = 'charge' Then X := Atoms^[I].Charge
  Else X := Atoms^[I].Temp;
  Res := False;
  If Op = '=' Then Res := X = V
  Else If Op = '<' Then Res := X < V
  Else If Op = '>' Then Res := X > V
  Else If Op = '<=' Then Res := X <= V
  Else If Op = '>=' Then Res := X >= V
  Else If Op = '!=' Then Res := X <> V;
  SStk[D]^[I] := Ord(Res);
 End;
End;

Procedure ParsePrim(D: Integer);
Var
 I, J, K, N1, N2 : Integer;
 W, Ch : String;
 Dist, Dx, Dy, Dz : Real;
 Found : Boolean;
 Z : Byte;
 Hit : Boolean;
Begin
 If SErr <> '' Then Exit;
 If STokK = 0 Then Begin SErr := 'expression incomplete'; Exit; End;
 If (STok = '(') Or (STok = '{') Then Begin
  SNext;
  ParseOr(D);
  If (STok = ')') Or (STok = '}') Then SNext
  Else If SErr = '' Then SErr := 'parenthese fermante attendue';
  Exit;
 End;
 If STok = '*' Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := 1;
  SNext;
  If STok = ':' Then Begin
   SNext; Ch := STok; SNext;
   For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(UpStr(Atoms^[I].Chain) = UpStr(Ch));
  End
  Else If STok = '.' Then Begin
   SNext; W := STok; SNext;
   For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(UpStr(TrimStr(Atoms^[I].Name)) = UpStr(W));
  End;
  Exit;
 End;
 If STok = ':' Then Begin
  SNext; Ch := STok; SNext;
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(UpStr(Atoms^[I].Chain) = UpStr(Ch));
  Exit;
 End;
 If STok = '.' Then Begin
  SNext; W := STok; SNext;
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(UpStr(TrimStr(Atoms^[I].Name)) = UpStr(W));
  Exit;
 End;
 If STok = '[' Then Begin
  W := '';
  While (SPos <= Length(SExpr)) And (SExpr[SPos] <> ']') Do Begin W := W + SExpr[SPos]; Inc(SPos); End;
  Inc(SPos);
  SNext;
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(UpStr(TrimStr(Atoms^[I].Res)) = UpStr(TrimStr(W)));
  Exit;
 End;
 If STokK = 2 Then Begin
  N1 := StrToIntDef(STok, 0); N2 := N1;
  SNext;
  If STok = '-' Then Begin
   SNext;
   If STokK = 2 Then Begin N2 := StrToIntDef(STok, N1); SNext; End;
  End;
  Ch := '';
  If STok = ':' Then Begin SNext; Ch := STok; SNext; End;
  MatchResRange(Ch, N1, N2, D);
  Exit;
 End;
 If STokK <> 1 Then Begin SErr := 'terme inattendu "' + STok + '"'; Exit; End;
 W := STok;
 SNext;
 { comparaisons }
 If ((W = 'atomno') Or (W = 'resno') Or (W = 'elemno') Or (W = 'model') Or
     (W = 'charge') Or (W = 'temperature') Or (W = 'chain') Or (W = 'atomname') Or
     (W = 'resname') Or (W = 'element')) And
    ((STok = '=') Or (STok = '<') Or (STok = '>') Or (STok = '<=') Or
     (STok = '>=') Or (STok = '!=')) Then Begin
  { ParseCompare relit l'operateur : reculer d'un jeton }
  SPos := SPos - Length(STok);
  If ParseCompare(W, D) Then SNext;
  Exit;
 End;
 If W = 'within' Then Begin
  If STok <> '(' Then Begin SErr := '"(" attendu apres within'; Exit; End;
  SNext;
  If Not StrToReal(STok, Dist) Then Begin SErr := 'distance attendue'; Exit; End;
  SNext;
  If STok = ',' Then SNext;
  ParseOr(D + 1);
  If STok = ')' Then SNext Else If SErr = '' Then SErr := '")" attendu';
  If SErr <> '' Then Exit;
  For I := 0 To NAtoms - 1 Do Begin
   Hit := False;
   If InCur(I) Then
    For J := 0 To NAtoms - 1 Do Begin
     If SStk[D + 1]^[J] = 0 Then Continue;
     If Atoms^[I].Model <> Atoms^[J].Model Then Continue;
     Dx := Atoms^[I].X - Atoms^[J].X; If AbsR(Dx) > Dist Then Continue;
     Dy := Atoms^[I].Y - Atoms^[J].Y; If AbsR(Dy) > Dist Then Continue;
     Dz := Atoms^[I].Z - Atoms^[J].Z; If AbsR(Dz) > Dist Then Continue;
     If Dx * Dx + Dy * Dy + Dz * Dz <= Dist * Dist Then Begin Hit := True; Break; End;
    End;
   SStk[D]^[I] := Ord(Hit);
  End;
  Exit;
 End;
 For I := 0 To NAtoms - 1 Do SStk[D]^[I] := 0;
 If (W = 'all') Or (W = 'tout') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := 1;
 End
 Else If (W = 'none') Or (W = 'rien') Or (W = 'aucun') Then Begin
 End
 Else If (W = 'protein') Or (W = 'proteine') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Kind = ATOMK_PROT);
 End
 Else If (W = 'nucleic') Or (W = 'nucleique') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Kind = ATOMK_NUCL);
 End
 Else If (W = 'dna') Or (W = 'rna') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Kind = ATOMK_NUCL);
 End
 Else If (W = 'water') Or (W = 'solvent') Or (W = 'eau') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Kind = ATOMK_WATER);
 End
 Else If (W = 'hetero') Or (W = 'heteroatomes') Then Begin
  For I := 0 To NAtoms - 1 Do
   SStk[D]^[I] := Ord((Atoms^[I].Kind = ATOMK_LIGAND) Or (Atoms^[I].Kind = ATOMK_WATER) Or
                      (Atoms^[I].Kind = ATOMK_HETERO));
 End
 Else If W = 'ligand' Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Kind = ATOMK_LIGAND);
 End
 Else If (W = 'hydrogen') Or (W = 'hydrogens') Or (W = 'hydrogene') Or (W = 'hydrogenes') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].ElemZ = 1);
 End
 Else If (W = 'helix') Or (W = 'helice') Or (W = 'helices') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(SSOf(I) = 'H');
 End
 Else If (W = 'sheet') Or (W = 'feuillet') Or (W = 'feuillets') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(SSOf(I) = 'E');
 End
 Else If (W = 'turn') Or (W = 'coude') Or (W = 'coudes') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(SSOf(I) = 'T');
 End
 Else If (W = 'backbone') Or (W = 'squelette') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(IsBackboneAtom(I));
 End
 Else If (W = 'sidechain') Or (W = 'chaine_laterale') Then Begin
  For I := 0 To NAtoms - 1 Do
   SStk[D]^[I] := Ord((Atoms^[I].Kind = ATOMK_PROT) And Not IsBackboneAtom(I));
 End
 Else If (W = 'selected') Or (W = 'selection') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Sel);
 End
 Else If (W = 'hidden') Or (W = 'cache') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].Hid);
 End
 Else If (W = 'visible') Or (W = 'displayed') Or (W = 'affiche') Then Begin
  For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Not Atoms^[I].Hid);
 End
 Else Begin
  Found := False;
  { symbole d'element prefixe par _ }
  If (Length(W) >= 2) And (W[1] = '_') Then Begin
   Z := ElemBySymbol(Copy(W, 2, 2));
   If Z > 0 Then Begin
    For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].ElemZ = Z);
    Found := True;
   End;
  End;
  { nom d'element (anglais / francais) }
  If Not Found Then Begin
   Z := 0;
   For K := 1 To NUM_ENAMES Do If EnNames[K] = W Then Z := EnZ[K];
   For K := 1 To 103 Do If LoStr(ElemDB[K].Name) = W Then Z := K;
   If Z > 0 Then Begin
    For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].ElemZ = Z);
    Found := True;
   End;
  End;
  { nom de residu }
  If (Not Found) And (Length(W) <= 4) Then Begin
   For I := 0 To NAtoms - 1 Do
    If LoStr(TrimStr(Atoms^[I].Res)) = W Then Begin SStk[D]^[I] := 1; Found := True; End;
  End;
  { symbole d'element }
  If (Not Found) And (Length(W) <= 2) Then Begin
   Z := ElemBySymbol(W);
   If (Z > 0) And (LoStr(ElemDB[Z].Sym) = W) Then Begin
    For I := 0 To NAtoms - 1 Do SStk[D]^[I] := Ord(Atoms^[I].ElemZ = Z);
    Found := True;
   End;
  End;
  If Not Found Then SErr := 'terme inconnu "' + W + '"';
 End;
End;

Procedure ParseNot(D: Integer);
Var I: Integer;
Begin
 If (STok = 'not') Or (STok = '!') Or (STok = 'non') Then Begin
  SNext;
  ParseNot(D);
  If SErr = '' Then
   For I := 0 To NAtoms - 1 Do SStk[D]^[I] := 1 - SStk[D]^[I];
 End
 Else ParsePrim(D);
End;

Procedure ParseAnd(D: Integer);
Var I: Integer;
Begin
 ParseNot(D);
 While (SErr = '') And ((STok = 'and') Or (STok = '&') Or (STok = 'et')) Do Begin
  SNext;
  ParseNot(D + 1);
  If SErr = '' Then
   For I := 0 To NAtoms - 1 Do SStk[D]^[I] := SStk[D]^[I] And SStk[D + 1]^[I];
 End;
End;

Procedure ParseOr(D: Integer);
Var I: Integer;
Begin
 ParseAnd(D);
 While (SErr = '') And ((STok = 'or') Or (STok = '|') Or (STok = ',') Or (STok = 'ou')) Do Begin
  SNext;
  ParseAnd(D + 1);
  If SErr = '' Then
   For I := 0 To NAtoms - 1 Do SStk[D]^[I] := SStk[D]^[I] Or SStk[D + 1]^[I];
 End;
End;

{ Evalue une expression : resultat dans SStk[0] (tous les modeles) }
Function EvalSel(Const Expr: String): Boolean;
Begin
 SExpr := Expr; SPos := 1; SErr := '';
 SNext;
 If STokK = 0 Then Begin SErr := 'expression vide'; EvalSel := False; Exit; End;
 ParseOr(0);
 If (SErr = '') And (STokK <> 0) Then SErr := 'terme inattendu "' + STok + '"';
 If SErr <> '' Then Begin EvalSel := False; Exit; End;
 EvalSel := True;
End;

{ Selectionne les atomes de l'expression ; renvoie le nombre d'atomes }
Function SelectExprFn(Const Expr: String): Integer;
Var I, N: Integer;
Begin
 SelectExprFn := -1;
 If Not Loaded Then Begin SErr := 'aucune molecule chargee'; Exit; End;
 If Not EvalSel(Expr) Then Exit;
 N := 0;
 For I := 0 To NAtoms - 1 Do Begin
  Atoms^[I].Sel := SStk[0]^[I] <> 0;
  If Atoms^[I].Sel Then Inc(N);
 End;
 NSelCount := N;
 NeedScene := True;
 SelectExprFn := N;
End;

{ --------------------------------------------------------------------- }
{ Application de styles a la cible (selection ou toute l'image)         }
{ --------------------------------------------------------------------- }

Var
 TgtSel : Boolean;

Procedure BeginTarget;
Begin
 TgtSel := CountSelected > 0;
End;

Function InTarget(I: Integer): Boolean;
Begin
 If TgtSel Then InTarget := Atoms^[I].Sel
 Else InTarget := True;
End;

Function BondInTarget(B: Integer): Boolean;
Begin
 BondInTarget := InTarget(Bonds^[B].A) And InTarget(Bonds^[B].B);
End;

Procedure ApplyScheme(Cm: Integer);
Var I: Integer;
Begin
 BeginTarget;
 For I := 0 To NAtoms - 1 Do
  If InTarget(I) Then
   Case Cm Of
    CM_SCH_CPK  : Atoms^[I].Sf := 100;
    CM_SCH_BALL : Atoms^[I].Sf := 23;
    CM_SCH_STICK: Atoms^[I].Sf := -15;
    CM_SCH_WIRE : Atoms^[I].Sf := 0;
   End;
 For I := 0 To NBonds - 1 Do
  If BondInTarget(I) Then
   Case Cm Of
    CM_SCH_CPK  : Bonds^[I].BR := 0;
    CM_SCH_BALL : Bonds^[I].BR := 15;
    CM_SCH_STICK: Bonds^[I].BR := 15;
    CM_SCH_WIRE : Bonds^[I].BR := 1;
   End;
 NeedScene := True;
End;

Procedure ApplyAtomPct(Pct: Integer);
Var I: Integer;
Begin
 BeginTarget;
 For I := 0 To NAtoms - 1 Do
  If InTarget(I) Then Atoms^[I].Sf := Pct;
 NeedScene := True;
End;

Procedure ApplyBondRad(R: Integer);
Var I: Integer;
Begin
 BeginTarget;
 For I := 0 To NBonds - 1 Do
  If BondInTarget(I) Then Bonds^[I].BR := R;
 NeedScene := True;
End;

Procedure ApplyStruct(Rib: Integer);
Var R, Ca: Integer;
Begin
 BeginTarget;
 For R := 0 To NRes - 1 Do Begin
  Ca := Resi^[R].CA;
  If Ca < 0 Then Continue;
  If InTarget(Ca) Then Resi^[R].Rib := Rib;
 End;
 NeedScene := True;
End;

Function LabelIndexFor(Const Fmt: String): Integer;
Var I: Integer;
Begin
 LabelIndexFor := 0;
 If Fmt = '' Then Exit;
 For I := 1 To NLblFmt Do
  If LblFmt[I] = Fmt Then Begin LabelIndexFor := I; Exit; End;
 If NLblFmt < MAX_LBLFMT Then Begin
  Inc(NLblFmt);
  LblFmt[NLblFmt] := Copy(Fmt, 1, 24);
  LabelIndexFor := NLblFmt;
 End
 Else Begin
  LblFmt[MAX_LBLFMT] := Copy(Fmt, 1, 24);
  LabelIndexFor := MAX_LBLFMT;
 End;
End;

Procedure ApplyLabel(Const Fmt: String);
Var I, K: Integer;
Begin
 BeginTarget;
 K := LabelIndexFor(Fmt);
 For I := 0 To NAtoms - 1 Do
  If InTarget(I) Then Atoms^[I].Lbl := K;
 NeedScene := True;
End;

Procedure ApplyDots(On: Boolean; AllAtoms: Boolean);
Var I: Integer;
Begin
 BeginTarget;
 If AllAtoms Then TgtSel := False;
 For I := 0 To NAtoms - 1 Do
  If InTarget(I) Then Atoms^[I].Dots := On;
 DotsDirty := True;
 NeedScene := True;
End;

Procedure ApplyColor(C: Byte);
Var I: Integer;
Begin
 BeginTarget;
 For I := 0 To NAtoms - 1 Do
  If InTarget(I) Then Atoms^[I].Col := C;
 NeedScene := True;
End;

Function AminoColor(Const R: String; Kind: Byte): Byte;
Var U: String;
Begin
 U := UpStr(TrimStr(R));
 AminoColor := LightGray;
 If Kind = ATOMK_NUCL Then Begin
  If (U = 'A') Or (U = 'DA') Then AminoColor := LightRed
  Else If (U = 'G') Or (U = 'DG') Then AminoColor := LightGreen
  Else If (U = 'C') Or (U = 'DC') Then AminoColor := LightBlue
  Else AminoColor := Yellow;
  Exit;
 End;
 If (U = 'ALA') Or (U = 'VAL') Or (U = 'LEU') Or (U = 'ILE') Or (U = 'MET') Or
    (U = 'PRO') Then AminoColor := LightGreen
 Else If (U = 'SER') Or (U = 'THR') Or (U = 'ASN') Or (U = 'GLN') Then AminoColor := LightCyan
 Else If (U = 'ASP') Or (U = 'GLU') Then AminoColor := LightRed
 Else If (U = 'LYS') Or (U = 'ARG') Or (U = 'HIS') Then AminoColor := LightBlue
 Else If (U = 'PHE') Or (U = 'TRP') Or (U = 'TYR') Then AminoColor := LightMagenta
 Else If U = 'CYS' Then AminoColor := Yellow
 Else If U = 'GLY' Then AminoColor := White;
End;

Function RainbowOf(T: Real): Byte;
Var K: Integer;
Begin
 K := Trunc(T * (NUM_RAINBOW - 1) + 0.5);
 If K < 0 Then K := 0;
 If K > NUM_RAINBOW - 1 Then K := NUM_RAINBOW - 1;
 RainbowOf := Rainbow[K];
End;

{ Position relative (0..1) d'un residu dans sa chaine }
Procedure ChainFractions(Var Fr: TRealArr);
Var A, B, R, Cnt: Integer;
Begin
 R := 0;
 While R < NRes Do Begin
  A := R; B := R;
  While (B + 1 < NRes) And (Resi^[B + 1].Chain = Resi^[A].Chain) And
        (Resi^[B + 1].Model = Resi^[A].Model) Do Inc(B);
  Cnt := B - A;
  For R := A To B Do
   If Cnt > 0 Then Fr[R] := (R - A) / Cnt Else Fr[R] := 0;
  R := B + 1;
 End;
End;

Var
 ChainAll  : Array[0..25] Of Char;
 NChainAll : Integer;

Procedure BuildChainAll;
Var I, K: Integer; F: Boolean;
Begin
 NChainAll := 0;
 For I := 0 To NRes - 1 Do Begin
  F := False;
  For K := 0 To NChainAll - 1 Do If ChainAll[K] = Resi^[I].Chain Then F := True;
  If (Not F) And (NChainAll < 26) Then Begin
   ChainAll[NChainAll] := Resi^[I].Chain;
   Inc(NChainAll);
  End;
 End;
End;

Function ChainIndexOf(Ch: Char): Integer;
Var K: Integer;
Begin
 ChainIndexOf := 0;
 For K := 0 To NChainAll - 1 Do If ChainAll[K] = Ch Then ChainIndexOf := K;
End;

Procedure ApplyColorScheme(Sch: Integer);
Var
 I, R : Integer;
 Fr : PRealArr;
 Mn, Mx : Real;
 C : Byte;
Begin
 BeginTarget;
 If Sch = 0 Then Begin
  For I := 0 To NAtoms - 1 Do If InTarget(I) Then Atoms^[I].Col := 255;
  NeedScene := True;
  Exit;
 End;
 New(Fr);
 BuildChainAll;
 ChainFractions(Fr^);
 Mn := 1E9; Mx := -1E9;
 For I := 0 To NAtoms - 1 Do Begin
  If Atoms^[I].Temp < Mn Then Mn := Atoms^[I].Temp;
  If Atoms^[I].Temp > Mx Then Mx := Atoms^[I].Temp;
 End;
 If Mx - Mn < 1E-6 Then Mx := Mn + 1;
 For I := 0 To NAtoms - 1 Do Begin
  If Not InTarget(I) Then Continue;
  R := Atoms^[I].ResIdx;
  Case Sch Of
   1: Begin
    If R >= 0 Then Atoms^[I].Col := ChainCols[ChainIndexOf(Atoms^[I].Chain) Mod 8];
   End;
   2: Atoms^[I].Col := AminoColor(Atoms^[I].Res, Atoms^[I].Kind);
   3: If R >= 0 Then Begin
    Case Resi^[R].SS Of
     'H': C := Magenta;
     'E': C := Yellow;
     'T': C := LightBlue;
    Else
     C := White;
    End;
    Atoms^[I].Col := C;
   End;
   4: If (R >= 0) And (Resi^[R].IsAA Or Resi^[R].IsNA) Then Atoms^[I].Col := RainbowOf(Fr^[R]);
   5: Atoms^[I].Col := RainbowOf((Atoms^[I].Temp - Mn) / (Mx - Mn));
  End;
 End;
 Dispose(Fr);
 NeedScene := True;
End;

Procedure ApplyRibColor(Sch: Integer; C: Integer);
Var
 R, Ca : Integer;
 Fr : PRealArr;
Begin
 BeginTarget;
 New(Fr);
 BuildChainAll;
 ChainFractions(Fr^);
 For R := 0 To NRes - 1 Do Begin
  Ca := Resi^[R].CA;
  If (Ca < 0) Or Not InTarget(Ca) Then Continue;
  If Sch = 0 Then Resi^[R].RCol := 255
  Else If Sch = 1 Then Resi^[R].RCol := RainbowOf(Fr^[R])
  Else If Sch = 2 Then Resi^[R].RCol := ChainCols[ChainIndexOf(Resi^[R].Chain) Mod 8]
  Else Resi^[R].RCol := C;
 End;
 Dispose(Fr);
 NeedScene := True;
End;

{ --------------------------------------------------------------------- }
{ Modeles, vue, mesures                                                 }
{ --------------------------------------------------------------------- }

Procedure SetModel(M: Integer);
Begin
 If M < 0 Then M := 0;
 If M > NModels Then M := NModels;
 CurModel := M;
 DotsDirty := True;
 NeedScene := True;
 NeedChrome := True;
End;

Function FindSerial(S: LongInt): Integer;
Var I: Integer;
Begin
 FindSerial := -1;
 For I := 0 To NAtoms - 1 Do
  If (Atoms^[I].Serial = S) And InCur(I) Then Begin FindSerial := I; Exit; End;
End;

Procedure AddMeasure(N: Integer; A0, A1, A2, A3: Integer);
Var I: Integer;
Begin
 If NMeas >= MAX_MEAS Then Begin ConsoleErr('trop de mesures'); Exit; End;
 Inc(NMeas);
 Meas[NMeas].N := N;
 Meas[NMeas].A[0] := A0; Meas[NMeas].A[1] := A1;
 Meas[NMeas].A[2] := A2; Meas[NMeas].A[3] := A3;
 { ne pas dupliquer }
 For I := 1 To NMeas - 1 Do
  If (Meas[I].N = N) And (Meas[I].A[0] = A0) And (Meas[I].A[1] = A1) And
     (Meas[I].A[2] = A2) And (Meas[I].A[3] = A3) Then Begin Dec(NMeas); Exit; End;
 NeedScene := True;
End;

Procedure DeleteMeasures(N: Integer);
Var I, J: Integer;
Begin
 J := 0;
 For I := 1 To NMeas Do
  If (N <> 0) And (Meas[I].N <> N) Then Begin Inc(J); Meas[J] := Meas[I]; End;
 NMeas := J;
 NeedScene := True;
End;

{ Description d'un atome pour la barre d'etat }
Function AtomInfo(I: Integer): String;
Var S: String;
Begin
 If (I < 0) Or (I >= NAtoms) Then Begin AtomInfo := ''; Exit; End;
 If (Atoms^[I].ResIdx >= 0) And (Atoms^[I].Res <> 'UNL') Then Begin
  S := '[' + TrimStr(Atoms^[I].Res) + ']' + IntToStr(Atoms^[I].ResSeq);
  If Atoms^[I].Chain <> ' ' Then S := S + ':' + Atoms^[I].Chain;
  S := S + '.' + TrimStr(Atoms^[I].Name);
 End
 Else S := TrimStr(Atoms^[I].Name) + IntToStr(I + 1);
 S := S + ' #' + IntToStr(Atoms^[I].Serial) + ' ' + RealToStr(Atoms^[I].X, 3) + ' ' +
      RealToStr(Atoms^[I].Y, 3) + ' ' + RealToStr(Atoms^[I].Z, 3);
 AtomInfo := S;
End;

Procedure CenterOnTarget;
Var
 I, N : Integer;
 Sx, Sy, Sz : Real;
Begin
 N := 0; Sx := 0; Sy := 0; Sz := 0;
 BeginTarget;
 For I := 0 To NAtoms - 1 Do
  If InTarget(I) And InCur(I) Then Begin
   Sx := Sx + WX^[I]; Sy := Sy + WY^[I]; Sz := Sz + WZ^[I]; Inc(N);
  End;
 If N = 0 Then Exit;
 { translater : le nouveau centre de rotation est le centre des atomes cibles }
 CenX := Sx / N; CenY := Sy / N; CenZ := Sz / N;
 PanX := 0; PanY := 0;
 NeedScene := True;
End;

{ Image : positions de travail pour les calculs hors rendu }
Procedure UpdateWork;
Var I: Integer;
Begin
 For I := 0 To NAtoms - 1 Do Begin
  WX^[I] := Atoms^[I].X; WY^[I] := Atoms^[I].Y; WZ^[I] := Atoms^[I].Z;
 End;
End;

{ --------------------------------------------------------------------- }
{ Enregistrement : PDB, XYZ, MOL, image BMP                             }
{ --------------------------------------------------------------------- }

Function PdbName(I: Integer): String;
Var N: String; Sy: String;
Begin
 N := TrimStr(Atoms^[I].Name);
 Sy := ElemDB[Atoms^[I].ElemZ].Sym;
 If (Length(Sy) = 1) And (Length(N) < 4) Then N := ' ' + N;
 PdbName := FixedStr(N, 4);
End;

Function SavePDBFile(Const FN: String): Boolean;
Var
 F : Text;
 I, M : Integer;
 Rec, Ln : String;
 Multi : Boolean;
Begin
 SavePDBFile := False;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 WriteLn(F, 'HEADER    ' + Copy(FixedStr(UpStr(MolTitle), 40), 1, 40));
 WriteLn(F, 'TITLE     ' + Copy(MolTitle, 1, 60));
 If HasCell Then
  WriteLn(F, 'CRYST1' + RealFmt(CellP[0], 9, 3) + RealFmt(CellP[1], 9, 3) +
             RealFmt(CellP[2], 9, 3) + RealFmt(CellP[3], 7, 2) +
             RealFmt(CellP[4], 7, 2) + RealFmt(CellP[5], 7, 2) + ' P 1           1');
 Multi := NModels > 1;
 For M := 1 To NModels Do Begin
  If Multi Then WriteLn(F, 'MODEL     ' + PadL(IntToStr(M), 4));
  For I := 0 To NAtoms - 1 Do Begin
   If Atoms^[I].Model <> M Then Continue;
   If (Atoms^[I].Kind = ATOMK_PROT) Or (Atoms^[I].Kind = ATOMK_NUCL) Then Rec := 'ATOM  '
   Else Rec := 'HETATM';
   Ln := Rec + PadL(IntToStr(Atoms^[I].Serial), 5) + ' ' + PdbName(I) + ' ' +
         PadL(TrimStr(Atoms^[I].Res), 3) + ' ' + Atoms^[I].Chain +
         PadL(IntToStr(Atoms^[I].ResSeq), 4) + Atoms^[I].ICode + '   ' +
         RealFmt(Atoms^[I].X, 8, 3) + RealFmt(Atoms^[I].Y, 8, 3) +
         RealFmt(Atoms^[I].Z, 8, 3) + '  1.00' + RealFmt(Atoms^[I].Temp, 6, 2) +
         '          ' + PadL(ElemDB[Atoms^[I].ElemZ].Sym, 2);
   WriteLn(F, Ln);
  End;
  If Multi Then WriteLn(F, 'ENDMDL');
 End;
 WriteLn(F, 'END');
 Close(F);
 SavePDBFile := True;
End;

Function SaveXYZFile(Const FN: String): Boolean;
Var
 F : Text;
 I, M, N : Integer;
Begin
 SaveXYZFile := False;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 For M := 1 To NModels Do Begin
  N := 0;
  For I := 0 To NAtoms - 1 Do If Atoms^[I].Model = M Then Inc(N);
  WriteLn(F, N);
  If Models^[M].Name <> '' Then WriteLn(F, Models^[M].Name) Else WriteLn(F, MolTitle);
  For I := 0 To NAtoms - 1 Do
   If Atoms^[I].Model = M Then
    WriteLn(F, FixedStr(ElemDB[Atoms^[I].ElemZ].Sym, 3), ' ', RealFmt(Atoms^[I].X, 12, 6),
            ' ', RealFmt(Atoms^[I].Y, 12, 6), ' ', RealFmt(Atoms^[I].Z, 12, 6));
 End;
 Close(F);
 SaveXYZFile := True;
End;

Function SaveMOLFileOut(Const FN: String): Boolean;
Var
 F : Text;
 I, B, M, First, Cnt, NB, A1, A2 : Integer;
Begin
 SaveMOLFileOut := False;
 M := MaxI(CurModel, 1);
 First := Models^[M].First;
 If First < 0 Then Exit;
 Cnt := Models^[M].Last - First + 1;
 NB := 0;
 For B := 0 To NBonds - 1 Do
  If (Bonds^[B].A >= First) And (Bonds^[B].A <= Models^[M].Last) And (Bonds^[B].Ord <= 4) Then Inc(NB);
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 WriteLn(F, Copy(MolTitle, 1, 70));
 WriteLn(F, '  JMOL.PAS  3D');
 WriteLn(F, '');
 WriteLn(F, PadL(IntToStr(Cnt), 3), PadL(IntToStr(NB), 3), '  0  0  0  0  0  0  0  0999 V2000');
 For I := First To Models^[M].Last Do
  WriteLn(F, RealFmt(Atoms^[I].X, 10, 4), RealFmt(Atoms^[I].Y, 10, 4), RealFmt(Atoms^[I].Z, 10, 4),
          ' ', FixedStr(ElemDB[Atoms^[I].ElemZ].Sym, 3), ' 0  0  0  0  0  0  0  0  0  0  0  0');
 For B := 0 To NBonds - 1 Do
  If (Bonds^[B].A >= First) And (Bonds^[B].A <= Models^[M].Last) And (Bonds^[B].Ord <= 4) Then Begin
   A1 := Bonds^[B].A - First + 1; A2 := Bonds^[B].B - First + 1;
   WriteLn(F, PadL(IntToStr(A1), 3), PadL(IntToStr(A2), 3), PadL(IntToStr(Bonds^[B].Ord), 3), '  0');
  End;
 WriteLn(F, 'M  END');
 Close(F);
 SaveMOLFileOut := True;
End;

{ Image 4 bits par pixel (palette VGA standard) }
Function SaveBMPFile(Const FN: String): Boolean;
Const
 Pal : Array[0..15, 0..2] Of Byte = (
  (0, 0, 0), (0, 0, 170), (0, 170, 0), (0, 170, 170), (170, 0, 0), (170, 0, 170),
  (170, 85, 0), (170, 170, 170), (85, 85, 85), (85, 85, 255), (85, 255, 85),
  (85, 255, 255), (255, 85, 85), (255, 85, 255), (255, 255, 85), (255, 255, 255));
Var
 F : File;
 W, H, RowB, X, Y, I : Integer;
 Hdr : Array[0..53] Of Byte;
 Row : Array[0..511] Of Byte;
 Sz, Off : LongInt;
 P : Array[0..3] Of Byte;
 C1, C2 : Byte;
Begin
 SaveBMPFile := False;
 W := CvX2 - CvX1 + 1; H := CvY2 - CvY1 + 1;
 If W > 1000 Then W := 1000;
 RowB := (((W + 1) Div 2) + 3) And (Not 3);
 If RowB > 512 Then Exit;
 Off := 54 + 64;
 Sz := Off + LongInt(RowB) * H;
 FillChar(Hdr, SizeOf(Hdr), 0);
 Hdr[0] := Ord('B'); Hdr[1] := Ord('M');
 Hdr[2] := Sz And 255; Hdr[3] := (Sz Shr 8) And 255; Hdr[4] := (Sz Shr 16) And 255; Hdr[5] := (Sz Shr 24) And 255;
 Hdr[10] := Off And 255; Hdr[11] := (Off Shr 8) And 255;
 Hdr[14] := 40;
 Hdr[18] := W And 255; Hdr[19] := (W Shr 8) And 255;
 Hdr[22] := H And 255; Hdr[23] := (H Shr 8) And 255;
 Hdr[26] := 1; Hdr[28] := 4;
 Hdr[46] := 16;
 {$I-}
 Assign(F, FN); Rewrite(F, 1);
 {$I+}
 If IOResult <> 0 Then Exit;
 BlockWrite(F, Hdr, 54);
 For I := 0 To 15 Do Begin
  P[0] := Pal[I, 2]; P[1] := Pal[I, 1]; P[2] := Pal[I, 0]; P[3] := 0;
  BlockWrite(F, P, 4);
 End;
 HideMouse;
 For Y := H - 1 DownTo 0 Do Begin
  FillChar(Row, SizeOf(Row), 0);
  X := 0;
  While X < W Do Begin
   C1 := GetPixel(CvX1 + X, CvY1 + Y) And 15;
   If X + 1 < W Then C2 := GetPixel(CvX1 + X + 1, CvY1 + Y) And 15 Else C2 := 0;
   Row[X Div 2] := (C1 Shl 4) Or C2;
   Inc(X, 2);
  End;
  BlockWrite(F, Row, RowB);
 End;
 ShowMouse;
 Close(F);
 SaveBMPFile := True;
End;

{ --------------------------------------------------------------------- }
{ Calculs : formule, masse, centre de masse, giration, dimensions       }
{ --------------------------------------------------------------------- }

Procedure FrameRange(Var F, L: Integer);
Var M: Integer;
Begin
 F := 0; L := NAtoms - 1;
 M := CurModel;
 If M = 0 Then M := 1;
 If (M >= 1) And (M <= NModels) And (Models^[M].First >= 0) Then Begin
  F := Models^[M].First; L := Models^[M].Last;
 End;
End;

Function FormulaStr: String;
Var
 Cnt : Array[0..103] Of Integer;
 I, F, L, Z : Integer;
 S : String;

 Procedure Put(Zz: Integer);
 Begin
  If Cnt[Zz] > 0 Then Begin
   S := S + ElemDB[Zz].Sym;
   If Cnt[Zz] > 1 Then S := S + IntToStr(Cnt[Zz]);
   Cnt[Zz] := 0;
  End;
 End;

 Procedure PutAlpha;
 Var Best, K : Integer;
 Begin
  Repeat
   Best := -1;
   For K := 1 To 103 Do
    If Cnt[K] > 0 Then
     If (Best < 0) Or (ElemDB[K].Sym < ElemDB[Best].Sym) Then Best := K;
   If Best > 0 Then Put(Best);
  Until Best < 0;
 End;

Begin
 FrameRange(F, L);
 FillChar(Cnt, SizeOf(Cnt), 0);
 For I := F To L Do Begin
  Z := Atoms^[I].ElemZ; If Z > 103 Then Z := 0;
  Inc(Cnt[Z]);
 End;
 S := '';
 If Cnt[6] > 0 Then Begin Put(6); Put(1); End;
 PutAlpha;
 FormulaStr := S;
End;

Function TotalMass: Real;
Var I, F, L: Integer; M: Real;
Begin
 FrameRange(F, L);
 M := 0;
 For I := F To L Do M := M + ElemDB[Atoms^[I].ElemZ].Mass;
 TotalMass := M;
End;

Procedure CenterOfMass(Var X, Y, Z: Real);
Var I, F, L: Integer; M, T: Real;
Begin
 FrameRange(F, L);
 X := 0; Y := 0; Z := 0; T := 0;
 For I := F To L Do Begin
  M := ElemDB[Atoms^[I].ElemZ].Mass;
  X := X + M * Atoms^[I].X; Y := Y + M * Atoms^[I].Y; Z := Z + M * Atoms^[I].Z;
  T := T + M;
 End;
 If T > 0 Then Begin X := X / T; Y := Y / T; Z := Z / T; End;
End;

Function GyrationRadius: Real;
Var I, F, L: Integer; X, Y, Z, S: Real;
Begin
 FrameRange(F, L);
 CenterOfMass(X, Y, Z);
 S := 0;
 For I := F To L Do
  S := S + Sqr(Atoms^[I].X - X) + Sqr(Atoms^[I].Y - Y) + Sqr(Atoms^[I].Z - Z);
 If L >= F Then GyrationRadius := Sqrt(S / (L - F + 1)) Else GyrationRadius := 0;
End;

{ Remplit TextBuf avec une ligne }
Procedure TB(Const S: String);
Begin
 If TextN < MAX_TEXT Then Begin TextBuf[TextN] := Copy(S, 1, 78); Inc(TextN); End;
End;

Procedure DoCalc(Cm: Integer);
Var
 X, Y, Z, Mnx, Mxx, Mny, Mxy, Mnz, Mxz : Real;
 I, F, L : Integer;
Begin
 If Not Loaded Then Begin ConsoleErr('aucune molecule chargee'); Exit; End;
 Case Cm Of
  CM_CALC_FORM: Begin
   ConsoleSay('Formule moleculaire : ' + FormulaStr);
   ShowMessage('Formule moleculaire', FormulaStr);
  End;
  CM_CALC_MASS: Begin
   ConsoleSay('Masse molaire : ' + RealToStr(TotalMass, 3) + ' g/mol');
   ShowMessage('Masse molaire', FormulaStr + '|' + RealToStr(TotalMass, 3) + ' g/mol');
  End;
  CM_CALC_COM: Begin
   CenterOfMass(X, Y, Z);
   ConsoleSay('Centre de masse : ' + RealToStr(X, 3) + ' ' + RealToStr(Y, 3) + ' ' + RealToStr(Z, 3));
   ShowMessage('Centre de masse', 'X = ' + RealToStr(X, 4) + ' ' + #143 +
               '|Y = ' + RealToStr(Y, 4) + ' ' + #143 +
               '|Z = ' + RealToStr(Z, 4) + ' ' + #143);
  End;
  CM_CALC_GYR: Begin
   ConsoleSay('Rayon de giration : ' + RealToStr(GyrationRadius, 3) + ' ' + #143);
   ShowMessage('Rayon de giration', RealToStr(GyrationRadius, 3) + ' ' + #143);
  End;
  CM_CALC_DIM: Begin
   FrameRange(F, L);
   Mnx := 1E9; Mny := 1E9; Mnz := 1E9; Mxx := -1E9; Mxy := -1E9; Mxz := -1E9;
   For I := F To L Do Begin
    If Atoms^[I].X < Mnx Then Mnx := Atoms^[I].X; If Atoms^[I].X > Mxx Then Mxx := Atoms^[I].X;
    If Atoms^[I].Y < Mny Then Mny := Atoms^[I].Y; If Atoms^[I].Y > Mxy Then Mxy := Atoms^[I].Y;
    If Atoms^[I].Z < Mnz Then Mnz := Atoms^[I].Z; If Atoms^[I].Z > Mxz Then Mxz := Atoms^[I].Z;
   End;
   ConsoleSay('Dimensions : ' + RealToStr(Mxx - Mnx, 2) + ' x ' +
              RealToStr(Mxy - Mny, 2) + ' x ' + RealToStr(Mxz - Mnz, 2) + ' ' + #143);
   ShowMessage('Dimensions', 'Largeur : ' + RealToStr(Mxx - Mnx, 3) + ' ' + #143 +
               '|Hauteur : ' + RealToStr(Mxy - Mny, 3) + ' ' + #143 +
               '|Profondeur : ' + RealToStr(Mxz - Mnz, 3) + ' ' + #143);
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Minimisation d'energie (champ de force simplifie)                     }
{ --------------------------------------------------------------------- }

Type
 TPairRec = Record A, B : Integer; D : Real; End;
 TPairArr = Array[0..MAX_PAIRS-1] Of TPairRec;

Procedure Minimize;
Var
 Gx, Gy, Gz : PRealArr;
 Pairs : ^TPairArr;
 NP, I, J, K, F, L, It, A, B, C : Integer;
 D, D0, Dx, Dy, Dz, Fc, E, E0, Th, Da, Db, Step, Mx, Dmin, Rc : Real;
 Hyb : Real;
 Nonb : Boolean;
 BondD0 : ^TRealArr;
Begin
 If Not Loaded Then Begin ConsoleErr('aucune molecule chargee'); Exit; End;
 FrameRange(F, L);
 If L - F + 1 > 600 Then Begin
  ConsoleErr('minimisation limitee a 600 atomes');
  Exit;
 End;
 New(Gx); New(Gy); New(Gz); New(Pairs); New(BondD0);
 { longueurs de liaison ideales }
 For K := 0 To NBonds - 1 Do Begin
  A := Bonds^[K].A; B := Bonds^[K].B;
  If (A < F) Or (A > L) Then Begin BondD0^[K] := 0; Continue; End;
  D0 := ElemDB[Atoms^[A].ElemZ].Radius + ElemDB[Atoms^[B].ElemZ].Radius;
  If Bonds^[K].Ord = 2 Then D0 := D0 - 0.12
  Else If Bonds^[K].Ord = 3 Then D0 := D0 - 0.2
  Else If Bonds^[K].Ord = 4 Then D0 := D0 - 0.06;
  BondD0^[K] := D0;
 End;
 { paires 1-3 }
 NP := 0;
 For I := 0 To NBonds - 1 Do Begin
  If (Bonds^[I].A < F) Or (Bonds^[I].A > L) Then Continue;
  For J := I + 1 To NBonds - 1 Do Begin
   If (Bonds^[J].A < F) Or (Bonds^[J].A > L) Then Continue;
   C := -1; A := -1; B := -1;
   If Bonds^[I].A = Bonds^[J].A Then Begin C := Bonds^[I].A; A := Bonds^[I].B; B := Bonds^[J].B; End
   Else If Bonds^[I].A = Bonds^[J].B Then Begin C := Bonds^[I].A; A := Bonds^[I].B; B := Bonds^[J].A; End
   Else If Bonds^[I].B = Bonds^[J].A Then Begin C := Bonds^[I].B; A := Bonds^[I].A; B := Bonds^[J].B; End
   Else If Bonds^[I].B = Bonds^[J].B Then Begin C := Bonds^[I].B; A := Bonds^[I].A; B := Bonds^[J].A; End;
   If (C < 0) Or (NP >= MAX_PAIRS) Then Continue;
   { hybridation du centre }
   Hyb := 109.5;
   If (Bonds^[I].Ord >= 2) Or (Bonds^[J].Ord >= 2) Then Hyb := 120.0;
   If (Bonds^[I].Ord = 3) Or (Bonds^[J].Ord = 3) Then Hyb := 180.0;
   If (Bonds^[I].Ord = 4) Or (Bonds^[J].Ord = 4) Then Hyb := 120.0;
   If (Atoms^[C].ElemZ = 8) And (Hyb < 110) Then Hyb := 104.5;
   Da := BondD0^[I]; Db := BondD0^[J];
   Th := Hyb * DEG2RAD;
   Pairs^[NP].A := A; Pairs^[NP].B := B;
   Pairs^[NP].D := Sqrt(Da * Da + Db * Db - 2 * Da * Db * Cos(Th));
   Inc(NP);
  End;
 End;
 E0 := 0; E := 0; Step := 0.04;
 For It := 1 To 250 Do Begin
  For I := F To L Do Begin Gx^[I] := 0; Gy^[I] := 0; Gz^[I] := 0; End;
  E := 0;
  For K := 0 To NBonds - 1 Do Begin
   A := Bonds^[K].A; B := Bonds^[K].B;
   If (A < F) Or (A > L) Or (BondD0^[K] <= 0) Then Continue;
   Dx := Atoms^[B].X - Atoms^[A].X; Dy := Atoms^[B].Y - Atoms^[A].Y; Dz := Atoms^[B].Z - Atoms^[A].Z;
   D := Sqrt(Dx * Dx + Dy * Dy + Dz * Dz); If D < 0.01 Then D := 0.01;
   E := E + Sqr(D - BondD0^[K]);
   Fc := 2.0 * (D - BondD0^[K]) / D;
   Gx^[A] := Gx^[A] - Fc * Dx; Gy^[A] := Gy^[A] - Fc * Dy; Gz^[A] := Gz^[A] - Fc * Dz;
   Gx^[B] := Gx^[B] + Fc * Dx; Gy^[B] := Gy^[B] + Fc * Dy; Gz^[B] := Gz^[B] + Fc * Dz;
  End;
  For K := 0 To NP - 1 Do Begin
   A := Pairs^[K].A; B := Pairs^[K].B;
   Dx := Atoms^[B].X - Atoms^[A].X; Dy := Atoms^[B].Y - Atoms^[A].Y; Dz := Atoms^[B].Z - Atoms^[A].Z;
   D := Sqrt(Dx * Dx + Dy * Dy + Dz * Dz); If D < 0.01 Then D := 0.01;
   E := E + 0.5 * Sqr(D - Pairs^[K].D);
   Fc := 0.5 * (D - Pairs^[K].D) / D;
   Gx^[A] := Gx^[A] - Fc * Dx; Gy^[A] := Gy^[A] - Fc * Dy; Gz^[A] := Gz^[A] - Fc * Dz;
   Gx^[B] := Gx^[B] + Fc * Dx; Gy^[B] := Gy^[B] + Fc * Dy; Gz^[B] := Gz^[B] + Fc * Dz;
  End;
  { repulsion entre atomes non lies }
  If L - F < 400 Then
   For A := F To L - 1 Do
    For B := A + 1 To L Do Begin
     Dx := Atoms^[B].X - Atoms^[A].X; If AbsR(Dx) > 3.5 Then Continue;
     Dy := Atoms^[B].Y - Atoms^[A].Y; If AbsR(Dy) > 3.5 Then Continue;
     Dz := Atoms^[B].Z - Atoms^[A].Z; If AbsR(Dz) > 3.5 Then Continue;
     D := Sqrt(Dx * Dx + Dy * Dy + Dz * Dz); If D < 0.01 Then D := 0.01;
     Dmin := 0.62 * (ElemDB[Atoms^[A].ElemZ].VdW + ElemDB[Atoms^[B].ElemZ].VdW);
     If D >= Dmin Then Continue;
     Nonb := True;
     For K := 0 To NBonds - 1 Do
      If ((Bonds^[K].A = A) And (Bonds^[K].B = B)) Or ((Bonds^[K].A = B) And (Bonds^[K].B = A)) Then Nonb := False;
     If Nonb Then
      For K := 0 To NP - 1 Do
       If ((Pairs^[K].A = A) And (Pairs^[K].B = B)) Or ((Pairs^[K].A = B) And (Pairs^[K].B = A)) Then Nonb := False;
     If Not Nonb Then Continue;
     E := E + 0.2 * Sqr(Dmin - D);
     Fc := -0.4 * (Dmin - D) / D;
     Gx^[A] := Gx^[A] - Fc * Dx; Gy^[A] := Gy^[A] - Fc * Dy; Gz^[A] := Gz^[A] - Fc * Dz;
     Gx^[B] := Gx^[B] + Fc * Dx; Gy^[B] := Gy^[B] + Fc * Dy; Gz^[B] := Gz^[B] + Fc * Dz;
    End;
  If It = 1 Then E0 := E;
  Mx := 0;
  For I := F To L Do Begin
   D := Sqrt(Gx^[I] * Gx^[I] + Gy^[I] * Gy^[I] + Gz^[I] * Gz^[I]);
   If D > Mx Then Mx := D;
  End;
  If Mx < 1E-4 Then Break;
  Rc := Step;
  If Mx * Rc > 0.15 Then Rc := 0.15 / Mx;
  For I := F To L Do Begin
   Atoms^[I].X := Atoms^[I].X - Rc * Gx^[I];
   Atoms^[I].Y := Atoms^[I].Y - Rc * Gy^[I];
   Atoms^[I].Z := Atoms^[I].Z - Rc * Gz^[I];
  End;
 End;
 Dispose(Gx); Dispose(Gy); Dispose(Gz); Dispose(Pairs); Dispose(BondD0);
 DotsDirty := True;
 NeedScene := True;
 ConsoleSay('Minimisation terminee : energie ' + RealToStr(E0, 3) + ' -> ' + RealToStr(E, 3));
End;

{ --------------------------------------------------------------------- }
{ Textes d'aide                                                         }
{ --------------------------------------------------------------------- }

Function OneLetter(Const R: String): Char;
Var U: String;
Begin
 U := UpStr(TrimStr(R));
 OneLetter := 'X';
 If U = 'ALA' Then OneLetter := 'A' Else If U = 'ARG' Then OneLetter := 'R'
 Else If U = 'ASN' Then OneLetter := 'N' Else If U = 'ASP' Then OneLetter := 'D'
 Else If U = 'CYS' Then OneLetter := 'C' Else If U = 'GLN' Then OneLetter := 'Q'
 Else If U = 'GLU' Then OneLetter := 'E' Else If U = 'GLY' Then OneLetter := 'G'
 Else If U = 'HIS' Then OneLetter := 'H' Else If U = 'ILE' Then OneLetter := 'I'
 Else If U = 'LEU' Then OneLetter := 'L' Else If U = 'LYS' Then OneLetter := 'K'
 Else If U = 'MET' Then OneLetter := 'M' Else If U = 'PHE' Then OneLetter := 'F'
 Else If U = 'PRO' Then OneLetter := 'P' Else If U = 'SER' Then OneLetter := 'S'
 Else If U = 'THR' Then OneLetter := 'T' Else If U = 'TRP' Then OneLetter := 'W'
 Else If U = 'TYR' Then OneLetter := 'Y' Else If U = 'VAL' Then OneLetter := 'V'
 Else If (U = 'A') Or (U = 'DA') Then OneLetter := 'A'
 Else If (U = 'G') Or (U = 'DG') Then OneLetter := 'G'
 Else If (U = 'C') Or (U = 'DC') Then OneLetter := 'C'
 Else If (U = 'T') Or (U = 'DT') Then OneLetter := 'T'
 Else If (U = 'U') Or (U = 'DU') Then OneLetter := 'U';
End;

Procedure HelpText(Kind: Integer);
Begin
 TextN := 0;
 Case Kind Of
  0: Begin
   TB('JMOL - Visualiseur moleculaire 3D (clone en Pascal)');
   TB('');
   TB('Fichier > Ouvrir pour charger une molecule : PDB, XYZ, MOL, SDF, MOL2, CIF.');
   TB('');
   TB('Bouton droit de la souris : menu contextuel (selection, style, couleur...).');
   TB('Bouton gauche + glisser    : rotation de la molecule.');
   TB('Bouton milieu + glisser    : zoom (ou + / - au clavier).');
   TB('Bouton droit + glisser     : deplacement (translation).');
   TB('Clic gauche sur un atome   : selon le mode "Choisir" (identifier, mesurer, ...).');
   TB('');
   TB('La barre d''etat affiche l''atome sous le curseur :');
   TB('   [RES]numero:chaine.atome #serie x y z');
   TB('');
   TB('Les menus Style et Couleur agissent sur la selection courante ;');
   TB('si rien n''est selectionne, ils agissent sur tous les atomes du modele.');
   TB('');
   TB('La console (F4) accepte les commandes du langage de script de Jmol :');
   TB('   select, color, spacefill, wireframe, cartoon, label, rotate, zoom...');
   TB('Tapez "aide" dans la console pour la liste complete.');
  End;
  1: Begin
   TB('SOURIS');
   TB('  Gauche (glisser)   Rotation');
   TB('  Gauche (clic)      Action du mode Choisir (identifier / mesurer / selectionner)');
   TB('  Milieu (glisser)   Zoom');
   TB('  Droit (glisser)    Translation');
   TB('  Droit (clic)       Menu contextuel');
   TB('');
   TB('CLAVIER');
   TB('  Fleches            Rotation de 5 degres');
   TB('  + / -  PgHaut/Bas  Zoom avant / arriere');
   TB('  Home               Reinitialiser la vue');
   TB('  Espace             Rotation automatique');
   TB('  < / >              Modele (image) precedent / suivant');
   TB('  H  A  B            Hydrogenes, Axes, Boite englobante');
   TB('  V                  Vibration');
   TB('  Ctrl+O  Ctrl+S     Ouvrir / enregistrer en PDB');
   TB('  Ctrl+A             Tout selectionner');
   TB('  F1                 Aide          F3 Recharger');
   TB('  F4                 Console       F10 Barre de menus');
   TB('  Echap              Ferme les menus et la console');
   TB('  Ctrl+Q             Quitter');
  End;
  2: Begin
   TB('COMMANDES DE LA CONSOLE (syntaxe de Jmol ; alias francais entre parentheses)');
   TB('');
   TB('load fichier                (charger)      Charge une molecule');
   TB('select expression           (selectionner) Selectionne des atomes');
   TB('   all, none, protein, nucleic, water, hetero, hydrogen, helix, sheet, turn,');
   TB('   backbone, sidechain, carbon, ala, 10-20:A, *:A, .CA, _Fe, atomno=5,');
   TB('   within(3.0, expr), not, and, or, ( )');
   TB('display / hide / restrict expression   (afficher / masquer / restreindre)');
   TB('color [cible] couleur|schema (couleur) Cibles : atoms, cartoon, background');
   TB('   Schemas : cpk, chain, group, structure, amino, temperature');
   TB('spacefill [n% | off | on]   (remplissage)  wireframe [on | off | rayon] (filaire)');
   TB('cartoon | ribbons | trace | backbone | strands [on|off]');
   TB('dots [on|off]               (points)       label [format|off] (etiquette)');
   TB('hbonds [on|off]             axes / boundbox / unitcell [on|off]');
   TB('rotate x|y|z degres         (tourner)      spin on|off|x|y|z [vitesse] (rotation)');
   TB('zoom [n|in|out|reset]       translate x|y pourcent');
   TB('center [expression]         (centrer)      reset (reinitialiser)');
   TB('background couleur          (fond)         measure a b [c [d]] (mesurer)');
   TB('measure delete|list         model n|next|prev|first|last|all (modele)');
   TB('animation on|off|loop|once|palindrome      vibration on|off|scale n');
   TB('vector on|off|scale n       set cle valeur   echo texte');
   TB('write pdb|xyz|mol|bmp fichier (ecrire)     show file|model|atoms|sequence|...');
   TB('script fichier              (executer)     minimize (minimiser)');
   TB('calculate hbonds|formula|mass              console on|off');
   TB('help (aide)                 quit (quitter) clear (effacer)');
   TB('');
   TB('Plusieurs commandes peuvent etre separees par un point-virgule.');
  End;
 Else
  Begin
   TB(PROG_TITLE);
   TB('');
   TB('Clone du visualiseur moleculaire Jmol, ecrit en Turbo Pascal 7 / Free Pascal');
   TB('(mode graphique, sans dependance externe).');
   TB('');
   TB('Auteur : Sylvain Maltais (support@gladir.com)');
   TB('Site   : https://www.gladir.com/qchimie');
   TB('');
   TB('Jmol est un projet libre (http://jmol.sourceforge.net).');
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Chargement d'un fichier (interface)                                   }
{ --------------------------------------------------------------------- }

Procedure AfterLoad;
Begin
 IdentityRM;
 ZoomF := 1.0; PanX := 0; PanY := 0;
 SpinOn := False; AnimOn := False; VibOn := False; VibPanel := False;
 ShowHB := False; ShowCell := False; ShowVec := False;
 PickMode := PICK_IDENT; PickN := 0;
 NSelCount := 0;
 NHB := 0;
 HoverAtom := -1;
 VibPhase := 0;
 SetupView;
 If AtomLimit Then
  ConsoleErr('limite de ' + IntToStr(MAX_ATOMS) + ' atomes atteinte : fichier tronque');
 If Loaded Then
  ConsoleSay('Fichier charge : ' + ExtractName(FileNameCur) + ' (' + IntToStr(NAtoms) +
             ' atomes, ' + IntToStr(NBonds) + ' liaisons, ' + IntToStr(NRes) +
             ' residus, ' + IntToStr(NModels) + ' modele(s))');
 NeedScene := True; NeedChrome := True;
End;

Function LoadFileUI(Const FN: String; Quiet: Boolean): Boolean;
Var Full: String;
Begin
 Full := FExpand(FN);
 LoadFileUI := False;
 If Not LoadMolecule(Full) Then Begin
  If Not Quiet Then ShowMessage('Erreur', 'Impossible de charger le fichier :|' + FN);
  ConsoleErr('impossible de charger ' + FN);
  NeedScene := True; NeedChrome := True;
  Exit;
 End;
 AfterLoad;
 LoadFileUI := True;
End;

Procedure EnsureHBonds;
Begin
 If (NHB = 0) And Loaded Then ComputeHBonds;
End;

Function WindowTitle: String;
Var S: String;
Begin
 If Not Loaded Then Begin WindowTitle := PROG_NAME; Exit; End;
 S := ExtractName(FileNameCur);
 If MolHeader <> '' Then S := S + ' - ' + MolHeader
 Else If MolTitle <> '' Then S := S + ' - ' + Copy(MolTitle, 1, 30);
 WindowTitle := S;
End;

{ --------------------------------------------------------------------- }
{ Animation, vibration, rotation                                        }
{ --------------------------------------------------------------------- }

Procedure FrameStep(Dir: Integer);
Var M: Integer;
Begin
 If NModels < 2 Then Exit;
 M := CurModel + Dir;
 If M < 1 Then M := NModels;
 If M > NModels Then M := 1;
 SetModel(M);
End;

Procedure AnimTick;
Var M: Integer;
Begin
 If NModels < 2 Then Begin AnimOn := False; Exit; End;
 M := CurModel;
 If M < 1 Then M := 1;
 If AnimDir = 0 Then AnimDir := 1;
 M := M + AnimDir;
 If M > NModels Then Begin
  If AnimMode = 0 Then Begin M := NModels; AnimOn := False; NeedChrome := True; End
  Else If AnimMode = 1 Then M := 1
  Else Begin AnimDir := -1; M := NModels - 1; End;
 End
 Else If M < 1 Then Begin
  If AnimMode = 2 Then Begin AnimDir := 1; M := 2; End
  Else M := 1;
 End;
 If M < 1 Then M := 1;
 If M > NModels Then M := NModels;
 SetModel(M);
End;

Procedure VibNextMode(Dir: Integer);
Var M, I: Integer;
Begin
 If Not HasVibData Then Exit;
 M := MaxI(CurModel, 1);
 For I := 1 To NModels Do Begin
  M := M + Dir;
  If M > NModels Then M := 1;
  If M < 1 Then M := NModels;
  If Models^[M].Vib Then Begin SetModel(M); Exit; End;
 End;
End;

Procedure SetVib(On: Boolean);
Begin
 If On And Not HasVibData Then Begin
  ConsoleErr('aucune donnee de vibration dans ce fichier');
  Exit;
 End;
 VibOn := On;
 If On And (CurModel = 0) Then SetModel(1);
 If On And (CurModel >= 1) And Not Models^[CurModel].Vib Then VibNextMode(1);
 NeedScene := True; NeedChrome := True;
End;

Procedure ApplyView(Cm: Integer);
Begin
 Case Cm Of
  CM_VIEW_FRONT : SetOrientation(1, 0);
  CM_VIEW_BACK  : SetOrientation(1, 180);
  CM_VIEW_LEFT  : SetOrientation(1, 90);
  CM_VIEW_RIGHT : SetOrientation(1, -90);
  CM_VIEW_TOP   : SetOrientation(0, 90);
  CM_VIEW_BOTTOM: SetOrientation(0, -90);
 End;
End;

Procedure SetZoom(Z: Real);
Begin
 ZoomF := ClampR(Z, 0.02, 60.0);
 SetupView;
 NeedScene := True;
End;

Procedure SetPick(M: Integer);
Begin
 PickMode := M;
 PickN := 0;
 Case M Of
  PICK_IDENT   : ConsoleSay('Mode : identification des atomes.');
  PICK_DIST    : ConsoleSay('Mesure de distance : cliquez sur 2 atomes.');
  PICK_ANGLE   : ConsoleSay('Mesure d''angle : cliquez sur 3 atomes.');
  PICK_TORSION : ConsoleSay('Mesure de diedre : cliquez sur 4 atomes.');
  PICK_SELATOM : Begin ConsoleSay('Selection d''atomes : cliquez sur un atome.'); ShowHalo := True; End;
  PICK_SELRES  : Begin ConsoleSay('Selection de residus : cliquez sur un atome.'); ShowHalo := True; End;
  PICK_SELCHAIN: Begin ConsoleSay('Selection de chaines : cliquez sur un atome.'); ShowHalo := True; End;
  PICK_CENTER  : ConsoleSay('Centrage : cliquez sur un atome.');
  PICK_LABEL   : ConsoleSay('Etiquettes : cliquez sur un atome.');
 End;
 NeedScene := True; NeedChrome := True;
End;

{ Clic gauche sur la zone 3D (atome sous le curseur ou -1) }
Procedure DoPick(A: Integer);
Var
 I, R, K : Integer;
 Nn : Integer;
Begin
 Case PickMode Of
  PICK_IDENT: If A >= 0 Then ConsoleSay(AtomInfo(A));
  PICK_DIST, PICK_ANGLE, PICK_TORSION: Begin
   If A < 0 Then Exit;
   Nn := PickMode + 1;
   PickQ[PickN] := A; Inc(PickN);
   If PickN >= Nn Then Begin
    Case Nn Of
     2: Begin
      AddMeasure(2, PickQ[0], PickQ[1], 0, 0);
      UpdateWork;
      ConsoleSay('Distance : ' + DistStr(DistW(PickQ[0], PickQ[1])));
     End;
     3: Begin
      AddMeasure(3, PickQ[0], PickQ[1], PickQ[2], 0);
      UpdateWork;
      ConsoleSay('Angle : ' + RealToStr(AngleW(PickQ[0], PickQ[1], PickQ[2]), 2) + #248);
     End;
    Else
     Begin
      AddMeasure(4, PickQ[0], PickQ[1], PickQ[2], PickQ[3]);
      UpdateWork;
      ConsoleSay('Diedre : ' + RealToStr(TorsionW(PickQ[0], PickQ[1], PickQ[2], PickQ[3]), 2) + #248);
     End;
    End;
    PickN := 0;
   End
   Else ConsoleSay('Atome ' + IntToStr(PickN) + '/' + IntToStr(Nn) + ' : ' + AtomInfo(A));
  End;
  PICK_SELATOM: Begin
   If A < 0 Then Begin
    For I := 0 To NAtoms - 1 Do Atoms^[I].Sel := False;
   End
   Else Atoms^[A].Sel := Not Atoms^[A].Sel;
  End;
  PICK_SELRES: If A >= 0 Then Begin
   R := Atoms^[A].ResIdx;
   If R >= 0 Then Begin
    K := Ord(Not Atoms^[A].Sel);
    For I := Resi^[R].First To Resi^[R].Last Do Atoms^[I].Sel := K <> 0;
   End
   Else Atoms^[A].Sel := Not Atoms^[A].Sel;
  End;
  PICK_SELCHAIN: If A >= 0 Then Begin
   K := Ord(Not Atoms^[A].Sel);
   For I := 0 To NAtoms - 1 Do
    If (Atoms^[I].Chain = Atoms^[A].Chain) And (Atoms^[I].Model = Atoms^[A].Model) Then
     Atoms^[I].Sel := K <> 0;
  End;
  PICK_CENTER: If A >= 0 Then Begin
   CenX := WX^[A]; CenY := WY^[A]; CenZ := WZ^[A];
   PanX := 0; PanY := 0;
  End;
  PICK_LABEL: If A >= 0 Then Begin
   If Atoms^[A].Lbl = 0 Then Atoms^[A].Lbl := LabelIndexFor('%a') Else Atoms^[A].Lbl := 0;
  End;
 End;
 NSelCount := CountSelected;
 NeedScene := True; NeedChrome := True;
End;

{ Atome sous le curseur (coordonnees ecran absolues) ou -1 }
Function AtomAt(Mx, My: Integer): Integer;
Var
 K, I, Best, R : Integer;
 Bz, Dx, Dy : Real;
 Rx, Ry : Integer;
Begin
 AtomAt := -1;
 If Not Loaded Then Exit;
 Rx := Mx - CvX1; Ry := My - CvY1;
 Best := -1; Bz := -1E30;
 For K := 0 To NAct - 1 Do Begin
  I := Act^[K];
  If PVis^[I] = 0 Then Continue;
  R := MaxI(PR^[I], 5);
  Dx := PX^[I] - Rx; Dy := PY^[I] - Ry;
  If Dx * Dx + Dy * Dy <= R * R Then
   If PZ^[I] > Bz Then Begin Bz := PZ^[I]; Best := I; End;
 End;
 AtomAt := Best;
End;

{ --------------------------------------------------------------------- }
{ Affichage d'informations                                              }
{ --------------------------------------------------------------------- }

Procedure ShowInfo(Cm: Integer);
Var
 I, K, F, L, N, M, C : Integer;
 S : String;
 Cnt : Array[0..103] Of Integer;
Begin
 If (Not Loaded) And (Cm <> CM_SHOW_ORIENT) Then Begin
  ConsoleErr('aucune molecule chargee');
  Exit;
 End;
 TextN := 0;
 FrameRange(F, L);
 Case Cm Of
  CM_SHOW_FILE, CM_FILEINFO: Begin
   TB('Fichier       : ' + FileNameCur);
   TB('Titre         : ' + MolTitle);
   TB('Format        : ' + MolFormat);
   If MolHeader <> '' Then TB('Code          : ' + MolHeader);
   TB('Atomes        : ' + IntToStr(NAtoms));
   TB('Liaisons      : ' + IntToStr(NBonds));
   TB('Residus       : ' + IntToStr(NRes));
   TB('Modeles       : ' + IntToStr(NModels));
   If HasVibData Then TB('Vibrations    : oui');
   If HasCell Then
    TB('Cellule       : a=' + RealToStr(CellP[0], 3) + ' b=' + RealToStr(CellP[1], 3) +
       ' c=' + RealToStr(CellP[2], 3));
   TB('Formule       : ' + FormulaStr);
   TB('Masse molaire : ' + RealToStr(TotalMass, 3) + ' g/mol');
   ShowTextDlg('Informations sur le fichier');
  End;
  CM_SHOW_MODEL: Begin
   M := MaxI(CurModel, 1);
   TB('Modele        : ' + IntToStr(M) + ' / ' + IntToStr(NModels));
   TB('Nom           : ' + Models^[M].Name);
   TB('Atomes        : ' + IntToStr(L - F + 1));
   If Models^[M].Vib Then TB('Frequence     : ' + RealToStr(Models^[M].Freq, 2) + ' cm^-1');
   TB('Formule       : ' + FormulaStr);
   TB('Masse molaire : ' + RealToStr(TotalMass, 3) + ' g/mol');
   ShowTextDlg('Informations sur le modele');
  End;
  CM_SHOW_ATOMS: Begin
   TB('Serie Nom  Res  Chn Num        X         Y         Z');
   For I := F To L Do Begin
    If TextN >= MAX_TEXT - 6 Then Begin TB('... (liste tronquee)'); Break; End;
    TB(PadL(IntToStr(Atoms^[I].Serial), 5) + ' ' + FixedStr(TrimStr(Atoms^[I].Name), 4) + ' ' +
       FixedStr(TrimStr(Atoms^[I].Res), 4) + ' ' + Atoms^[I].Chain + '  ' +
       PadL(IntToStr(Atoms^[I].ResSeq), 4) + RealFmt(Atoms^[I].X, 10, 3) +
       RealFmt(Atoms^[I].Y, 10, 3) + RealFmt(Atoms^[I].Z, 10, 3));
   End;
   ShowTextDlg('Liste des atomes');
  End;
  CM_SHOW_SEQ: Begin
   BuildChainAll;
   N := 0;
   For K := 0 To NChainAll - 1 Do Begin
    S := ''; C := 0;
    For I := 0 To NRes - 1 Do
     If (Resi^[I].Chain = ChainAll[K]) And (Resi^[I].IsAA Or Resi^[I].IsNA) And
        (Resi^[I].Model = MaxI(CurModel, 1)) Then Begin
      S := S + OneLetter(Resi^[I].Name); Inc(C);
     End;
    If C = 0 Then Continue;
    Inc(N);
    TB('Chaine ' + ChainAll[K] + ' (' + IntToStr(C) + ' residus) :');
    While S <> '' Do Begin
     TB('  ' + Copy(S, 1, 60));
     Delete(S, 1, 60);
    End;
   End;
   If N = 0 Then TB('Aucune sequence (pas de residus proteiques ou nucleiques).');
   ShowTextDlg('Sequence');
  End;
  CM_SHOW_CELL: Begin
   If Not HasCell Then TB('Pas de cellule unitaire dans ce fichier.')
   Else Begin
    TB('a = ' + RealToStr(CellP[0], 4) + ' ' + #143);
    TB('b = ' + RealToStr(CellP[1], 4) + ' ' + #143);
    TB('c = ' + RealToStr(CellP[2], 4) + ' ' + #143);
    TB('alpha = ' + RealToStr(CellP[3], 3) + #248);
    TB('beta  = ' + RealToStr(CellP[4], 3) + #248);
    TB('gamma = ' + RealToStr(CellP[5], 3) + #248);
   End;
   ShowTextDlg('Cellule unitaire');
  End;
  CM_SHOW_ORIENT: Begin
   TB('Matrice de rotation :');
   For I := 0 To 2 Do
    TB('  ' + RealFmt(RM[I, 0], 9, 5) + RealFmt(RM[I, 1], 9, 5) + RealFmt(RM[I, 2], 9, 5));
   TB('Zoom        : ' + RealToStr(ZoomF * 100, 0) + '%');
   TB('Translation : ' + IntToStr(PanX) + ', ' + IntToStr(PanY) + ' pixels');
   TB('Centre      : ' + RealToStr(CenX, 3) + ' ' + RealToStr(CenY, 3) + ' ' + RealToStr(CenZ, 3));
   ShowTextDlg('Orientation');
  End;
  CM_SHOW_SEL: Begin
   FillChar(Cnt, SizeOf(Cnt), 0);
   N := 0;
   For I := 0 To NAtoms - 1 Do
    If Atoms^[I].Sel Then Begin Inc(N); If Atoms^[I].ElemZ <= 103 Then Inc(Cnt[Atoms^[I].ElemZ]); End;
   TB('Atomes selectionnes : ' + IntToStr(N));
   For I := 0 To 103 Do
    If Cnt[I] > 0 Then TB('   ' + ElemDB[I].Name + ' (' + ElemDB[I].Sym + ') : ' + IntToStr(Cnt[I]));
   ShowTextDlg('Selection');
  End;
 End;
End;

Procedure ShowMeasureList;
Var I: Integer;
Begin
 UpdateWork;
 TextN := 0;
 If NMeas = 0 Then TB('Aucune mesure.')
 Else
  For I := 1 To NMeas Do Begin
   Case Meas[I].N Of
    2: TB('Distance : #' + IntToStr(Atoms^[Meas[I].A[0]].Serial) + ' - #' +
          IntToStr(Atoms^[Meas[I].A[1]].Serial) + '   ' + MeasValueStr(Meas[I]));
    3: TB('Angle    : #' + IntToStr(Atoms^[Meas[I].A[0]].Serial) + ' - #' +
          IntToStr(Atoms^[Meas[I].A[1]].Serial) + ' - #' +
          IntToStr(Atoms^[Meas[I].A[2]].Serial) + '   ' + MeasValueStr(Meas[I]));
   Else
    TB('Diedre   : #' + IntToStr(Atoms^[Meas[I].A[0]].Serial) + ' - #' +
       IntToStr(Atoms^[Meas[I].A[1]].Serial) + ' - #' +
       IntToStr(Atoms^[Meas[I].A[2]].Serial) + ' - #' +
       IntToStr(Atoms^[Meas[I].A[3]].Serial) + '   ' + MeasValueStr(Meas[I]));
   End;
  End;
 ShowTextDlg('Liste des mesures');
End;

{ --------------------------------------------------------------------- }
{ Macros                                                                }
{ --------------------------------------------------------------------- }

Procedure ResetAll;
Var I: Integer;
Begin
 For I := 0 To NAtoms - 1 Do Begin
  Atoms^[I].Sel := False; Atoms^[I].Hid := False; Atoms^[I].Dots := False;
  Atoms^[I].Col := 255; Atoms^[I].Lbl := 0;
 End;
 ApplyDefaultStyle;
 NMeas := 0; PickN := 0;
 ShowHydro := True; ShowWater := True; ShowHB := False; ShowHalo := False;
 ShowAxes := False; ShowBBox := False; ShowCell := False; ShowVec := False;
 Persp := True; SpinOn := False;
 BgCol := Black;
 DotsDirty := True;
 ResetView;
End;

Procedure SelectAllQuiet;
Var I: Integer;
Begin
 For I := 0 To NAtoms - 1 Do Atoms^[I].Sel := False;
 NSelCount := 0;
End;

Procedure RunMacro(K: Integer);
Var I: Integer;
Begin
 If Not Loaded Then Begin ConsoleErr('aucune molecule chargee'); Exit; End;
 SelectAllQuiet;
 Case K Of
  0: Begin ApplyScheme(CM_SCH_BALL); ApplyStruct(RIB_NONE); End;
  1: Begin ApplyScheme(CM_SCH_CPK); ApplyStruct(RIB_NONE); End;
  2: Begin ApplyScheme(CM_SCH_WIRE); ApplyStruct(RIB_NONE); End;
  3: Begin
   ApplyScheme(CM_SCH_WIRE); ApplyAtomPct(0); ApplyBondRad(0);
   ApplyStruct(RIB_CARTOON); ApplyRibColor(0, 0);
   BeginTarget;
  End;
  4: Begin
   ApplyScheme(CM_SCH_WIRE); ApplyAtomPct(0); ApplyBondRad(0);
   ApplyStruct(RIB_CARTOON); ApplyRibColor(1, 0);
  End;
  5: Begin
   ApplyScheme(CM_SCH_BALL);
   ApplyStruct(RIB_CARTOON); ApplyRibColor(0, 0);
   For I := 0 To NAtoms - 1 Do
    If (Atoms^[I].Kind = ATOMK_PROT) And IsBackboneAtom(I) Then Begin
     Atoms^[I].Sf := 0;
    End;
  End;
  6: Begin
   ApplyScheme(CM_SCH_WIRE);
   ApplyDots(True, True);
  End;
  7: Begin
   ApplyScheme(CM_SCH_WIRE); ApplyAtomPct(0); ApplyBondRad(0);
   ApplyStruct(RIB_CARTOON); ApplyRibColor(2, 0);
  End;
  8: Begin
   ComputeHBonds; ShowHB := True;
   ConsoleSay(IntToStr(NHB) + ' liaisons hydrogene.');
  End;
  9: Begin SpinOn := True; SpinAxis := 1; End;
  10: ResetAll;
 End;
 NeedScene := True; NeedChrome := True;
End;

{ --------------------------------------------------------------------- }
{ Dispatcheur de commandes (menus, barre d'outils, clavier)             }
{ --------------------------------------------------------------------- }

Procedure ExecLine(Const Line: String); Forward;
Procedure WaitMs(Ms: LongInt); Forward;

Procedure DoOpen;
Var FN: String;
Begin
 FN := '';
 If FileOpenDialog('Ouvrir un fichier', FN, '*.*') Then LoadFileUI(FN, False);
End;

Procedure DoSave(Cm: Integer);
Var
 FN, Def : String;
 Ok : Boolean;
Begin
 If Not Loaded Then Begin ConsoleErr('aucune molecule chargee'); Exit; End;
 Def := ExtractName(FileNameCur);
 Case Cm Of
  CM_SAVEPDB: Def := 'SORTIE.PDB';
  CM_SAVEXYZ: Def := 'SORTIE.XYZ';
  CM_SAVEMOL: Def := 'SORTIE.MOL';
  CM_SAVEBMP: Def := 'JMOL.BMP';
 End;
 FN := Def;
 If Not FileSaveDialog('Enregistrer', FN) Then Exit;
 Case Cm Of
  CM_SAVEPDB: Ok := SavePDBFile(FN);
  CM_SAVEXYZ: Ok := SaveXYZFile(FN);
  CM_SAVEMOL: Ok := SaveMOLFileOut(FN);
 Else
  Begin
   RenderCanvas;
   Ok := SaveBMPFile(FN);
  End;
 End;
 If Ok Then ConsoleSay('Fichier enregistre : ' + FN)
 Else ShowMessage('Erreur', 'Impossible d''ecrire le fichier :|' + FN);
 NeedScene := True;
End;

Procedure DoCmd(Cm: Integer);
Var
 S : String;
 I, N : Integer;
Begin
 If Cm = 0 Then Exit;
 Case Cm Of
  CM_OPEN      : DoOpen;
  CM_RELOAD    : If Loaded Then LoadFileUI(FileNameCur, False);
  CM_SAVEPDB, CM_SAVEXYZ, CM_SAVEMOL, CM_SAVEBMP: DoSave(Cm);
  CM_RUNSCRIPT : Begin
   S := '';
   If FileOpenDialog('Executer un script', S, '*.*') Then ExecLine('script ' + S);
  End;
  CM_FILEINFO, CM_SHOW_FILE, CM_SHOW_MODEL, CM_SHOW_ATOMS, CM_SHOW_SEQ,
  CM_SHOW_CELL, CM_SHOW_ORIENT, CM_SHOW_SEL: ShowInfo(Cm);
  CM_QUIT      : Running := False;
  CM_SELALL    : Begin I := SelectExprFn('all'); ConsoleSay(IntToStr(I) + ' atomes selectionnes.'); End;
  CM_SELNONE   : Begin SelectAllQuiet; NeedScene := True; ConsoleSay('Selection vide.'); End;
  CM_SELINV    : Begin
   N := 0;
   For I := 0 To NAtoms - 1 Do Begin
    Atoms^[I].Sel := Not Atoms^[I].Sel;
    If Atoms^[I].Sel Then Inc(N);
   End;
   NSelCount := N; NeedScene := True;
   ConsoleSay(IntToStr(N) + ' atomes selectionnes.');
  End;
  CM_SELEXPR   : Begin
   S := '';
   If InputString('Selection', 'Expression (ex: protein and helix) :', S) Then ExecLine('select ' + S);
  End;
  CM_HIDESEL   : Begin
   For I := 0 To NAtoms - 1 Do If Atoms^[I].Sel Then Atoms^[I].Hid := True;
   NeedScene := True;
  End;
  CM_SHOWALL   : Begin
   For I := 0 To NAtoms - 1 Do Atoms^[I].Hid := False;
   ShowHydro := True; ShowWater := True;
   NeedScene := True;
  End;
  CM_CENTERSEL : CenterOnTarget;
  CM_CLEARMEAS : Begin DeleteMeasures(0); NMeas := 0; PickN := 0; NeedScene := True; End;
  CM_DEL_DIST  : DeleteMeasures(2);
  CM_CLEARLABELS: Begin
   For I := 0 To NAtoms - 1 Do Atoms^[I].Lbl := 0;
   NeedScene := True;
  End;
  CM_SCH_CPK, CM_SCH_BALL, CM_SCH_STICK, CM_SCH_WIRE: ApplyScheme(Cm);
  CM_TOG_HYDRO : Begin ShowHydro := Not ShowHydro; NeedScene := True; End;
  CM_TOG_WATER : Begin ShowWater := Not ShowWater; NeedScene := True; End;
  CM_TOG_HBOND : Begin
   ShowHB := Not ShowHB;
   If ShowHB Then Begin EnsureHBonds; ConsoleSay(IntToStr(NHB) + ' liaisons hydrogene.'); End;
   NeedScene := True;
  End;
  CM_TOG_HALO  : Begin ShowHalo := Not ShowHalo; NeedScene := True; End;
  CM_TOG_PERSP : Begin Persp := Not Persp; NeedScene := True; End;
  CM_TOG_AXES  : Begin ShowAxes := Not ShowAxes; NeedScene := True; End;
  CM_TOG_BBOX  : Begin ShowBBox := Not ShowBBox; NeedScene := True; End;
  CM_TOG_CELL  : If HasCell Then Begin ShowCell := Not ShowCell; NeedScene := True; End
                 Else ConsoleErr('pas de cellule unitaire dans ce fichier');
  CM_TOG_VECTOR: Begin ShowVec := Not ShowVec; NeedScene := True; End;
  CM_TOG_CONSOLE: Begin ConsoleOn := Not ConsoleOn; NeedLayout := True; End;
  CM_TOG_SPIN  : Begin SpinOn := Not SpinOn; NeedChrome := True; End;
  CM_TOG_VIB   : SetVib(Not VibOn);
  CM_TOG_ANIM  : Begin
   If NModels > 1 Then Begin
    AnimOn := Not AnimOn;
    If AnimOn And (CurModel < 1) Then SetModel(1);
    NeedChrome := True;
   End;
  End;
  CM_VIEW_FRONT..CM_VIEW_BOTTOM: ApplyView(Cm);
  CM_VIEW_RESET: ResetView;
  CM_ZOOM_IN   : SetZoom(ZoomF * 1.25);
  CM_ZOOM_OUT  : SetZoom(ZoomF / 1.25);
  CM_ZOOM_FIT  : Begin PanX := 0; PanY := 0; SetZoom(1.0); End;
  CM_SPIN_X    : SpinAxis := 0;
  CM_SPIN_Y    : SpinAxis := 1;
  CM_SPIN_Z    : SpinAxis := 2;
  CM_SPIN_SLOW : SpinSpeed := 1.0;
  CM_SPIN_MED  : SpinSpeed := 3.0;
  CM_SPIN_FAST : SpinSpeed := 7.0;
  CM_FR_FIRST  : SetModel(1);
  CM_FR_PREV   : FrameStep(-1);
  CM_FR_NEXT   : FrameStep(1);
  CM_FR_LAST   : SetModel(NModels);
  CM_FR_ALL    : If CurModel = 0 Then SetModel(1) Else SetModel(0);
  CM_AN_ONCE   : AnimMode := 0;
  CM_AN_LOOP   : AnimMode := 1;
  CM_AN_PALIN  : AnimMode := 2;
  CM_AN_SLOWER : AnimDelay := MinI(2000, AnimDelay + AnimDelay Div 2 + 10);
  CM_AN_FASTER : AnimDelay := MaxI(20, AnimDelay - AnimDelay Div 3);
  CM_VIB_PANEL : Begin
   VibPanel := Not VibPanel;
   If VibPanel Then SetVib(True);
   NeedChrome := True; NeedScene := True;
  End;
  CM_VIB_BIGGER : Begin VibScale := MinR(4.0, VibScale * 1.4); NeedScene := True; End;
  CM_VIB_SMALLER: Begin VibScale := MaxR(0.05, VibScale / 1.4); NeedScene := True; End;
  CM_VIB_NEXT  : VibNextMode(1);
  CM_VIB_PREV  : VibNextMode(-1);
  CM_MEAS_DIST : SetPick(PICK_DIST);
  CM_MEAS_ANGLE: SetPick(PICK_ANGLE);
  CM_MEAS_TORS : SetPick(PICK_TORSION);
  CM_DEL_ANGLE : DeleteMeasures(3);
  CM_DEL_TORS  : DeleteMeasures(4);
  CM_MEAS_LIST : ShowMeasureList;
  CM_UNIT_ANG  : Begin MeasUnit := 0; NeedScene := True; End;
  CM_UNIT_NM   : Begin MeasUnit := 1; NeedScene := True; End;
  CM_UNIT_PM   : Begin MeasUnit := 2; NeedScene := True; End;
  CM_SURF_SEL  : ApplyDots(True, False);
  CM_SURF_ALL  : ApplyDots(True, True);
  CM_SURF_OFF  : Begin
   For I := 0 To NAtoms - 1 Do Atoms^[I].Dots := False;
   NeedScene := True;
  End;
  CM_CALC_FORM, CM_CALC_MASS, CM_CALC_COM, CM_CALC_GYR, CM_CALC_DIM: DoCalc(Cm);
  CM_CALC_MIN  : Minimize;
  CM_CALC_HBOND: Begin
   If Loaded Then Begin
    ComputeHBonds; ShowHB := True;
    ConsoleSay(IntToStr(NHB) + ' liaisons hydrogene detectees.');
    NeedScene := True;
   End;
  End;
  CM_HELP      : Begin HelpText(0); ShowTextDlg('Aide de Jmol'); End;
  CM_HELP_MOUSE: Begin HelpText(1); ShowTextDlg('Souris et clavier'); End;
  CM_HELP_CMDS : Begin HelpText(2); ShowTextDlg('Commandes de la console'); End;
  CM_ABOUT     : Begin HelpText(3); ShowTextDlg('A propos de Jmol'); End;
  CM_LANG_FR   : ;
 Else
  Begin
   If (Cm >= CM_STRUCT_BASE) And (Cm <= CM_STRUCT_BASE + 9) Then ApplyStruct(Cm - CM_STRUCT_BASE)
   Else If (Cm >= CM_PICK_BASE) And (Cm <= CM_PICK_BASE + 9) Then SetPick(Cm - CM_PICK_BASE)
   Else If (Cm >= CM_SEL_BASE) And (Cm < CM_SEL_BASE + NUM_SELEXPR) Then Begin
    N := SelectExprFn(SelExpr[Cm - CM_SEL_BASE]);
    If N >= 0 Then ConsoleSay(IntToStr(N) + ' atomes selectionnes.');
   End
   Else If (Cm >= CM_COLSCH_BASE) And (Cm < CM_COLSCH_BASE + 10) Then ApplyColorScheme(Cm - CM_COLSCH_BASE)
   Else If (Cm >= CM_COLATOM_BASE) And (Cm < CM_COLATOM_BASE + NUM_COLORS) Then
    ApplyColor(ColIdx[Cm - CM_COLATOM_BASE])
   Else If (Cm >= CM_COLRIBSCH_BASE) And (Cm < CM_COLRIBSCH_BASE + 5) Then ApplyRibColor(Cm - CM_COLRIBSCH_BASE, 0)
   Else If (Cm >= CM_COLRIB_BASE) And (Cm < CM_COLRIB_BASE + NUM_COLORS) Then
    ApplyRibColor(9, ColIdx[Cm - CM_COLRIB_BASE])
   Else If (Cm >= CM_COLBG_BASE) And (Cm < CM_COLBG_BASE + NUM_COLORS) Then Begin
    BgCol := ColIdx[Cm - CM_COLBG_BASE]; NeedScene := True;
   End
   Else If (Cm >= CM_ATOMPCT_BASE) And (Cm < CM_ATOMPCT_BASE + 10) Then ApplyAtomPct(AtomPct[Cm - CM_ATOMPCT_BASE])
   Else If (Cm >= CM_BONDRAD_BASE) And (Cm < CM_BONDRAD_BASE + 7) Then ApplyBondRad(BondRad[Cm - CM_BONDRAD_BASE])
   Else If (Cm >= CM_ZOOM_BASE) And (Cm < CM_ZOOM_BASE + 7) Then SetZoom(ZoomPct[Cm - CM_ZOOM_BASE] / 100.0)
   Else If (Cm >= CM_LABEL_BASE) And (Cm < CM_LABEL_BASE + NUM_LBLPRE) Then ApplyLabel(LblPre[Cm - CM_LABEL_BASE])
   Else If (Cm >= CM_MACRO_BASE) And (Cm < CM_MACRO_BASE + 11) Then RunMacro(Cm - CM_MACRO_BASE)
   Else If (Cm >= CM_ELEM_BASE) And (Cm < CM_ELEM_BASE + 104) Then Begin
    N := SelectExprFn('_' + ElemDB[Cm - CM_ELEM_BASE].Sym);
    If N >= 0 Then ConsoleSay(IntToStr(N) + ' atomes selectionnes.');
   End
   Else If (Cm >= CM_CHAIN_BASE) And (Cm < CM_CHAIN_BASE + 26) Then Begin
    I := Cm - CM_CHAIN_BASE;
    If I < NChainList Then Begin
     N := 0;
     For N := 0 To NAtoms - 1 Do
      Atoms^[N].Sel := Atoms^[N].Chain = ChainList[I];
     NSelCount := CountSelected;
     ConsoleSay(IntToStr(NSelCount) + ' atomes selectionnes (chaine ' + ChainList[I] + ').');
     NeedScene := True;
    End;
   End;
  End;
 End;
 NeedChrome := True;
End;

{ --------------------------------------------------------------------- }
{ Console : interpreteur de commandes (langage de script de type Jmol)  }
{ --------------------------------------------------------------------- }

Function NextWord(Var S: String): String;
Var I: Integer;
Begin
 S := TrimStr(S);
 I := 1;
 While (I <= Length(S)) And (S[I] <> ' ') And (S[I] <> #9) Do Inc(I);
 NextWord := Copy(S, 1, I - 1);
 Delete(S, 1, I);
 S := TrimStr(S);
End;

{ Alias francais -> commande Jmol }
Function CanonCmd(Const W: String): String;
Begin
 CanonCmd := W;
 If (W = 'charger') Then CanonCmd := 'load'
 Else If (W = 'selectionner') Or (W = 'choisir') Then CanonCmd := 'select'
 Else If W = 'afficher' Then CanonCmd := 'display'
 Else If W = 'masquer' Then CanonCmd := 'hide'
 Else If W = 'restreindre' Then CanonCmd := 'restrict'
 Else If W = 'couleur' Then CanonCmd := 'color'
 Else If (W = 'remplissage') Or (W = 'cpk') Then CanonCmd := 'spacefill'
 Else If (W = 'filaire') Or (W = 'fil') Then CanonCmd := 'wireframe'
 Else If W = 'etiquette' Then CanonCmd := 'label'
 Else If W = 'tourner' Then CanonCmd := 'rotate'
 Else If W = 'rotation' Then CanonCmd := 'spin'
 Else If W = 'translater' Then CanonCmd := 'translate'
 Else If W = 'centrer' Then CanonCmd := 'center'
 Else If W = 'reinitialiser' Then CanonCmd := 'reset'
 Else If W = 'fond' Then CanonCmd := 'background'
 Else If W = 'mesurer' Then CanonCmd := 'measure'
 Else If (W = 'modele') Or (W = 'frame') Then CanonCmd := 'model'
 Else If (W = 'animer') Or (W = 'anim') Then CanonCmd := 'animation'
 Else If W = 'vecteur' Then CanonCmd := 'vector'
 Else If (W = 'ecrire') Or (W = 'enregistrer') Or (W = 'save') Then CanonCmd := 'write'
 Else If (W = 'montrer') Or (W = 'voir') Then CanonCmd := 'show'
 Else If W = 'aide' Then CanonCmd := 'help'
 Else If (W = 'quitter') Or (W = 'exit') Then CanonCmd := 'quit'
 Else If W = 'executer' Then CanonCmd := 'script'
 Else If W = 'minimiser' Then CanonCmd := 'minimize'
 Else If W = 'calculer' Then CanonCmd := 'calculate'
 Else If W = 'boite' Then CanonCmd := 'boundbox'
 Else If W = 'cellule' Then CanonCmd := 'unitcell'
 Else If (W = 'points') Then CanonCmd := 'dots'
 Else If W = 'effacer' Then CanonCmd := 'clear'
 Else If W = 'rubans' Then CanonCmd := 'ribbons'
 Else If W = 'squelette' Then CanonCmd := 'backbone'
 Else If W = 'brins' Then CanonCmd := 'strands'
 Else If W = 'liaisonsh' Then CanonCmd := 'hbonds'
 Else If W = 'vue' Then CanonCmd := 'view'
 Else If W = 'attendre' Then CanonCmd := 'delay';
End;

Function ParseOnOff(Const A: String; Cur: Boolean; Var Ok: Boolean): Boolean;
Var U: String;
Begin
 U := LoStr(TrimStr(A));
 Ok := True;
 If (U = '') Or (U = 'on') Or (U = 'oui') Or (U = 'true') Or (U = '1') Then ParseOnOff := True
 Else If (U = 'off') Or (U = 'non') Or (U = 'false') Or (U = '0') Then ParseOnOff := False
 Else If (U = 'toggle') Or (U = 'inverser') Then ParseOnOff := Not Cur
 Else Begin Ok := False; ParseOnOff := Cur; End;
End;

Function ColorByName(Const N: String; Var C: Byte): Boolean;
Var I: Integer; U: String;
Begin
 U := LoStr(TrimStr(N));
 ColorByName := True;
 For I := 0 To NUM_COLORS - 1 Do
  If (U = ColNameEn[I]) Or (U = LoStr(ColNameFr[I])) Then Begin C := ColIdx[I]; Exit; End;
 If (U = 'grey') Or (U = 'lightgray') Or (U = 'argent') Then C := LightGray
 Else If (U = 'violet') Or (U = 'magenta') Then C := Magenta
 Else If (U = 'bleuclair') Or (U = 'lightblue') Then C := LightBlue
 Else If (U = 'marron') Or (U = 'brown') Then C := Brown
 Else If U = 'darkgreen' Then C := Green
 Else If U = 'darkblue' Then C := Blue
 Else If U = 'darkred' Then C := Red
 Else ColorByName := False;
End;

Function NumArg(Const A: String; Var V: Real): Boolean;
Begin
 NumArg := StrToReal(A, V);
End;

Procedure SetSpacefill(Const Arg: String);
Var A: String; V: Real; Ok: Boolean;
Begin
 A := LoStr(TrimStr(Arg));
 If (A = '') Or (A = 'on') Or (A = 'oui') Then Begin ApplyAtomPct(100); Exit; End;
 If (A = 'off') Or (A = 'non') Then Begin ApplyAtomPct(0); Exit; End;
 If A = 'reset' Then Begin ApplyAtomPct(23); Exit; End;
 If A[Length(A)] = '%' Then Begin
  Delete(A, Length(A), 1);
  If NumArg(A, V) Then ApplyAtomPct(RoundI(V)) Else ConsoleErr('valeur invalide : ' + Arg);
  Exit;
 End;
 Ok := NumArg(A, V);
 If Not Ok Then Begin ConsoleErr('valeur invalide : ' + Arg); Exit; End;
 If Pos('.', A) > 0 Then ApplyAtomPct(-RoundI(V * 100)) Else ApplyAtomPct(RoundI(V));
End;

Procedure SetWireframe(Const Arg: String);
Var A: String; V: Real;
Begin
 A := LoStr(TrimStr(Arg));
 If (A = '') Or (A = 'on') Or (A = 'oui') Then Begin ApplyBondRad(1); Exit; End;
 If (A = 'off') Or (A = 'non') Then Begin ApplyBondRad(0); Exit; End;
 If NumArg(A, V) Then Begin
  If Pos('.', A) > 0 Then ApplyBondRad(MaxI(1, RoundI(V * 100)))
  Else ApplyBondRad(MaxI(1, RoundI(V)));
 End
 Else ConsoleErr('valeur invalide : ' + Arg);
End;

Procedure CmdColor(Const Arg: String);
Var
 A, W : String;
 C : Byte;
 Target : Integer;
Begin
 A := Arg;
 W := LoStr(NextWord(A));
 Target := 0;   { 0 atomes, 1 ruban, 2 fond }
 If (W = 'atoms') Or (W = 'atomes') Then Begin W := LoStr(NextWord(A)); End
 Else If (W = 'cartoon') Or (W = 'cartoons') Or (W = 'ribbon') Or (W = 'ribbons') Or
         (W = 'trace') Or (W = 'backbone') Or (W = 'strands') Or (W = 'rockets') Or
         (W = 'ruban') Or (W = 'rubans') Then Begin Target := 1; W := LoStr(NextWord(A)); End
 Else If (W = 'background') Or (W = 'fond') Then Begin Target := 2; W := LoStr(NextWord(A)); End
;
 If W = '' Then Begin ConsoleErr('couleur ou schema attendu'); Exit; End;
 If Target = 2 Then Begin
  If ColorByName(W, C) Then Begin BgCol := C; NeedScene := True; End
  Else ConsoleErr('couleur inconnue : ' + W);
  Exit;
 End;
 If (W = 'cpk') Or (W = 'element') Or (W = 'none') Or (W = 'aucun') Or (W = 'name') Then Begin
  If Target = 0 Then ApplyColorScheme(0) Else ApplyRibColor(0, 0);
 End
 Else If W = 'chain' Then Begin If Target = 0 Then ApplyColorScheme(1) Else ApplyRibColor(2, 0); End
 Else If (W = 'amino') Or (W = 'residue') Or (W = 'residu') Then Begin
  If Target = 0 Then ApplyColorScheme(2) Else ApplyRibColor(0, 0);
 End
 Else If (W = 'structure') Or (W = 'secondary') Then Begin
  If Target = 0 Then ApplyColorScheme(3) Else ApplyRibColor(0, 0);
 End
 Else If (W = 'group') Or (W = 'rainbow') Or (W = 'monomer') Or (W = 'groupe') Then Begin
  If Target = 0 Then ApplyColorScheme(4) Else ApplyRibColor(1, 0);
 End
 Else If (W = 'temperature') Or (W = 'bfactor') Then Begin
  If Target = 0 Then ApplyColorScheme(5) Else ApplyRibColor(1, 0);
 End
 Else If ColorByName(W, C) Then Begin
  If Target = 0 Then ApplyColor(C) Else ApplyRibColor(9, C);
 End
 Else ConsoleErr('couleur ou schema inconnu : ' + W);
End;

Procedure CmdRotate(Const Arg: String);
Var
 A, W : String;
 V : Real;
 Ax : Integer;
Begin
 A := Arg;
 While A <> '' Do Begin
  W := LoStr(NextWord(A));
  Ax := -1;
  If W = 'x' Then Ax := 0 Else If W = 'y' Then Ax := 1 Else If W = 'z' Then Ax := 2;
  If Ax < 0 Then Begin ConsoleErr('axe attendu (x, y ou z)'); Exit; End;
  W := NextWord(A);
  If Not StrToReal(W, V) Then Begin ConsoleErr('angle attendu'); Exit; End;
  RotateView(Ax, V);
 End;
 OrthoRM;
 NeedScene := True;
End;

Procedure CmdMeasure(Const Arg: String);
Var
 A, W : String;
 Idx : Array[0..3] Of Integer;
 N, S : Integer;
Begin
 A := Arg;
 If A = '' Then Begin SetPick(PICK_DIST); Exit; End;
 W := LoStr(NextWord(A));
 If (W = 'delete') Or (W = 'clear') Or (W = 'off') Or (W = 'effacer') Then Begin
  NMeas := 0; PickN := 0; NeedScene := True; ConsoleSay('Mesures supprimees.'); Exit;
 End;
 If (W = 'list') Or (W = 'liste') Then Begin ShowMeasureList; Exit; End;
 N := 0;
 While (W <> '') And (N < 4) Do Begin
  S := StrToIntDef(W, -1);
  Idx[N] := FindSerial(S);
  If Idx[N] < 0 Then Begin ConsoleErr('atome introuvable : ' + W); Exit; End;
  Inc(N);
  W := NextWord(A);
 End;
 If N < 2 Then Begin ConsoleErr('au moins 2 atomes sont necessaires'); Exit; End;
 UpdateWork;
 If N = 2 Then Begin
  AddMeasure(2, Idx[0], Idx[1], 0, 0);
  ConsoleSay('Distance : ' + DistStr(DistW(Idx[0], Idx[1])));
 End
 Else If N = 3 Then Begin
  AddMeasure(3, Idx[0], Idx[1], Idx[2], 0);
  ConsoleSay('Angle : ' + RealToStr(AngleW(Idx[0], Idx[1], Idx[2]), 2) + #248);
 End
 Else Begin
  AddMeasure(4, Idx[0], Idx[1], Idx[2], Idx[3]);
  ConsoleSay('Diedre : ' + RealToStr(TorsionW(Idx[0], Idx[1], Idx[2], Idx[3]), 2) + #248);
 End;
End;

Procedure CmdModel(Const Arg: String);
Var A: String; V: LongInt;
Begin
 A := LoStr(TrimStr(Arg));
 If A = 'next' Then FrameStep(1)
 Else If (A = 'prev') Or (A = 'previous') Then FrameStep(-1)
 Else If A = 'first' Then SetModel(1)
 Else If A = 'last' Then SetModel(NModels)
 Else If (A = 'all') Or (A = '0') Then SetModel(0)
 Else Begin
  V := StrToIntDef(A, -1);
  If (V >= 1) And (V <= NModels) Then SetModel(V)
  Else ConsoleErr('modele invalide (1 a ' + IntToStr(NModels) + ')');
 End;
 ConsoleSay('Modele ' + IntToStr(CurModel) + '/' + IntToStr(NModels));
End;

Procedure CmdSet(Const Arg: String);
Var
 A, K, V : String;
 Ok : Boolean;
 R : Real;
Begin
 A := Arg;
 K := LoStr(NextWord(A));
 V := A;
 If (K = 'selectionhalos') Or (K = 'halos') Then ShowHalo := ParseOnOff(V, ShowHalo, Ok)
 Else If (K = 'perspectivedepth') Or (K = 'perspective') Then Persp := ParseOnOff(V, Persp, Ok)
 Else If K = 'showhydrogens' Then ShowHydro := ParseOnOff(V, ShowHydro, Ok)
 Else If K = 'showwater' Then ShowWater := ParseOnOff(V, ShowWater, Ok)
 Else If K = 'showaxes' Then ShowAxes := ParseOnOff(V, ShowAxes, Ok)
 Else If K = 'showboundbox' Then ShowBBox := ParseOnOff(V, ShowBBox, Ok)
 Else If K = 'showunitcell' Then ShowCell := ParseOnOff(V, ShowCell, Ok)
 Else If K = 'measurementunits' Then Begin
  V := LoStr(V);
  If (V = 'nm') Or (V = 'nanometers') Then MeasUnit := 1
  Else If (V = 'pm') Or (V = 'picometers') Then MeasUnit := 2
  Else MeasUnit := 0;
 End
 Else If K = 'vibrationscale' Then Begin If NumArg(V, R) Then VibScale := ClampR(R, 0.01, 10); End
 Else If K = 'vectorscale' Then Begin If NumArg(V, R) Then VecScale := ClampR(R, 0.1, 20); End
 Else If (K = 'animationfps') Or (K = 'fps') Then Begin
  If NumArg(V, R) And (R > 0) Then AnimDelay := RoundI(1000.0 / R);
 End
 Else If (K = 'spinfps') Or (K = 'spinspeed') Then Begin If NumArg(V, R) Then SpinSpeed := ClampR(R, 0.1, 30); End
 Else If K = 'picking' Then Begin
  V := LoStr(V);
  If (V = 'distance') Or (V = 'measure') Then SetPick(PICK_DIST)
  Else If V = 'angle' Then SetPick(PICK_ANGLE)
  Else If (V = 'torsion') Or (V = 'dihedral') Then SetPick(PICK_TORSION)
  Else If (V = 'select') Or (V = 'atom') Then SetPick(PICK_SELATOM)
  Else If V = 'group' Then SetPick(PICK_SELRES)
  Else If V = 'chain' Then SetPick(PICK_SELCHAIN)
  Else If V = 'center' Then SetPick(PICK_CENTER)
  Else If V = 'label' Then SetPick(PICK_LABEL)
  Else SetPick(PICK_IDENT);
 End
 Else Begin
  ConsoleSay('Parametre "' + K + '" accepte (sans effet).');
 End;
 NeedScene := True; NeedChrome := True;
End;

Procedure CmdWrite(Const Arg: String);
Var
 A, W, FN : String;
 Ok : Boolean;
 E : String;
Begin
 If Not Loaded Then Begin ConsoleErr('aucune molecule chargee'); Exit; End;
 A := Arg;
 W := LoStr(NextWord(A));
 FN := A;
 If (W <> 'pdb') And (W <> 'xyz') And (W <> 'mol') And (W <> 'bmp') And (W <> 'image') Then Begin
  FN := W + ' ' + A; FN := TrimStr(FN);
  E := LoStr(ExtractExt(FN));
  W := E;
 End;
 If FN = '' Then Begin ConsoleErr('nom de fichier attendu'); Exit; End;
 If W = 'pdb' Then Ok := SavePDBFile(FN)
 Else If W = 'xyz' Then Ok := SaveXYZFile(FN)
 Else If W = 'mol' Then Ok := SaveMOLFileOut(FN)
 Else If (W = 'bmp') Or (W = 'image') Then Begin RenderCanvas; Ok := SaveBMPFile(FN); End
 Else Begin ConsoleErr('format inconnu (pdb, xyz, mol, bmp)'); Exit; End;
 If Ok Then ConsoleSay('Fichier ecrit : ' + FN) Else ConsoleErr('impossible d''ecrire ' + FN);
End;

Procedure CmdScript(Const FN: String);
Var
 F : Text;
 Ln, Old, Dir, Full : String;
 I : Integer;
Begin
 If ScriptDepth >= 4 Then Begin ConsoleErr('scripts imbriques trop profondement'); Exit; End;
 Full := FExpand(FN);
 If Not OpenText(Full, F) Then Begin ConsoleErr('script introuvable : ' + FN); Exit; End;
 { les chemins relatifs du script sont resolus depuis son dossier }
 GetDir(0, Old);
 Dir := '';
 For I := Length(Full) DownTo 1 Do
  If (Full[I] = '\') Or (Full[I] = '/') Then Begin Dir := Copy(Full, 1, I - 1); Break; End;
 If Dir <> '' Then Begin
  {$I-} ChDir(Dir); {$I+}
  If IOResult <> 0 Then ;
 End;
 If ScriptDepth = 0 Then ScriptAbort := False;
 Inc(ScriptDepth);
 While Not Eof(F) And Running And Not ScriptAbort Do Begin
  ReadLn(F, Ln);
  ExecLine(Ln);
  If KeyPressed Then Begin
   If ReadKey = #27 Then ScriptAbort := True;
  End;
 End;
 If ScriptAbort And (ScriptDepth = 1) Then ConsoleSay('Script interrompu.');
 Close(F);
 Dec(ScriptDepth);
 {$I-} ChDir(Old); {$I+}
 If IOResult <> 0 Then ;
End;

Procedure CenterExpr(Const Arg: String);
Var
 I, N : Integer;
 Sx, Sy, Sz : Real;
Begin
 If Arg = '' Then Begin
  For I := 0 To NAtoms - 1 Do SStk[0]^[I] := Ord(InCur(I));
 End
 Else If Not EvalSel(Arg) Then Begin ConsoleErr(SErr); Exit; End;
 UpdateWork;
 N := 0; Sx := 0; Sy := 0; Sz := 0;
 For I := 0 To NAtoms - 1 Do
  If (SStk[0]^[I] <> 0) And InCur(I) Then Begin
   Sx := Sx + WX^[I]; Sy := Sy + WY^[I]; Sz := Sz + WZ^[I]; Inc(N);
  End;
 If N = 0 Then Begin ConsoleErr('aucun atome'); Exit; End;
 CenX := Sx / N; CenY := Sy / N; CenZ := Sz / N;
 PanX := 0; PanY := 0;
 NeedScene := True;
End;

Procedure ExecLine(Const Line: String);
Var
 L, Arg, W, W2, Rest : String;
 P, I, N : Integer;
 Ok : Boolean;
 V : Real;
 C : Byte;
Begin
 L := TrimStr(Line);
 If L = '' Then Exit;
 If (L[1] = '#') Or ((Length(L) > 1) And (L[1] = '/') And (L[2] = '/')) Then Exit;
 P := Pos(';', L);
 If (P > 0) And (Not StartsWithI(L, 'echo')) Then Begin
  ExecLine(Copy(L, 1, P - 1));
  ExecLine(Copy(L, P + 1, 250));
  Exit;
 End;
 Arg := L;
 W := LoStr(NextWord(Arg));
 W := CanonCmd(W);
 If Not Loaded Then Begin
  If (W <> 'load') And (W <> 'help') And (W <> 'quit') And (W <> 'clear') And
     (W <> 'script') And (W <> 'console') And (W <> 'echo') And (W <> 'background') And
     (W <> 'set') Then Begin
   ConsoleErr('aucune molecule chargee (utilisez : load fichier)');
   Exit;
  End;
 End;
 If W = 'load' Then Begin
  If Arg = '' Then DoOpen
  Else LoadFileUI(Arg, True);
 End
 Else If W = 'select' Then Begin
  N := SelectExprFn(Arg);
  If N < 0 Then ConsoleErr(SErr) Else ConsoleSay(IntToStr(N) + ' atomes selectionnes.');
 End
 Else If (W = 'display') Or (W = 'restrict') Then Begin
  If Arg = '' Then Arg := 'all';
  If EvalSel(Arg) Then Begin
   N := 0;
   For I := 0 To NAtoms - 1 Do Begin
    Atoms^[I].Hid := SStk[0]^[I] = 0;
    If Not Atoms^[I].Hid Then Inc(N);
    If W = 'restrict' Then Atoms^[I].Sel := SStk[0]^[I] <> 0;
   End;
   NeedScene := True;
   ConsoleSay(IntToStr(N) + ' atomes affiches.');
  End
  Else ConsoleErr(SErr);
 End
 Else If W = 'hide' Then Begin
  If Arg = '' Then Arg := 'selected';
  If EvalSel(Arg) Then Begin
   For I := 0 To NAtoms - 1 Do If SStk[0]^[I] <> 0 Then Atoms^[I].Hid := True;
   NeedScene := True;
  End
  Else ConsoleErr(SErr);
 End
 Else If W = 'color' Then CmdColor(Arg)
 Else If W = 'spacefill' Then SetSpacefill(Arg)
 Else If W = 'wireframe' Then SetWireframe(Arg)
 Else If (W = 'cartoon') Or (W = 'ribbons') Or (W = 'trace') Or (W = 'backbone') Or
         (W = 'strands') Or (W = 'rockets') Then Begin
  W2 := LoStr(Arg);
  I := RIB_CARTOON;
  If W = 'ribbons' Then I := RIB_RIBBON
  Else If W = 'trace' Then I := RIB_TRACE
  Else If W = 'backbone' Then I := RIB_BACKBONE
  Else If W = 'strands' Then I := RIB_STRAND;
  If (W2 = 'off') Or (W2 = 'non') Then I := RIB_NONE;
  ApplyStruct(I);
 End
 Else If W = 'dots' Then Begin
  Ok := True;
  ApplyDots(ParseOnOff(Arg, True, Ok), False);
 End
 Else If W = 'label' Then Begin
  If (LoStr(Arg) = 'off') Or (LoStr(Arg) = 'non') Then ApplyLabel('')
  Else If Arg = '' Then ApplyLabel('%a')
  Else Begin
   If (Arg[1] = '"') Or (Arg[1] = '''') Then Arg := Copy(Arg, 2, Length(Arg) - 2);
   ApplyLabel(Arg);
  End;
 End
 Else If W = 'hbonds' Then Begin
  ShowHB := ParseOnOff(Arg, ShowHB, Ok);
  If ShowHB Then Begin EnsureHBonds; ConsoleSay(IntToStr(NHB) + ' liaisons hydrogene.'); End;
  NeedScene := True;
 End
 Else If W = 'ssbonds' Then Begin ConsoleSay('Les ponts disulfure sont affiches comme liaisons normales.'); End
 Else If W = 'axes' Then Begin ShowAxes := ParseOnOff(Arg, ShowAxes, Ok); NeedScene := True; End
 Else If W = 'boundbox' Then Begin ShowBBox := ParseOnOff(Arg, ShowBBox, Ok); NeedScene := True; End
 Else If W = 'unitcell' Then Begin
  If HasCell Then Begin ShowCell := ParseOnOff(Arg, ShowCell, Ok); NeedScene := True; End
  Else ConsoleErr('pas de cellule unitaire dans ce fichier');
 End
 Else If W = 'rotate' Then CmdRotate(Arg)
 Else If W = 'spin' Then Begin
  W2 := LoStr(NextWord(Arg));
  If (W2 = 'x') Or (W2 = 'y') Or (W2 = 'z') Then Begin
   If W2 = 'x' Then SpinAxis := 0 Else If W2 = 'y' Then SpinAxis := 1 Else SpinAxis := 2;
   SpinOn := True;
   If StrToReal(NextWord(Arg), V) Then SpinSpeed := ClampR(V / 10.0, 0.1, 30);
  End
  Else SpinOn := ParseOnOff(W2, SpinOn, Ok);
  NeedChrome := True;
 End
 Else If W = 'zoom' Then Begin
  W2 := LoStr(Arg);
  If W2 = 'in' Then SetZoom(ZoomF * 1.25)
  Else If W2 = 'out' Then SetZoom(ZoomF / 1.25)
  Else If (W2 = 'reset') Or (W2 = 'fit') Then Begin PanX := 0; PanY := 0; SetZoom(1.0); End
  Else If StrToReal(W2, V) Then SetZoom(V / 100.0)
  Else ConsoleErr('valeur de zoom invalide');
 End
 Else If W = 'translate' Then Begin
  W2 := LoStr(NextWord(Arg));
  If StrToReal(NextWord(Arg), V) Then Begin
   If W2 = 'x' Then PanX := PanX + RoundI(V * CvW / 100.0)
   Else If W2 = 'y' Then PanY := PanY + RoundI(V * CvH / 100.0);
   NeedScene := True;
  End
  Else ConsoleErr('syntaxe : translate x|y pourcentage');
 End
 Else If W = 'center' Then CenterExpr(Arg)
 Else If W = 'reset' Then ResetView
 Else If W = 'view' Then Begin
  W2 := LoStr(Arg);
  If (W2 = 'front') Or (W2 = 'avant') Then ApplyView(CM_VIEW_FRONT)
  Else If (W2 = 'back') Or (W2 = 'arriere') Then ApplyView(CM_VIEW_BACK)
  Else If (W2 = 'left') Or (W2 = 'gauche') Then ApplyView(CM_VIEW_LEFT)
  Else If (W2 = 'right') Or (W2 = 'droite') Then ApplyView(CM_VIEW_RIGHT)
  Else If (W2 = 'top') Or (W2 = 'dessus') Then ApplyView(CM_VIEW_TOP)
  Else If (W2 = 'bottom') Or (W2 = 'dessous') Then ApplyView(CM_VIEW_BOTTOM)
  Else ConsoleErr('orientation inconnue : ' + Arg);
 End
 Else If W = 'background' Then Begin
  If ColorByName(Arg, C) Then Begin BgCol := C; NeedScene := True; End
  Else ConsoleErr('couleur inconnue : ' + Arg);
 End
 Else If W = 'measure' Then CmdMeasure(Arg)
 Else If W = 'model' Then CmdModel(Arg)
 Else If W = 'animation' Then Begin
  W2 := LoStr(NextWord(Arg));
  If W2 = 'loop' Then Begin AnimMode := 1; AnimOn := True; End
  Else If W2 = 'once' Then Begin AnimMode := 0; AnimOn := True; End
  Else If W2 = 'palindrome' Then Begin AnimMode := 2; AnimOn := True; End
  Else If W2 = 'fps' Then Begin
   If StrToReal(NextWord(Arg), V) And (V > 0) Then AnimDelay := RoundI(1000.0 / V);
  End
  Else If NModels > 1 Then AnimOn := ParseOnOff(W2, AnimOn, Ok)
  Else ConsoleErr('un seul modele dans ce fichier');
  NeedChrome := True;
 End
 Else If W = 'vibration' Then Begin
  W2 := LoStr(NextWord(Arg));
  If W2 = 'scale' Then Begin If StrToReal(NextWord(Arg), V) Then VibScale := ClampR(V, 0.01, 10); End
  Else If (W2 = 'panel') Or (W2 = 'panneau') Then DoCmd(CM_VIB_PANEL)
  Else SetVib(ParseOnOff(W2, VibOn, Ok));
 End
 Else If W = 'vector' Then Begin
  W2 := LoStr(NextWord(Arg));
  If W2 = 'scale' Then Begin If StrToReal(NextWord(Arg), V) Then VecScale := ClampR(V, 0.1, 20); End
  Else Begin ShowVec := ParseOnOff(W2, ShowVec, Ok); End;
  NeedScene := True;
 End
 Else If W = 'set' Then CmdSet(Arg)
 Else If W = 'echo' Then ConsoleSay(Arg)
 Else If W = 'write' Then CmdWrite(Arg)
 Else If W = 'show' Then Begin
  W2 := LoStr(Arg);
  If (W2 = 'file') Or (W2 = 'fichier') Or (W2 = 'info') Then ShowInfo(CM_SHOW_FILE)
  Else If W2 = 'model' Then ShowInfo(CM_SHOW_MODEL)
  Else If (W2 = 'atoms') Or (W2 = 'atomes') Then ShowInfo(CM_SHOW_ATOMS)
  Else If (W2 = 'sequence') Or (W2 = 'chain') Then ShowInfo(CM_SHOW_SEQ)
  Else If (W2 = 'cell') Or (W2 = 'unitcell') Then ShowInfo(CM_SHOW_CELL)
  Else If W2 = 'orientation' Then ShowInfo(CM_SHOW_ORIENT)
  Else If (W2 = 'selection') Or (W2 = 'selected') Then ShowInfo(CM_SHOW_SEL)
  Else If (W2 = 'measures') Or (W2 = 'mesures') Then ShowMeasureList
  Else ConsoleErr('show : file, model, atoms, sequence, cell, orientation, selection, measures');
 End
 Else If W = 'script' Then CmdScript(Arg)
 Else If W = 'delay' Then Begin
  If StrToReal(Arg, V) Then Begin
   If V < 20 Then V := V * 1000;
   WaitMs(Round(MinR(V, 60000)));
  End;
 End
 Else If W = 'minimize' Then Minimize
 Else If W = 'calculate' Then Begin
  W2 := LoStr(Arg);
  If W2 = 'hbonds' Then DoCmd(CM_CALC_HBOND)
  Else If W2 = 'formula' Then DoCalc(CM_CALC_FORM)
  Else If W2 = 'mass' Then DoCalc(CM_CALC_MASS)
  Else If (W2 = 'com') Or (W2 = 'centermass') Then DoCalc(CM_CALC_COM)
  Else If W2 = 'gyration' Then DoCalc(CM_CALC_GYR)
  Else ConsoleErr('calculate : hbonds, formula, mass, com, gyration');
 End
 Else If W = 'console' Then Begin
  ConsoleOn := ParseOnOff(Arg, ConsoleOn, Ok);
  NeedLayout := True;
 End
 Else If W = 'help' Then Begin HelpText(2); ShowTextDlg('Commandes de la console'); End
 Else If W = 'quit' Then Running := False
 Else If W = 'clear' Then Begin ConCount := 0; NeedChrome := True; End
 Else If W = 'hydrogens' Then Begin ShowHydro := ParseOnOff(Arg, ShowHydro, Ok); NeedScene := True; End
 Else If W = 'halos' Then Begin ShowHalo := ParseOnOff(Arg, ShowHalo, Ok); NeedScene := True; End
 Else If W = 'perspective' Then Begin Persp := ParseOnOff(Arg, Persp, Ok); NeedScene := True; End
 Else If W = 'ballandstick' Then ApplyScheme(CM_SCH_BALL)
 Else If W = 'sticks' Then ApplyScheme(CM_SCH_STICK)
 Else Begin
  Rest := W;
  ConsoleErr('commande inconnue : ' + Rest + ' (tapez "aide")');
 End;
End;

{ --------------------------------------------------------------------- }
{ Interface : titre, menus, barre d'outils, barre d'etat, console       }
{ --------------------------------------------------------------------- }

Const
 NUM_TB    = 17;
 TB_SEP_AFTER = 5;

 { Icones }
 IC_OPEN = 0;  IC_CAMERA = 1; IC_PALETTE = 2; IC_ZOOM = 3;  IC_INFO = 4;
 IC_SPIN = 5;  IC_ARROW = 6;  IC_RULER = 7;   IC_SELECT = 8; IC_CENTER = 9;
 IC_HOME = 10; IC_FIRST = 11; IC_PREV = 12;   IC_NEXT = 13;  IC_LAST = 14;
 IC_PLAY = 15; IC_VIB = 16;

 CM_TB_COLOR = 900;

Type
 TToolBtn = Record
  Cmd  : Integer;
  Icon : Integer;
  Tip  : String[40];
  X    : Integer;
  Sep  : Boolean;
 End;

Var
 TBtn      : Array[0..NUM_TB-1] Of TToolBtn;
 TBY       : Integer;
 TBDown    : Integer;
 VpX, VpY  : Integer;          { panneau de vibration }
 VibDrag   : Integer;
 LastSpin, LastAnim, LastVib : LongInt;

Const
 VP_W = 170;
 VP_H = 258;

Procedure InitToolbar;
Var I, X: Integer;
 Procedure T(K, C, Ic: Integer; Const Tip: String; Sp: Boolean);
 Begin
  TBtn[K].Cmd := C; TBtn[K].Icon := Ic; TBtn[K].Tip := Tip; TBtn[K].Sep := Sp;
 End;
Begin
 T(0, CM_OPEN, IC_OPEN, 'Ouvrir un fichier', False);
 T(1, CM_SAVEBMP, IC_CAMERA, 'Exporter l''image (BMP)', False);
 T(2, CM_TB_COLOR, IC_PALETTE, 'Couleurs', False);
 T(3, CM_ZOOM_IN, IC_ZOOM, 'Zoom avant', False);
 T(4, CM_FILEINFO, IC_INFO, 'Informations sur le fichier', False);
 T(5, CM_TOG_SPIN, IC_SPIN, 'Rotation automatique', True);
 T(6, CM_PICK_BASE + PICK_IDENT, IC_ARROW, 'Identifier les atomes', False);
 T(7, CM_MEAS_DIST, IC_RULER, 'Mesurer une distance', False);
 T(8, CM_PICK_BASE + PICK_SELATOM, IC_SELECT, 'Selectionner des atomes', False);
 T(9, CM_PICK_BASE + PICK_CENTER, IC_CENTER, 'Centrer sur un atome', False);
 T(10, CM_VIEW_RESET, IC_HOME, 'Reinitialiser la vue', False);
 T(11, CM_FR_FIRST, IC_FIRST, 'Premier modele', True);
 T(12, CM_FR_PREV, IC_PREV, 'Modele precedent', False);
 T(13, CM_FR_NEXT, IC_NEXT, 'Modele suivant', False);
 T(14, CM_FR_LAST, IC_LAST, 'Dernier modele', False);
 T(15, CM_TOG_ANIM, IC_PLAY, 'Animation des modeles', False);
 T(16, CM_TOG_VIB, IC_VIB, 'Vibration', False);
 X := 6;
 For I := 0 To NUM_TB - 1 Do Begin
  If TBtn[I].Sep Then X := X + 8;
  TBtn[I].X := X;
  X := X + 26;
 End;
End;

Procedure Layout;
Begin
 CvX1 := 3; CvX2 := SW - 4;
 CvY1 := TITLE_H + MENU_H + TOOL_H + 1;
 If ConsoleOn Then CvY2 := SH - STAT_H - CONS_H - 2
 Else CvY2 := SH - STAT_H - 2;
 SetupView;
 VpX := CvX2 - VP_W - 6;
 VpY := CvY1 + 6;
End;

Function TBActive(Cmd: Integer): Boolean;
Begin
 Case Cmd Of
  CM_TOG_SPIN : TBActive := SpinOn;
  CM_TOG_ANIM : TBActive := AnimOn;
  CM_TOG_VIB  : TBActive := VibOn;
 Else
  If (Cmd >= CM_PICK_BASE) And (Cmd <= CM_PICK_BASE + 9) Then TBActive := PickMode = Cmd - CM_PICK_BASE
  Else If Cmd = CM_MEAS_DIST Then TBActive := PickMode = PICK_DIST
  Else TBActive := False;
 End;
End;

Function TBEnabled(Cmd: Integer): Boolean;
Begin
 TBEnabled := True;
 Case Cmd Of
  CM_FR_FIRST, CM_FR_PREV, CM_FR_NEXT, CM_FR_LAST, CM_TOG_ANIM: TBEnabled := NModels > 1;
  CM_TOG_VIB: TBEnabled := HasVibData;
  CM_SAVEBMP, CM_FILEINFO, CM_ZOOM_IN, CM_VIEW_RESET, CM_TOG_SPIN, CM_MEAS_DIST: TBEnabled := Loaded;
 End;
 If (Cmd >= CM_PICK_BASE) And (Cmd <= CM_PICK_BASE + 9) Then TBEnabled := Loaded;
End;

{ Dessine l'icone d'un bouton dans un carre 20x20 d'origine (X, Y) }
Procedure DrawToolIcon(Ic, X, Y: Integer; En: Boolean);
Var
 P : Array[0..6] Of PointType;
 C1, C2, C3 : Byte;
Begin
 SetLineStyle(SolidLn, 0, NormWidth);
 If Not En Then Begin C1 := DarkGray; C2 := DarkGray; C3 := LightGray; End
 Else Begin C1 := Black; C2 := Blue; C3 := Yellow; End;
 Case Ic Of
  IC_OPEN: Begin
   If En Then C3 := Yellow Else C3 := White;
   SetColor(C1); SetFillStyle(SolidFill, C3);
   Bar(X + 2, Y + 5, X + 17, Y + 16); Rectangle(X + 2, Y + 5, X + 17, Y + 16);
   Bar(X + 2, Y + 3, X + 8, Y + 5); Rectangle(X + 2, Y + 3, X + 8, Y + 5);
   Line(X + 5, Y + 8, X + 19, Y + 8);
  End;
  IC_CAMERA: Begin
   SetColor(C1); SetFillStyle(SolidFill, DarkGray);
   Bar(X + 2, Y + 6, X + 17, Y + 16); Rectangle(X + 2, Y + 6, X + 17, Y + 16);
   Bar(X + 6, Y + 4, X + 11, Y + 6);
   SetFillStyle(SolidFill, LightCyan); SetColor(C1);
   FillEllipse(X + 9, Y + 11, 3, 3);
  End;
  IC_PALETTE: Begin
   SetColor(C1); SetFillStyle(SolidFill, White);
   FillEllipse(X + 10, Y + 10, 8, 6);
   If En Then Begin
    SetFillStyle(SolidFill, LightRed); FillEllipse(X + 6, Y + 8, 1, 1);
    SetFillStyle(SolidFill, LightBlue); FillEllipse(X + 10, Y + 6, 1, 1);
    SetFillStyle(SolidFill, LightGreen); FillEllipse(X + 14, Y + 8, 1, 1);
    SetFillStyle(SolidFill, Yellow); FillEllipse(X + 13, Y + 12, 1, 1);
   End;
  End;
  IC_ZOOM: Begin
   SetColor(C1); SetFillStyle(SolidFill, LightCyan);
   FillEllipse(X + 8, Y + 8, 5, 5);
   SetLineStyle(SolidLn, 0, ThickWidth);
   Line(X + 12, Y + 12, X + 17, Y + 17);
   SetLineStyle(SolidLn, 0, NormWidth);
   SetColor(C2); Line(X + 6, Y + 8, X + 10, Y + 8); Line(X + 8, Y + 6, X + 8, Y + 10);
  End;
  IC_INFO: Begin
   SetColor(C1); SetFillStyle(SolidFill, LightBlue);
   FillEllipse(X + 10, Y + 10, 7, 7);
   SetColor(White); Line(X + 10, Y + 9, X + 10, Y + 14); PutPixel(X + 10, Y + 6, White);
   PutPixel(X + 9, Y + 9, White);
  End;
  IC_SPIN: Begin
   SetColor(Green);
   SetLineStyle(SolidLn, 0, ThickWidth);
   Arc(X + 10, Y + 10, 20, 320, 6);
   SetLineStyle(SolidLn, 0, NormWidth);
   P[0].X := X + 15; P[0].Y := Y + 2;
   P[1].X := X + 17; P[1].Y := Y + 9;
   P[2].X := X + 10; P[2].Y := Y + 7;
   SetFillStyle(SolidFill, Green);
   If En Then FillPoly(3, P);
  End;
  IC_ARROW: Begin
   P[0].X := X + 5; P[0].Y := Y + 3;
   P[1].X := X + 5; P[1].Y := Y + 16;
   P[2].X := X + 8; P[2].Y := Y + 13;
   P[3].X := X + 11; P[3].Y := Y + 18;
   P[4].X := X + 13; P[4].Y := Y + 17;
   P[5].X := X + 10; P[5].Y := Y + 12;
   P[6].X := X + 14; P[6].Y := Y + 12;
   SetColor(C1); SetFillStyle(SolidFill, White);
   FillPoly(7, P);
  End;
  IC_RULER: Begin
   P[0].X := X + 2; P[0].Y := Y + 14;
   P[1].X := X + 14; P[1].Y := Y + 2;
   P[2].X := X + 18; P[2].Y := Y + 6;
   P[3].X := X + 6; P[3].Y := Y + 18;
   SetColor(C1); SetFillStyle(SolidFill, C3);
   FillPoly(4, P);
   Line(X + 5, Y + 11, X + 7, Y + 13); Line(X + 8, Y + 8, X + 10, Y + 10);
   Line(X + 11, Y + 5, X + 13, Y + 7);
  End;
  IC_SELECT: Begin
   SetColor(C1); SetFillStyle(SolidFill, White);
   FillEllipse(X + 10, Y + 10, 6, 6);
   SetColor(C2);
   Line(X + 2, Y + 10, X + 18, Y + 10); Line(X + 10, Y + 2, X + 10, Y + 18);
  End;
  IC_CENTER: Begin
   SetColor(C2);
   SetLineStyle(SolidLn, 0, ThickWidth);
   Line(X + 4, Y + 4, X + 16, Y + 16); Line(X + 16, Y + 4, X + 4, Y + 16);
   SetLineStyle(SolidLn, 0, NormWidth);
  End;
  IC_HOME: Begin
   P[0].X := X + 10; P[0].Y := Y + 3;
   P[1].X := X + 18; P[1].Y := Y + 10;
   P[2].X := X + 2; P[2].Y := Y + 10;
   SetColor(C1);
   If En Then SetFillStyle(SolidFill, LightRed) Else SetFillStyle(SolidFill, White);
   FillPoly(3, P);
   If En Then SetFillStyle(SolidFill, Brown) Else SetFillStyle(SolidFill, White);
   Bar(X + 5, Y + 10, X + 15, Y + 17); Rectangle(X + 5, Y + 10, X + 15, Y + 17);
  End;
  IC_FIRST, IC_PREV, IC_NEXT, IC_LAST: Begin
   SetColor(C2); SetFillStyle(SolidFill, C2);
   If (Ic = IC_FIRST) Or (Ic = IC_PREV) Then Begin
    P[0].X := X + 14; P[0].Y := Y + 4; P[1].X := X + 14; P[1].Y := Y + 16; P[2].X := X + 6; P[2].Y := Y + 10;
    FillPoly(3, P);
    If Ic = IC_FIRST Then Bar(X + 4, Y + 4, X + 5, Y + 16);
   End
   Else Begin
    P[0].X := X + 6; P[0].Y := Y + 4; P[1].X := X + 6; P[1].Y := Y + 16; P[2].X := X + 14; P[2].Y := Y + 10;
    FillPoly(3, P);
    If Ic = IC_LAST Then Bar(X + 15, Y + 4, X + 16, Y + 16);
   End;
  End;
  IC_PLAY: Begin
   If AnimOn Then Begin
    SetColor(C2); SetFillStyle(SolidFill, C2);
    Bar(X + 5, Y + 4, X + 8, Y + 16); Bar(X + 12, Y + 4, X + 15, Y + 16);
   End
   Else Begin
    SetColor(Green); SetFillStyle(SolidFill, Green);
    If Not En Then Begin SetColor(DarkGray); SetFillStyle(SolidFill, DarkGray); End;
    P[0].X := X + 6; P[0].Y := Y + 3; P[1].X := X + 6; P[1].Y := Y + 17; P[2].X := X + 16; P[2].Y := Y + 10;
    FillPoly(3, P);
   End;
  End;
  IC_VIB: Begin
   SetColor(C2);
   SetLineStyle(SolidLn, 0, ThickWidth);
   Line(X + 2, Y + 10, X + 5, Y + 4); Line(X + 5, Y + 4, X + 8, Y + 16);
   Line(X + 8, Y + 16, X + 11, Y + 4); Line(X + 11, Y + 4, X + 14, Y + 16);
   Line(X + 14, Y + 16, X + 17, Y + 10);
   SetLineStyle(SolidLn, 0, NormWidth);
  End;
 End;
End;

Procedure DrawToolbar;
Var
 I, X, Y : Integer;
 Down, En : Boolean;
Begin
 Y := TITLE_H + MENU_H;
 PanelBar(0, Y, SW - 1, Y + TOOL_H - 1, LightGray);
 SetColor(White); Line(0, Y, SW - 1, Y);
 SetColor(DarkGray); Line(0, Y + TOOL_H - 1, SW - 1, Y + TOOL_H - 1);
 TBY := Y + 2;
 For I := 0 To NUM_TB - 1 Do Begin
  X := TBtn[I].X;
  En := TBEnabled(TBtn[I].Cmd);
  Down := TBActive(TBtn[I].Cmd) Or (TBDown = I);
  If TBtn[I].Sep Then Begin
   SetColor(DarkGray); Line(X - 5, TBY + 2, X - 5, TBY + 21);
   SetColor(White); Line(X - 4, TBY + 2, X - 4, TBY + 21);
  End;
  If Down Then Begin
   PanelBar(X, TBY, X + 23, TBY + 23, LightGray);
   Frame3D(X, TBY, X + 23, TBY + 23, False);
  End
  Else If (ToolHover = I) And En Then Begin
   PanelBar(X, TBY, X + 23, TBY + 23, LightGray);
   Frame3D(X, TBY, X + 23, TBY + 23, True);
  End
  Else PanelBar(X, TBY, X + 23, TBY + 23, LightGray);
  DrawToolIcon(TBtn[I].Icon, X + 2 + Ord(Down), TBY + 2 + Ord(Down), En);
 End;
End;

Procedure DrawTitleBar;
Var X: Integer;
Begin
 PanelBar(0, 0, SW - 1, TITLE_H - 1, Blue);
 { icone }
 SetColor(White); SetFillStyle(SolidFill, LightRed);
 FillEllipse(8, 9, 4, 4);
 SetFillStyle(SolidFill, LightCyan); FillEllipse(13, 6, 3, 3);
 Txt(24, 5, WindowTitle, White);
 { boutons de fenetre }
 X := SW - 19;
 PanelBar(X, 2, X + 16, TITLE_H - 3, LightRed);
 Frame3D(X, 2, X + 16, TITLE_H - 3, True);
 SetColor(White); SetLineStyle(SolidLn, 0, NormWidth);
 Line(X + 5, 5, X + 11, TITLE_H - 6); Line(X + 11, 5, X + 5, TITLE_H - 6);
 Line(X + 6, 5, X + 12, TITLE_H - 6); Line(X + 12, 5, X + 6, TITLE_H - 6);
 X := SW - 37;
 Button3D(X, 2, X + 16, TITLE_H - 3, False);
 SetColor(Black); Rectangle(X + 4, 5, X + 12, TITLE_H - 6); Line(X + 4, 6, X + 12, 6);
 X := SW - 54;
 Button3D(X, 2, X + 16, TITLE_H - 3, False);
 SetColor(Black); Line(X + 4, TITLE_H - 6, X + 11, TITLE_H - 6); Line(X + 4, TITLE_H - 5, X + 11, TITLE_H - 5);
End;

Procedure DrawMenuBar(Active: Integer);
Var I: Integer;
Begin
 BarLayout;
 PanelBar(0, TITLE_H, SW - 1, TITLE_H + MENU_H - 1, LightGray);
 For I := 0 To NUM_BAR - 1 Do Begin
  If I = Active Then Begin
   PanelBar(BarX[I], TITLE_H + 1, BarX[I] + BarW[I] - 1, TITLE_H + MENU_H - 2, Blue);
   Txt(BarX[I] + 7, TITLE_H + 5, PM^[PM_B_FILE + I].Title, White);
  End
  Else Txt(BarX[I] + 7, TITLE_H + 5, PM^[PM_B_FILE + I].Title, Black);
 End;
End;

Procedure DrawStatusBar;
Var
 Y, W1, W2, X2, X3 : Integer;
 S : String;
Begin
 Y := SH - STAT_H;
 PanelBar(0, Y, SW - 1, SH - 1, LightGray);
 W1 := (SW * 60) Div 100;
 W2 := (SW * 14) Div 100;
 X2 := W1 + 2; X3 := X2 + W2 + 2;
 Frame3D(2, Y + 2, W1, SH - 3, False);
 Frame3D(X2, Y + 2, X2 + W2, SH - 3, False);
 Frame3D(X3, Y + 2, SW - 3, SH - 3, False);
 If (HoverAtom >= 0) And (HoverAtom < NAtoms) Then S := AtomInfo(HoverAtom)
 Else If ToolHover >= 0 Then S := TBtn[ToolHover].Tip
 Else S := StatMsg;
 SetViewPort(4, Y + 3, W1 - 2, SH - 4, True);
 Txt(0, 4, Copy(S, 1, (W1 - 6) Div 8), Black);
 SetViewPort(0, 0, SW - 1, SH - 1, True);
 S := IntToStr(CvW) + ' x ' + IntToStr(CvH);
 Txt(X2 + (W2 - TextWidth(S)) Div 2, Y + 6, S, Black);
 If NModels > 1 Then
  S := 'Modele ' + IntToStr(MaxI(CurModel, 0)) + '/' + IntToStr(NModels) + '  ' +
       IntToStr(RenderMs) + ' ms'
 Else S := IntToStr(NAtoms) + ' atomes  ' + IntToStr(RenderMs) + ' ms';
 SetViewPort(X3 + 2, Y + 3, SW - 4, SH - 4, True);
 Txt(0, 4, Copy(S, 1, (SW - X3 - 8) Div 8), Black);
 SetViewPort(0, 0, SW - 1, SH - 1, True);
End;

Procedure DrawConsolePanel;
Var
 Y1, Y2, I, N, Top, Y : Integer;
 S : String;
 C : Byte;
Begin
 If Not ConsoleOn Then Exit;
 Y1 := CvY2 + 2; Y2 := SH - STAT_H - 1;
 PanelBar(0, Y1, SW - 1, Y2, LightGray);
 PanelBar(2, Y1 + 1, SW - 3, Y1 + 11, Blue);
 Txt(6, Y1 + 2, 'Console Jmol', White);
 PanelBar(2, Y1 + 12, SW - 3, Y2 - 14, White);
 Frame3D(2, Y1 + 12, SW - 3, Y2 - 14, False);
 N := (Y2 - 14 - Y1 - 14) Div 10;
 Top := ConCount - N;
 If Top < 0 Then Top := 0;
 SetViewPort(4, Y1 + 14, SW - 5, Y2 - 15, True);
 Y := 0;
 For I := Top To ConCount - 1 Do Begin
  S := Copy(ConLines[I], 2, 79);
  Case ConLines[I][1] Of
   '!': C := Red;
   '>': C := Blue;
  Else
   C := Black;
  End;
  Txt(2, Y, S, C);
  Inc(Y, 10);
 End;
 SetViewPort(0, 0, SW - 1, SH - 1, True);
 PanelBar(2, Y2 - 13, SW - 3, Y2 - 1, White);
 Frame3D(2, Y2 - 13, SW - 3, Y2 - 1, False);
 S := '$ ' + ConInput + '_';
 SetViewPort(4, Y2 - 12, SW - 5, Y2 - 2, True);
 If Length(S) > 76 Then S := Copy(S, Length(S) - 75, 76);
 Txt(2, 2, S, Black);
 SetViewPort(0, 0, SW - 1, SH - 1, True);
End;

Procedure DrawChrome;
Begin
 HideMouse;
 SetLineStyle(SolidLn, 0, NormWidth);
 SetTextStyle(DefaultFont, HorizDir, 1);
 DrawTitleBar;
 DrawMenuBar(-1);
 DrawToolbar;
 { bordures laterales }
 PanelBar(0, CvY1 - 1, 2, SH - STAT_H - 1, LightGray);
 PanelBar(SW - 3, CvY1 - 1, SW - 1, SH - STAT_H - 1, LightGray);
 PanelBar(0, CvY2 + 1, SW - 1, CvY2 + 1, LightGray);
 If ConsoleOn Then DrawConsolePanel
 Else PanelBar(0, CvY2 + 1, SW - 1, SH - STAT_H - 1, LightGray);
 Frame3D(2, CvY1 - 1, SW - 3, CvY2 + 1, False);
 DrawStatusBar;
 ShowMouse;
 NeedChrome := False; NeedStatus := False;
End;

{ --------------------------------------------------------------------- }
{ Panneau de vibration                                                  }
{ --------------------------------------------------------------------- }

Procedure DrawVibPanel;
Var
 X, Y, I, Tx : Integer;
 S : String;
Begin
 If Not (VibPanel And HasVibData) Then Exit;
 X := VpX; Y := VpY;
 HideMouse;
 SetTextStyle(DefaultFont, HorizDir, 1);
 PanelBar(X, Y, X + VP_W, Y + VP_H, LightGray);
 Frame3D(X, Y, X + VP_W, Y + VP_H, True);
 Frame3D(X + 1, Y + 1, X + VP_W - 1, Y + VP_H - 1, True);
 PanelBar(X + 3, Y + 3, X + VP_W - 3, Y + 17, Blue);
 Txt(X + 8, Y + 7, 'Vibration', White);
 PanelBar(X + VP_W - 19, Y + 4, X + VP_W - 4, Y + 16, LightRed);
 Frame3D(X + VP_W - 19, Y + 4, X + VP_W - 4, Y + 16, True);
 SetColor(White); Line(X + VP_W - 15, Y + 7, X + VP_W - 8, Y + 13); Line(X + VP_W - 8, Y + 7, X + VP_W - 15, Y + 13);
 Txt(X + 8, Y + 22, 'Image d''evaluation', Black);
 PanelBar(X + 8, Y + 33, X + VP_W - 8, Y + 47, LightGray);
 Frame3D(X + 8, Y + 33, X + VP_W - 8, Y + 47, False);
 Txt(X + 12, Y + 37, 'Image ' + IntToStr(MaxI(CurModel, 1)) + '/' + IntToStr(NModels), DarkGray);
 Txt(X + 8, Y + 54, 'Frequence du mode', Black);
 PanelBar(X + 8, Y + 65, X + VP_W - 8, Y + 81, White);
 Frame3D(X + 8, Y + 65, X + VP_W - 8, Y + 81, False);
 If (CurModel >= 1) And (Models^[CurModel].Freq <> 0) Then S := RealToStr(Models^[CurModel].Freq, 2) + ' cm^-1'
 Else S := '(sans frequence)';
 Txt(X + 12, Y + 70, S, Black);
 Button3D(X + VP_W - 40, Y + 66, X + VP_W - 25, Y + 80, False);
 Txt(X + VP_W - 36, Y + 70, '<', Black);
 Button3D(X + VP_W - 24, Y + 66, X + VP_W - 9, Y + 80, False);
 Txt(X + VP_W - 20, Y + 70, '>', Black);
 Txt(X + 8, Y + 90, 'Progression', Black);
 PanelBar(X + 8, Y + 103, X + VP_W - 8, Y + 106, DarkGray);
 Frame3D(X + 8, Y + 103, X + VP_W - 8, Y + 106, False);
 Tx := X + 12 + RoundI((Sin(VibPhase) * 0.5 + 0.5) * (VP_W - 36));
 Button3D(Tx, Y + 97, Tx + 10, Y + 112, False);
 Txt(X + 8, Y + 120, 'Controles', Black);
 For I := 0 To 5 Do Begin
  Button3D(X + 8 + I * 26, Y + 133, X + 8 + I * 26 + 24, Y + 155, False);
  SetColor(Blue); SetFillStyle(SolidFill, Blue);
  Case I Of
   0: Begin
    Bar(X + 13, Y + 139, X + 14, Y + 149);
    Line(X + 21, Y + 139, X + 15, Y + 144); Line(X + 15, Y + 144, X + 21, Y + 149);
   End;
   1: Begin
    Line(X + 36, Y + 139, X + 30, Y + 144); Line(X + 30, Y + 144, X + 36, Y + 149);
    Line(X + 37, Y + 139, X + 31, Y + 144); Line(X + 31, Y + 144, X + 37, Y + 149);
   End;
   2: Begin
    Line(X + 58, Y + 139, X + 64, Y + 144); Line(X + 64, Y + 144, X + 58, Y + 149);
    Line(X + 57, Y + 139, X + 63, Y + 144); Line(X + 63, Y + 144, X + 57, Y + 149);
   End;
   3: Begin
    Bar(X + 84, Y + 139, X + 86, Y + 149); Bar(X + 89, Y + 139, X + 91, Y + 149);
   End;
   4: Begin
    Line(X + 110, Y + 139, X + 116, Y + 144); Line(X + 116, Y + 144, X + 110, Y + 149);
    Line(X + 111, Y + 139, X + 117, Y + 144); Line(X + 117, Y + 144, X + 111, Y + 149);
   End;
   5: Begin
    Bar(X + 143, Y + 139, X + 144, Y + 149);
    Line(X + 135, Y + 139, X + 141, Y + 144); Line(X + 141, Y + 144, X + 135, Y + 149);
   End;
  End;
 End;
 Txt(X + 8, Y + 164, 'Delai (secondes)', Black);
 PanelBar(X + 8, Y + 177, X + VP_W - 8, Y + 180, DarkGray);
 Frame3D(X + 8, Y + 177, X + VP_W - 8, Y + 180, False);
 Tx := X + 10 + RoundI(VibDelay * (VP_W - 30));
 Button3D(Tx, Y + 171, Tx + 10, Y + 186, False);
 Txt(X + 8, Y + 190, '0.0', Black);
 Txt(X + VP_W Div 2 - 12, Y + 190, '0.5', Black);
 Txt(X + VP_W - 32, Y + 190, '1.0', Black);
 SetColor(DarkGray); Rectangle(X + 6, Y + 204, X + VP_W - 6, Y + VP_H - 6);
 PanelBar(X + 12, Y + 200, X + 12 + TextWidth('Sauver le mode') + 6, Y + 208, LightGray);
 Txt(X + 15, Y + 201, 'Sauver le mode', Black);
 Button3D(X + 12, Y + 214, X + VP_W - 12, Y + 232, False);
 Txt(X + 18, Y + 219, 'Enregistrer (XYZ)', Black);
 Button3D(X + 12, Y + 236, X + VP_W - 12, Y + 254, False);
 Txt(X + 18, Y + 241, 'Exporter un film', Black);
 ShowMouse;
End;

Procedure SaveVibFile(Movie: Boolean);
Var
 FN : String;
 F : Text;
 I, M, K, NF : Integer;
 Nrm, Mg, Mx, Ph : Real;
Begin
 M := MaxI(CurModel, 1);
 If Not Models^[M].Vib Then Begin ConsoleErr('ce modele n''a pas de vibration'); Exit; End;
 If Movie Then FN := 'VIBRATION.XYZ' Else FN := 'MODE.XYZ';
 If Not FileSaveDialog('Enregistrer le mode normal', FN) Then Exit;
 {$I-} Assign(F, FN); Rewrite(F); {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire :|' + FN); Exit; End;
 Mx := 0;
 For I := Models^[M].First To Models^[M].Last Do Begin
  Mg := Sqrt(Sqr(Vib^[I].X) + Sqr(Vib^[I].Y) + Sqr(Vib^[I].Z));
  If Mg > Mx Then Mx := Mg;
 End;
 If Mx > 1E-6 Then Nrm := 1.0 / Mx Else Nrm := 0;
 If Movie Then NF := 20 Else NF := 1;
 For K := 0 To NF - 1 Do Begin
  Ph := Sin(2 * PI_VAL * K / NF) * VibScale * Nrm;
  WriteLn(F, Models^[M].Last - Models^[M].First + 1);
  WriteLn(F, Models^[M].Name, '  ', RealToStr(Models^[M].Freq, 2), ' cm-1');
  For I := Models^[M].First To Models^[M].Last Do Begin
   If Movie Then
    WriteLn(F, FixedStr(ElemDB[Atoms^[I].ElemZ].Sym, 3), ' ',
            RealFmt(Atoms^[I].X + Vib^[I].X * Ph, 12, 6), ' ',
            RealFmt(Atoms^[I].Y + Vib^[I].Y * Ph, 12, 6), ' ',
            RealFmt(Atoms^[I].Z + Vib^[I].Z * Ph, 12, 6))
   Else
    WriteLn(F, FixedStr(ElemDB[Atoms^[I].ElemZ].Sym, 3), ' ',
            RealFmt(Atoms^[I].X, 12, 6), ' ', RealFmt(Atoms^[I].Y, 12, 6), ' ',
            RealFmt(Atoms^[I].Z, 12, 6), ' ', RealFmt(Vib^[I].X, 10, 5), ' ',
            RealFmt(Vib^[I].Y, 10, 5), ' ', RealFmt(Vib^[I].Z, 10, 5));
  End;
 End;
 Close(F);
 ConsoleSay('Mode enregistre : ' + FN);
End;

{ Clic dans le panneau de vibration. Renvoie True si le clic est dans le panneau }
Function VibPanelClick(Mx, My: Integer; Press: Boolean): Boolean;
Var
 X, Y, I : Integer;
 Fr : Real;
Begin
 VibPanelClick := False;
 If Not (VibPanel And HasVibData) Then Exit;
 X := VpX; Y := VpY;
 If Not PtIn(Mx, My, X, Y, X + VP_W, Y + VP_H) Then Exit;
 VibPanelClick := True;
 If Not Press Then Begin
  { glissement du curseur }
  If VibDrag = 1 Then Begin
   Fr := (Mx - X - 12) / (VP_W - 36);
   VibPhase := ArcTan(0) + (ClampR(Fr, 0, 1) * 2 - 1) * (PI_VAL / 2);
   VibPaused := True;
   NeedScene := True;
  End
  Else If VibDrag = 2 Then Begin
   VibDelay := ClampR((Mx - X - 10) / (VP_W - 30), 0, 1);
   NeedScene := True;
  End;
  Exit;
 End;
 VibDrag := 0;
 If PtIn(Mx, My, X + VP_W - 19, Y + 4, X + VP_W - 4, Y + 16) Then Begin
  VibPanel := False; NeedScene := True; NeedChrome := True; Exit;
 End;
 If PtIn(Mx, My, X + VP_W - 40, Y + 66, X + VP_W - 25, Y + 80) Then Begin VibNextMode(-1); Exit; End;
 If PtIn(Mx, My, X + VP_W - 24, Y + 66, X + VP_W - 9, Y + 80) Then Begin VibNextMode(1); Exit; End;
 If PtIn(Mx, My, X + 8, Y + 95, X + VP_W - 8, Y + 114) Then Begin
  VibDrag := 1; VibPanelClick(Mx, My, False); Exit;
 End;
 If PtIn(Mx, My, X + 8, Y + 169, X + VP_W - 8, Y + 188) Then Begin
  VibDrag := 2; VibPanelClick(Mx, My, False); Exit;
 End;
 If PtIn(Mx, My, X + 8, Y + 133, X + VP_W - 8, Y + 155) Then Begin
  I := (Mx - X - 8) Div 26;
  Case I Of
   0: Begin VibPhase := 0; VibPaused := True; End;
   1: Begin VibPhase := VibPhase - 0.4; VibPaused := True; End;
   2: VibPaused := False;
   3: VibPaused := True;
   4: Begin VibPhase := VibPhase + 0.4; VibPaused := True; End;
   5: Begin VibPhase := PI_VAL / 2; VibPaused := True; End;
  End;
  NeedScene := True;
  Exit;
 End;
 If PtIn(Mx, My, X + 12, Y + 214, X + VP_W - 12, Y + 232) Then Begin
  SaveVibFile(False); NeedScene := True; NeedChrome := True; Exit;
 End;
 If PtIn(Mx, My, X + 12, Y + 236, X + VP_W - 12, Y + 254) Then Begin
  SaveVibFile(True); NeedScene := True; NeedChrome := True; Exit;
 End;
End;

{ --------------------------------------------------------------------- }
{ Redessin                                                              }
{ --------------------------------------------------------------------- }

Procedure FullRepaint; Forward;

Procedure RedrawAll;
Begin
 If NeedLayout Then Begin
  NeedLayout := False;
  FullRepaint;
  Exit;
 End;
 If NeedChrome Then DrawChrome;
 If NeedScene Then Begin
  NeedScene := False;
  RenderCanvas;
  DrawVibPanel;
  NeedStatus := True;
 End;
 If NeedStatus Then Begin
  HideMouse;
  DrawStatusBar;
  ShowMouse;
  NeedStatus := False;
 End;
End;

Procedure FullRepaint;
Begin
 NeedLayout := False;
 HideMouse;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, 0, SW - 1, SH - 1);
 ShowMouse;
 Layout;
 NeedChrome := True; NeedScene := True;
 RedrawAll;
End;

{ --------------------------------------------------------------------- }
{ Entrees clavier / souris de la boucle principale                      }
{ --------------------------------------------------------------------- }

Procedure OpenBarMenu(Idx: Integer);
Var Cmd: Integer;
Begin
 BarLayout;
 HideMouse;
 DrawMenuBar(Idx);
 ShowMouse;
 Cmd := RunMenu(PM_B_FILE + Idx, BarX[Idx], TITLE_H + MENU_H, Idx);
 NeedChrome := True;
 If Cmd <> 0 Then Begin
  DoCmd(Cmd);
  NeedScene := True;
 End;
 RedrawAll;
 PollMouse;
 OldMB := MouseB;
End;

Procedure OpenPopup(X, Y: Integer);
Var Cmd: Integer;
Begin
 Cmd := RunMenu(PM_ROOT, X, Y, -1);
 If Cmd <> 0 Then DoCmd(Cmd);
 NeedScene := True; NeedChrome := True;
 RedrawAll;
 PollMouse;
 OldMB := MouseB;
End;

Procedure OpenColorPopup(X, Y: Integer);
Var Cmd: Integer;
Begin
 Cmd := RunMenu(PM_COLOR, X, Y, -1);
 If Cmd <> 0 Then DoCmd(Cmd);
 NeedScene := True; NeedChrome := True;
 RedrawAll;
 PollMouse;
 OldMB := MouseB;
End;

Procedure RunToolButton(I: Integer);
Begin
 If Not TBEnabled(TBtn[I].Cmd) Then Exit;
 If TBtn[I].Cmd = CM_TB_COLOR Then OpenColorPopup(TBtn[I].X, TBY + 24)
 Else Begin
  If TBActive(TBtn[I].Cmd) And (TBtn[I].Cmd >= CM_PICK_BASE) And (TBtn[I].Cmd <= CM_PICK_BASE + 9) Then
   DoCmd(CM_PICK_BASE + PICK_IDENT)
  Else If (TBtn[I].Cmd = CM_MEAS_DIST) And (PickMode = PICK_DIST) Then DoCmd(CM_PICK_BASE + PICK_IDENT)
  Else DoCmd(TBtn[I].Cmd);
 End;
 NeedChrome := True;
End;

Procedure ConsoleKey(Ch: Char; Ext: Boolean);
Var I: Integer;
Begin
 If Ext Then Begin
  Case Ch Of
   #72: If ConHistN > 0 Then Begin
    If ConHistPos > 0 Then Dec(ConHistPos);
    ConInput := ConHist[ConHistPos];
   End;
   #80: If ConHistN > 0 Then Begin
    If ConHistPos < ConHistN - 1 Then Begin Inc(ConHistPos); ConInput := ConHist[ConHistPos]; End
    Else Begin ConHistPos := ConHistN; ConInput := ''; End;
   End;
  End;
  NeedChrome := True;
  Exit;
 End;
 If Ch = #13 Then Begin
  If ConInput <> '' Then Begin
   ConsoleEcho('$ ' + ConInput);
   If (ConHistN = 0) Or (ConHist[ConHistN - 1] <> ConInput) Then Begin
    If ConHistN >= 16 Then Begin
     For I := 0 To 14 Do ConHist[I] := ConHist[I + 1];
     ConHistN := 15;
    End;
    ConHist[ConHistN] := Copy(ConInput, 1, 80);
    Inc(ConHistN);
   End;
   ConHistPos := ConHistN;
   ExecLine(ConInput);
   ConInput := '';
  End;
 End
 Else If Ch = #27 Then Begin
  ConsoleOn := False;
  NeedLayout := True;
 End
 Else If Ch = #8 Then Begin
  If Length(ConInput) > 0 Then Delete(ConInput, Length(ConInput), 1);
 End
 Else If (Ord(Ch) >= 32) And (Length(ConInput) < 78) Then ConInput := ConInput + Ch;
 NeedChrome := True;
End;

Procedure HandleKey;
Var
 Ch : Char;
 Ext : Boolean;
Begin
 Ch := GetKey(Ext);
 If Ext Then Begin
  Case Ch Of
   #59: DoCmd(CM_HELP);
   #61: DoCmd(CM_RELOAD);
   #62: DoCmd(CM_TOG_CONSOLE);
   #68: OpenBarMenu(0);
  Else
   If ConsoleOn And ((Ch = #72) Or (Ch = #80)) Then ConsoleKey(Ch, True)
   Else If Loaded Then Begin
    Case Ch Of
     #72: Begin RotateView(0, -5); NeedScene := True; End;
     #80: Begin RotateView(0, 5); NeedScene := True; End;
     #75: Begin RotateView(1, -5); NeedScene := True; End;
     #77: Begin RotateView(1, 5); NeedScene := True; End;
     #73: SetZoom(ZoomF * 1.15);
     #81: SetZoom(ZoomF / 1.15);
     #71: ResetView;
    End;
   End;
  End;
  Exit;
 End;
 { touches de controle }
 Case Ch Of
  #15: Begin DoCmd(CM_OPEN); Exit; End;
  #19: Begin DoCmd(CM_SAVEPDB); Exit; End;
  #17: Begin Running := False; Exit; End;
  #1 : Begin DoCmd(CM_SELALL); Exit; End;
 End;
 If ConsoleOn Then Begin ConsoleKey(Ch, False); Exit; End;
 Case Ch Of
  ' ': DoCmd(CM_TOG_SPIN);
  '+', '=': DoCmd(CM_ZOOM_IN);
  '-', '_': DoCmd(CM_ZOOM_OUT);
  'h', 'H': DoCmd(CM_TOG_HYDRO);
  'a', 'A': DoCmd(CM_TOG_AXES);
  'b', 'B': DoCmd(CM_TOG_BBOX);
  'v', 'V': If HasVibData Then SetVib(Not VibOn);
  '<', ',': DoCmd(CM_FR_PREV);
  '>', '.': DoCmd(CM_FR_NEXT);
  #27: Begin
   If VibPanel Then Begin VibPanel := False; NeedScene := True; End;
  End;
 End;
End;

Function ToolAt(Mx, My: Integer): Integer;
Var A: Integer;
Begin
 ToolAt := -1;
 If (My < TBY) Or (My > TBY + 23) Then Exit;
 For A := 0 To NUM_TB - 1 Do
  If (Mx >= TBtn[A].X) And (Mx <= TBtn[A].X + 23) Then ToolAt := A;
End;

Procedure HandleMouse;
Var
 Dx, Dy, I, A : Integer;
 L, R, Mi : Boolean;
 OL, OR_, OM : Boolean;
 InCv : Boolean;
Begin
 PollMouse;
 L := (MouseB And 1) <> 0; R := (MouseB And 2) <> 0; Mi := (MouseB And 4) <> 0;
 OL := (OldMB And 1) <> 0; OR_ := (OldMB And 2) <> 0; OM := (OldMB And 4) <> 0;
 InCv := PtIn(MouseX, MouseY, CvX1, CvY1, CvX2, CvY2);
 { survol de la barre d'outils }
 If (DragMode = 0) And Not L And Not R Then Begin
  I := ToolAt(MouseX, MouseY);
  If I <> ToolHover Then Begin ToolHover := I; NeedChrome := True; End;
  { atome sous le curseur }
  If InCv And Loaded And Not (VibPanel And PtIn(MouseX, MouseY, VpX, VpY, VpX + VP_W, VpY + VP_H)) Then Begin
   A := AtomAt(MouseX, MouseY);
   If A <> HoverAtom Then Begin HoverAtom := A; NeedStatus := True; End;
  End
  Else If HoverAtom <> -1 Then Begin HoverAtom := -1; NeedStatus := True; End;
 End;
 { appui bouton gauche }
 If L And Not OL Then Begin
  If PtIn(MouseX, MouseY, SW - 19, 2, SW - 3, TITLE_H - 3) Then Begin Running := False; End
  Else If (MouseY >= TITLE_H) And (MouseY < TITLE_H + MENU_H) Then Begin
   For I := 0 To NUM_BAR - 1 Do
    If (MouseX >= BarX[I]) And (MouseX < BarX[I] + BarW[I]) Then Begin OpenBarMenu(I); Exit; End;
  End
  Else If ToolAt(MouseX, MouseY) >= 0 Then Begin
   ToolHover := ToolAt(MouseX, MouseY);
   TBDown := ToolHover;
   NeedChrome := True;
   RedrawAll;
  End
  Else If VibPanelClick(MouseX, MouseY, True) Then Begin
   DragMode := 9;
  End
  Else If InCv Then Begin
   DragMode := 1; DragMoved := False;
   DragX := MouseX; DragY := MouseY; PressX := MouseX; PressY := MouseY;
  End;
 End;
 { appui bouton droit }
 If R And Not OR_ And InCv And (DragMode = 0) Then Begin
  DragMode := 2; DragMoved := False;
  DragX := MouseX; DragY := MouseY; PressX := MouseX; PressY := MouseY;
 End;
 { appui bouton milieu }
 If Mi And Not OM And InCv And (DragMode = 0) Then Begin
  DragMode := 3; DragMoved := True;
  DragX := MouseX; DragY := MouseY;
 End;
 { glissement }
 If (DragMode <> 0) And (L Or R Or Mi) Then Begin
  Dx := MouseX - DragX; Dy := MouseY - DragY;
  If (Abs(MouseX - PressX) > 3) Or (Abs(MouseY - PressY) > 3) Then DragMoved := True;
  Case DragMode Of
   1: If DragMoved And Loaded And ((Dx <> 0) Or (Dy <> 0)) Then Begin
    RotateView(1, Dx * 0.6); RotateView(0, Dy * 0.6);
    DragX := MouseX; DragY := MouseY;
    NeedScene := True;
   End;
   2: If DragMoved And ((Dx <> 0) Or (Dy <> 0)) Then Begin
    PanX := PanX + Dx; PanY := PanY + Dy;
    DragX := MouseX; DragY := MouseY;
    NeedScene := True;
   End;
   3: If Dy <> 0 Then Begin
    SetZoom(ZoomF * (1.0 - Dy * 0.01));
    DragX := MouseX; DragY := MouseY;
   End;
   9: VibPanelClick(MouseX, MouseY, False);
  End;
 End;
 { relachement }
 If (MouseB = 0) And (OldMB <> 0) Then Begin
  If TBDown >= 0 Then Begin
   I := TBDown; TBDown := -1;
   NeedChrome := True; RedrawAll;
   If ToolAt(MouseX, MouseY) = I Then RunToolButton(I);
  End;
  If (DragMode = 1) And Not DragMoved Then Begin
   OrthoRM;
   If Loaded Then Begin
    A := AtomAt(PressX, PressY);
    DoPick(A);
   End;
  End
  Else If (DragMode = 2) And Not DragMoved Then Begin
   DragMode := 0;
   OldMB := 0;
   OpenPopup(MouseX, MouseY);
  End;
  If DragMode = 9 Then VibDrag := 0;
  DragMode := 0;
  OrthoRM;
 End;
 OldMB := MouseB;
End;

{ --------------------------------------------------------------------- }
{ Boucle principale                                                     }
{ --------------------------------------------------------------------- }

Procedure Animate;
Begin
 If SpinOn And Loaded Then Begin
  If ElapsedMs(LastSpin) >= 30 Then Begin
   RotateView(SpinAxis, SpinSpeed);
   LastSpin := TickMs;
   NeedScene := True;
  End;
 End;
 If AnimOn And (NModels > 1) Then Begin
  If ElapsedMs(LastAnim) >= AnimDelay Then Begin
   AnimTick;
   LastAnim := TickMs;
  End;
 End;
 If VibOn And HasVibData And Not VibPaused Then Begin
  If ElapsedMs(LastVib) >= 25 + Round(VibDelay * 400) Then Begin
   VibPhase := VibPhase + 0.35;
   If VibPhase > 100 * PI_VAL Then VibPhase := 0;
   LastVib := TickMs;
   NeedScene := True;
  End;
 End;
End;

{ Attente avec animation (commande delay des scripts) ; Echap interrompt }
Procedure WaitMs(Ms: LongInt);
Var
 T0 : LongInt;
 K : Char;
Begin
 T0 := TickMs;
 While (ElapsedMs(T0) < Ms) And Running And Not ScriptAbort Do Begin
  If KeyPressed Then Begin
   K := ReadKey;
   If K = #0 Then Begin If KeyPressed Then ReadKey; End
   Else If K = #27 Then ScriptAbort := True;
  End;
  PollMouse;
  Animate;
  If NeedScene Or NeedChrome Or NeedStatus Or NeedLayout Then RedrawAll;
  Delay(10);
 End;
End;

Procedure MainLoop;
Begin
 FullRepaint;
 LastSpin := TickMs; LastAnim := TickMs; LastVib := TickMs;
 If PendingScript <> '' Then Begin
  ExecLine('script ' + PendingScript);
  PendingScript := '';
  RedrawAll;
 End;
 While Running Do Begin
  While KeyPressed And Running Do HandleKey;
  HandleMouse;
  Animate;
  If NeedScene Or NeedChrome Or NeedStatus Then RedrawAll;
  Delay(8);
 End;
End;

Procedure InitializeAll;
Begin
 ClearMol;
 IdentityRM;
 ZoomF := 1.0; PanX := 0; PanY := 0; MolRad := 5.0;
 Persp := True; BgCol := Black;
 ShowHydro := True; ShowWater := True; ShowHB := False; ShowHalo := False;
 ShowAxes := False; ShowBBox := False; ShowCell := False; ShowVec := False;
 SpinOn := False; SpinAxis := 1; SpinSpeed := 1.5;
 AnimOn := False; AnimMode := 1; AnimDir := 1; AnimDelay := 150;
 VibOn := False; VibPhase := 0; VibScale := 0.5; VibDelay := 0.25; VibPanel := False;
 VibPaused := False; VecScale := 1.0; VibDrag := 0;
 PickMode := PICK_IDENT; PickN := 0; MeasUnit := 0;
 NMeas := 0; NLblFmt := 0; NSelCount := 0;
 ConsoleOn := False; ConCount := 0; ConInput := ''; ConHistN := 0; ConHistPos := 0;
 DragMode := 0; OldMB := 0; ToolHover := -1; TBDown := -1; HoverAtom := -1;
 Running := True; NeedScene := True; NeedChrome := True; NeedStatus := True;
 NeedLayout := False;
 StatMsg := 'Pret. Fichier > Ouvrir... pour charger une molecule.';
 RenderMs := 0; ScriptDepth := 0; ScriptAbort := False; PendingScript := '';
 TextN := 0;
 BarCur := 0;
 InitUnitSphere;
 InitMenus;
 InitToolbar;
 BarLayout;
 Layout;
 ConsoleSay(PROG_TITLE);
 ConsoleSay('Tapez "aide" pour la liste des commandes. F4 = console, F10 = menus.');
End;

{ ===================================================================== }
{ Programme principal                                                   }
{ ===================================================================== }

Var
 InitFile, InitExt : String;
 K : Integer;

BEGIN
 If (ParamStr(1) = '/?') Or (ParamStr(1) = '--help') Or (ParamStr(1) = '-h') Then Begin
  WriteLn('JMOL - Cette commande permet de visualiser et manipuler des structures ',
          'moleculaires tridimensionnelles (proteines, acides nucleiques et ',
          'composes chimiques), a la maniere du logiciel Jmol.');
  WriteLn;
  WriteLn('Syntaxe : JMOL [fichier.PDB|.XYZ|.MOL|.SDF|.MOL2|.CIF|.SPT]');
  WriteLn;
  WriteLn(' fichier  Fichier de molecule ou script de commandes');
  Halt;
 End;
 AllocAll;
 For K := 0 To 7 Do New(SStk[K]);
 InitializeGraphics;
 MouseOk := MouseDriverFound;
 InitializeAll;
 ShowMouse;
 If ParamCount >= 1 Then Begin
  InitFile := ParamStr(1);
  InitExt := ExtractExt(InitFile);
  If (InitExt = 'SPT') Or (InitExt = 'TXT') Or (InitExt = 'JMOL') Then PendingScript := InitFile
  Else LoadFileUI(InitFile, True);
 End;
 MainLoop;
 CloseGraph;
END.
