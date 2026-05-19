{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2026
  @version: 1.00
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal 7, Free Pascal)
  @description: Clone de l'application PyMOL - visualiseur moleculaire 3D
                avec support PDB, MOL, MOL2, SDF, CIF et scripts PML.

  Ligne de commande :
    PYMOL [fichier.PDB|.MOL|.MOL2|.SDF|.CIF|.PML]

  Fonctionnalites :
    - Chargement et visualisation de molecules en 3D
    - Modes de rendu : lignes, batonnets, spheres CPK, rubans
    - Rotation / zoom / panoramique a la souris et au clavier
    - Console PyMOL> pour les commandes interactives
    - Execution de scripts .PML
    - Menus deroulants traduits en francais

  Inspirations :
    - Menus / dialogues : CLONEDRAW\PAINT.PAS
    - Boite d'ouverture de fichier : CLONEDRAW\PSP.PAS
    - Simulation / boucle / particules : AVIONIX\SOUFFLE.PAS
}

Program PYMOL;

{$IFDEF FPC}
 {$MODE TP}
 {$Q-}
 {$R-}
{$ENDIF}

Uses
 {$IFDEF FPC}
  Crt, PtcGraph, PtcCrt, PtcMouse, Dos
 {$ELSE}
  Crt, Graph, Dos
 {$ENDIF};

{ --------------------------------------------------------------------- }
{ Constantes generales                                                  }
{ --------------------------------------------------------------------- }

Const
 PROG_TITLE         = 'PyMOL Clone 1.00 - QChimie';
 PROG_PROMPT        = 'PyMOL> ';
 PROG_PROMPT_CONT   = '   ... ';
 PROG_VERSION       = '1.00';

 MENUBAR_HEIGHT     = 20;
 STATUSBAR_HEIGHT   = 20;
 CONSOLE_HEIGHT     = 96;     { hauteur de la zone console quand visible }

 MAX_ATOMS          = 4096;
 MAX_BONDS          = 8192;
 MAX_RESIDUES       = 1024;
 MAX_CHAINS         = 64;

 MAX_HISTORY        = 32;     { historique des commandes }
 MAX_CONSOLE_LINES  = 200;    { lignes conservees dans le tampon ecran }
 MAX_CMDLINE_LEN    = 120;
 CONSOLE_VISIBLE_LN = 7;

 NUM_MENUS          = 11;
 MENU_FILE          = 0;
 MENU_EDIT          = 1;
 MENU_BUILD         = 2;
 MENU_SELECT        = 3;
 MENU_DISPLAY       = 4;
 MENU_SETTING       = 5;
 MENU_SCENE         = 6;
 MENU_MOUSE         = 7;
 MENU_WIZARD        = 8;
 MENU_PLUGIN        = 9;
 MENU_HELP          = 10;

 MAX_MENU_ITEMS     = 18;

 { Modes de rendu }
 RENDER_LINES       = 0;
 RENDER_STICKS      = 1;
 RENDER_SPHERES     = 2;
 RENDER_BALLSTICK   = 3;
 RENDER_RIBBON      = 4;
 RENDER_DOTS        = 5;

 { Drapeaux d'affichage }
 SHOW_AXES          = 1;
 SHOW_GRID          = 2;
 SHOW_LABELS        = 4;
 SHOW_HYDROGENS     = 8;
 SHOW_BBOX          = 16;
 SHOW_INFO          = 32;
 SHOW_CONSOLE       = 64;

 { Type d'atome simplifie }
 ATOMK_HETERO       = 0;
 ATOMK_PROT         = 1;     { residu d'acide amine standard }
 ATOMK_NUCL         = 2;     { residu d'acide nucleique }
 ATOMK_WATER        = 3;
 ATOMK_LIGAND       = 4;

 PI_VAL             = 3.14159265358979;

{ --------------------------------------------------------------------- }
{ Types                                                                 }
{ --------------------------------------------------------------------- }

Type
 TVec3 = Record
  X, Y, Z : Real;
 End;

 TAtom = Record
  Serial   : LongInt;        { numero (PDB) }
  Name     : String[6];      { CA, N, O, FE etc. }
  AltLoc   : Char;
  ResName  : String[4];
  ChainID  : Char;
  ResSeq   : Integer;
  ICode    : Char;
  X, Y, Z  : Real;
  Occup    : Real;
  TempF    : Real;
  Element  : String[2];      { H, C, N, O, FE, etc. }
  ElemZ    : Byte;           { numero atomique }
  Charge   : Integer;        { charge formelle }
  Kind     : Byte;           { ATOMK_xxx }
  Selected : Boolean;
  Visible  : Boolean;
  Color    : Byte;           { couleur explicite (0 = couleur CPK auto) }
 End;

 TBond = Record
  A, B     : Integer;        { indices d'atomes (0..NAtoms-1) }
  Order    : Byte;           { 1, 2, 3 ou 4 (aromatique) }
  Visible  : Boolean;
 End;

 TResidue = Record
  ChainID  : Char;
  ResSeq   : Integer;
  ICode    : Char;
  ResName  : String[4];
  FirstAt  : Integer;        { premier indice d'atome }
  LastAt   : Integer;
  IsAA     : Boolean;
  IsNucl   : Boolean;
  SSType   : Char;           { 'H'=helice, 'E'=brin, 'L'=boucle }
 End;

 TChainInfo = Record
  ChainID  : Char;
  NRes     : Integer;
 End;

 TMolecule = Record
  Title    : String;
  Source   : String;         { en-tete (PROTEIN DATA BANK, HEADER...) }
  FileName : String;
  NAtoms   : Integer;
  NBonds   : Integer;
  NRes     : Integer;
  NChains  : Integer;
  Atoms    : Array[0..MAX_ATOMS-1] Of TAtom;
  Bonds    : Array[0..MAX_BONDS-1] Of TBond;
  Resi     : Array[0..MAX_RESIDUES-1] Of TResidue;
  Chains   : Array[0..MAX_CHAINS-1] Of TChainInfo;
  MinX, MaxX, MinY, MaxY, MinZ, MaxZ : Real;
  CenterX, CenterY, CenterZ : Real;
  Radius   : Real;
  Loaded   : Boolean;
 End;

 TMenuItem = Record
  Caption   : String[34];
  Shortcut  : String[12];
  Enabled   : Boolean;
  Checked   : Boolean;
  Separator : Boolean;
 End;

 TMenu = Record
  Title     : String[16];
  Items     : Array[0..MAX_MENU_ITEMS-1] Of TMenuItem;
  ItemCount : Byte;
  X, Y, Width : Integer;
 End;

 TElementInfo = Record
  Z       : Byte;
  Sym     : String[2];
  Name    : String[16];
  Color   : Byte;          { couleur CPK approchee (palette VGA 16) }
  Radius  : Real;          { rayon covalent (A) }
  VdW     : Real;          { rayon de van der Waals (A) }
 End;

 TElemDB = Array[0..103] Of TElementInfo;

 TConsoleLine = Record
  S    : String;
  Kind : Byte;            { 0=normal, 1=invite, 2=erreur, 3=echo }
 End;

 TCmdHistory = Array[0..MAX_HISTORY-1] Of String;

{ --------------------------------------------------------------------- }
{ Variables globales                                                    }
{ --------------------------------------------------------------------- }

Var
 { Etat graphique }
 ScreenWidth, ScreenHeight : Integer;
 CanvasTop, CanvasBottom   : Integer;

 { Molecule courante }
 Mol             : TMolecule;

 { Vue : rotation, zoom, translation }
 RotX, RotY, RotZ : Real;        { angles en degres }
 ViewScale        : Real;        { px par Angstrom }
 ViewPanX         : Integer;     { decalage en pixels }
 ViewPanY         : Integer;
 ViewDistance     : Real;        { pseudo-perspective }
 BgColor          : Byte;        { couleur de fond }

 { Affichage }
 RenderMode  : Byte;
 ViewFlags   : Word;
 SimSpinning : Boolean;          { rotation automatique }
 SpinSpeed   : Real;

 { Parametres PyMOL (`set NAME, VAL`) }
 SetCartoonDiscreteColors : Boolean;
 SetCartoonHelixColor     : Byte;
 SetCartoonSheetColor     : Byte;
 SetCartoonLoopColor      : Byte;
 SetRayShadows            : Boolean;
 SetOrthoscopic           : Boolean;
 SetSphereScale           : Real;
 SetStickRadius           : Real;
 SetCartoonTransparency   : Real;
 SetTransparency          : Real;
 SetAmbient               : Real;
 SetSpecular              : Real;
 SetAntialias             : Integer;
 SetRayTraceMode          : Integer;
 SetDashWidth             : Integer;
 SetLabelSize             : Integer;
 SetHashLimit             : LongInt;

 { Menus }
 Menus       : Array[0..NUM_MENUS-1] Of TMenu;
 MenuBarHeight : Integer;
 MenuActive  : Boolean;
 ActiveMenu  : Integer;
 ActiveMenuItem : Integer;

 { Souris }
 MouseX, MouseY       : Integer;
 MouseButtons         : Integer;
 OldMouseX, OldMouseY : Integer;
 OldMouseButtons      : Integer;
 MouseAvailable       : Boolean;
 DragMode             : Byte;    { 0=rien, 1=rotation, 2=pan, 3=zoom }
 DragStartX, DragStartY : Integer;
 DragStartRotX, DragStartRotY, DragStartRotZ : Real;
 DragStartPanX, DragStartPanY : Integer;
 DragStartScale : Real;

 { Boucle principale }
 Running    : Boolean;
 NeedRedraw : Boolean;

 { Tampon de message en bas }
 StatusMsg    : String;

 { Console PyMOL> }
 ConsoleLines : Array[0..MAX_CONSOLE_LINES-1] Of TConsoleLine;
 ConsoleCount : Integer;       { nombre total ajoute }
 ConsoleView  : Integer;       { decalage de scrolling (0 = bas) }
 ConsoleInput : String;
 ConsoleCursor: Integer;
 ConsoleVisible : Boolean;
 CmdHist      : TCmdHistory;
 CmdHistCount : Integer;
 CmdHistPos   : Integer;

 { Selection courante (par nom) }
 SelName      : String;        { nom de la selection active : 'all', 'sele', ... }

{ --------------------------------------------------------------------- }
{ Table d'elements chimiques (subset large)                             }
{ --------------------------------------------------------------------- }

Const
 ElemDB : TElemDB = (
  (Z:0;   Sym:'?';  Name:'Inconnu';     Color:LightGray;  Radius:0.70; VdW:1.50),
  (Z:1;   Sym:'H';  Name:'Hydrogene';   Color:White;      Radius:0.31; VdW:1.20),
  (Z:2;   Sym:'He'; Name:'Helium';      Color:LightCyan;  Radius:0.28; VdW:1.40),
  (Z:3;   Sym:'Li'; Name:'Lithium';     Color:LightMagenta;Radius:1.28; VdW:1.82),
  (Z:4;   Sym:'Be'; Name:'Beryllium';   Color:LightGreen; Radius:0.96; VdW:1.53),
  (Z:5;   Sym:'B';  Name:'Bore';        Color:Brown;      Radius:0.84; VdW:1.92),
  (Z:6;   Sym:'C';  Name:'Carbone';     Color:LightGray;  Radius:0.76; VdW:1.70),
  (Z:7;   Sym:'N';  Name:'Azote';       Color:LightBlue;  Radius:0.71; VdW:1.55),
  (Z:8;   Sym:'O';  Name:'Oxygene';     Color:LightRed;   Radius:0.66; VdW:1.52),
  (Z:9;   Sym:'F';  Name:'Fluor';       Color:Green;      Radius:0.57; VdW:1.47),
  (Z:10;  Sym:'Ne'; Name:'Neon';        Color:LightCyan;  Radius:0.58; VdW:1.54),
  (Z:11;  Sym:'Na'; Name:'Sodium';      Color:Magenta;    Radius:1.66; VdW:2.27),
  (Z:12;  Sym:'Mg'; Name:'Magnesium';   Color:Green;      Radius:1.41; VdW:1.73),
  (Z:13;  Sym:'Al'; Name:'Aluminium';   Color:LightGray;  Radius:1.21; VdW:1.84),
  (Z:14;  Sym:'Si'; Name:'Silicium';    Color:Brown;      Radius:1.11; VdW:2.10),
  (Z:15;  Sym:'P';  Name:'Phosphore';   Color:Brown;      Radius:1.07; VdW:1.80),
  (Z:16;  Sym:'S';  Name:'Soufre';      Color:Yellow;     Radius:1.05; VdW:1.80),
  (Z:17;  Sym:'Cl'; Name:'Chlore';      Color:LightGreen; Radius:1.02; VdW:1.75),
  (Z:18;  Sym:'Ar'; Name:'Argon';       Color:LightCyan;  Radius:1.06; VdW:1.88),
  (Z:19;  Sym:'K';  Name:'Potassium';   Color:Magenta;    Radius:2.03; VdW:2.75),
  (Z:20;  Sym:'Ca'; Name:'Calcium';     Color:LightGreen; Radius:1.76; VdW:2.31),
  (Z:21;  Sym:'Sc'; Name:'Scandium';    Color:LightGray;  Radius:1.70; VdW:2.11),
  (Z:22;  Sym:'Ti'; Name:'Titane';      Color:LightGray;  Radius:1.60; VdW:2.00),
  (Z:23;  Sym:'V';  Name:'Vanadium';    Color:LightGray;  Radius:1.53; VdW:2.00),
  (Z:24;  Sym:'Cr'; Name:'Chrome';      Color:LightGray;  Radius:1.39; VdW:2.00),
  (Z:25;  Sym:'Mn'; Name:'Manganese';   Color:Magenta;    Radius:1.39; VdW:2.00),
  (Z:26;  Sym:'Fe'; Name:'Fer';         Color:LightRed;   Radius:1.32; VdW:2.00),
  (Z:27;  Sym:'Co'; Name:'Cobalt';      Color:LightMagenta;Radius:1.26; VdW:2.00),
  (Z:28;  Sym:'Ni'; Name:'Nickel';      Color:Green;      Radius:1.24; VdW:1.63),
  (Z:29;  Sym:'Cu'; Name:'Cuivre';      Color:Brown;      Radius:1.32; VdW:1.40),
  (Z:30;  Sym:'Zn'; Name:'Zinc';        Color:Brown;      Radius:1.22; VdW:1.39),
  (Z:31;  Sym:'Ga'; Name:'Gallium';     Color:Brown;      Radius:1.22; VdW:1.87),
  (Z:32;  Sym:'Ge'; Name:'Germanium';   Color:LightGray;  Radius:1.20; VdW:2.11),
  (Z:33;  Sym:'As'; Name:'Arsenic';     Color:Magenta;    Radius:1.19; VdW:1.85),
  (Z:34;  Sym:'Se'; Name:'Selenium';    Color:Brown;      Radius:1.20; VdW:1.90),
  (Z:35;  Sym:'Br'; Name:'Brome';       Color:Brown;      Radius:1.20; VdW:1.85),
  (Z:36;  Sym:'Kr'; Name:'Krypton';     Color:LightCyan;  Radius:1.16; VdW:2.02),
  (Z:37;  Sym:'Rb'; Name:'Rubidium';    Color:Magenta;    Radius:2.20; VdW:3.03),
  (Z:38;  Sym:'Sr'; Name:'Strontium';   Color:LightGreen; Radius:1.95; VdW:2.49),
  (Z:39;  Sym:'Y';  Name:'Yttrium';     Color:LightCyan;  Radius:1.90; VdW:2.00),
  (Z:40;  Sym:'Zr'; Name:'Zirconium';   Color:LightCyan;  Radius:1.75; VdW:2.00),
  (Z:41;  Sym:'Nb'; Name:'Niobium';     Color:LightCyan;  Radius:1.64; VdW:2.00),
  (Z:42;  Sym:'Mo'; Name:'Molybdene';   Color:LightCyan;  Radius:1.54; VdW:2.00),
  (Z:43;  Sym:'Tc'; Name:'Technetium';  Color:LightCyan;  Radius:1.47; VdW:2.00),
  (Z:44;  Sym:'Ru'; Name:'Ruthenium';   Color:LightCyan;  Radius:1.46; VdW:2.00),
  (Z:45;  Sym:'Rh'; Name:'Rhodium';     Color:LightCyan;  Radius:1.42; VdW:2.00),
  (Z:46;  Sym:'Pd'; Name:'Palladium';   Color:LightGray;  Radius:1.39; VdW:1.63),
  (Z:47;  Sym:'Ag'; Name:'Argent';      Color:LightGray;  Radius:1.45; VdW:1.72),
  (Z:48;  Sym:'Cd'; Name:'Cadmium';     Color:LightGray;  Radius:1.44; VdW:1.58),
  (Z:49;  Sym:'In'; Name:'Indium';      Color:LightGray;  Radius:1.42; VdW:1.93),
  (Z:50;  Sym:'Sn'; Name:'Etain';       Color:LightGray;  Radius:1.39; VdW:2.17),
  (Z:51;  Sym:'Sb'; Name:'Antimoine';   Color:LightGray;  Radius:1.39; VdW:2.06),
  (Z:52;  Sym:'Te'; Name:'Tellure';     Color:Brown;      Radius:1.38; VdW:2.06),
  (Z:53;  Sym:'I';  Name:'Iode';        Color:Magenta;    Radius:1.39; VdW:1.98),
  (Z:54;  Sym:'Xe'; Name:'Xenon';       Color:LightCyan;  Radius:1.40; VdW:2.16),
  (Z:55;  Sym:'Cs'; Name:'Cesium';      Color:Magenta;    Radius:2.44; VdW:3.43),
  (Z:56;  Sym:'Ba'; Name:'Baryum';      Color:LightGreen; Radius:2.15; VdW:2.68),
  (Z:57;  Sym:'La'; Name:'Lanthane';    Color:LightCyan;  Radius:2.07; VdW:2.40),
  (Z:58;  Sym:'Ce'; Name:'Cerium';      Color:LightCyan;  Radius:2.04; VdW:2.35),
  (Z:59;  Sym:'Pr'; Name:'Praseodyme';  Color:LightCyan;  Radius:2.03; VdW:2.39),
  (Z:60;  Sym:'Nd'; Name:'Neodyme';     Color:LightCyan;  Radius:2.01; VdW:2.29),
  (Z:61;  Sym:'Pm'; Name:'Promethium';  Color:LightCyan;  Radius:1.99; VdW:2.36),
  (Z:62;  Sym:'Sm'; Name:'Samarium';    Color:LightCyan;  Radius:1.98; VdW:2.29),
  (Z:63;  Sym:'Eu'; Name:'Europium';    Color:LightCyan;  Radius:1.98; VdW:2.33),
  (Z:64;  Sym:'Gd'; Name:'Gadolinium';  Color:LightCyan;  Radius:1.96; VdW:2.37),
  (Z:65;  Sym:'Tb'; Name:'Terbium';     Color:LightCyan;  Radius:1.94; VdW:2.21),
  (Z:66;  Sym:'Dy'; Name:'Dysprosium';  Color:LightCyan;  Radius:1.92; VdW:2.29),
  (Z:67;  Sym:'Ho'; Name:'Holmium';     Color:LightCyan;  Radius:1.92; VdW:2.16),
  (Z:68;  Sym:'Er'; Name:'Erbium';      Color:LightCyan;  Radius:1.89; VdW:2.35),
  (Z:69;  Sym:'Tm'; Name:'Thulium';     Color:LightCyan;  Radius:1.90; VdW:2.27),
  (Z:70;  Sym:'Yb'; Name:'Ytterbium';   Color:LightCyan;  Radius:1.87; VdW:2.42),
  (Z:71;  Sym:'Lu'; Name:'Lutetium';    Color:LightCyan;  Radius:1.87; VdW:2.21),
  (Z:72;  Sym:'Hf'; Name:'Hafnium';     Color:LightCyan;  Radius:1.75; VdW:2.12),
  (Z:73;  Sym:'Ta'; Name:'Tantale';     Color:LightCyan;  Radius:1.70; VdW:2.17),
  (Z:74;  Sym:'W';  Name:'Tungstene';   Color:LightCyan;  Radius:1.62; VdW:2.10),
  (Z:75;  Sym:'Re'; Name:'Rhenium';     Color:LightCyan;  Radius:1.51; VdW:2.17),
  (Z:76;  Sym:'Os'; Name:'Osmium';      Color:LightCyan;  Radius:1.44; VdW:2.16),
  (Z:77;  Sym:'Ir'; Name:'Iridium';     Color:LightCyan;  Radius:1.41; VdW:2.02),
  (Z:78;  Sym:'Pt'; Name:'Platine';     Color:LightGray;  Radius:1.36; VdW:1.75),
  (Z:79;  Sym:'Au'; Name:'Or';          Color:Yellow;     Radius:1.36; VdW:1.66),
  (Z:80;  Sym:'Hg'; Name:'Mercure';     Color:LightGray;  Radius:1.32; VdW:1.55),
  (Z:81;  Sym:'Tl'; Name:'Thallium';    Color:LightGray;  Radius:1.45; VdW:1.96),
  (Z:82;  Sym:'Pb'; Name:'Plomb';       Color:DarkGray;   Radius:1.46; VdW:2.02),
  (Z:83;  Sym:'Bi'; Name:'Bismuth';     Color:Magenta;    Radius:1.48; VdW:2.07),
  (Z:84;  Sym:'Po'; Name:'Polonium';    Color:Brown;      Radius:1.40; VdW:1.97),
  (Z:85;  Sym:'At'; Name:'Astate';      Color:Brown;      Radius:1.50; VdW:2.02),
  (Z:86;  Sym:'Rn'; Name:'Radon';       Color:LightCyan;  Radius:1.50; VdW:2.20),
  (Z:87;  Sym:'Fr'; Name:'Francium';    Color:Magenta;    Radius:2.60; VdW:3.48),
  (Z:88;  Sym:'Ra'; Name:'Radium';      Color:LightGreen; Radius:2.21; VdW:2.83),
  (Z:89;  Sym:'Ac'; Name:'Actinium';    Color:LightCyan;  Radius:2.15; VdW:2.00),
  (Z:90;  Sym:'Th'; Name:'Thorium';     Color:LightCyan;  Radius:2.06; VdW:2.43),
  (Z:91;  Sym:'Pa'; Name:'Protactinium';Color:LightCyan;  Radius:2.00; VdW:2.43),
  (Z:92;  Sym:'U';  Name:'Uranium';     Color:LightCyan;  Radius:1.96; VdW:2.41),
  (Z:93;  Sym:'Np'; Name:'Neptunium';   Color:LightCyan;  Radius:1.90; VdW:2.39),
  (Z:94;  Sym:'Pu'; Name:'Plutonium';   Color:LightCyan;  Radius:1.87; VdW:2.43),
  (Z:95;  Sym:'Am'; Name:'Americium';   Color:LightCyan;  Radius:1.80; VdW:2.44),
  (Z:96;  Sym:'Cm'; Name:'Curium';      Color:LightCyan;  Radius:1.69; VdW:2.45),
  (Z:97;  Sym:'Bk'; Name:'Berkelium';   Color:LightCyan;  Radius:1.68; VdW:2.44),
  (Z:98;  Sym:'Cf'; Name:'Californium'; Color:LightCyan;  Radius:1.68; VdW:2.45),
  (Z:99;  Sym:'Es'; Name:'Einsteinium'; Color:LightCyan;  Radius:1.65; VdW:2.45),
  (Z:100; Sym:'Fm'; Name:'Fermium';     Color:LightCyan;  Radius:1.67; VdW:2.45),
  (Z:101; Sym:'Md'; Name:'Mendelevium'; Color:LightCyan;  Radius:1.73; VdW:2.46),
  (Z:102; Sym:'No'; Name:'Nobelium';    Color:LightCyan;  Radius:1.76; VdW:2.46),
  (Z:103; Sym:'Lr'; Name:'Lawrencium';  Color:LightCyan;  Radius:1.61; VdW:2.46)
 );

{ --------------------------------------------------------------------- }
{ Helpers : chaines, conversions, math                                  }
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
 If (C >= 'a') And (C <= 'z') Then UpCh := Chr(Ord(C) - 32)
 Else UpCh := C;
End;

Function LoCh(C: Char): Char;
Begin
 If (C >= 'A') And (C <= 'Z') Then LoCh := Chr(Ord(C) + 32)
 Else LoCh := C;
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
 If V < Lo Then ClampI := Lo
 Else If V > Hi Then ClampI := Hi
 Else ClampI := V;
End;

Function ClampR(V, Lo, Hi: Real): Real;
Begin
 If V < Lo Then ClampR := Lo
 Else If V > Hi Then ClampR := Hi
 Else ClampR := V;
End;

Function AbsR(R: Real): Real;
Begin If R < 0 Then AbsR := -R Else AbsR := R; End;

Function SgnR(R: Real): Integer;
Begin
 If R > 0 Then SgnR := 1
 Else If R < 0 Then SgnR := -1
 Else SgnR := 0;
End;

Function TrimStr(Const S: String): String;
Var I, J: Integer;
Begin
 I := 1; J := Length(S);
 While (I <= J) And ((S[I] = ' ') Or (S[I] = #9)) Do Inc(I);
 While (J >= I) And ((S[J] = ' ') Or (S[J] = #9) Or
                     (S[J] = #13) Or (S[J] = #10)) Do Dec(J);
 If J < I Then TrimStr := '' Else TrimStr := Copy(S, I, J - I + 1);
End;

Function StartsWith(Const S, Prefix: String): Boolean;
Begin
 If Length(Prefix) > Length(S) Then Begin StartsWith := False; Exit; End;
 StartsWith := Copy(S, 1, Length(Prefix)) = Prefix;
End;

Function StartsWithI(Const S, Prefix: String): Boolean;
Begin
 If Length(Prefix) > Length(S) Then Begin StartsWithI := False; Exit; End;
 StartsWithI := UpStr(Copy(S, 1, Length(Prefix))) = UpStr(Prefix);
End;

Function PosI(Const Needle, S: String): Integer;
Begin
 PosI := Pos(UpStr(Needle), UpStr(S));
End;

{ Convertit en entier ; renvoie False si echec. }
Function StrToIntDef(Const S: String; Def: LongInt): LongInt;
Var V: LongInt; Code: Integer; T: String;
Begin
 T := TrimStr(S);
 Val(T, V, Code);
 If Code = 0 Then StrToIntDef := V Else StrToIntDef := Def;
End;

Function StrToReal(Const S: String; Var V: Real): Boolean;
Var Code: Integer; T: String;
Begin
 T := TrimStr(S);
 If T = '' Then Begin V := 0; StrToReal := False; Exit; End;
 Val(T, V, Code);
 StrToReal := Code = 0;
End;

{ Extrait un champ entre les positions A et B (1-based, inclusive).
  Tronque silencieusement. }
Function FieldS(Const S: String; A, B: Integer): String;
Var L: Integer;
Begin
 L := Length(S);
 If A > L Then Begin FieldS := ''; Exit; End;
 If B > L Then B := L;
 FieldS := Copy(S, A, B - A + 1);
End;

Function FieldT(Const S: String; A, B: Integer): String;
Begin FieldT := TrimStr(FieldS(S, A, B)); End;

{ Tokenize une chaine selon les espaces et tabulations. Retourne dans
  Tokens[0..N-1] avec N renvoye dans Count. }
Procedure Tokenize(Const S: String; Var Tokens: Array Of String;
                   Var Count: Integer; Max: Integer);
Var I, L: Integer; T: String;
Begin
 Count := 0;
 L := Length(S);
 I := 1;
 While I <= L Do Begin
  While (I <= L) And ((S[I] = ' ') Or (S[I] = #9) Or (S[I] = ',')) Do Inc(I);
  If I > L Then Break;
  T := '';
  While (I <= L) And (S[I] <> ' ') And (S[I] <> #9) And (S[I] <> ',') Do Begin
   T := T + S[I]; Inc(I);
  End;
  If Count < Max Then Begin
   Tokens[Count] := T; Inc(Count);
  End;
 End;
End;

{ Pour un trace texte rapide a l'ecran graphique. }
Function PadL(Const S: String; W: Integer): String;
Var R: String;
Begin
 R := S;
 While Length(R) < W Do R := ' ' + R;
 PadL := R;
End;

Function PadR(Const S: String; W: Integer): String;
Var R: String;
Begin
 R := S;
 While Length(R) < W Do R := R + ' ';
 PadR := R;
End;

{ --------------------------------------------------------------------- }
{ Gestion de la souris                                                  }
{ --------------------------------------------------------------------- }

{$IFNDEF FPC}
Procedure InitMouseDriver;
Var Regs: Registers;
Begin
 Regs.AX := $0000;
 Intr($33, Regs);
 MouseAvailable := Regs.AX = $FFFF;
End;

Procedure ShowMouseCursor;
Var Regs: Registers;
Begin
 If Not MouseAvailable Then Exit;
 Regs.AX := $0001;
 Intr($33, Regs);
End;

Procedure HideMouseCursor;
Var Regs: Registers;
Begin
 If Not MouseAvailable Then Exit;
 Regs.AX := $0002;
 Intr($33, Regs);
End;

Procedure GetMouseLocal(Var X, Y, Buttons: Integer);
Var Regs: Registers;
Begin
 If Not MouseAvailable Then Begin
  X := 0; Y := 0; Buttons := 0;
  Exit;
 End;
 Regs.AX := $0003;
 Intr($33, Regs);
 X := Regs.CX;
 Y := Regs.DX;
 Buttons := Regs.BX;
End;
{$ELSE}
Procedure InitMouseDriver;
Begin
 PtcMouse.InitMouse;
 MouseAvailable := True;
End;

Procedure ShowMouseCursor;
Begin
 PtcMouse.ShowMouse;
End;

Procedure HideMouseCursor;
Begin
 PtcMouse.HideMouse;
End;

Procedure GetMouseLocal(Var X, Y, Buttons: Integer);
Var LX, LY, LB: LongInt;
Begin
 PtcMouse.GetMouseState(LX, LY, LB);
 X := LX; Y := LY; Buttons := LB;
End;
{$ENDIF}

{ --------------------------------------------------------------------- }
{ Recherche d'un element                                                }
{ --------------------------------------------------------------------- }

{ Renvoie l'indice (Z) dans ElemDB d'apres un symbole. Renvoie 0 si
  inconnu. Tolere casse et chiffres parasites. }
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

{ Tente de deduire l'element a partir du nom d'atome PDB (4 caracteres).
  Si la chaine commence par un chiffre, prendre le caractere suivant.
  Les atomes calpha, carbone, etc., sont prefixes "CA","C","N","O",etc. }
Function ElemFromAtomName(Const Nm: String): Byte;
Var
 S : String;
 I : Integer;
 Two : String[2];
 One : String[1];
 R : Byte;
Begin
 S := TrimStr(Nm);
 If S = '' Then Begin ElemFromAtomName := 0; Exit; End;
 I := 1;
 While (I <= Length(S)) And (S[I] >= '0') And (S[I] <= '9') Do Inc(I);
 If I > Length(S) Then Begin ElemFromAtomName := 0; Exit; End;
 If I + 1 <= Length(S) Then Begin
  Two := UpCh(S[I]) + LoCh(S[I + 1]);
  { Mais beaucoup d'atomes proteines : "CA" = C alpha = Carbone, pas
    Calcium. Si le residu est un AA et que le nom est CA/CB/CG/etc.,
    on doit retourner Carbone. Cette decision est faite par l'appelant
    en regardant le residu : ici on essaie d'abord deux lettres. }
  R := 0;
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
 { Sinon, premier caractere uniquement. }
 One := UpCh(S[I]);
 R := ElemBySymbol(One);
 ElemFromAtomName := R;
End;

{ Reconnaissance du type d'atome (hetero / proteique / etc.) selon le
  nom de residu standard. }
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
  (U = 'ASX') Or (U = 'GLX') Or (U = 'CSE') Or (U = 'CSO');
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
                (U = 'TIP3') Or (U = 'SOL') Or (U = 'H2O');
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

Procedure ClearMolecule(Var M: TMolecule);
Begin
 M.Title := '';
 M.Source := '';
 M.FileName := '';
 M.NAtoms := 0;
 M.NBonds := 0;
 M.NRes := 0;
 M.NChains := 0;
 M.MinX := 0; M.MaxX := 0;
 M.MinY := 0; M.MaxY := 0;
 M.MinZ := 0; M.MaxZ := 0;
 M.CenterX := 0; M.CenterY := 0; M.CenterZ := 0;
 M.Radius := 1.0;
 M.Loaded := False;
End;

{ Ajoute un atome a la molecule. Renvoie son indice ou -1. }
Function AddAtom(Var M: TMolecule; Const A: TAtom): Integer;
Begin
 If M.NAtoms >= MAX_ATOMS Then Begin AddAtom := -1; Exit; End;
 M.Atoms[M.NAtoms] := A;
 AddAtom := M.NAtoms;
 Inc(M.NAtoms);
End;

Function AddBond(Var M: TMolecule; A, B: Integer; Order: Byte): Integer;
Var
 I : Integer;
 BR : TBond;
Begin
 If (A < 0) Or (B < 0) Or (A >= M.NAtoms) Or (B >= M.NAtoms) Or
    (A = B) Then Begin AddBond := -1; Exit; End;
 { Evite les doublons. }
 For I := 0 To M.NBonds - 1 Do
  If ((M.Bonds[I].A = A) And (M.Bonds[I].B = B)) Or
     ((M.Bonds[I].A = B) And (M.Bonds[I].B = A)) Then Begin
   AddBond := I; Exit;
  End;
 If M.NBonds >= MAX_BONDS Then Begin AddBond := -1; Exit; End;
 BR.A := A; BR.B := B; BR.Order := Order; BR.Visible := True;
 M.Bonds[M.NBonds] := BR;
 AddBond := M.NBonds;
 Inc(M.NBonds);
End;

Procedure ComputeBounds(Var M: TMolecule);
Var
 I : Integer;
 D : Real;
Begin
 If M.NAtoms = 0 Then Begin
  M.MinX := 0; M.MaxX := 0;
  M.MinY := 0; M.MaxY := 0;
  M.MinZ := 0; M.MaxZ := 0;
  M.CenterX := 0; M.CenterY := 0; M.CenterZ := 0;
  M.Radius := 1.0;
  Exit;
 End;
 M.MinX := M.Atoms[0].X; M.MaxX := M.MinX;
 M.MinY := M.Atoms[0].Y; M.MaxY := M.MinY;
 M.MinZ := M.Atoms[0].Z; M.MaxZ := M.MinZ;
 For I := 1 To M.NAtoms - 1 Do Begin
  If M.Atoms[I].X < M.MinX Then M.MinX := M.Atoms[I].X;
  If M.Atoms[I].X > M.MaxX Then M.MaxX := M.Atoms[I].X;
  If M.Atoms[I].Y < M.MinY Then M.MinY := M.Atoms[I].Y;
  If M.Atoms[I].Y > M.MaxY Then M.MaxY := M.Atoms[I].Y;
  If M.Atoms[I].Z < M.MinZ Then M.MinZ := M.Atoms[I].Z;
  If M.Atoms[I].Z > M.MaxZ Then M.MaxZ := M.Atoms[I].Z;
 End;
 M.CenterX := (M.MinX + M.MaxX) * 0.5;
 M.CenterY := (M.MinY + M.MaxY) * 0.5;
 M.CenterZ := (M.MinZ + M.MaxZ) * 0.5;
 D := 0.0;
 For I := 0 To M.NAtoms - 1 Do Begin
  D := MaxR(D, Sqrt(Sqr(M.Atoms[I].X - M.CenterX) +
                    Sqr(M.Atoms[I].Y - M.CenterY) +
                    Sqr(M.Atoms[I].Z - M.CenterZ)));
 End;
 If D < 1.0 Then D := 1.0;
 M.Radius := D;
End;

{ Construit la table des residus et chaines a partir des atomes deja
  charges. Suppose que les atomes d'un meme residu sont contigus
  (ce qui est le cas dans un PDB bien forme). }
Procedure BuildResidueIndex(Var M: TMolecule);
Var
 I, R, Cnt : Integer;
 LastCh : Char;
 LastSeq : Integer;
 LastIC  : Char;
 LastName : String;
Begin
 M.NRes := 0;
 M.NChains := 0;
 If M.NAtoms = 0 Then Exit;
 LastCh := #0; LastSeq := -32768; LastIC := #0; LastName := '';
 R := -1;
 For I := 0 To M.NAtoms - 1 Do Begin
  If (M.Atoms[I].ChainID <> LastCh) Or
     (M.Atoms[I].ResSeq <> LastSeq) Or
     (M.Atoms[I].ICode <> LastIC) Or
     (M.Atoms[I].ResName <> LastName) Then Begin
   If R >= 0 Then M.Resi[R].LastAt := I - 1;
   If M.NRes < MAX_RESIDUES Then Begin
    R := M.NRes;
    M.Resi[R].ChainID := M.Atoms[I].ChainID;
    M.Resi[R].ResSeq := M.Atoms[I].ResSeq;
    M.Resi[R].ICode := M.Atoms[I].ICode;
    M.Resi[R].ResName := M.Atoms[I].ResName;
    M.Resi[R].FirstAt := I;
    M.Resi[R].LastAt := I;
    M.Resi[R].IsAA := IsAminoAcidName(M.Atoms[I].ResName);
    M.Resi[R].IsNucl := IsNucleicName(M.Atoms[I].ResName);
    M.Resi[R].SSType := 'L';
    Inc(M.NRes);
   End;
   LastCh := M.Atoms[I].ChainID;
   LastSeq := M.Atoms[I].ResSeq;
   LastIC := M.Atoms[I].ICode;
   LastName := M.Atoms[I].ResName;
  End;
 End;
 If R >= 0 Then M.Resi[R].LastAt := M.NAtoms - 1;
 { Chaines : compactage par ChainID consecutif }
 LastCh := #1;
 Cnt := 0;
 For I := 0 To M.NRes - 1 Do Begin
  If M.Resi[I].ChainID <> LastCh Then Begin
   If M.NChains < MAX_CHAINS Then Begin
    M.Chains[M.NChains].ChainID := M.Resi[I].ChainID;
    M.Chains[M.NChains].NRes := 1;
    Inc(M.NChains);
    Cnt := 1;
   End;
   LastCh := M.Resi[I].ChainID;
  End
  Else Begin
   If M.NChains > 0 Then
    Inc(M.Chains[M.NChains - 1].NRes);
   Inc(Cnt);
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Detection automatique des liaisons (covalentes)                       }
{ --------------------------------------------------------------------- }

{ Critere simple : distance(i,j) < radius(i) + radius(j) + tolerance.
  Tolerance = 0.45 A. Exclut H-H. Limite l'explosion combinatoire en
  utilisant le bounding box et en sautant les atomes trop eloignes
  sur X. }
Procedure AutoDetectBonds(Var M: TMolecule);
Var
 I, J : Integer;
 Dx, Dy, Dz, D2, R, Lim : Real;
 ZI, ZJ : Byte;
Begin
 M.NBonds := 0;
 If M.NAtoms < 2 Then Exit;
 For I := 0 To M.NAtoms - 2 Do Begin
  ZI := M.Atoms[I].ElemZ;
  If ZI > 103 Then ZI := 0;
  For J := I + 1 To M.NAtoms - 1 Do Begin
   ZJ := M.Atoms[J].ElemZ;
   If ZJ > 103 Then ZJ := 0;
   { Pas de liaison entre deux atomes de chaines/residus differents
     si l'un est de l'eau. }
   If (M.Atoms[I].Kind = ATOMK_WATER) Or (M.Atoms[J].Kind = ATOMK_WATER) Then
    If (M.Atoms[I].ResSeq <> M.Atoms[J].ResSeq) Or
       (M.Atoms[I].ChainID <> M.Atoms[J].ChainID) Then Continue;
   Dx := M.Atoms[I].X - M.Atoms[J].X;
   If AbsR(Dx) > 2.5 Then Continue;
   Dy := M.Atoms[I].Y - M.Atoms[J].Y;
   If AbsR(Dy) > 2.5 Then Continue;
   Dz := M.Atoms[I].Z - M.Atoms[J].Z;
   If AbsR(Dz) > 2.5 Then Continue;
   D2 := Dx*Dx + Dy*Dy + Dz*Dz;
   If D2 < 0.16 Then Continue; { atomes superposes }
   R := ElemDB[ZI].Radius + ElemDB[ZJ].Radius + 0.45;
   Lim := R * R;
   If D2 <= Lim Then
    If AddBond(M, I, J, 1) < 0 Then Exit;
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Chargeur PDB                                                          }
{ --------------------------------------------------------------------- }

{ Lit une ligne ATOM ou HETATM standard.
  Champs (1-based) :
   1- 6   Record name "ATOM  " ou "HETATM"
   7-11   Serial
  13-16   Atom name
  17      AltLoc
  18-20   ResName
  22      ChainID
  23-26   ResSeq
  27      ICode
  31-38   X
  39-46   Y
  47-54   Z
  55-60   Occupancy
  61-66   Temp factor
  77-78   Element
  79-80   Charge
}
Function ParseAtomLine(Const Line: String; HetAtm: Boolean;
                      Var A: TAtom): Boolean;
Var
 SX, SY, SZ : String;
 X, Y, Z, Oc, Tf : Real;
 Sn : LongInt;
 Sym : String;
 Z2 : Byte;
Begin
 ParseAtomLine := False;
 If Length(Line) < 54 Then Exit;
 Sn := StrToIntDef(FieldT(Line, 7, 11), 0);
 A.Serial := Sn;
 A.Name := FieldT(Line, 13, 16);
 A.AltLoc := ' ';
 If Length(Line) >= 17 Then A.AltLoc := Line[17];
 A.ResName := FieldT(Line, 18, 20);
 A.ChainID := ' ';
 If Length(Line) >= 22 Then A.ChainID := Line[22];
 A.ResSeq := StrToIntDef(FieldT(Line, 23, 26), 0);
 A.ICode := ' ';
 If Length(Line) >= 27 Then A.ICode := Line[27];
 SX := FieldT(Line, 31, 38);
 SY := FieldT(Line, 39, 46);
 SZ := FieldT(Line, 47, 54);
 If Not StrToReal(SX, X) Then Exit;
 If Not StrToReal(SY, Y) Then Exit;
 If Not StrToReal(SZ, Z) Then Exit;
 A.X := X; A.Y := Y; A.Z := Z;
 Oc := 1.0; Tf := 0.0;
 If Length(Line) >= 60 Then StrToReal(FieldT(Line, 55, 60), Oc);
 If Length(Line) >= 66 Then StrToReal(FieldT(Line, 61, 66), Tf);
 A.Occup := Oc; A.TempF := Tf;
 Sym := '';
 If Length(Line) >= 78 Then Sym := FieldT(Line, 77, 78);
 If Sym = '' Then Sym := A.Name;
 Z2 := ElemBySymbol(Sym);
 If Z2 = 0 Then Z2 := ElemFromAtomName(A.Name);
 { Cas special : 'CA' dans un AA est Carbone alpha (Z=6), pas Calcium (Z=20) }
 If (Z2 = 20) And (UpStr(A.Name) = 'CA') And IsAminoAcidName(A.ResName) Then
  Z2 := 6;
 A.ElemZ := Z2;
 If Z2 = 0 Then A.Element := '?'
 Else A.Element := ElemDB[Z2].Sym;
 A.Charge := 0;
 If Length(Line) >= 80 Then Begin
  Sym := FieldT(Line, 79, 80);
  If Sym <> '' Then A.Charge := StrToIntDef(Sym, 0);
 End;
 A.Kind := AtomKindFor(A.ResName, HetAtm);
 A.Selected := False;
 A.Visible := True;
 A.Color := 0;
 ParseAtomLine := True;
End;

{ Lit un fichier PDB. Reconnait HEADER, TITLE, ATOM, HETATM, CONECT,
  END, ENDMDL. Sur ENDMDL le premier modele est conserve. }
Function LoadPDB(Const FileName: String; Var M: TMolecule): Boolean;
Const
 MAX_SSE = 256;
Type
 TSSERng = Record
  Kind     : Char;     { 'H' helix, 'E' sheet }
  ChainID  : Char;
  StartSeq : Integer;
  EndSeq   : Integer;
 End;
Var
 F : Text;
 Line, Rec : String;
 A : TAtom;
 Idx, Sn : LongInt;
 IxA, IxB : LongInt;
 I, K, J : Integer;
 ModelSeen : Boolean;
 SerialMap : Array[0..MAX_ATOMS-1] Of LongInt;
 NSerial : Integer;
 SSEs : Array[0..MAX_SSE - 1] Of TSSERng;
 NSSE : Integer;
 Ch   : Char;
 S1, S2 : Integer;
Begin
 LoadPDB := False;
 ClearMolecule(M);
 NSSE := 0;
 M.FileName := FileName;
 {$I-} Assign(F, FileName); Reset(F); {$I+}
 If IOResult <> 0 Then Exit;
 NSerial := 0;
 ModelSeen := False;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  If Line = '' Then Continue;
  Rec := UpStr(Copy(Line, 1, 6));
  If Rec = 'HEADER' Then Begin
   M.Source := TrimStr(FieldS(Line, 11, 50));
   If M.Title = '' Then M.Title := M.Source;
  End
  Else If Rec = 'TITLE ' Then Begin
   If M.Title <> '' Then M.Title := M.Title + ' ';
   M.Title := M.Title + TrimStr(FieldS(Line, 11, 80));
  End
  Else If Rec = 'HELIX ' Then Begin
   If NSSE < MAX_SSE Then Begin
    SSEs[NSSE].Kind := 'H';
    If Length(Line) >= 22 Then SSEs[NSSE].ChainID := Line[20] Else SSEs[NSSE].ChainID := ' ';
    SSEs[NSSE].StartSeq := StrToIntDef(FieldT(Line, 22, 25), 0);
    SSEs[NSSE].EndSeq   := StrToIntDef(FieldT(Line, 34, 37), 0);
    Inc(NSSE);
   End;
  End
  Else If Rec = 'SHEET ' Then Begin
   If NSSE < MAX_SSE Then Begin
    SSEs[NSSE].Kind := 'E';
    If Length(Line) >= 22 Then SSEs[NSSE].ChainID := Line[22] Else SSEs[NSSE].ChainID := ' ';
    SSEs[NSSE].StartSeq := StrToIntDef(FieldT(Line, 23, 26), 0);
    SSEs[NSSE].EndSeq   := StrToIntDef(FieldT(Line, 34, 37), 0);
    Inc(NSSE);
   End;
  End
  Else If (Rec = 'ATOM  ') Or (Rec = 'HETATM') Then Begin
   If ParseAtomLine(Line, Rec = 'HETATM', A) Then Begin
    { Filtrer altloc autre que ' ' ou 'A' pour eviter doublons. }
    If (A.AltLoc = ' ') Or (A.AltLoc = 'A') Then Begin
     Idx := AddAtom(M, A);
     If Idx >= 0 Then Begin
      If NSerial < MAX_ATOMS Then Begin
       SerialMap[NSerial] := A.Serial;
       Inc(NSerial);
      End;
     End;
    End;
   End;
  End
  Else If Rec = 'CONECT' Then Begin
   { CONECT atom_a atom_b1 atom_b2 atom_b3 atom_b4 }
   Sn := StrToIntDef(FieldT(Line, 7, 11), 0);
   IxA := -1;
   For I := 0 To M.NAtoms - 1 Do
    If M.Atoms[I].Serial = Sn Then Begin IxA := I; Break; End;
   If IxA < 0 Then Continue;
   For K := 0 To 3 Do Begin
    Sn := StrToIntDef(FieldT(Line, 12 + K * 5, 16 + K * 5), 0);
    If Sn = 0 Then Continue;
    IxB := -1;
    For I := 0 To M.NAtoms - 1 Do
     If M.Atoms[I].Serial = Sn Then Begin IxB := I; Break; End;
    If IxB >= 0 Then AddBond(M, IxA, IxB, 1);
   End;
  End
  Else If Rec = 'ENDMDL' Then Begin
   ModelSeen := True;
   Break;
  End
  Else If Rec = 'END   ' Then Break;
 End;
 Close(F);
 If M.Title = '' Then M.Title := FileName;
 M.Source := 'PROTEIN DATA BANK';
 If M.NBonds = 0 Then AutoDetectBonds(M);
 ComputeBounds(M);
 BuildResidueIndex(M);
 { Application des annotations HELIX / SHEET aux residus }
 For J := 0 To NSSE - 1 Do Begin
  Ch := SSEs[J].ChainID;
  S1 := SSEs[J].StartSeq;
  S2 := SSEs[J].EndSeq;
  For I := 0 To M.NRes - 1 Do
   If ((Ch = ' ') Or (M.Resi[I].ChainID = Ch)) And
      (M.Resi[I].ResSeq >= S1) And (M.Resi[I].ResSeq <= S2) Then
    M.Resi[I].SSType := SSEs[J].Kind;
 End;
 M.Loaded := M.NAtoms > 0;
 LoadPDB := M.Loaded;
 If ModelSeen Then ; { eviter avertissement }
End;

{ --------------------------------------------------------------------- }
{ Chargeur MOL / SDF (V2000)                                            }
{ --------------------------------------------------------------------- }

{ Format MDL :
   Ligne 1 : nom
   Ligne 2 : ligne d'en-tete
   Ligne 3 : commentaire
   Ligne 4 : counts "nnn mmm  ... V2000" (atomes, liaisons)
   N x : "  X.XXXX  Y.YYYY  Z.ZZZZ Sym ..."
   M x : "aaa bbb ttt ..." (a,b,type)
  SDF : enchaine plusieurs entrees separees par "$$$$".
}
Function LoadMOLStream(Var F: Text; Var M: TMolecule;
                      ReadHeader: Boolean): Boolean;
Var
 Line  : String;
 NA, NB, I, A1, A2, Tp : Integer;
 Counts : String;
 Sym : String;
 SX, SY, SZ : String;
 X, Y, Z : Real;
 Atm : TAtom;
 Z2  : Byte;
Begin
 LoadMOLStream := False;
 If ReadHeader Then Begin
  If Eof(F) Then Exit;
  ReadLn(F, Line); { ligne 1 : nom }
  M.Title := TrimStr(Line);
  If Eof(F) Then Exit;
  ReadLn(F, Line); { ligne 2 : en-tete }
  If M.Source = '' Then M.Source := TrimStr(Line);
  If Eof(F) Then Exit;
  ReadLn(F, Line); { ligne 3 : commentaire }
 End;
 If Eof(F) Then Exit;
 ReadLn(F, Counts);
 NA := StrToIntDef(TrimStr(FieldS(Counts, 1, 3)), 0);
 NB := StrToIntDef(TrimStr(FieldS(Counts, 4, 6)), 0);
 If (NA <= 0) Or (NA > MAX_ATOMS) Then Exit;
 For I := 1 To NA Do Begin
  If Eof(F) Then Exit;
  ReadLn(F, Line);
  SX := TrimStr(FieldS(Line, 1, 10));
  SY := TrimStr(FieldS(Line, 11, 20));
  SZ := TrimStr(FieldS(Line, 21, 30));
  Sym := TrimStr(FieldS(Line, 32, 34));
  If Sym = '' Then Sym := TrimStr(FieldS(Line, 31, 34));
  If Not StrToReal(SX, X) Then Continue;
  If Not StrToReal(SY, Y) Then Continue;
  If Not StrToReal(SZ, Z) Then Continue;
  Z2 := ElemBySymbol(Sym);
  Atm.Serial := I;
  Atm.Name := Sym + IntToStr(I);
  Atm.AltLoc := ' ';
  Atm.ResName := 'UNL';
  Atm.ChainID := ' ';
  Atm.ResSeq := 1;
  Atm.ICode := ' ';
  Atm.X := X; Atm.Y := Y; Atm.Z := Z;
  Atm.Occup := 1.0; Atm.TempF := 0.0;
  Atm.Element := Sym;
  Atm.ElemZ := Z2;
  Atm.Charge := 0;
  Atm.Kind := ATOMK_LIGAND;
  Atm.Selected := False;
  Atm.Visible := True;
  Atm.Color := 0;
  AddAtom(M, Atm);
 End;
 For I := 1 To NB Do Begin
  If Eof(F) Then Exit;
  ReadLn(F, Line);
  A1 := StrToIntDef(TrimStr(FieldS(Line, 1, 3)), 0);
  A2 := StrToIntDef(TrimStr(FieldS(Line, 4, 6)), 0);
  Tp := StrToIntDef(TrimStr(FieldS(Line, 7, 9)), 1);
  If (A1 >= 1) And (A2 >= 1) And (A1 <= M.NAtoms) And (A2 <= M.NAtoms) Then
   AddBond(M, A1 - 1, A2 - 1, Tp);
 End;
 LoadMOLStream := True;
End;

Function LoadMOL(Const FileName: String; Var M: TMolecule): Boolean;
Var F: Text;
Begin
 LoadMOL := False;
 ClearMolecule(M);
 M.FileName := FileName;
 {$I-} Assign(F, FileName); Reset(F); {$I+}
 If IOResult <> 0 Then Exit;
 If LoadMOLStream(F, M, True) Then Begin
  Close(F);
  If M.Title = '' Then M.Title := FileName;
  ComputeBounds(M);
  BuildResidueIndex(M);
  M.Loaded := M.NAtoms > 0;
  LoadMOL := M.Loaded;
 End
 Else Begin
  Close(F);
 End;
End;

Function LoadSDF(Const FileName: String; Var M: TMolecule): Boolean;
Var F: Text;
Begin
 { SDF = MOL + "$$$$" + suivant. On charge la 1re entree uniquement. }
 LoadSDF := False;
 ClearMolecule(M);
 M.FileName := FileName;
 {$I-} Assign(F, FileName); Reset(F); {$I+}
 If IOResult <> 0 Then Exit;
 If LoadMOLStream(F, M, True) Then Begin
  Close(F);
  If M.Title = '' Then M.Title := FileName;
  ComputeBounds(M);
  BuildResidueIndex(M);
  M.Loaded := M.NAtoms > 0;
  LoadSDF := M.Loaded;
 End
 Else Begin
  Close(F);
 End;
End;

{ --------------------------------------------------------------------- }
{ Chargeur MOL2 (Tripos)                                                }
{ --------------------------------------------------------------------- }

{ Format MOL2 (extrait) :
   @<TRIPOS>MOLECULE
   <nom>
   <n_atoms n_bonds ...>
   <SMALL>
   @<TRIPOS>ATOM
   id name x y z type subst_id subst_name charge ...
   @<TRIPOS>BOND
   id a b type
}
Function LoadMOL2(Const FileName: String; Var M: TMolecule): Boolean;
Var
 F : Text;
 Line : String;
 Section : Byte;     { 0=rien, 1=ATOM, 2=BOND, 3=MOLECULE }
 NA, NB : Integer;
 MolLineCnt : Integer;
 Tokens : Array[0..15] Of String;
 NTok : Integer;
 Sym  : String;
 X, Y, Z, Q : Real;
 Atm : TAtom;
 A1, A2 : Integer;
 BondType : String;
 OrderB : Byte;
 Idx : Integer;
 ParsedCounts : Boolean;
Begin
 LoadMOL2 := False;
 ClearMolecule(M);
 M.FileName := FileName;
 {$I-} Assign(F, FileName); Reset(F); {$I+}
 If IOResult <> 0 Then Exit;
 Section := 0;
 NA := 0; NB := 0;
 MolLineCnt := 0;
 ParsedCounts := False;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  If Line = '' Then Continue;
  If StartsWith(Line, '@<TRIPOS>') Then Begin
   If UpStr(Line) = '@<TRIPOS>MOLECULE' Then Begin
    Section := 3; MolLineCnt := 0;
   End
   Else If UpStr(Line) = '@<TRIPOS>ATOM' Then Section := 1
   Else If UpStr(Line) = '@<TRIPOS>BOND' Then Section := 2
   Else Section := 0;
   Continue;
  End;
  Case Section Of
   3: Begin
    Inc(MolLineCnt);
    If MolLineCnt = 1 Then M.Title := Line
    Else If MolLineCnt = 2 Then Begin
     Tokenize(Line, Tokens, NTok, 16);
     If NTok >= 1 Then NA := StrToIntDef(Tokens[0], 0);
     If NTok >= 2 Then NB := StrToIntDef(Tokens[1], 0);
     ParsedCounts := True;
    End;
   End;
   1: Begin
    Tokenize(Line, Tokens, NTok, 16);
    If NTok < 6 Then Continue;
    Atm.Serial := StrToIntDef(Tokens[0], M.NAtoms + 1);
    Atm.Name := Copy(Tokens[1], 1, 6);
    If Not StrToReal(Tokens[2], X) Then Continue;
    If Not StrToReal(Tokens[3], Y) Then Continue;
    If Not StrToReal(Tokens[4], Z) Then Continue;
    { Le sybyl type est sous la forme C.3, C.ar, N.am etc. }
    Sym := Tokens[5];
    Idx := Pos('.', Sym);
    If Idx > 0 Then Sym := Copy(Sym, 1, Idx - 1);
    Atm.Element := Sym;
    Atm.ElemZ := ElemBySymbol(Sym);
    Atm.X := X; Atm.Y := Y; Atm.Z := Z;
    Atm.AltLoc := ' ';
    Atm.ResName := 'UNL';
    Atm.ChainID := ' ';
    Atm.ResSeq := 1;
    Atm.ICode := ' ';
    Atm.Occup := 1.0; Atm.TempF := 0.0;
    Atm.Charge := 0;
    If NTok >= 9 Then If StrToReal(Tokens[8], Q) Then Atm.Charge := Round(Q);
    If NTok >= 7 Then Atm.ResSeq := StrToIntDef(Tokens[6], 1);
    If NTok >= 8 Then Atm.ResName := Copy(Tokens[7], 1, 4);
    Atm.Kind := AtomKindFor(Atm.ResName, True);
    Atm.Selected := False;
    Atm.Visible := True;
    Atm.Color := 0;
    AddAtom(M, Atm);
   End;
   2: Begin
    Tokenize(Line, Tokens, NTok, 16);
    If NTok < 4 Then Continue;
    A1 := StrToIntDef(Tokens[1], 0) - 1;
    A2 := StrToIntDef(Tokens[2], 0) - 1;
    BondType := UpStr(Tokens[3]);
    If BondType = '1' Then OrderB := 1
    Else If BondType = '2' Then OrderB := 2
    Else If BondType = '3' Then OrderB := 3
    Else If (BondType = 'AR') Or (BondType = 'AM') Then OrderB := 4
    Else OrderB := 1;
    AddBond(M, A1, A2, OrderB);
   End;
  End;
 End;
 Close(F);
 If M.Title = '' Then M.Title := FileName;
 M.Source := 'Tripos MOL2';
 If M.NBonds = 0 Then AutoDetectBonds(M);
 ComputeBounds(M);
 BuildResidueIndex(M);
 M.Loaded := M.NAtoms > 0;
 LoadMOL2 := M.Loaded;
 If ParsedCounts Then ;
 If NA = 0 Then ;
 If NB = 0 Then ;
End;

{ --------------------------------------------------------------------- }
{ Chargeur CIF (mmCIF/PDBx simplifie)                                   }
{ --------------------------------------------------------------------- }

{ Lecture minimaliste de la boucle _atom_site. On lit l'ordre des colonnes
  declare dans la loop_ avant les valeurs, puis chaque ligne de valeurs
  est tokenisee. }
Function LoadCIF(Const FileName: String; Var M: TMolecule): Boolean;
Var
 F : Text;
 Line : String;
 InLoop : Boolean;
 LoopCols : Array[0..31] Of String;
 NCols : Integer;
 ColAtom, ColLabel, ColComp, ColSeq, ColAsym, ColAlt : Integer;
 ColX, ColY, ColZ, ColOcc, ColB, ColEl, ColIns : Integer;
 ColModel, ColAuthAsym, ColAuthSeq : Integer;
 Tokens : Array[0..31] Of String;
 NTok, I : Integer;
 Atm : TAtom;
 X, Y, Z, V : Real;
 IsAtomSiteLoop : Boolean;
 ColTagsCollected : Boolean;
 ModelOk : Integer;
 ModelStr : String;
Begin
 LoadCIF := False;
 ClearMolecule(M);
 M.FileName := FileName;
 {$I-} Assign(F, FileName); Reset(F); {$I+}
 If IOResult <> 0 Then Exit;
 InLoop := False;
 IsAtomSiteLoop := False;
 ColTagsCollected := False;
 NCols := 0;
 ColAtom := -1; ColLabel := -1; ColComp := -1; ColSeq := -1;
 ColAsym := -1; ColAlt := -1; ColX := -1; ColY := -1; ColZ := -1;
 ColOcc := -1; ColB := -1; ColEl := -1; ColIns := -1;
 ColModel := -1; ColAuthAsym := -1; ColAuthSeq := -1;
 ModelOk := 1;
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  If Line = '' Then Continue;
  If Line[1] = '#' Then Continue;
  If StartsWithI(Line, 'data_') Then Begin
   M.Title := Copy(Line, 6, Length(Line) - 5);
   Continue;
  End;
  If StartsWithI(Line, '_struct.title') Then Begin
   If Pos(' ', Line) > 0 Then
    M.Title := TrimStr(Copy(Line, Pos(' ', Line), Length(Line)));
   Continue;
  End;
  If LoStr(Line) = 'loop_' Then Begin
   InLoop := True;
   IsAtomSiteLoop := False;
   ColTagsCollected := False;
   NCols := 0;
   ColAtom := -1; ColLabel := -1; ColComp := -1; ColSeq := -1;
   ColAsym := -1; ColAlt := -1; ColX := -1; ColY := -1; ColZ := -1;
   ColOcc := -1; ColB := -1; ColEl := -1; ColIns := -1;
   ColModel := -1; ColAuthAsym := -1; ColAuthSeq := -1;
   Continue;
  End;
  If InLoop And (Length(Line) > 0) And (Line[1] = '_') Then Begin
   If NCols < 32 Then Begin
    LoopCols[NCols] := LoStr(Line);
    If StartsWith(LoStr(Line), '_atom_site.') Then Begin
     IsAtomSiteLoop := True;
     If LoStr(Line) = '_atom_site.group_pdb' Then ColAtom := NCols
     Else If LoStr(Line) = '_atom_site.label_atom_id' Then ColLabel := NCols
     Else If LoStr(Line) = '_atom_site.label_comp_id' Then ColComp := NCols
     Else If LoStr(Line) = '_atom_site.label_seq_id' Then ColSeq := NCols
     Else If LoStr(Line) = '_atom_site.label_asym_id' Then ColAsym := NCols
     Else If LoStr(Line) = '_atom_site.label_alt_id' Then ColAlt := NCols
     Else If LoStr(Line) = '_atom_site.cartn_x' Then ColX := NCols
     Else If LoStr(Line) = '_atom_site.cartn_y' Then ColY := NCols
     Else If LoStr(Line) = '_atom_site.cartn_z' Then ColZ := NCols
     Else If LoStr(Line) = '_atom_site.occupancy' Then ColOcc := NCols
     Else If LoStr(Line) = '_atom_site.b_iso_or_equiv' Then ColB := NCols
     Else If LoStr(Line) = '_atom_site.type_symbol' Then ColEl := NCols
     Else If LoStr(Line) = '_atom_site.pdbx_pdb_ins_code' Then ColIns := NCols
     Else If LoStr(Line) = '_atom_site.pdbx_pdb_model_num' Then ColModel := NCols
     Else If LoStr(Line) = '_atom_site.auth_asym_id' Then ColAuthAsym := NCols
     Else If LoStr(Line) = '_atom_site.auth_seq_id' Then ColAuthSeq := NCols;
    End;
    Inc(NCols);
   End;
   Continue;
  End;
  If InLoop And IsAtomSiteLoop Then Begin
   If (Length(Line) > 0) And (Line[1] = '_') Then Begin
    InLoop := False; Continue;
   End;
   If (Length(Line) > 0) And (Line[1] = '#') Then Begin
    InLoop := False; Continue;
   End;
   If Not ColTagsCollected Then ColTagsCollected := True;
   Tokenize(Line, Tokens, NTok, 32);
   If NTok < NCols Then Continue;
   { Filtre : conserver uniquement le 1er modele. }
   If ColModel >= 0 Then Begin
    ModelStr := Tokens[ColModel];
    If Atm.Serial = 0 Then Atm.Serial := 0; { no-op }
    If ModelOk = 1 Then ModelOk := StrToIntDef(ModelStr, 1);
    If StrToIntDef(ModelStr, 1) <> ModelOk Then Continue;
   End;
   { Construction de l'atome }
   Atm.Serial := M.NAtoms + 1;
   If ColLabel >= 0 Then Atm.Name := Copy(Tokens[ColLabel], 1, 6) Else Atm.Name := 'X';
   { Enleve les guillemets eventuels }
   If (Length(Atm.Name) > 0) And (Atm.Name[1] = '"') Then
    Atm.Name := Copy(Atm.Name, 2, Length(Atm.Name) - 2);
   Atm.AltLoc := ' ';
   If (ColAlt >= 0) And (Length(Tokens[ColAlt]) > 0) And
      (Tokens[ColAlt] <> '.') And (Tokens[ColAlt] <> '?') Then
    Atm.AltLoc := Tokens[ColAlt][1];
   If ColComp >= 0 Then Atm.ResName := Copy(Tokens[ColComp], 1, 4)
   Else Atm.ResName := 'UNL';
   Atm.ChainID := ' ';
   If (ColAuthAsym >= 0) And (Length(Tokens[ColAuthAsym]) > 0) Then
    Atm.ChainID := Tokens[ColAuthAsym][1]
   Else If (ColAsym >= 0) And (Length(Tokens[ColAsym]) > 0) Then
    Atm.ChainID := Tokens[ColAsym][1];
   If ColAuthSeq >= 0 Then Atm.ResSeq := StrToIntDef(Tokens[ColAuthSeq], 0)
   Else If ColSeq >= 0 Then Atm.ResSeq := StrToIntDef(Tokens[ColSeq], 0)
   Else Atm.ResSeq := M.NRes + 1;
   Atm.ICode := ' ';
   If (ColIns >= 0) And (Length(Tokens[ColIns]) > 0) And
      (Tokens[ColIns] <> '?') And (Tokens[ColIns] <> '.') Then
    Atm.ICode := Tokens[ColIns][1];
   If ColX < 0 Then Continue;
   If Not StrToReal(Tokens[ColX], X) Then Continue;
   If Not StrToReal(Tokens[ColY], Y) Then Continue;
   If Not StrToReal(Tokens[ColZ], Z) Then Continue;
   Atm.X := X; Atm.Y := Y; Atm.Z := Z;
   Atm.Occup := 1.0; Atm.TempF := 0.0;
   If ColOcc >= 0 Then If StrToReal(Tokens[ColOcc], V) Then Atm.Occup := V;
   If ColB >= 0 Then If StrToReal(Tokens[ColB], V) Then Atm.TempF := V;
   If ColEl >= 0 Then Begin
    Atm.Element := Tokens[ColEl];
    Atm.ElemZ := ElemBySymbol(Atm.Element);
   End
   Else Begin
    Atm.ElemZ := ElemFromAtomName(Atm.Name);
    Atm.Element := ElemDB[Atm.ElemZ].Sym;
   End;
   If (Atm.ElemZ = 20) And (UpStr(Atm.Name) = 'CA') And
      IsAminoAcidName(Atm.ResName) Then Atm.ElemZ := 6;
   Atm.Charge := 0;
   Atm.Kind := AtomKindFor(Atm.ResName,
                           (ColAtom >= 0) And (UpStr(Tokens[ColAtom]) = 'HETATM'));
   Atm.Selected := False;
   Atm.Visible := True;
   Atm.Color := 0;
   If (Atm.AltLoc = ' ') Or (Atm.AltLoc = 'A') Then AddAtom(M, Atm);
  End;
 End;
 Close(F);
 If M.Title = '' Then M.Title := FileName;
 M.Source := 'mmCIF';
 AutoDetectBonds(M);
 ComputeBounds(M);
 BuildResidueIndex(M);
 M.Loaded := M.NAtoms > 0;
 LoadCIF := M.Loaded;
End;

{ --------------------------------------------------------------------- }
{ Dispatcher selon extension                                            }
{ --------------------------------------------------------------------- }

Function ExtractExt(Const FN: String): String;
Var I: Integer; E: String;
Begin
 E := '';
 For I := Length(FN) DownTo 1 Do
  If FN[I] = '.' Then Begin
   E := UpStr(Copy(FN, I + 1, Length(FN) - I));
   Break;
  End
  Else If (FN[I] = '\') Or (FN[I] = '/') Then Break;
 ExtractExt := E;
End;

Function LoadMolecule(Const FileName: String; Var M: TMolecule): Boolean;
Var
 E : String;
 R : Boolean;
Begin
 E := ExtractExt(FileName);
 R := False;
 If (E = 'PDB') Or (E = 'ENT') Then R := LoadPDB(FileName, M)
 Else If E = 'MOL' Then R := LoadMOL(FileName, M)
 Else If E = 'SDF' Then R := LoadSDF(FileName, M)
 Else If E = 'MOL2' Then R := LoadMOL2(FileName, M)
 Else If E = 'CIF' Then R := LoadCIF(FileName, M)
 Else Begin
  { Detection par contenu : ouvrir et lire les premieres lignes. }
  R := LoadPDB(FileName, M);
  If Not R Then R := LoadMOL(FileName, M);
 End;
 LoadMolecule := R;
End;

{ --------------------------------------------------------------------- }
{ Sauvegarde PDB                                                        }
{ --------------------------------------------------------------------- }

Procedure WritePDBField(Var F: Text; Const A: TAtom; Idx: LongInt);
Var
 Line : String;
 S : String;
 I : Integer;
 RecName : String;
Begin
 If A.Kind = ATOMK_LIGAND Then RecName := 'HETATM'
 Else If A.Kind = ATOMK_WATER Then RecName := 'HETATM'
 Else RecName := 'ATOM  ';
 Line := RecName;
 Str(Idx:5, S); Line := Line + S + ' ';
 S := PadR(Copy(A.Name, 1, 4), 4); Line := Line + S + ' ';
 S := PadL(Copy(A.ResName, 1, 3), 3); Line := Line + S + ' ';
 If A.ChainID = #0 Then Line := Line + ' ' Else Line := Line + A.ChainID;
 Str(A.ResSeq:4, S); Line := Line + S;
 Line := Line + '    ';
 Str(A.X:8:3, S); Line := Line + S;
 Str(A.Y:8:3, S); Line := Line + S;
 Str(A.Z:8:3, S); Line := Line + S;
 Str(A.Occup:6:2, S); Line := Line + S;
 Str(A.TempF:6:2, S); Line := Line + S;
 Line := Line + '          ';
 If A.ElemZ > 0 Then Begin
  S := PadL(ElemDB[A.ElemZ].Sym, 2); Line := Line + S;
 End
 Else Line := Line + '  ';
 WriteLn(F, Line);
 I := 0; If I <> 0 Then ;
End;

Function SavePDB(Const FileName: String; Const M: TMolecule): Boolean;
Var
 F : Text;
 I : Integer;
Begin
 SavePDB := False;
 {$I-} Assign(F, FileName); Rewrite(F); {$I+}
 If IOResult <> 0 Then Exit;
 WriteLn(F, 'HEADER    ', M.Title);
 WriteLn(F, 'REMARK   1 Sauvegarde par PYMOL Clone ', PROG_VERSION);
 For I := 0 To M.NAtoms - 1 Do
  WritePDBField(F, M.Atoms[I], I + 1);
 { Liaisons CONECT explicites (utile pour HETATM/ligands) }
 For I := 0 To M.NBonds - 1 Do
  WriteLn(F, 'CONECT', M.Atoms[M.Bonds[I].A].Serial:5,
                       M.Atoms[M.Bonds[I].B].Serial:5);
 WriteLn(F, 'END');
 Close(F);
 SavePDB := True;
End;

{ Sauvegarde au format MOL (V2000) }
Function SaveMOL(Const FileName: String; Const M: TMolecule): Boolean;
Var
 F : Text;
 I : Integer;
 S : String;
Begin
 SaveMOL := False;
 {$I-} Assign(F, FileName); Rewrite(F); {$I+}
 If IOResult <> 0 Then Exit;
 WriteLn(F, M.Title);
 WriteLn(F, '  PyMOL Clone   ');
 WriteLn(F, '');
 Str(M.NAtoms:3, S);
 WriteLn(F, S, M.NBonds:3, '  0  0  0  0  0  0  0  0999 V2000');
 For I := 0 To M.NAtoms - 1 Do Begin
  Str(M.Atoms[I].X:10:4, S); Write(F, S);
  Str(M.Atoms[I].Y:10:4, S); Write(F, S);
  Str(M.Atoms[I].Z:10:4, S); Write(F, S);
  Write(F, ' ');
  If M.Atoms[I].ElemZ > 0 Then
   Write(F, PadR(ElemDB[M.Atoms[I].ElemZ].Sym, 3))
  Else Write(F, 'X  ');
  WriteLn(F, ' 0  0  0  0  0  0  0  0  0  0  0  0');
 End;
 For I := 0 To M.NBonds - 1 Do Begin
  Str((M.Bonds[I].A + 1):3, S); Write(F, S);
  Str((M.Bonds[I].B + 1):3, S); Write(F, S);
  Str(M.Bonds[I].Order:3, S); Write(F, S);
  WriteLn(F, '  0  0  0  0');
 End;
 WriteLn(F, 'M  END');
 Close(F);
 SaveMOL := True;
End;

{ --------------------------------------------------------------------- }
{ Console PyMOL>                                                        }
{ --------------------------------------------------------------------- }

Procedure ConsolePush(Const S: String; Kind: Byte);
Var Idx, I: Integer;
Begin
 If ConsoleCount < MAX_CONSOLE_LINES Then Begin
  ConsoleLines[ConsoleCount].S := S;
  ConsoleLines[ConsoleCount].Kind := Kind;
  Inc(ConsoleCount);
 End
 Else Begin
  For I := 0 To MAX_CONSOLE_LINES - 2 Do
   ConsoleLines[I] := ConsoleLines[I + 1];
  Idx := MAX_CONSOLE_LINES - 1;
  ConsoleLines[Idx].S := S;
  ConsoleLines[Idx].Kind := Kind;
 End;
 ConsoleView := 0;
End;

Procedure ConsoleSay(Const S: String);
Begin ConsolePush(S, 0); End;

Procedure ConsoleErr(Const S: String);
Begin ConsolePush('Erreur: ' + S, 2); End;

Procedure ConsoleEcho(Const S: String);
Begin ConsolePush(PROG_PROMPT + S, 1); End;

{ --------------------------------------------------------------------- }
{ Selections (langage tres reduit)                                      }
{ --------------------------------------------------------------------- }

{ Marque les atomes selectionnes selon une expression. Renvoie le nombre
  d'atomes selectionnes. Reconnait :
    all                : tous les atomes
    none               : aucun
    hetatm | hetero    : atomes HETATM
    water | waters     : eau
    protein            : residus d'acides amines
    nucleic            : acides nucleiques
    hydrogens          : H
    backbone           : N, CA, C, O des AA
    sidechain          : non-backbone des AA
    elem X             : par symbole
    chain X            : par chaine
    resi N             : par numero
    resn XYZ           : par nom de residu
    name XYZ           : par nom d'atome
    visible            : atome visible
    sele               : atome actuellement Selected
    not <expr>         : negation simple sur token initial
}
Function SelectAtoms(Var M: TMolecule; Const Expr: String): Integer;
Var
 Tokens : Array[0..15] Of String;
 NTok : Integer;
 I, K, Count : Integer;
 Negate : Boolean;
 Crit : String;
 Arg : String;
 SymZ : Byte;
 Sn : LongInt;
 R : Boolean;
 Up : String;
Begin
 Tokenize(Expr, Tokens, NTok, 16);
 Negate := False;
 K := 0;
 If (NTok > 0) And (UpStr(Tokens[0]) = 'NOT') Then Begin
  Negate := True; K := 1;
 End;
 If NTok = K Then Begin
  For I := 0 To M.NAtoms - 1 Do M.Atoms[I].Selected := True;
  SelectAtoms := M.NAtoms;
  Exit;
 End;
 Crit := UpStr(Tokens[K]);
 Arg := '';
 If NTok > K + 1 Then Arg := Tokens[K + 1];
 Count := 0;
 For I := 0 To M.NAtoms - 1 Do Begin
  R := False;
  If Crit = 'ALL' Then R := True
  Else If Crit = 'NONE' Then R := False
  Else If (Crit = 'HETATM') Or (Crit = 'HETERO') Then
   R := M.Atoms[I].Kind = ATOMK_LIGAND
  Else If (Crit = 'WATER') Or (Crit = 'WATERS') Or (Crit = 'RESN') Then Begin
   If Crit = 'RESN' Then R := UpStr(TrimStr(M.Atoms[I].ResName)) = UpStr(Arg)
   Else R := M.Atoms[I].Kind = ATOMK_WATER;
  End
  Else If Crit = 'PROTEIN' Then R := M.Atoms[I].Kind = ATOMK_PROT
  Else If Crit = 'NUCLEIC' Then R := M.Atoms[I].Kind = ATOMK_NUCL
  Else If (Crit = 'HYDROGENS') Or (Crit = 'HYDRO') Or (Crit = 'H') Then
   R := M.Atoms[I].ElemZ = 1
  Else If Crit = 'BACKBONE' Then Begin
   Up := UpStr(TrimStr(M.Atoms[I].Name));
   R := M.Atoms[I].Kind = ATOMK_PROT;
   R := R And ((Up = 'N') Or (Up = 'CA') Or (Up = 'C') Or (Up = 'O'));
  End
  Else If Crit = 'SIDECHAIN' Then Begin
   Up := UpStr(TrimStr(M.Atoms[I].Name));
   R := M.Atoms[I].Kind = ATOMK_PROT;
   R := R And Not ((Up = 'N') Or (Up = 'CA') Or (Up = 'C') Or (Up = 'O'));
  End
  Else If (Crit = 'ELEM') Or (Crit = 'ELEMENT') Then Begin
   SymZ := ElemBySymbol(Arg);
   R := (SymZ <> 0) And (M.Atoms[I].ElemZ = SymZ);
  End
  Else If Crit = 'CHAIN' Then Begin
   If Length(Arg) > 0 Then R := UpCh(M.Atoms[I].ChainID) = UpCh(Arg[1]);
  End
  Else If Crit = 'RESI' Then Begin
   Sn := StrToIntDef(Arg, -32768);
   R := M.Atoms[I].ResSeq = Sn;
  End
  Else If Crit = 'NAME' Then
   R := UpStr(TrimStr(M.Atoms[I].Name)) = UpStr(Arg)
  Else If Crit = 'VISIBLE' Then R := M.Atoms[I].Visible
  Else If (Crit = 'SELE') Or (Crit = 'SELECTED') Then
   R := M.Atoms[I].Selected;
  If Negate Then R := Not R;
  M.Atoms[I].Selected := R;
  If R Then Inc(Count);
 End;
 SelectAtoms := Count;
End;

{ --------------------------------------------------------------------- }
{ Application de couleurs                                               }
{ --------------------------------------------------------------------- }

Function ParseColorName(Const N: String): Byte;
Var U: String;
Begin
 U := UpStr(N);
 If (U = 'WHITE') Or (U = 'BLANC') Then ParseColorName := White
 Else If U = 'BLACK' Then ParseColorName := Black
 Else If (U = 'RED') Or (U = 'ROUGE') Then ParseColorName := LightRed
 Else If (U = 'GREEN') Or (U = 'VERT') Then ParseColorName := LightGreen
 Else If (U = 'BLUE') Or (U = 'BLEU') Then ParseColorName := LightBlue
 Else If (U = 'YELLOW') Or (U = 'JAUNE') Then ParseColorName := Yellow
 Else If U = 'CYAN' Then ParseColorName := LightCyan
 Else If (U = 'MAGENTA') Or (U = 'PINK') Then ParseColorName := LightMagenta
 Else If (U = 'BROWN') Or (U = 'BRUN') Then ParseColorName := Brown
 Else If (U = 'GRAY') Or (U = 'GREY') Then ParseColorName := LightGray
 Else If U = 'DARKGRAY' Then ParseColorName := DarkGray
 Else If U = 'ORANGE' Then ParseColorName := LightRed
 Else If U = 'PURPLE' Then ParseColorName := LightMagenta
 Else ParseColorName := 0;
End;

Procedure ColorSelected(Var M: TMolecule; ColIdx: Byte);
Var I, Cnt: Integer;
Begin
 Cnt := 0;
 For I := 0 To M.NAtoms - 1 Do
  If M.Atoms[I].Selected Then Begin
   M.Atoms[I].Color := ColIdx;
   Inc(Cnt);
  End;
 ConsoleSay(' Coloration appliquee a ' + IntToStr(Cnt) + ' atomes.');
End;

Procedure ShowSelected(Var M: TMolecule; Vis: Boolean);
Var I, Cnt: Integer;
Begin
 Cnt := 0;
 For I := 0 To M.NAtoms - 1 Do
  If M.Atoms[I].Selected Then Begin
   M.Atoms[I].Visible := Vis;
   Inc(Cnt);
  End;
 If Vis Then ConsoleSay(' Affichage active pour ' + IntToStr(Cnt) + ' atomes.')
 Else ConsoleSay(' Masque ' + IntToStr(Cnt) + ' atomes.');
End;

{ --------------------------------------------------------------------- }
{ Math 3D (rotations, projection)                                       }
{ --------------------------------------------------------------------- }

Procedure ApplyRotation(Var X, Y, Z: Real;
                        AxF, AyF, AzF: Real);
{ Applique en place Rx, Ry, Rz dans l'ordre Z-Y-X (typique PyMOL like). }
Var
 C, S, T : Real;
Begin
 { Rx }
 C := Cos(AxF * PI_VAL / 180.0);
 S := Sin(AxF * PI_VAL / 180.0);
 T := Y * C - Z * S;
 Z := Y * S + Z * C;
 Y := T;
 { Ry }
 C := Cos(AyF * PI_VAL / 180.0);
 S := Sin(AyF * PI_VAL / 180.0);
 T := Z * C - X * S;
 X := Z * S + X * C;
 Z := T;
 { Rz }
 C := Cos(AzF * PI_VAL / 180.0);
 S := Sin(AzF * PI_VAL / 180.0);
 T := X * C - Y * S;
 Y := X * S + Y * C;
 X := T;
End;

{ Projette les coordonnees 3D (en A) vers les pixels ecran.
  Retourne aussi Zp pour le tri Z (depth) }
Procedure ProjectAtom(Const M: TMolecule; AtIdx: Integer;
                      Var Sx, Sy: Integer; Var Zp: Real);
Var
 X, Y, Z, Persp : Real;
Begin
 X := M.Atoms[AtIdx].X - M.CenterX;
 Y := M.Atoms[AtIdx].Y - M.CenterY;
 Z := M.Atoms[AtIdx].Z - M.CenterZ;
 ApplyRotation(X, Y, Z, RotX, RotY, RotZ);
 Zp := Z;
 If ViewDistance > 0.0 Then Begin
  Persp := ViewDistance / (ViewDistance - Z);
  If Persp < 0.1 Then Persp := 0.1;
  If Persp > 8.0 Then Persp := 8.0;
 End
 Else Persp := 1.0;
 Sx := (ScreenWidth Div 2) + ViewPanX + Round(X * ViewScale * Persp);
 Sy := ((CanvasTop + CanvasBottom) Div 2) + ViewPanY - Round(Y * ViewScale * Persp);
End;

{ --------------------------------------------------------------------- }
{ Parseur de scripts .PML (commandes PyMOL)                             }
{ --------------------------------------------------------------------- }

Procedure ExecPyMOLCommand(Const Cmd: String); Forward;

Function LoadPML(Const FileName: String): Boolean;
Var
 F : Text;
 Line : String;
Begin
 LoadPML := False;
 {$I-} Assign(F, FileName); Reset(F); {$I+}
 If IOResult <> 0 Then Exit;
 ConsoleSay(' Execution du script ' + FileName);
 While Not Eof(F) Do Begin
  ReadLn(F, Line);
  Line := TrimStr(Line);
  If Line = '' Then Continue;
  If Line[1] = '#' Then Continue;
  ConsoleEcho(Line);
  ExecPyMOLCommand(Line);
 End;
 Close(F);
 LoadPML := True;
End;

{ --------------------------------------------------------------------- }
{ Cadrage automatique de la vue                                         }
{ --------------------------------------------------------------------- }

Procedure ZoomToFit;
Var
 Avail : Integer;
Begin
 Avail := ScreenHeight - MenuBarHeight - STATUSBAR_HEIGHT;
 If ConsoleVisible Then Dec(Avail, CONSOLE_HEIGHT);
 If Avail < 100 Then Avail := 100;
 If Mol.Radius > 0.001 Then
  ViewScale := (Avail Div 2) / (Mol.Radius * 1.4)
 Else
  ViewScale := 30.0;
 If ViewScale < 0.5 Then ViewScale := 0.5;
 ViewPanX := 0;
 ViewPanY := 0;
 RotX := 0; RotY := 0; RotZ := 0;
 ViewDistance := Mol.Radius * 6.0;
 If ViewDistance < 10.0 Then ViewDistance := 30.0;
End;

{ --------------------------------------------------------------------- }
{ Couleurs CPK et utilitaires                                           }
{ --------------------------------------------------------------------- }

Function AtomColor(Const M: TMolecule; I: Integer): Byte;
Begin
 If M.Atoms[I].Color <> 0 Then AtomColor := M.Atoms[I].Color
 Else If M.Atoms[I].ElemZ > 0 Then AtomColor := ElemDB[M.Atoms[I].ElemZ].Color
 Else AtomColor := LightGray;
End;

Function AtomVdWPx(Const M: TMolecule; I: Integer): Integer;
Var R: Real;
Begin
 If M.Atoms[I].ElemZ > 0 Then R := ElemDB[M.Atoms[I].ElemZ].VdW
 Else R := 1.5;
 AtomVdWPx := Round(R * ViewScale * 0.5);
End;

Function AtomCovPx(Const M: TMolecule; I: Integer): Integer;
Var R: Real;
Begin
 If M.Atoms[I].ElemZ > 0 Then R := ElemDB[M.Atoms[I].ElemZ].Radius
 Else R := 0.7;
 AtomCovPx := Round(R * ViewScale * 0.3);
End;

{ --------------------------------------------------------------------- }
{ Sphere remplie (approximation par cercles concentriques + ombrage)    }
{ --------------------------------------------------------------------- }

Procedure DrawSphere(Cx, Cy, R: Integer; Col: Byte);
Var
 I, RR : Integer;
 ShadeCol : Byte;
Begin
 If R < 1 Then Begin
  PutPixel(Cx, Cy, Col); Exit;
 End;
 SetColor(Col);
 SetFillStyle(SolidFill, Col);
 FillEllipse(Cx, Cy, R, R);
 { Reflexion speculaire : petit cercle plus clair en haut a gauche }
 If R >= 4 Then Begin
  Case Col Of
   Black     : ShadeCol := DarkGray;
   Blue      : ShadeCol := LightBlue;
   Green     : ShadeCol := LightGreen;
   Cyan      : ShadeCol := LightCyan;
   Red       : ShadeCol := LightRed;
   Magenta   : ShadeCol := LightMagenta;
   Brown     : ShadeCol := Yellow;
   LightGray : ShadeCol := White;
   DarkGray  : ShadeCol := LightGray;
  Else
   ShadeCol := White;
  End;
  SetFillStyle(SolidFill, ShadeCol);
  RR := R Div 3;
  If RR < 1 Then RR := 1;
  FillEllipse(Cx - R Div 3, Cy - R Div 3, RR, RR);
 End;
 SetColor(Black);
 { Contour }
 SetLineStyle(SolidLn, 0, NormWidth);
 Circle(Cx, Cy, R);
 I := 0; If I <> 0 Then ;
End;

{ --------------------------------------------------------------------- }
{ Tri Z (Painter's algorithm) sur les indices                           }
{ --------------------------------------------------------------------- }

Var
 ZOrder : Array[0..MAX_ATOMS-1] Of Integer;
 ZDepth : Array[0..MAX_ATOMS-1] Of Real;
 ScreenX, ScreenY : Array[0..MAX_ATOMS-1] Of Integer;
 ScreenDepth : Array[0..MAX_ATOMS-1] Of Real;

Procedure ComputeProjection;
Var I, Sx, Sy: Integer; Zp: Real;
Begin
 For I := 0 To Mol.NAtoms - 1 Do Begin
  ProjectAtom(Mol, I, Sx, Sy, Zp);
  ScreenX[I] := Sx;
  ScreenY[I] := Sy;
  ScreenDepth[I] := Zp;
 End;
End;

Procedure SortByDepth;
Var
 I, J, T : Integer;
 Tmp : Real;
 Left, Right : Integer;
 Stack : Array[0..63] Of Integer;
 SP : Integer;
 Pivot : Real;
Begin
 { Initialisation index identite }
 For I := 0 To Mol.NAtoms - 1 Do Begin
  ZOrder[I] := I;
  ZDepth[I] := ScreenDepth[I];
 End;
 If Mol.NAtoms < 2 Then Exit;
 { Quicksort iteratif }
 SP := 0;
 Stack[SP] := 0; Inc(SP);
 Stack[SP] := Mol.NAtoms - 1; Inc(SP);
 While SP > 0 Do Begin
  Dec(SP); Right := Stack[SP];
  Dec(SP); Left := Stack[SP];
  If Left >= Right Then Continue;
  Pivot := ZDepth[(Left + Right) Div 2];
  I := Left; J := Right;
  While I <= J Do Begin
   While ZDepth[I] < Pivot Do Inc(I);
   While ZDepth[J] > Pivot Do Dec(J);
   If I <= J Then Begin
    Tmp := ZDepth[I]; ZDepth[I] := ZDepth[J]; ZDepth[J] := Tmp;
    T := ZOrder[I]; ZOrder[I] := ZOrder[J]; ZOrder[J] := T;
    Inc(I); Dec(J);
   End;
  End;
  If Left < J Then Begin
   If SP + 2 <= 64 Then Begin
    Stack[SP] := Left; Inc(SP);
    Stack[SP] := J; Inc(SP);
   End;
  End;
  If I < Right Then Begin
   If SP + 2 <= 64 Then Begin
    Stack[SP] := I; Inc(SP);
    Stack[SP] := Right; Inc(SP);
   End;
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Trace 3D : lignes, batonnets, spheres, balls-and-sticks, ruban        }
{ --------------------------------------------------------------------- }

Procedure DrawAxes;
Const Len = 30;
Var
 Cx, Cy : Integer;
 X0, Y0, Z0 : Real;
 X1, Y1, Z1 : Real;
 Sx, Sy : Integer;
Begin
 Cx := 30;
 Cy := CanvasBottom - 40;
 SetTextStyle(DefaultFont, HorizDir, 1);
 { Axe X (rouge) }
 X0 := 0; Y0 := 0; Z0 := 0;
 X1 := 1; Y1 := 0; Z1 := 0;
 ApplyRotation(X0, Y0, Z0, RotX, RotY, RotZ);
 ApplyRotation(X1, Y1, Z1, RotX, RotY, RotZ);
 Sx := Cx + Round(X1 * Len);
 Sy := Cy - Round(Y1 * Len);
 SetColor(LightRed);
 Line(Cx, Cy, Sx, Sy);
 OutTextXY(Sx + 2, Sy - 4, 'X');
 { Axe Y (vert) }
 X1 := 0; Y1 := 1; Z1 := 0;
 ApplyRotation(X1, Y1, Z1, RotX, RotY, RotZ);
 Sx := Cx + Round(X1 * Len);
 Sy := Cy - Round(Y1 * Len);
 SetColor(LightGreen);
 Line(Cx, Cy, Sx, Sy);
 OutTextXY(Sx + 2, Sy - 4, 'Y');
 { Axe Z (bleu) }
 X1 := 0; Y1 := 0; Z1 := 1;
 ApplyRotation(X1, Y1, Z1, RotX, RotY, RotZ);
 Sx := Cx + Round(X1 * Len);
 Sy := Cy - Round(Y1 * Len);
 SetColor(LightBlue);
 Line(Cx, Cy, Sx, Sy);
 OutTextXY(Sx + 2, Sy - 4, 'Z');
End;

Procedure DrawBBox;
Var
 P : Array[0..7] Of Record X, Y, Z: Real; Sx, Sy: Integer; End;
 I : Integer;
 Bx0, Bx1, By0, By1, Bz0, Bz1 : Real;
 Persp : Real;
 Procedure ProjectCorner(K: Integer);
 Begin
  ApplyRotation(P[K].X, P[K].Y, P[K].Z, RotX, RotY, RotZ);
  If ViewDistance > 0.0 Then
   Persp := ViewDistance / (ViewDistance - P[K].Z)
  Else Persp := 1.0;
  P[K].Sx := (ScreenWidth Div 2) + ViewPanX + Round(P[K].X * ViewScale * Persp);
  P[K].Sy := ((CanvasTop + CanvasBottom) Div 2) + ViewPanY - Round(P[K].Y * ViewScale * Persp);
 End;
Begin
 Bx0 := Mol.MinX - Mol.CenterX; Bx1 := Mol.MaxX - Mol.CenterX;
 By0 := Mol.MinY - Mol.CenterY; By1 := Mol.MaxY - Mol.CenterY;
 Bz0 := Mol.MinZ - Mol.CenterZ; Bz1 := Mol.MaxZ - Mol.CenterZ;
 P[0].X := Bx0; P[0].Y := By0; P[0].Z := Bz0;
 P[1].X := Bx1; P[1].Y := By0; P[1].Z := Bz0;
 P[2].X := Bx1; P[2].Y := By1; P[2].Z := Bz0;
 P[3].X := Bx0; P[3].Y := By1; P[3].Z := Bz0;
 P[4].X := Bx0; P[4].Y := By0; P[4].Z := Bz1;
 P[5].X := Bx1; P[5].Y := By0; P[5].Z := Bz1;
 P[6].X := Bx1; P[6].Y := By1; P[6].Z := Bz1;
 P[7].X := Bx0; P[7].Y := By1; P[7].Z := Bz1;
 For I := 0 To 7 Do ProjectCorner(I);
 SetColor(DarkGray);
 SetLineStyle(DottedLn, 0, NormWidth);
 For I := 0 To 3 Do Begin
  Line(P[I].Sx, P[I].Sy, P[(I+1) Mod 4].Sx, P[(I+1) Mod 4].Sy);
  Line(P[I+4].Sx, P[I+4].Sy, P[((I+1) Mod 4)+4].Sx, P[((I+1) Mod 4)+4].Sy);
  Line(P[I].Sx, P[I].Sy, P[I+4].Sx, P[I+4].Sy);
 End;
 SetLineStyle(SolidLn, 0, NormWidth);
End;

Procedure DrawGridXY;
Var
 X, Y, Step : Integer;
Begin
 Step := 32;
 SetColor(DarkGray);
 X := 0;
 While X < ScreenWidth Do Begin
  Line(X, CanvasTop + 1, X, CanvasBottom - 1);
  Inc(X, Step);
 End;
 Y := CanvasTop + 1;
 While Y < CanvasBottom Do Begin
  Line(0, Y, ScreenWidth - 1, Y);
  Inc(Y, Step);
 End;
End;

Procedure DrawBondAtomic(I, J: Integer; Col1, Col2: Byte; Thick: Integer);
Var
 Mx, My : Integer;
 K : Integer;
Begin
 Mx := (ScreenX[I] + ScreenX[J]) Div 2;
 My := (ScreenY[I] + ScreenY[J]) Div 2;
 If Thick <= 1 Then Begin
  SetColor(Col1);
  Line(ScreenX[I], ScreenY[I], Mx, My);
  SetColor(Col2);
  Line(Mx, My, ScreenX[J], ScreenY[J]);
 End
 Else Begin
  For K := -Thick Div 2 To Thick Div 2 Do Begin
   SetColor(Col1);
   Line(ScreenX[I], ScreenY[I] + K, Mx, My + K);
   SetColor(Col2);
   Line(Mx, My + K, ScreenX[J], ScreenY[J] + K);
  End;
 End;
End;

Procedure DrawMoleculeLines;
Var B, I, J: Integer;
Begin
 For B := 0 To Mol.NBonds - 1 Do Begin
  If Not Mol.Bonds[B].Visible Then Continue;
  I := Mol.Bonds[B].A;
  J := Mol.Bonds[B].B;
  If Not Mol.Atoms[I].Visible Then Continue;
  If Not Mol.Atoms[J].Visible Then Continue;
  If ((ViewFlags And SHOW_HYDROGENS) = 0) And
     ((Mol.Atoms[I].ElemZ = 1) Or (Mol.Atoms[J].ElemZ = 1)) Then Continue;
  DrawBondAtomic(I, J, AtomColor(Mol, I), AtomColor(Mol, J), 1);
 End;
 { Atomes seuls / non lies : point }
 For I := 0 To Mol.NAtoms - 1 Do
  If Mol.Atoms[I].Visible Then Begin
   If ((ViewFlags And SHOW_HYDROGENS) = 0) And (Mol.Atoms[I].ElemZ = 1) Then Continue;
   PutPixel(ScreenX[I], ScreenY[I], AtomColor(Mol, I));
  End;
End;

Procedure DrawMoleculeSticks;
Var B, I, J: Integer;
Begin
 For B := 0 To Mol.NBonds - 1 Do Begin
  If Not Mol.Bonds[B].Visible Then Continue;
  I := Mol.Bonds[B].A;
  J := Mol.Bonds[B].B;
  If Not Mol.Atoms[I].Visible Then Continue;
  If Not Mol.Atoms[J].Visible Then Continue;
  If ((ViewFlags And SHOW_HYDROGENS) = 0) And
     ((Mol.Atoms[I].ElemZ = 1) Or (Mol.Atoms[J].ElemZ = 1)) Then Continue;
  DrawBondAtomic(I, J, AtomColor(Mol, I), AtomColor(Mol, J), 3);
  { Double / triple : trait parallele }
  If Mol.Bonds[B].Order >= 2 Then
   DrawBondAtomic(I, J, AtomColor(Mol, I), AtomColor(Mol, J), 1);
 End;
End;

Procedure DrawMoleculeSpheres;
Var K, I, R: Integer;
Begin
 For K := 0 To Mol.NAtoms - 1 Do Begin
  I := ZOrder[K];
  If Not Mol.Atoms[I].Visible Then Continue;
  If ((ViewFlags And SHOW_HYDROGENS) = 0) And (Mol.Atoms[I].ElemZ = 1) Then Continue;
  R := AtomVdWPx(Mol, I);
  DrawSphere(ScreenX[I], ScreenY[I], R, AtomColor(Mol, I));
 End;
End;

Procedure DrawMoleculeBallStick;
Var K, I, R, B: Integer;
Begin
 { batons d'abord }
 For B := 0 To Mol.NBonds - 1 Do Begin
  If Not Mol.Bonds[B].Visible Then Continue;
  I := Mol.Bonds[B].A;
  If Not Mol.Atoms[I].Visible Then Continue;
  If Not Mol.Atoms[Mol.Bonds[B].B].Visible Then Continue;
  If ((ViewFlags And SHOW_HYDROGENS) = 0) And
     ((Mol.Atoms[I].ElemZ = 1) Or (Mol.Atoms[Mol.Bonds[B].B].ElemZ = 1)) Then Continue;
  DrawBondAtomic(I, Mol.Bonds[B].B,
                 AtomColor(Mol, I),
                 AtomColor(Mol, Mol.Bonds[B].B), 3);
 End;
 { puis spheres (couvre les jonctions) }
 For K := 0 To Mol.NAtoms - 1 Do Begin
  I := ZOrder[K];
  If Not Mol.Atoms[I].Visible Then Continue;
  If ((ViewFlags And SHOW_HYDROGENS) = 0) And (Mol.Atoms[I].ElemZ = 1) Then Continue;
  R := AtomCovPx(Mol, I);
  If R < 2 Then R := 2;
  DrawSphere(ScreenX[I], ScreenY[I], R, AtomColor(Mol, I));
 End;
End;

Procedure DrawMoleculeRibbon;
Var
 R, FirstCA, PrevCA, CurCA : Integer;
 I : Integer;
 SegCol, EdgeCol : Byte;
 SS : Char;
Begin
 PrevCA := -1;
 SetColor(LightCyan);
 For R := 0 To Mol.NRes - 1 Do Begin
  If Not Mol.Resi[R].IsAA Then Begin PrevCA := -1; Continue; End;
  CurCA := -1;
  For I := Mol.Resi[R].FirstAt To Mol.Resi[R].LastAt Do
   If UpStr(TrimStr(Mol.Atoms[I].Name)) = 'CA' Then Begin
    CurCA := I; Break;
   End;
  If CurCA < 0 Then Begin PrevCA := -1; Continue; End;
  If PrevCA >= 0 Then Begin
   { Couleur du segment selon le SSE si cartoon_discrete_colors = on }
   If SetCartoonDiscreteColors Then Begin
    SS := Mol.Resi[R].SSType;
    Case SS Of
     'H': SegCol := SetCartoonHelixColor;
     'E': SegCol := SetCartoonSheetColor;
    Else   SegCol := SetCartoonLoopColor;
    End;
    EdgeCol := SegCol;
   End Else Begin
    SegCol  := LightCyan;
    EdgeCol := Cyan;
   End;
   SetColor(SegCol);
   Line(ScreenX[PrevCA], ScreenY[PrevCA] - 1,
        ScreenX[CurCA],  ScreenY[CurCA] - 1);
   Line(ScreenX[PrevCA], ScreenY[PrevCA],
        ScreenX[CurCA],  ScreenY[CurCA]);
   Line(ScreenX[PrevCA], ScreenY[PrevCA] + 1,
        ScreenX[CurCA],  ScreenY[CurCA] + 1);
   SetColor(EdgeCol);
   Line(ScreenX[PrevCA], ScreenY[PrevCA] + 2,
        ScreenX[CurCA],  ScreenY[CurCA] + 2);
  End;
  PrevCA := CurCA;
 End;
 FirstCA := 0; If FirstCA <> 0 Then ;
End;

Procedure DrawMoleculeDots;
Var I: Integer;
Begin
 For I := 0 To Mol.NAtoms - 1 Do Begin
  If Not Mol.Atoms[I].Visible Then Continue;
  If ((ViewFlags And SHOW_HYDROGENS) = 0) And (Mol.Atoms[I].ElemZ = 1) Then Continue;
  PutPixel(ScreenX[I], ScreenY[I], AtomColor(Mol, I));
  PutPixel(ScreenX[I] + 1, ScreenY[I], AtomColor(Mol, I));
  PutPixel(ScreenX[I], ScreenY[I] + 1, AtomColor(Mol, I));
 End;
End;

Procedure DrawSelectionMarkers;
Var I: Integer;
Begin
 For I := 0 To Mol.NAtoms - 1 Do
  If Mol.Atoms[I].Selected And Mol.Atoms[I].Visible Then Begin
   SetColor(Yellow);
   Circle(ScreenX[I], ScreenY[I], 5);
  End;
End;

Procedure DrawAtomLabels;
Var I: Integer; S: String;
Begin
 SetColor(White);
 For I := 0 To Mol.NAtoms - 1 Do
  If Mol.Atoms[I].Selected And Mol.Atoms[I].Visible Then Begin
   S := Mol.Atoms[I].Name;
   OutTextXY(ScreenX[I] + 6, ScreenY[I] - 4, S);
  End;
End;

Procedure DrawCanvas;
Begin
 SetFillStyle(SolidFill, BgColor);
 Bar(0, CanvasTop, ScreenWidth - 1, CanvasBottom);
 If Not Mol.Loaded Then Exit;
 ComputeProjection;
 SortByDepth;
 If (ViewFlags And SHOW_GRID) <> 0 Then DrawGridXY;
 If (ViewFlags And SHOW_BBOX) <> 0 Then DrawBBox;
 Case RenderMode Of
  RENDER_LINES     : DrawMoleculeLines;
  RENDER_STICKS    : DrawMoleculeSticks;
  RENDER_SPHERES   : DrawMoleculeSpheres;
  RENDER_BALLSTICK : DrawMoleculeBallStick;
  RENDER_RIBBON    : Begin DrawMoleculeRibbon; DrawMoleculeSticks; End;
  RENDER_DOTS      : DrawMoleculeDots;
 End;
 DrawSelectionMarkers;
 If (ViewFlags And SHOW_LABELS) <> 0 Then DrawAtomLabels;
 If (ViewFlags And SHOW_AXES) <> 0 Then DrawAxes;
End;

Procedure DrawStatusBar;
Var
 Y : Integer;
 S : String;
Begin
 Y := ScreenHeight - STATUSBAR_HEIGHT;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, Y, ScreenWidth - 1, ScreenHeight - 1);
 SetColor(DarkGray);
 Line(0, Y, ScreenWidth - 1, Y);
 SetColor(Black);
 If Mol.Loaded Then Begin
  S := Mol.Title;
  If Length(S) > 38 Then S := Copy(S, 1, 38);
  OutTextXY(5, Y + 5, S);
  Str(Mol.NAtoms, S);
  OutTextXY(290, Y + 5, S + ' atomes');
  Str(Mol.NBonds, S);
  OutTextXY(380, Y + 5, S + ' liaisons');
  Str(Mol.NRes, S);
  OutTextXY(470, Y + 5, S + ' residus');
 End
 Else
  OutTextXY(5, Y + 5, 'Aucune molecule chargee.');
 Case RenderMode Of
  RENDER_LINES     : S := '[Lignes]';
  RENDER_STICKS    : S := '[Batonnets]';
  RENDER_SPHERES   : S := '[CPK]';
  RENDER_BALLSTICK : S := '[Boules+batonnets]';
  RENDER_RIBBON    : S := '[Ruban]';
  RENDER_DOTS      : S := '[Points]';
 Else
  S := '[?]';
 End;
 OutTextXY(ScreenWidth - 8 * Length(S) - 6, Y + 5, S);
End;

Procedure DrawConsoleArea;
Var
 Y, Yt, I, Line0, K : Integer;
 S : String;
 Col : Byte;
Begin
 If Not ConsoleVisible Then Exit;
 Y := ScreenHeight - STATUSBAR_HEIGHT - CONSOLE_HEIGHT;
 SetFillStyle(SolidFill, Black);
 Bar(0, Y, ScreenWidth - 1, Y + CONSOLE_HEIGHT - 1);
 SetColor(DarkGray);
 Line(0, Y, ScreenWidth - 1, Y);
 { Affiche CONSOLE_VISIBLE_LN dernieres lignes }
 Line0 := ConsoleCount - CONSOLE_VISIBLE_LN - ConsoleView;
 If Line0 < 0 Then Line0 := 0;
 Yt := Y + 4;
 For K := 0 To CONSOLE_VISIBLE_LN - 1 Do Begin
  I := Line0 + K;
  If I >= ConsoleCount Then Break;
  Case ConsoleLines[I].Kind Of
   1 : Col := LightGreen;
   2 : Col := LightRed;
   3 : Col := Yellow;
  Else
   Col := LightGray;
  End;
  SetColor(Col);
  S := ConsoleLines[I].S;
  If Length(S) > 78 Then S := Copy(S, 1, 78);
  OutTextXY(4, Yt, S);
  Inc(Yt, 10);
 End;
 { Ligne d'invite (saisie) }
 SetColor(LightGreen);
 S := PROG_PROMPT + ConsoleInput;
 If Length(S) > 78 Then S := Copy(S, Length(S) - 78, 78);
 OutTextXY(4, Y + CONSOLE_HEIGHT - 14, S + '_');
End;

{ --------------------------------------------------------------------- }
{ Menus deroulants                                                      }
{ --------------------------------------------------------------------- }

Procedure ZeroItem(Var M: TMenu);
Var I: Integer;
Begin
 For I := 0 To MAX_MENU_ITEMS - 1 Do Begin
  M.Items[I].Caption   := '';
  M.Items[I].Shortcut  := '';
  M.Items[I].Enabled   := True;
  M.Items[I].Checked   := False;
  M.Items[I].Separator := False;
 End;
End;

Procedure MI(Var Mn: TMenu; Idx: Integer; Const Cap, Sc: String);
Begin
 Mn.Items[Idx].Caption := Cap;
 Mn.Items[Idx].Shortcut := Sc;
End;

Procedure MS(Var Mn: TMenu; Idx: Integer);
Begin Mn.Items[Idx].Separator := True; End;

Procedure InitMenus;
Begin
 { Fichier }
 ZeroItem(Menus[MENU_FILE]);
 Menus[MENU_FILE].Title := 'Fichier';
 Menus[MENU_FILE].ItemCount := 14;
 MI(Menus[MENU_FILE], 0, 'Ouvrir...',            'Ctrl+O');
 MI(Menus[MENU_FILE], 1, 'Recharger',            'F3');
 MI(Menus[MENU_FILE], 2, 'Executer script (.PML)...', '');
 MS(Menus[MENU_FILE], 3);
 MI(Menus[MENU_FILE], 4, 'Enregistrer .PDB...',  'Ctrl+S');
 MI(Menus[MENU_FILE], 5, 'Enregistrer .MOL...',  '');
 MI(Menus[MENU_FILE], 6, 'Enregistrer session', '');
 MI(Menus[MENU_FILE], 7, 'Enregistrer image',   '');
 MS(Menus[MENU_FILE], 8);
 MI(Menus[MENU_FILE], 9, 'Importer fragment...',  '');
 MI(Menus[MENU_FILE], 10,'Exporter image PPM...', '');
 MS(Menus[MENU_FILE], 11);
 MI(Menus[MENU_FILE], 12,'Reinitialiser',        'Ctrl+R');
 MI(Menus[MENU_FILE], 13,'Quitter',              'Alt+F4');

 { Edition }
 ZeroItem(Menus[MENU_EDIT]);
 Menus[MENU_EDIT].Title := 'Edition';
 Menus[MENU_EDIT].ItemCount := 8;
 MI(Menus[MENU_EDIT], 0, 'Annuler',            'Ctrl+Z');
 MI(Menus[MENU_EDIT], 1, 'Refaire',            'Ctrl+Y');
 MS(Menus[MENU_EDIT], 2);
 MI(Menus[MENU_EDIT], 3, 'Copier',             'Ctrl+C');
 MI(Menus[MENU_EDIT], 4, 'Coller',             'Ctrl+V');
 MS(Menus[MENU_EDIT], 5);
 MI(Menus[MENU_EDIT], 6, 'Effacer molecule',   '');
 MI(Menus[MENU_EDIT], 7, 'Effacer selection',  'Del');

 { Construire }
 ZeroItem(Menus[MENU_BUILD]);
 Menus[MENU_BUILD].Title := 'Construire';
 Menus[MENU_BUILD].ItemCount := 10;
 MI(Menus[MENU_BUILD], 0, 'Fragment H2O',         '');
 MI(Menus[MENU_BUILD], 1, 'Fragment CH4',         '');
 MI(Menus[MENU_BUILD], 2, 'Fragment NH3',         '');
 MI(Menus[MENU_BUILD], 3, 'Fragment Benzene',     '');
 MS(Menus[MENU_BUILD], 4);
 MI(Menus[MENU_BUILD], 5, 'Residu Glycine',       '');
 MI(Menus[MENU_BUILD], 6, 'Residu Alanine',       '');
 MS(Menus[MENU_BUILD], 7);
 MI(Menus[MENU_BUILD], 8, 'Detecter liaisons',    '');
 MI(Menus[MENU_BUILD], 9, 'Effacer hydrogenes',   '');

 { Selection }
 ZeroItem(Menus[MENU_SELECT]);
 Menus[MENU_SELECT].Title := 'Selection';
 Menus[MENU_SELECT].ItemCount := 12;
 MI(Menus[MENU_SELECT], 0,  'Tous',                  'Ctrl+A');
 MI(Menus[MENU_SELECT], 1,  'Aucun',                 '');
 MI(Menus[MENU_SELECT], 2,  'Inverser',              '');
 MS(Menus[MENU_SELECT], 3);
 MI(Menus[MENU_SELECT], 4,  'Hetero atomes',         '');
 MI(Menus[MENU_SELECT], 5,  'Eau',                   '');
 MI(Menus[MENU_SELECT], 6,  'Proteine',              '');
 MI(Menus[MENU_SELECT], 7,  'Acides nucleiques',     '');
 MS(Menus[MENU_SELECT], 8);
 MI(Menus[MENU_SELECT], 9,  'Squelette (backbone)',  '');
 MI(Menus[MENU_SELECT], 10, 'Chaine laterale',       '');
 MI(Menus[MENU_SELECT], 11, 'Hydrogenes',            '');

 { Affichage }
 ZeroItem(Menus[MENU_DISPLAY]);
 Menus[MENU_DISPLAY].Title := 'Affichage';
 Menus[MENU_DISPLAY].ItemCount := 14;
 MI(Menus[MENU_DISPLAY], 0,  'Lignes',                'L');
 MI(Menus[MENU_DISPLAY], 1,  'Batonnets',             'S');
 MI(Menus[MENU_DISPLAY], 2,  'Spheres (CPK)',         'C');
 MI(Menus[MENU_DISPLAY], 3,  'Boules + batonnets',    'B');
 MI(Menus[MENU_DISPLAY], 4,  'Ruban (cartoon)',       'R');
 MI(Menus[MENU_DISPLAY], 5,  'Points',                '');
 MS(Menus[MENU_DISPLAY], 6);
 MI(Menus[MENU_DISPLAY], 7,  'Axes',                  'A');
 Menus[MENU_DISPLAY].Items[7].Checked := True;
 MI(Menus[MENU_DISPLAY], 8,  'Grille',                'G');
 MI(Menus[MENU_DISPLAY], 9,  'Boite englobante',      '');
 MI(Menus[MENU_DISPLAY], 10, 'Etiquettes selection',  '');
 MI(Menus[MENU_DISPLAY], 11, 'Hydrogenes',            'H');
 MS(Menus[MENU_DISPLAY], 12);
 MI(Menus[MENU_DISPLAY], 13, 'Console',               'Tab');

 { Configuration (Setting) }
 ZeroItem(Menus[MENU_SETTING]);
 Menus[MENU_SETTING].Title := 'Configuration';
 Menus[MENU_SETTING].ItemCount := 10;
 MI(Menus[MENU_SETTING], 0,  'Couleur CPK',           '');
 MI(Menus[MENU_SETTING], 1,  'Coloration par chaine', '');
 MI(Menus[MENU_SETTING], 2,  'Coloration par residu', '');
 MI(Menus[MENU_SETTING], 3,  'Coloration unique...',  '');
 MS(Menus[MENU_SETTING], 4);
 MI(Menus[MENU_SETTING], 5,  'Fond noir',             '');
 MI(Menus[MENU_SETTING], 6,  'Fond blanc',            '');
 MI(Menus[MENU_SETTING], 7,  'Fond bleu',             '');
 MS(Menus[MENU_SETTING], 8);
 MI(Menus[MENU_SETTING], 9,  'Reinitialiser vue',     'Ctrl+Home');

 { Scene }
 ZeroItem(Menus[MENU_SCENE]);
 Menus[MENU_SCENE].Title := 'Scene';
 Menus[MENU_SCENE].ItemCount := 8;
 MI(Menus[MENU_SCENE], 0, 'Memoriser scene 1',  'F5');
 MI(Menus[MENU_SCENE], 1, 'Memoriser scene 2',  'F6');
 MI(Menus[MENU_SCENE], 2, 'Memoriser scene 3',  'F7');
 MI(Menus[MENU_SCENE], 3, 'Memoriser scene 4',  'F8');
 MS(Menus[MENU_SCENE], 4);
 MI(Menus[MENU_SCENE], 5, 'Rappeler scene 1',   'Shift+F5');
 MI(Menus[MENU_SCENE], 6, 'Rappeler scene 2',   'Shift+F6');
 MI(Menus[MENU_SCENE], 7, 'Rappeler precedente','');

 { Mouvement / Souris }
 ZeroItem(Menus[MENU_MOUSE]);
 Menus[MENU_MOUSE].Title := 'Mouvement';
 Menus[MENU_MOUSE].ItemCount := 8;
 MI(Menus[MENU_MOUSE], 0, 'Rotation auto',       'Espace');
 MI(Menus[MENU_MOUSE], 1, 'Stopper rotation',    '');
 MS(Menus[MENU_MOUSE], 2);
 MI(Menus[MENU_MOUSE], 3, 'Tourner +X',          'Up');
 MI(Menus[MENU_MOUSE], 4, 'Tourner -X',          'Down');
 MI(Menus[MENU_MOUSE], 5, 'Tourner +Y',          'Right');
 MI(Menus[MENU_MOUSE], 6, 'Tourner -Y',          'Left');
 MI(Menus[MENU_MOUSE], 7, 'Centrer',             '0');

 { Wizard / Assistant }
 ZeroItem(Menus[MENU_WIZARD]);
 Menus[MENU_WIZARD].Title := 'Assistant';
 Menus[MENU_WIZARD].ItemCount := 6;
 MI(Menus[MENU_WIZARD], 0, 'Mesure distance...',    '');
 MI(Menus[MENU_WIZARD], 1, 'Mesure angle...',       '');
 MI(Menus[MENU_WIZARD], 2, 'Mutagenese (residu)...','');
 MI(Menus[MENU_WIZARD], 3, 'Densite (carte)...',    '');
 MS(Menus[MENU_WIZARD], 4);
 MI(Menus[MENU_WIZARD], 5, 'Statistiques',          '');

 { Plugin }
 ZeroItem(Menus[MENU_PLUGIN]);
 Menus[MENU_PLUGIN].Title := 'Plugin';
 Menus[MENU_PLUGIN].ItemCount := 4;
 MI(Menus[MENU_PLUGIN], 0, 'Gestionnaire de plugins...', '');
 MS(Menus[MENU_PLUGIN], 1);
 MI(Menus[MENU_PLUGIN], 2, 'Convertir vers .MOL',  '');
 MI(Menus[MENU_PLUGIN], 3, 'Convertir vers .PDB',  '');

 { Aide }
 ZeroItem(Menus[MENU_HELP]);
 Menus[MENU_HELP].Title := 'Aide';
 Menus[MENU_HELP].ItemCount := 6;
 MI(Menus[MENU_HELP], 0, 'Documentation...',  'F1');
 MI(Menus[MENU_HELP], 1, 'Commandes PyMOL...',  '');
 MI(Menus[MENU_HELP], 2, 'Raccourcis clavier...', '');
 MI(Menus[MENU_HELP], 3, 'Souris...',  '');
 MS(Menus[MENU_HELP], 4);
 MI(Menus[MENU_HELP], 5, 'A propos...',     '');
End;

Procedure LayoutMenuBar;
Var I, X, W, MaxX, Row, TitleW : Integer;
Begin
 TitleW := 8 * Length(PROG_TITLE) + 8;
 MaxX := ScreenWidth - TitleW;
 X := 5; Row := 0;
 For I := 0 To NUM_MENUS - 1 Do Begin
  W := TextWidth(Menus[I].Title);
  If (Row = 0) And ((X + W + 12) > MaxX) Then Begin
   X := 5; Row := 1;
   MaxX := ScreenWidth - 5;
  End;
  Menus[I].X := X;
  Menus[I].Y := Row * MENUBAR_HEIGHT;
  X := X + W + 12;
 End;
 If Row = 0 Then MenuBarHeight := MENUBAR_HEIGHT
 Else MenuBarHeight := 2 * MENUBAR_HEIGHT;
End;

Procedure DrawMenuBar;
Var
 I, W : Integer;
Begin
 LayoutMenuBar;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, 0, ScreenWidth - 1, MenuBarHeight - 1);
 SetColor(DarkGray);
 Line(0, MenuBarHeight - 1, ScreenWidth - 1, MenuBarHeight - 1);
 For I := 0 To NUM_MENUS - 1 Do Begin
  W := TextWidth(Menus[I].Title);
  If MenuActive And (ActiveMenu = I) Then Begin
   SetFillStyle(SolidFill, Blue);
   Bar(Menus[I].X - 2, Menus[I].Y + 2,
       Menus[I].X + W + 2, Menus[I].Y + MENUBAR_HEIGHT - 3);
   SetColor(White);
  End
  Else SetColor(Black);
  OutTextXY(Menus[I].X, Menus[I].Y + 5, Menus[I].Title);
 End;
 SetColor(DarkGray);
 OutTextXY(ScreenWidth - 8 * Length(PROG_TITLE) - 4, 5, PROG_TITLE);
End;

Procedure DrawMenuPopup(MenuIdx: Integer);
Var
 I, ItemHeight, X, Y, W, H, Wmax, WSc : Integer;
 Cap, Sc : String;
Begin
 If (MenuIdx < 0) Or (MenuIdx >= NUM_MENUS) Then Exit;
 ItemHeight := 16;
 Wmax := 0;
 For I := 0 To Menus[MenuIdx].ItemCount - 1 Do Begin
  Cap := Menus[MenuIdx].Items[I].Caption;
  Sc  := Menus[MenuIdx].Items[I].Shortcut;
  W := TextWidth(Cap) + 30 + TextWidth(Sc) + 12;
  If W > Wmax Then Wmax := W;
 End;
 If Wmax < 180 Then Wmax := 180;
 Menus[MenuIdx].Width := Wmax;
 X := Menus[MenuIdx].X - 2;
 Y := Menus[MenuIdx].Y + MENUBAR_HEIGHT;
 W := Wmax;
 H := Menus[MenuIdx].ItemCount * ItemHeight + 4;
 { Decale a gauche si depasse }
 If X + W >= ScreenWidth Then X := ScreenWidth - W - 2;
 SetFillStyle(SolidFill, LightGray);
 Bar(X, Y, X + W, Y + H);
 SetColor(White);
 Line(X, Y, X + W, Y); Line(X, Y, X, Y + H);
 SetColor(DarkGray);
 Line(X + W, Y, X + W, Y + H); Line(X, Y + H, X + W, Y + H);
 SetColor(Black);
 Line(X + W + 1, Y + 1, X + W + 1, Y + H + 1);
 Line(X + 1, Y + H + 1, X + W + 1, Y + H + 1);
 For I := 0 To Menus[MenuIdx].ItemCount - 1 Do Begin
  If Menus[MenuIdx].Items[I].Separator Then Begin
   SetColor(DarkGray);
   Line(X + 4, Y + 2 + I * ItemHeight + ItemHeight Div 2,
        X + W - 4, Y + 2 + I * ItemHeight + ItemHeight Div 2);
  End
  Else Begin
   If ActiveMenuItem = I Then Begin
    SetFillStyle(SolidFill, Blue);
    Bar(X + 2, Y + 2 + I * ItemHeight,
        X + W - 2, Y + 2 + (I + 1) * ItemHeight);
    SetColor(White);
   End
   Else Begin
    If Menus[MenuIdx].Items[I].Enabled Then SetColor(Black)
    Else SetColor(DarkGray);
   End;
   If Menus[MenuIdx].Items[I].Checked Then
    OutTextXY(X + 5, Y + 4 + I * ItemHeight, #251);
   OutTextXY(X + 22, Y + 4 + I * ItemHeight,
             Menus[MenuIdx].Items[I].Caption);
   Sc := Menus[MenuIdx].Items[I].Shortcut;
   If Sc <> '' Then Begin
    WSc := TextWidth(Sc);
    OutTextXY(X + W - WSc - 8, Y + 4 + I * ItemHeight, Sc);
   End;
  End;
 End;
End;

Function MenuBarHitTest(Mx, My: Integer): Integer;
Var I, X1, X2: Integer;
Begin
 MenuBarHitTest := -1;
 If My >= MenuBarHeight Then Exit;
 For I := 0 To NUM_MENUS - 1 Do Begin
  X1 := Menus[I].X - 2;
  X2 := Menus[I].X + TextWidth(Menus[I].Title) + 2;
  If (Mx >= X1) And (Mx <= X2) And
     (My >= Menus[I].Y) And (My < Menus[I].Y + MENUBAR_HEIGHT) Then Begin
   MenuBarHitTest := I; Exit;
  End;
 End;
End;

Function MenuItemHitTest(Mx, My, MenuIdx: Integer): Integer;
Var Idx, Top, Bot, X, W: Integer;
Begin
 MenuItemHitTest := -1;
 If (MenuIdx < 0) Or (MenuIdx >= NUM_MENUS) Then Exit;
 Top := Menus[MenuIdx].Y + MENUBAR_HEIGHT;
 Bot := Top + Menus[MenuIdx].ItemCount * 16 + 4;
 X := Menus[MenuIdx].X - 2;
 W := Menus[MenuIdx].Width;
 If X + W >= ScreenWidth Then X := ScreenWidth - W - 2;
 If (Mx < X) Or (Mx > X + W) Or (My < Top) Or (My > Bot) Then Exit;
 Idx := (My - Top - 2) Div 16;
 If (Idx < 0) Or (Idx >= Menus[MenuIdx].ItemCount) Then Exit;
 MenuItemHitTest := Idx;
End;

Procedure MenuNextItem(Step: Integer);
Var N, I: Integer;
Begin
 If (ActiveMenu < 0) Or (ActiveMenu >= NUM_MENUS) Then Exit;
 N := Menus[ActiveMenu].ItemCount;
 If N <= 0 Then Exit;
 If ActiveMenuItem < 0 Then Begin
  If Step > 0 Then ActiveMenuItem := -1 Else ActiveMenuItem := N;
 End;
 For I := 1 To N Do Begin
  ActiveMenuItem := ActiveMenuItem + Step;
  If ActiveMenuItem < 0 Then ActiveMenuItem := N - 1;
  If ActiveMenuItem >= N Then ActiveMenuItem := 0;
  If Not Menus[ActiveMenu].Items[ActiveMenuItem].Separator Then Exit;
 End;
End;

{ --------------------------------------------------------------------- }
{ Dialogues : message, saisie                                           }
{ --------------------------------------------------------------------- }

Function PtInRect(MX, MY, X1, Y1, X2, Y2: Integer): Boolean;
Begin
 PtInRect := (MX >= X1) And (MX <= X2) And (MY >= Y1) And (MY <= Y2);
End;

Procedure DrawButton(X1, Y1, X2, Y2: Integer;
                     Const Caption: String; Focused: Boolean);
Var TX, TY: Integer;
Begin
 SetFillStyle(SolidFill, LightGray);
 Bar(X1+1, Y1+1, X2-1, Y2-1);
 { Bord 3D releve }
 SetColor(White);
 Line(X1, Y1, X2, Y1);
 Line(X1, Y1, X1, Y2);
 SetColor(DarkGray);
 Line(X1+1, Y2-1, X2-1, Y2-1);
 Line(X2-1, Y1+1, X2-1, Y2-1);
 SetColor(Black);
 Line(X1, Y2, X2, Y2);
 Line(X2, Y1, X2, Y2);
 TX := (X1 + X2 - TextWidth(Caption)) Div 2;
 TY := (Y1 + Y2 - TextHeight(Caption)) Div 2;
 SetColor(Black);
 OutTextXY(TX, TY, Caption);
 If Focused Then Begin
  SetColor(Black);
  Rectangle(X1+3, Y1+3, X2-3, Y2-3);
 End;
End;

Procedure DrawDialogBox(X1, Y1, X2, Y2: Integer; Const Title: String);
Begin
 SetFillStyle(SolidFill, LightGray);
 Bar(X1, Y1, X2, Y2);
 SetColor(White);
 Rectangle(X1, Y1, X2, Y2);
 SetColor(DarkGray);
 Rectangle(X1+1, Y1+1, X2-1, Y2-1);
 SetFillStyle(SolidFill, Blue);
 Bar(X1+2, Y1+2, X2-2, Y1+16);
 SetColor(White);
 OutTextXY(X1+6, Y1+5, Title);
 SetColor(Black);
End;

Procedure ShowMessage(Const Title, Msg: String);
Const BW = 70; BH = 22;
Var
 X1, Y1, X2, Y2 : Integer;
 BX1, BY1, BX2, BY2 : Integer;
 W : Integer;
 K : Char;
 Done : Boolean;
 LeftDown, OldLeft : Boolean;
Begin
 W := TextWidth(Msg) + 40;
 If W < TextWidth(Title) + 40 Then W := TextWidth(Title) + 40;
 If W < 240 Then W := 240;
 X1 := (ScreenWidth - W) Div 2;
 Y1 := (ScreenHeight - 110) Div 2;
 X2 := X1 + W;
 Y2 := Y1 + 110;
 BX1 := (X1 + X2 - BW) Div 2; BY1 := Y2 - BH - 10;
 BX2 := BX1 + BW;             BY2 := BY1 + BH;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 SetColor(Black);
 OutTextXY(X1 + 20, Y1 + 35, Msg);
 DrawButton(BX1, BY1, BX2, BY2, 'OK', True);
 ShowMouseCursor;
 { Drainer clavier et boutons souris en attente }
 While KeyPressed Do Begin
  K := ReadKey; If K = #0 Then If KeyPressed Then ReadKey;
 End;
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 OldLeft := (MouseButtons And 1) <> 0;
 Done := False;
 While Not Done Do Begin
  If KeyPressed Then Begin
   K := ReadKey;
   If K = #0 Then Begin
    If KeyPressed Then ReadKey;
   End
   Else If (K = #13) Or (K = #27) Or (K = ' ') Then Done := True;
  End;
  GetMouseLocal(MouseX, MouseY, MouseButtons);
  LeftDown := (MouseButtons And 1) <> 0;
  If LeftDown And Not OldLeft And
     PtInRect(MouseX, MouseY, BX1, BY1, BX2, BY2) Then
   Done := True;
  OldLeft := LeftDown;
  Delay(15);
 End;
 OldMouseButtons := MouseButtons;
End;

Function InputString(Const Title, Prompt: String;
                     Var Value: String): Boolean;
Const BW = 80; BH = 22;
Var
 X1, Y1, X2, Y2 : Integer;
 BOk1, BOk2, BCa1, BCa2, BBY1, BBY2 : Integer;
 W              : Integer;
 K              : Char;
 Done, OK       : Boolean;
 LeftDown, OldLeft : Boolean;
Begin
 W := 400;
 X1 := (ScreenWidth - W) Div 2;
 Y1 := (ScreenHeight - 140) Div 2;
 X2 := X1 + W;
 Y2 := Y1 + 140;
 BBY1 := Y2 - BH - 10; BBY2 := BBY1 + BH;
 BOk1 := X1 + W Div 2 - BW - 10; BOk2 := BOk1 + BW;
 BCa1 := X1 + W Div 2 + 10;       BCa2 := BCa1 + BW;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 SetColor(Black);
 OutTextXY(X1 + 15, Y1 + 30, Prompt);
 Done := False;
 OK   := False;
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 OldLeft := (MouseButtons And 1) <> 0;
 ShowMouseCursor;
 While Not Done Do Begin
  HideMouseCursor;
  SetFillStyle(SolidFill, White);
  Bar(X1 + 15, Y1 + 50, X2 - 15, Y1 + 70);
  SetColor(Black);
  Rectangle(X1 + 14, Y1 + 49, X2 - 14, Y1 + 71);
  OutTextXY(X1 + 18, Y1 + 56, Value + '_');
  DrawButton(BOk1, BBY1, BOk2, BBY2, 'OK',      False);
  DrawButton(BCa1, BBY1, BCa2, BBY2, 'Annuler', False);
  ShowMouseCursor;
  { Poll clavier + souris }
  Repeat
   GetMouseLocal(MouseX, MouseY, MouseButtons);
   LeftDown := (MouseButtons And 1) <> 0;
   If LeftDown And Not OldLeft Then Begin
    If PtInRect(MouseX, MouseY, BOk1, BBY1, BOk2, BBY2) Then Begin
     Done := True; OK := True;
    End
    Else If PtInRect(MouseX, MouseY, BCa1, BBY1, BCa2, BBY2) Then Begin
     Done := True; OK := False;
    End;
   End;
   OldLeft := LeftDown;
   If Done Or KeyPressed Then Break;
   Delay(15);
  Until False;
  If Done Then Break;
  K := ReadKey;
  If K = #0 Then Begin
   ReadKey;
  End
  Else If K = #13 Then Begin Done := True; OK := True; End
  Else If K = #27 Then Begin Done := True; OK := False; End
  Else If K = #8 Then Begin
   If Length(Value) > 0 Then Delete(Value, Length(Value), 1);
  End
  Else If (Ord(K) >= 32) And (Length(Value) < 80) Then
   Value := Value + K;
 End;
 OldMouseButtons := MouseButtons;
 InputString := OK;
End;

{ --------------------------------------------------------------------- }
{ Boite d'ouverture de fichier (style PSP.PAS)                          }
{ --------------------------------------------------------------------- }

Const
 FD_MAX_FILES = 256;

Type
 TFDEntry = Record
  Name : String[14];
  IsDir: Boolean;
 End;

Var
 FDFiles      : Array[0..FD_MAX_FILES-1] Of TFDEntry;
 FDCount      : Integer;
 FDTop        : Integer;
 FDSel        : Integer;
 FDMask       : String;
 FDFileName   : String;

Procedure FD_ScanDir(Const Mask: String);
Var
 Sr   : SearchRec;
 K    : Integer;
 Tmp  : TFDEntry;
 I, J : Integer;
Begin
 FDCount := 0;
 { Repertoires }
 FindFirst('*.*', Directory, Sr);
 While (DosError = 0) And (FDCount < FD_MAX_FILES) Do Begin
  If ((Sr.Attr And Directory) <> 0) And (Sr.Name <> '.') Then Begin
   FDFiles[FDCount].Name := Sr.Name;
   FDFiles[FDCount].IsDir := True;
   Inc(FDCount);
  End;
  FindNext(Sr);
 End;
 {$IFDEF FPC}FindClose(Sr);{$ENDIF}
 { Fichiers correspondant au masque }
 FindFirst(Mask, AnyFile, Sr);
 While (DosError = 0) And (FDCount < FD_MAX_FILES) Do Begin
  If ((Sr.Attr And Directory) = 0) Then Begin
   FDFiles[FDCount].Name := Sr.Name;
   FDFiles[FDCount].IsDir := False;
   Inc(FDCount);
  End;
  FindNext(Sr);
 End;
 {$IFDEF FPC}FindClose(Sr);{$ENDIF}
 { Tri par nom (insertion) }
 For I := 1 To FDCount - 1 Do Begin
  Tmp := FDFiles[I];
  J := I - 1;
  While (J >= 0) And
        ((FDFiles[J].IsDir < Tmp.IsDir) Or
         ((FDFiles[J].IsDir = Tmp.IsDir) And
          (UpStr(FDFiles[J].Name) > UpStr(Tmp.Name)))) Do Begin
   FDFiles[J+1] := FDFiles[J];
   Dec(J);
  End;
  FDFiles[J+1] := Tmp;
 End;
 K := 0; If K <> 0 Then ;
End;

Procedure FD_DrawList(X1, Y1, X2, Y2: Integer);
Var
 I, Y, MaxLn : Integer;
 S : String;
Begin
 SetFillStyle(SolidFill, White);
 Bar(X1, Y1, X2, Y2);
 SetColor(Black);
 Rectangle(X1, Y1, X2, Y2);
 MaxLn := (Y2 - Y1 - 2) Div 10;
 For I := 0 To MaxLn - 1 Do Begin
  If FDTop + I >= FDCount Then Break;
  Y := Y1 + 2 + I * 10;
  If FDTop + I = FDSel Then Begin
   SetFillStyle(SolidFill, Blue);
   Bar(X1 + 1, Y, X2 - 1, Y + 9);
   SetColor(White);
  End
  Else SetColor(Black);
  If FDFiles[FDTop + I].IsDir Then S := '[' + FDFiles[FDTop + I].Name + ']'
  Else S := FDFiles[FDTop + I].Name;
  OutTextXY(X1 + 4, Y + 1, S);
 End;
End;

Function FileOpenDialog(Const Title: String; Var FileName: String;
                        Const InitialMask: String): Boolean;
Const BW = 80; BH = 22;
Var
 X1, Y1, X2, Y2 : Integer;
 LX1, LY1, LX2, LY2 : Integer;
 BOu1, BOu2, BAn1, BAn2, BMa1, BMa2, BBY1, BBY2 : Integer;
 K, SimK : Char;
 Done, Result : Boolean;
 OldDir, CurDir : String;
 NewMask : String;
 MaxLn : Integer;
 LeftDown, OldLeft : Boolean;
Begin
 FileOpenDialog := False;
 GetDir(0, OldDir);
 FDMask := InitialMask;
 FD_ScanDir(FDMask);
 FDTop := 0; FDSel := 0;
 X1 := 60; Y1 := 60;
 X2 := ScreenWidth - 60; Y2 := ScreenHeight - 60;
 LX1 := X1 + 10; LY1 := Y1 + 50; LX2 := X2 - 10; LY2 := Y2 - 90;
 MaxLn := (LY2 - LY1 - 2) Div 10;
 BBY1 := Y2 - BH - 8; BBY2 := BBY1 + BH;
 BOu1 := X2 - 3*(BW + 8) - 10; BOu2 := BOu1 + BW;
 BAn1 := X2 - 2*(BW + 8) - 10; BAn2 := BAn1 + BW;
 BMa1 := X2 -   (BW + 8) - 10; BMa2 := BMa1 + BW;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 GetDir(0, CurDir);
 SetColor(Black);
 OutTextXY(X1 + 10, Y1 + 22, 'Dossier : ' + CurDir);
 OutTextXY(X1 + 10, Y1 + 35, 'Masque : ' + FDMask);
 FD_DrawList(LX1, LY1, LX2, LY2);
 OutTextXY(X1 + 10, Y2 - 80, 'Fichier:');
 DrawButton(BOu1, BBY1, BOu2, BBY2, 'Ouvrir',  True);
 DrawButton(BAn1, BBY1, BAn2, BBY2, 'Annuler', False);
 DrawButton(BMa1, BBY1, BMa2, BBY2, 'Masque',  False);
 ShowMouseCursor;
 Done := False; Result := False;
 FDFileName := '';
 SimK := #0;
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 OldLeft := (MouseButtons And 1) <> 0;
 While Not Done Do Begin
  HideMouseCursor;
  SetFillStyle(SolidFill, White);
  Bar(X1 + 80, Y2 - 83, X2 - 10, Y2 - 70);
  SetColor(Black);
  Rectangle(X1 + 80, Y2 - 84, X2 - 9, Y2 - 69);
  OutTextXY(X1 + 84, Y2 - 81, FDFileName + '_');
  FD_DrawList(LX1, LY1, LX2, LY2);
  ShowMouseCursor;
  { Poll clavier + souris pour les boutons }
  Repeat
   GetMouseLocal(MouseX, MouseY, MouseButtons);
   LeftDown := (MouseButtons And 1) <> 0;
   If LeftDown And Not OldLeft Then Begin
    If PtInRect(MouseX, MouseY, BOu1, BBY1, BOu2, BBY2) Then SimK := #13
    Else If PtInRect(MouseX, MouseY, BAn1, BBY1, BAn2, BBY2) Then SimK := #27
    Else If PtInRect(MouseX, MouseY, BMa1, BBY1, BMa2, BBY2) Then SimK := #9;
   End;
   OldLeft := LeftDown;
   If KeyPressed Or (SimK <> #0) Then Break;
   Delay(15);
  Until False;
  If SimK <> #0 Then Begin
   K := SimK; SimK := #0;
  End
  Else K := ReadKey;
  If K = #0 Then Begin
   K := ReadKey;
   Case K Of
    #72: Begin  { Up }
     If FDSel > 0 Then Dec(FDSel);
     If FDSel < FDTop Then FDTop := FDSel;
    End;
    #80: Begin  { Down }
     If FDSel < FDCount - 1 Then Inc(FDSel);
     If FDSel >= FDTop + MaxLn Then FDTop := FDSel - MaxLn + 1;
    End;
    #73: Begin  { PgUp }
     Dec(FDSel, MaxLn); If FDSel < 0 Then FDSel := 0;
     FDTop := FDSel;
    End;
    #81: Begin  { PgDn }
     Inc(FDSel, MaxLn); If FDSel >= FDCount Then FDSel := FDCount - 1;
     FDTop := FDSel - MaxLn + 1; If FDTop < 0 Then FDTop := 0;
    End;
    #71: Begin FDSel := 0; FDTop := 0; End;
    #79: Begin FDSel := FDCount - 1;
               FDTop := FDSel - MaxLn + 1; If FDTop < 0 Then FDTop := 0;
         End;
   End;
  End
  Else If K = #13 Then Begin
   If FDFileName <> '' Then Begin
    FileName := FDFileName;
    Result := True; Done := True;
   End
   Else If (FDSel >= 0) And (FDSel < FDCount) Then Begin
    If FDFiles[FDSel].IsDir Then Begin
     {$I-}ChDir(FDFiles[FDSel].Name);{$I+}
     If IOResult <> 0 Then ;
     FD_ScanDir(FDMask);
     FDTop := 0; FDSel := 0;
     GetDir(0, CurDir);
     HideMouseCursor;
     SetFillStyle(SolidFill, LightGray);
     Bar(X1 + 60, Y1 + 18, X2 - 5, Y1 + 30);
     OutTextXY(X1 + 10, Y1 + 22, 'Dossier : ' + CurDir);
     ShowMouseCursor;
    End
    Else Begin
     FileName := FDFiles[FDSel].Name;
     Result := True; Done := True;
    End;
   End;
  End
  Else If K = #27 Then Begin
   Result := False; Done := True;
  End
  Else If K = #9 Then Begin
   NewMask := FDMask;
   ShowMouseCursor;
   If InputString('Masque', 'Entrer le masque (ex: *.PDB):', NewMask) Then Begin
    FDMask := NewMask;
    FD_ScanDir(FDMask);
    FDTop := 0; FDSel := 0;
   End;
   HideMouseCursor;
   DrawDialogBox(X1, Y1, X2, Y2, Title);
   SetColor(Black);
   OutTextXY(X1 + 10, Y1 + 22, 'Dossier : ' + CurDir);
   OutTextXY(X1 + 10, Y1 + 35, 'Masque : ' + FDMask);
   OutTextXY(X1 + 10, Y2 - 80, 'Fichier:');
   DrawButton(BOu1, BBY1, BOu2, BBY2, 'Ouvrir',  True);
   DrawButton(BAn1, BBY1, BAn2, BBY2, 'Annuler', False);
   DrawButton(BMa1, BBY1, BMa2, BBY2, 'Masque',  False);
   ShowMouseCursor;
  End
  Else If K = #8 Then Begin
   If Length(FDFileName) > 0 Then Delete(FDFileName, Length(FDFileName), 1);
  End
  Else If (Ord(K) >= 32) And (Length(FDFileName) < 60) Then
   FDFileName := FDFileName + K;
 End;
 If Not Result Then ChDir(OldDir);
 FileOpenDialog := Result;
End;

Function FileSaveDialog(Const Title: String; Var FileName: String;
                        Const Mask: String): Boolean;
Var
 X1, Y1, X2, Y2 : Integer;
 Cur : String;
Begin
 Cur := FileName;
 X1 := 80; Y1 := 130;
 X2 := ScreenWidth - 80; Y2 := Y1 + 140;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 OutTextXY(X1 + 10, Y1 + 24, 'Masque : ' + Mask);
 OutTextXY(X1 + 10, Y1 + 40, 'Nom du fichier a enregistrer :');
 ShowMouseCursor;
 If InputString(Title, 'Fichier :', Cur) And (Cur <> '') Then Begin
  FileName := Cur;
  FileSaveDialog := True;
 End
 Else FileSaveDialog := False;
End;

{ --------------------------------------------------------------------- }
{ Constructeurs de fragments simples                                    }
{ --------------------------------------------------------------------- }

Function MakeAtom(Const Nm, ResNm: String; ChId: Char; ResSq: Integer;
                  X, Y, Z: Real; Z2: Byte; Const Sym: String;
                  Het: Boolean): TAtom;
Var A: TAtom;
Begin
 FillChar(A, SizeOf(A), 0);
 A.Serial := 0;
 A.Name := Nm;
 A.AltLoc := ' ';
 A.ResName := ResNm;
 A.ChainID := ChId;
 A.ResSeq := ResSq;
 A.ICode := ' ';
 A.X := X; A.Y := Y; A.Z := Z;
 A.Occup := 1.0; A.TempF := 0.0;
 A.Element := Sym;
 A.ElemZ := Z2;
 A.Charge := 0;
 A.Kind := AtomKindFor(ResNm, Het);
 A.Selected := False;
 A.Visible := True;
 A.Color := 0;
 MakeAtom := A;
End;

Procedure BuildWater;
Var O, H1, H2: Integer;
Begin
 ClearMolecule(Mol);
 Mol.Title := 'H2O';
 Mol.Source := 'PYMOL CLONE FRAGMENT';
 O  := AddAtom(Mol, MakeAtom('O',  'HOH','A',1, 0.0,  0.0,  0.0, 8,'O',True));
 H1 := AddAtom(Mol, MakeAtom('H1', 'HOH','A',1, 0.96, 0.0,  0.0, 1,'H',True));
 H2 := AddAtom(Mol, MakeAtom('H2', 'HOH','A',1,-0.24, 0.93, 0.0, 1,'H',True));
 AddBond(Mol, O, H1, 1);
 AddBond(Mol, O, H2, 1);
 ComputeBounds(Mol);
 BuildResidueIndex(Mol);
 Mol.Loaded := True;
End;

Procedure BuildMethane;
Var C, H: Integer;
Begin
 ClearMolecule(Mol);
 Mol.Title := 'CH4';
 Mol.Source := 'PYMOL CLONE FRAGMENT';
 C := AddAtom(Mol, MakeAtom('C',  'LIG','A',1, 0,    0,    0,    6,'C',False));
 H := AddAtom(Mol, MakeAtom('H1', 'LIG','A',1, 0.63, 0.63, 0.63, 1,'H',False));
 AddBond(Mol, C, H, 1);
 H := AddAtom(Mol, MakeAtom('H2', 'LIG','A',1,-0.63,-0.63, 0.63, 1,'H',False));
 AddBond(Mol, C, H, 1);
 H := AddAtom(Mol, MakeAtom('H3', 'LIG','A',1,-0.63, 0.63,-0.63, 1,'H',False));
 AddBond(Mol, C, H, 1);
 H := AddAtom(Mol, MakeAtom('H4', 'LIG','A',1, 0.63,-0.63,-0.63, 1,'H',False));
 AddBond(Mol, C, H, 1);
 ComputeBounds(Mol);
 BuildResidueIndex(Mol);
 Mol.Loaded := True;
End;

Procedure BuildAmmonia;
Var N, H: Integer;
Begin
 ClearMolecule(Mol);
 Mol.Title := 'NH3';
 Mol.Source := 'PYMOL CLONE FRAGMENT';
 N := AddAtom(Mol, MakeAtom('N',  'LIG','A',1, 0,    0,    0,    7,'N',False));
 H := AddAtom(Mol, MakeAtom('H1', 'LIG','A',1, 0.94, 0.0,  0.34, 1,'H',False));
 AddBond(Mol, N, H, 1);
 H := AddAtom(Mol, MakeAtom('H2', 'LIG','A',1,-0.47, 0.81, 0.34, 1,'H',False));
 AddBond(Mol, N, H, 1);
 H := AddAtom(Mol, MakeAtom('H3', 'LIG','A',1,-0.47,-0.81, 0.34, 1,'H',False));
 AddBond(Mol, N, H, 1);
 ComputeBounds(Mol);
 BuildResidueIndex(Mol);
 Mol.Loaded := True;
End;

Procedure BuildBenzene;
Var
 I : Integer;
 Ang : Real;
 R : Real;
 IdxC : Array[0..5] Of Integer;
 IdxH : Array[0..5] Of Integer;
 NameS : String;
Begin
 ClearMolecule(Mol);
 Mol.Title := 'BENZENE C6H6';
 Mol.Source := 'PYMOL CLONE FRAGMENT';
 R := 1.40;
 For I := 0 To 5 Do Begin
  Ang := I * (PI_VAL / 3.0);
  NameS := 'C' + IntToStr(I+1);
  IdxC[I] := AddAtom(Mol, MakeAtom(NameS, 'BNZ', 'A', 1,
                     R * Cos(Ang), R * Sin(Ang), 0.0,
                     6, 'C', False));
 End;
 For I := 0 To 5 Do Begin
  Ang := I * (PI_VAL / 3.0);
  NameS := 'H' + IntToStr(I+1);
  IdxH[I] := AddAtom(Mol, MakeAtom(NameS, 'BNZ', 'A', 1,
                     2.49 * Cos(Ang), 2.49 * Sin(Ang), 0.0,
                     1, 'H', False));
 End;
 For I := 0 To 5 Do Begin
  AddBond(Mol, IdxC[I], IdxC[(I + 1) Mod 6], 1 + (I Mod 2));
  AddBond(Mol, IdxC[I], IdxH[I], 1);
 End;
 ComputeBounds(Mol);
 BuildResidueIndex(Mol);
 Mol.Loaded := True;
End;

{ --------------------------------------------------------------------- }
{ Mesures                                                               }
{ --------------------------------------------------------------------- }

Function AtomDistance(I, J: Integer): Real;
Var Dx, Dy, Dz: Real;
Begin
 Dx := Mol.Atoms[I].X - Mol.Atoms[J].X;
 Dy := Mol.Atoms[I].Y - Mol.Atoms[J].Y;
 Dz := Mol.Atoms[I].Z - Mol.Atoms[J].Z;
 AtomDistance := Sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
End;

Procedure SelectedAtomList(Var Idx: Array Of Integer; Var N: Integer;
                           Max: Integer);
Var I: Integer;
Begin
 N := 0;
 For I := 0 To Mol.NAtoms - 1 Do
  If Mol.Atoms[I].Selected And (N < Max) Then Begin
   Idx[N] := I; Inc(N);
  End;
End;

Procedure CmdDistance;
Var
 Idx: Array[0..3] Of Integer;
 N: Integer;
 D: Real;
Begin
 SelectedAtomList(Idx, N, 4);
 If N < 2 Then Begin
  ConsoleErr('Selectionner deux atomes pour mesurer la distance.');
  Exit;
 End;
 D := AtomDistance(Idx[0], Idx[1]);
 ConsoleSay('Distance: ' + RealToStr(D, 3) + ' A');
End;

{ --------------------------------------------------------------------- }
{ Coloration globale par chaine et par residu                           }
{ --------------------------------------------------------------------- }

Procedure ColorByChain;
Var I, K: Integer; C: Byte;
Begin
 For I := 0 To Mol.NAtoms - 1 Do Begin
  K := Ord(Mol.Atoms[I].ChainID);
  Case (K Mod 8) Of
   0: C := LightGreen;
   1: C := LightCyan;
   2: C := Yellow;
   3: C := LightRed;
   4: C := LightMagenta;
   5: C := LightBlue;
   6: C := Brown;
  Else
   C := White;
  End;
  Mol.Atoms[I].Color := C;
 End;
 ConsoleSay('Coloration par chaine appliquee.');
End;

Procedure ColorByResidue;
Var I: Integer; C: Byte;
Begin
 For I := 0 To Mol.NAtoms - 1 Do Begin
  Case (Mol.Atoms[I].ResSeq Mod 7) Of
   0: C := LightCyan;
   1: C := Yellow;
   2: C := LightGreen;
   3: C := LightRed;
   4: C := LightBlue;
   5: C := LightMagenta;
  Else
   C := White;
  End;
  Mol.Atoms[I].Color := C;
 End;
 ConsoleSay('Coloration par residu appliquee.');
End;

Procedure ColorByCPK;
Var I: Integer;
Begin
 For I := 0 To Mol.NAtoms - 1 Do Mol.Atoms[I].Color := 0;
 ConsoleSay('Coloration CPK retablie.');
End;

{ --------------------------------------------------------------------- }
{ Action de menu                                                        }
{ --------------------------------------------------------------------- }

Procedure DoOpenFile;
Var FN: String;
Begin
 FN := '';
 If FileOpenDialog('Ouvrir un fichier moleculaire', FN,
                   '*.PDB') Then Begin
  If LoadMolecule(FN, Mol) Then Begin
   ZoomToFit;
   ConsoleSay('Charge: ' + FN + '  (' + IntToStr(Mol.NAtoms) + ' atomes)');
   If PosI('PROTEIN DATA BANK', UpStr(Mol.Source)) > 0 Then
    ConsoleSay('Source: Protein Data Bank reconnue.');
  End
  Else ConsoleErr('Echec du chargement: ' + FN);
 End;
End;

Procedure DoOpenScript;
Var FN: String;
Begin
 FN := '';
 If FileOpenDialog('Executer un script .PML', FN, '*.PML') Then Begin
  ConsoleSay('@ ' + FN);
  If Not LoadPML(FN) Then ConsoleErr('Echec script: ' + FN);
 End;
End;

Procedure DoSavePDB;
Var FN: String;
Begin
 If Not Mol.Loaded Then Begin ConsoleErr('Aucune molecule chargee.'); Exit; End;
 FN := 'OUT.PDB';
 If FileSaveDialog('Enregistrer .PDB', FN, '*.PDB') Then Begin
  If SavePDB(FN, Mol) Then ConsoleSay('Enregistre: ' + FN)
  Else ConsoleErr('Erreur d''ecriture: ' + FN);
 End;
End;

Procedure DoSaveMOL;
Var FN: String;
Begin
 If Not Mol.Loaded Then Begin ConsoleErr('Aucune molecule chargee.'); Exit; End;
 FN := 'OUT.MOL';
 If FileSaveDialog('Enregistrer .MOL', FN, '*.MOL') Then Begin
  If SaveMOL(FN, Mol) Then ConsoleSay('Enregistre: ' + FN)
  Else ConsoleErr('Erreur d''ecriture: ' + FN);
 End;
End;

Procedure DoReinit;
Begin
 ClearMolecule(Mol);
 RotX := 0; RotY := 0; RotZ := 0;
 ViewPanX := 0; ViewPanY := 0;
 ViewScale := 30.0;
 RenderMode := RENDER_BALLSTICK;
 BgColor := Black;
 SimSpinning := False;
 ConsoleSay('reinitialize');
End;

Procedure DoReload;
Begin
 If Mol.FileName <> '' Then Begin
  If LoadMolecule(Mol.FileName, Mol) Then Begin
   ZoomToFit;
   ConsoleSay('Recharge: ' + Mol.FileName);
  End
  Else ConsoleErr('Echec rechargement.');
 End
 Else ConsoleErr('Aucun fichier a recharger.');
End;

Procedure DoSetRenderMode(M: Byte);
Begin
 RenderMode := M;
 Case M Of
  RENDER_LINES     : ConsoleSay('show lines');
  RENDER_STICKS    : ConsoleSay('show sticks');
  RENDER_SPHERES   : ConsoleSay('show spheres');
  RENDER_BALLSTICK : ConsoleSay('show ball-and-stick');
  RENDER_RIBBON    : ConsoleSay('show cartoon');
  RENDER_DOTS      : ConsoleSay('show dots');
 End;
End;

Procedure ToggleFlag(F: Word; Const N: String);
Begin
 If (ViewFlags And F) <> 0 Then Begin
  ViewFlags := ViewFlags And (Not F);
  ConsoleSay(N + ' : off');
 End
 Else Begin
  ViewFlags := ViewFlags Or F;
  ConsoleSay(N + ' : on');
 End;
End;

Procedure DispatchMenu(MenuIdx, ItemIdx: Integer);
Var S: String;
Begin
 Case MenuIdx Of
  MENU_FILE:
   Case ItemIdx Of
    0: DoOpenFile;
    1: DoReload;
    2: DoOpenScript;
    4: DoSavePDB;
    5: DoSaveMOL;
    6: ConsoleSay('save session (non implemente)');
    7: ConsoleSay('png (non implemente)');
    9: BuildBenzene;
    10: ConsoleSay('Export PPM (non implemente)');
    12: DoReinit;
    13: Running := False;
   End;
  MENU_EDIT:
   Case ItemIdx Of
    0,1: ConsoleSay('annuler/refaire (non implemente)');
    3,4: ConsoleSay('copier/coller (non implemente)');
    6: Begin ClearMolecule(Mol); ConsoleSay('Molecule effacee.'); End;
    7: ConsoleSay('delete selection (non implemente)');
   End;
  MENU_BUILD:
   Case ItemIdx Of
    0: Begin BuildWater;    ZoomToFit; ConsoleSay('build water'); End;
    1: Begin BuildMethane;  ZoomToFit; ConsoleSay('build methane'); End;
    2: Begin BuildAmmonia;  ZoomToFit; ConsoleSay('build ammonia'); End;
    3: Begin BuildBenzene;  ZoomToFit; ConsoleSay('build benzene'); End;
    5,6: ConsoleSay('Residu (non implemente).');
    8: Begin AutoDetectBonds(Mol);
            ConsoleSay('Liaisons detectees: ' + IntToStr(Mol.NBonds)); End;
    9: ConsoleSay('remove hydrogens (non implemente)');
   End;
  MENU_SELECT:
   Case ItemIdx Of
    0: ConsoleSay('select all (' +
                  IntToStr(SelectAtoms(Mol, 'all')) + ')');
    1: ConsoleSay('select none (' +
                  IntToStr(SelectAtoms(Mol, 'none')) + ')');
    2: ConsoleSay('invert (non implemente)');
    4: ConsoleSay('select hetatm (' +
                  IntToStr(SelectAtoms(Mol, 'hetatm')) + ')');
    5: ConsoleSay('select water (' +
                  IntToStr(SelectAtoms(Mol, 'water')) + ')');
    6: ConsoleSay('select protein (' +
                  IntToStr(SelectAtoms(Mol, 'protein')) + ')');
    7: ConsoleSay('select nucleic (' +
                  IntToStr(SelectAtoms(Mol, 'nucleic')) + ')');
    9: ConsoleSay('select backbone (' +
                  IntToStr(SelectAtoms(Mol, 'backbone')) + ')');
    10: ConsoleSay('select sidechain (' +
                  IntToStr(SelectAtoms(Mol, 'sidechain')) + ')');
    11: ConsoleSay('select hydrogens (' +
                  IntToStr(SelectAtoms(Mol, 'hydrogens')) + ')');
   End;
  MENU_DISPLAY:
   Case ItemIdx Of
    0: DoSetRenderMode(RENDER_LINES);
    1: DoSetRenderMode(RENDER_STICKS);
    2: DoSetRenderMode(RENDER_SPHERES);
    3: DoSetRenderMode(RENDER_BALLSTICK);
    4: DoSetRenderMode(RENDER_RIBBON);
    5: DoSetRenderMode(RENDER_DOTS);
    7: Begin ToggleFlag(SHOW_AXES, 'axes');
            Menus[MENU_DISPLAY].Items[7].Checked :=
              (ViewFlags And SHOW_AXES) <> 0; End;
    8: Begin ToggleFlag(SHOW_GRID, 'grille');
            Menus[MENU_DISPLAY].Items[8].Checked :=
              (ViewFlags And SHOW_GRID) <> 0; End;
    9: Begin ToggleFlag(SHOW_BBOX, 'boite');
            Menus[MENU_DISPLAY].Items[9].Checked :=
              (ViewFlags And SHOW_BBOX) <> 0; End;
    10: Begin ToggleFlag(SHOW_LABELS, 'etiquettes');
            Menus[MENU_DISPLAY].Items[10].Checked :=
              (ViewFlags And SHOW_LABELS) <> 0; End;
    11: Begin ToggleFlag(SHOW_HYDROGENS, 'hydrogenes');
            Menus[MENU_DISPLAY].Items[11].Checked :=
              (ViewFlags And SHOW_HYDROGENS) <> 0; End;
    13: Begin ConsoleVisible := Not ConsoleVisible;
              Menus[MENU_DISPLAY].Items[13].Checked := ConsoleVisible;
              If ConsoleVisible Then
               CanvasBottom := ScreenHeight - STATUSBAR_HEIGHT - CONSOLE_HEIGHT - 1
              Else
               CanvasBottom := ScreenHeight - STATUSBAR_HEIGHT - 1;
              ZoomToFit; End;
   End;
  MENU_SETTING:
   Case ItemIdx Of
    0: ColorByCPK;
    1: ColorByChain;
    2: ColorByResidue;
    3: ConsoleSay('color <couleur> sele (saisir en console)');
    5: Begin BgColor := Black; ConsoleSay('bg_color black'); End;
    6: Begin BgColor := White; ConsoleSay('bg_color white'); End;
    7: Begin BgColor := Blue;  ConsoleSay('bg_color blue'); End;
    9: ZoomToFit;
   End;
  MENU_SCENE:
   ConsoleSay('Scenes (non implementees) - item ' + IntToStr(ItemIdx));
  MENU_MOUSE:
   Case ItemIdx Of
    0: Begin SimSpinning := True; ConsoleSay('spin on'); End;
    1: Begin SimSpinning := False; ConsoleSay('spin off'); End;
    3: RotX := RotX + 10;
    4: RotX := RotX - 10;
    5: RotY := RotY + 10;
    6: RotY := RotY - 10;
    7: ZoomToFit;
   End;
  MENU_WIZARD:
   Case ItemIdx Of
    0: CmdDistance;
    1: ConsoleSay('Angle (non implemente)');
    2: ConsoleSay('Mutagenese (non implementee)');
    3: ConsoleSay('Densite (non implementee)');
    5: ShowMessage('Statistiques',
                   IntToStr(Mol.NAtoms) + ' atomes, ' +
                   IntToStr(Mol.NBonds) + ' liaisons, ' +
                   IntToStr(Mol.NRes) + ' residus');
   End;
  MENU_PLUGIN:
   Case ItemIdx Of
    0: ShowMessage('Plugins', 'Aucun plugin externe disponible.');
    2: DoSaveMOL;
    3: DoSavePDB;
   End;
  MENU_HELP:
   Case ItemIdx Of
    0: ShowMessage('Documentation',
                   'PyMOL Clone - QChimie - documentation : voir README.');
    1: ShowMessage('Commandes',
                   'load, save, select, color, show, hide, bg_color, ' +
                   'rotate, zoom, reset, quit, spin, ray, run.');
    2: ShowMessage('Raccourcis',
                   'F1=Aide F3=Recharger Tab=Console L/S/C/B/R=Mode ' +
                   'A=Axes G=Grille H=Hydrogenes Esc=Menu');
    3: ShowMessage('Souris',
                   'Bouton gauche=Rotation  droit=Translation  ' +
                   'central=Zoom');
    5: ShowMessage('A propos',
                   PROG_TITLE + #13 +
                   'Auteur: Sylvain Maltais 2026');
   End;
 End;
 NeedRedraw := True;
End;

{ --------------------------------------------------------------------- }
{ Mini-interpreteur Python embarque                                     }
{                                                                       }
{ Supporte expressions, comparaisons, booleens, chaines, listes,        }
{ variables persistantes, indexation, ternaire (x if c else y) et plus  }
{ de 30 fonctions natives (print, len, range, sqrt, sin, sum, ...).     }
{ --------------------------------------------------------------------- }

Const
 PY_MAX_VARS    = 64;
 PY_MAX_LIST    = 64;
 PY_MAX_LISTS   = 32;
 PY_MAX_TOKENS  = 256;
 PY_MAX_ARGS    = 16;

Type
 TPyTokKind = (
  PtkEOF, PtkNumber, PtkString, PtkIdent,
  PtkPlus, PtkMinus, PtkStar, PtkSlash, PtkDSlash, PtkPercent, PtkDStar,
  PtkEq, PtkEqEq, PtkNotEq, PtkLT, PtkLE, PtkGT, PtkGE,
  PtkLPar, PtkRPar, PtkLBrk, PtkRBrk, PtkComma, PtkColon, PtkDot, PtkSemi,
  PtkAnd, PtkOr, PtkNot, PtkIn, PtkIf, PtkElse,
  PtkTrue, PtkFalse, PtkNone, PtkPrint,
  PtkFor, PtkWhile, PtkImport, PtkBreak, PtkContinue, PtkPass
 );

 TPyToken = Record
  Kind  : TPyTokKind;
  S     : String[80];
  N     : LongInt;
  F     : Real;
  IsFlt : Boolean;
 End;

 TPyValueKind = (PvNone, PvBool, PvInt, PvFloat, PvStr, PvList);

 TPyValue = Record
  Kind : TPyValueKind;
  B    : Boolean;
  I    : LongInt;
  F    : Real;
  S    : String;
  L    : Integer;
 End;

 TPyList = Record
  InUse : Boolean;
  Count : Integer;
  Items : Array[0..PY_MAX_LIST - 1] Of TPyValue;
 End;

 TPyVarSlot = Record
  InUse : Boolean;
  Name  : String[24];
  Value : TPyValue;
 End;

Var
 PyTokens         : Array[0..PY_MAX_TOKENS - 1] Of TPyToken;
 PyTokCount       : Integer;
 PyTokPos         : Integer;
 PyVars           : Array[0..PY_MAX_VARS - 1] Of TPyVarSlot;
 PyLists          : Array[0..PY_MAX_LISTS - 1] Of TPyList;
 PyError          : Boolean;
 PyErrorMsg       : String;
 PyBreakSignal    : Boolean;
 PyContinueSignal : Boolean;
 PyContBuf        : String;
 PyContPending    : Boolean;
 PyLoopDepth      : Integer;

Procedure PyInit;
Var I: Integer;
Begin
 For I := 0 To PY_MAX_VARS - 1 Do PyVars[I].InUse := False;
 For I := 0 To PY_MAX_LISTS - 1 Do PyLists[I].InUse := False;
 PyError := False;
 PyErrorMsg := '';
 PyBreakSignal := False;
 PyContinueSignal := False;
 PyContBuf := '';
 PyContPending := False;
 PyLoopDepth := 0;
End;

{ Heure systeme (epoch Unix en secondes, comme time.time() de Python) }
Function PyTimeNow: Real;
Var
 Y, Mo, D, DOW    : Word;
 H, Mi, Se, S100  : Word;
 A, Y2, M2, JD    : LongInt;
Begin
 GetDate(Y, Mo, D, DOW);
 GetTime(H, Mi, Se, S100);
 A  := (14 - LongInt(Mo)) Div 12;
 Y2 := LongInt(Y) + 4800 - A;
 M2 := LongInt(Mo) + 12 * A - 3;
 JD := LongInt(D) + (153 * M2 + 2) Div 5 + 365 * Y2
       + Y2 Div 4 - Y2 Div 100 + Y2 Div 400 - 32045;
 { Jour julien du 1970-01-01 = 2440588 }
 PyTimeNow := (JD - 2440588) * 86400.0
              + H * 3600.0 + Mi * 60.0 + Se + S100 / 100.0;
End;

Procedure PySetErr(Const M: String);
Begin
 If Not PyError Then Begin PyError := True; PyErrorMsg := M; End;
End;

{ Constructeurs de valeurs }

Function PyMakeNone: TPyValue;
Var V: TPyValue;
Begin
 FillChar(V, SizeOf(V), 0);
 V.Kind := PvNone; V.L := -1; PyMakeNone := V;
End;

Function PyMakeBool(BV: Boolean): TPyValue;
Var V: TPyValue;
Begin
 FillChar(V, SizeOf(V), 0);
 V.Kind := PvBool; V.B := BV; V.L := -1; PyMakeBool := V;
End;

Function PyMakeInt(IV: LongInt): TPyValue;
Var V: TPyValue;
Begin
 FillChar(V, SizeOf(V), 0);
 V.Kind := PvInt; V.I := IV; V.L := -1; PyMakeInt := V;
End;

Function PyMakeFloat(FV: Real): TPyValue;
Var V: TPyValue;
Begin
 FillChar(V, SizeOf(V), 0);
 V.Kind := PvFloat; V.F := FV; V.L := -1; PyMakeFloat := V;
End;

Function PyMakeStr(Const SV: String): TPyValue;
Var V: TPyValue;
Begin
 FillChar(V, SizeOf(V), 0);
 V.Kind := PvStr; V.S := SV; V.L := -1; PyMakeStr := V;
End;

Function PyNewList: Integer;
Var I: Integer;
Begin
 PyNewList := -1;
 For I := 0 To PY_MAX_LISTS - 1 Do
  If Not PyLists[I].InUse Then Begin
   PyLists[I].InUse := True;
   PyLists[I].Count := 0;
   PyNewList := I; Exit;
  End;
 PySetErr('Memoire de listes epuisee');
End;

Function PyMakeList(L: Integer): TPyValue;
Var V: TPyValue;
Begin
 FillChar(V, SizeOf(V), 0);
 V.Kind := PvList; V.L := L; PyMakeList := V;
End;

Procedure PyListAppend(L: Integer; Const V: TPyValue);
Begin
 If (L < 0) Or (L >= PY_MAX_LISTS) Then Exit;
 If PyLists[L].Count >= PY_MAX_LIST Then Begin
  PySetErr('Liste pleine'); Exit;
 End;
 PyLists[L].Items[PyLists[L].Count] := V;
 Inc(PyLists[L].Count);
End;

{ Conversions et tests }

Function PyIsNumeric(Const V: TPyValue): Boolean;
Begin
 PyIsNumeric := (V.Kind = PvBool) Or (V.Kind = PvInt) Or (V.Kind = PvFloat);
End;

Function PyToFloat(Const V: TPyValue): Real;
Begin
 PyToFloat := 0.0;
 Case V.Kind Of
  PvBool : If V.B Then PyToFloat := 1.0 Else PyToFloat := 0.0;
  PvInt  : PyToFloat := V.I;
  PvFloat: PyToFloat := V.F;
 End;
End;

Function PyToInt(Const V: TPyValue): LongInt;
Begin
 PyToInt := 0;
 Case V.Kind Of
  PvBool : If V.B Then PyToInt := 1 Else PyToInt := 0;
  PvInt  : PyToInt := V.I;
  PvFloat: PyToInt := Trunc(V.F);
 End;
End;

Function PyTruth(Const V: TPyValue): Boolean;
Begin
 PyTruth := False;
 Case V.Kind Of
  PvNone : PyTruth := False;
  PvBool : PyTruth := V.B;
  PvInt  : PyTruth := V.I <> 0;
  PvFloat: PyTruth := V.F <> 0.0;
  PvStr  : PyTruth := Length(V.S) > 0;
  PvList : If V.L >= 0 Then PyTruth := PyLists[V.L].Count > 0;
 End;
End;

Function PyTypeName(Const V: TPyValue): String;
Begin
 Case V.Kind Of
  PvNone : PyTypeName := 'NoneType';
  PvBool : PyTypeName := 'bool';
  PvInt  : PyTypeName := 'int';
  PvFloat: PyTypeName := 'float';
  PvStr  : PyTypeName := 'str';
  PvList : PyTypeName := 'list';
 Else PyTypeName := '?';
 End;
End;

Function PyEqVal(Const A, B: TPyValue): Boolean;
Var I: Integer;
Begin
 PyEqVal := False;
 If PyIsNumeric(A) And PyIsNumeric(B) Then
  PyEqVal := PyToFloat(A) = PyToFloat(B)
 Else If (A.Kind = PvStr) And (B.Kind = PvStr) Then
  PyEqVal := A.S = B.S
 Else If (A.Kind = PvNone) And (B.Kind = PvNone) Then
  PyEqVal := True
 Else If (A.Kind = PvList) And (B.Kind = PvList) And
         (A.L >= 0) And (B.L >= 0) Then Begin
  If PyLists[A.L].Count <> PyLists[B.L].Count Then Exit;
  For I := 0 To PyLists[A.L].Count - 1 Do
   If Not PyEqVal(PyLists[A.L].Items[I], PyLists[B.L].Items[I]) Then Exit;
  PyEqVal := True;
 End;
End;

Function PyLtVal(Const A, B: TPyValue): Boolean;
Begin
 PyLtVal := False;
 If PyIsNumeric(A) And PyIsNumeric(B) Then
  PyLtVal := PyToFloat(A) < PyToFloat(B)
 Else If (A.Kind = PvStr) And (B.Kind = PvStr) Then
  PyLtVal := A.S < B.S
 Else PySetErr('Comparaison invalide: ' + PyTypeName(A) + ' / ' + PyTypeName(B));
End;

{ Representation textuelle }

Function PyFmtFloat(R: Real): String;
Var S: String; I: Integer;
Begin
 Str(R:0:6, S);
 I := Length(S);
 While (I > 1) And (S[I] = '0') Do Dec(I);
 If (I >= 1) And (S[I] = '.') Then Inc(I);
 PyFmtFloat := Copy(S, 1, I);
End;

Function PyReprVal(Const V: TPyValue): String; Forward;

Function PyStrVal(Const V: TPyValue): String;
Begin
 If V.Kind = PvStr Then PyStrVal := V.S
 Else PyStrVal := PyReprVal(V);
End;

Function PyReprVal(Const V: TPyValue): String;
Var I: Integer; S: String;
Begin
 PyReprVal := '';
 Case V.Kind Of
  PvNone : PyReprVal := 'None';
  PvBool : If V.B Then PyReprVal := 'True' Else PyReprVal := 'False';
  PvInt  : PyReprVal := IntToStr(V.I);
  PvFloat: PyReprVal := PyFmtFloat(V.F);
  PvStr  : PyReprVal := '''' + V.S + '''';
  PvList : Begin
   S := '[';
   If V.L >= 0 Then
    For I := 0 To PyLists[V.L].Count - 1 Do Begin
     If I > 0 Then S := S + ', ';
     S := S + PyReprVal(PyLists[V.L].Items[I]);
    End;
   S := S + ']';
   PyReprVal := S;
  End;
 End;
End;

{ Lexer }

Function PyDigitCh(C: Char): Boolean;
Begin PyDigitCh := (C >= '0') And (C <= '9'); End;

Function PyAlphaCh(C: Char): Boolean;
Begin
 PyAlphaCh := ((C >= 'a') And (C <= 'z')) Or
              ((C >= 'A') And (C <= 'Z')) Or (C = '_');
End;

Function PyAlnumCh(C: Char): Boolean;
Begin PyAlnumCh := PyAlphaCh(C) Or PyDigitCh(C); End;

Procedure PyTokenize(Const Src: String);
Var
 I, L : Integer;
 C    : Char;
 S    : String;
 N    : LongInt;
 F    : Real;
 IsFlt: Boolean;
 Code : Integer;
 Quote: Char;
Begin
 PyTokCount := 0;
 I := 1; L := Length(Src);
 While (I <= L) And (PyTokCount < PY_MAX_TOKENS - 1) Do Begin
  C := Src[I];
  If (C = ' ') Or (C = #9) Or (C = #13) Then Inc(I)
  Else If C = #10 Then Begin
   PyTokens[PyTokCount].Kind := PtkSemi;
   Inc(PyTokCount); Inc(I);
  End
  Else If C = '#' Then Begin
   While (I <= L) And (Src[I] <> #10) Do Inc(I);
  End
  Else If PyDigitCh(C) Then Begin
   S := ''; IsFlt := False;
   While (I <= L) And (PyDigitCh(Src[I]) Or (Src[I] = '.')) Do Begin
    If Src[I] = '.' Then Begin
     If IsFlt Then Break;
     IsFlt := True;
    End;
    S := S + Src[I]; Inc(I);
   End;
   If (I <= L) And ((Src[I] = 'e') Or (Src[I] = 'E')) Then Begin
    IsFlt := True;
    S := S + Src[I]; Inc(I);
    If (I <= L) And ((Src[I] = '+') Or (Src[I] = '-')) Then Begin
     S := S + Src[I]; Inc(I);
    End;
    While (I <= L) And PyDigitCh(Src[I]) Do Begin
     S := S + Src[I]; Inc(I);
    End;
   End;
   PyTokens[PyTokCount].Kind := PtkNumber;
   PyTokens[PyTokCount].S := Copy(S, 1, 80);
   PyTokens[PyTokCount].IsFlt := IsFlt;
   If IsFlt Then Begin
    Val(S, F, Code);
    PyTokens[PyTokCount].F := F;
    PyTokens[PyTokCount].N := Trunc(F);
   End Else Begin
    Val(S, N, Code);
    PyTokens[PyTokCount].N := N;
    PyTokens[PyTokCount].F := N;
   End;
   If Code <> 0 Then PySetErr('Nombre invalide: ' + S);
   Inc(PyTokCount);
  End
  Else If (C = '"') Or (C = '''') Then Begin
   Quote := C; Inc(I); S := '';
   While (I <= L) And (Src[I] <> Quote) Do Begin
    If (Src[I] = '\') And (I < L) Then Begin
     Inc(I);
     Case Src[I] Of
      'n': S := S + #10;
      't': S := S + #9;
      'r': S := S + #13;
      '\': S := S + '\';
      '''': S := S + '''';
      '"': S := S + '"';
     Else S := S + Src[I];
     End;
     Inc(I);
    End Else Begin
     S := S + Src[I]; Inc(I);
    End;
   End;
   If I <= L Then Inc(I);
   PyTokens[PyTokCount].Kind := PtkString;
   PyTokens[PyTokCount].S := Copy(S, 1, 80);
   Inc(PyTokCount);
  End
  Else If PyAlphaCh(C) Then Begin
   S := '';
   While (I <= L) And PyAlnumCh(Src[I]) Do Begin
    S := S + Src[I]; Inc(I);
   End;
   PyTokens[PyTokCount].S := Copy(S, 1, 80);
   If      S = 'and'   Then PyTokens[PyTokCount].Kind := PtkAnd
   Else If S = 'or'    Then PyTokens[PyTokCount].Kind := PtkOr
   Else If S = 'not'   Then PyTokens[PyTokCount].Kind := PtkNot
   Else If S = 'in'    Then PyTokens[PyTokCount].Kind := PtkIn
   Else If S = 'if'    Then PyTokens[PyTokCount].Kind := PtkIf
   Else If S = 'else'  Then PyTokens[PyTokCount].Kind := PtkElse
   Else If S = 'True'  Then PyTokens[PyTokCount].Kind := PtkTrue
   Else If S = 'False' Then PyTokens[PyTokCount].Kind := PtkFalse
   Else If S = 'None'  Then PyTokens[PyTokCount].Kind := PtkNone
   Else If S = 'print' Then PyTokens[PyTokCount].Kind := PtkPrint
   Else If S = 'for'      Then PyTokens[PyTokCount].Kind := PtkFor
   Else If S = 'while'    Then PyTokens[PyTokCount].Kind := PtkWhile
   Else If S = 'import'   Then PyTokens[PyTokCount].Kind := PtkImport
   Else If S = 'break'    Then PyTokens[PyTokCount].Kind := PtkBreak
   Else If S = 'continue' Then PyTokens[PyTokCount].Kind := PtkContinue
   Else If S = 'pass'     Then PyTokens[PyTokCount].Kind := PtkPass
   Else PyTokens[PyTokCount].Kind := PtkIdent;
   Inc(PyTokCount);
  End
  Else Begin
   Case C Of
    '+': Begin PyTokens[PyTokCount].Kind := PtkPlus;   Inc(I); Inc(PyTokCount); End;
    '-': Begin PyTokens[PyTokCount].Kind := PtkMinus;  Inc(I); Inc(PyTokCount); End;
    '*': Begin
     If (I < L) And (Src[I+1] = '*') Then Begin
      PyTokens[PyTokCount].Kind := PtkDStar; Inc(I, 2);
     End Else Begin
      PyTokens[PyTokCount].Kind := PtkStar; Inc(I);
     End;
     Inc(PyTokCount);
    End;
    '/': Begin
     If (I < L) And (Src[I+1] = '/') Then Begin
      PyTokens[PyTokCount].Kind := PtkDSlash; Inc(I, 2);
     End Else Begin
      PyTokens[PyTokCount].Kind := PtkSlash; Inc(I);
     End;
     Inc(PyTokCount);
    End;
    '%': Begin PyTokens[PyTokCount].Kind := PtkPercent; Inc(I); Inc(PyTokCount); End;
    '=': Begin
     If (I < L) And (Src[I+1] = '=') Then Begin
      PyTokens[PyTokCount].Kind := PtkEqEq; Inc(I, 2);
     End Else Begin
      PyTokens[PyTokCount].Kind := PtkEq; Inc(I);
     End;
     Inc(PyTokCount);
    End;
    '!': Begin
     If (I < L) And (Src[I+1] = '=') Then Begin
      PyTokens[PyTokCount].Kind := PtkNotEq; Inc(I, 2); Inc(PyTokCount);
     End Else Begin PySetErr('! doit etre suivi de ='); Inc(I); End;
    End;
    '<': Begin
     If (I < L) And (Src[I+1] = '=') Then Begin
      PyTokens[PyTokCount].Kind := PtkLE; Inc(I, 2);
     End Else Begin
      PyTokens[PyTokCount].Kind := PtkLT; Inc(I);
     End;
     Inc(PyTokCount);
    End;
    '>': Begin
     If (I < L) And (Src[I+1] = '=') Then Begin
      PyTokens[PyTokCount].Kind := PtkGE; Inc(I, 2);
     End Else Begin
      PyTokens[PyTokCount].Kind := PtkGT; Inc(I);
     End;
     Inc(PyTokCount);
    End;
    '(': Begin PyTokens[PyTokCount].Kind := PtkLPar;  Inc(I); Inc(PyTokCount); End;
    ')': Begin PyTokens[PyTokCount].Kind := PtkRPar;  Inc(I); Inc(PyTokCount); End;
    '[': Begin PyTokens[PyTokCount].Kind := PtkLBrk;  Inc(I); Inc(PyTokCount); End;
    ']': Begin PyTokens[PyTokCount].Kind := PtkRBrk;  Inc(I); Inc(PyTokCount); End;
    ',': Begin PyTokens[PyTokCount].Kind := PtkComma; Inc(I); Inc(PyTokCount); End;
    ':': Begin PyTokens[PyTokCount].Kind := PtkColon; Inc(I); Inc(PyTokCount); End;
    '.': Begin PyTokens[PyTokCount].Kind := PtkDot;   Inc(I); Inc(PyTokCount); End;
    ';': Begin PyTokens[PyTokCount].Kind := PtkSemi;  Inc(I); Inc(PyTokCount); End;
   Else
    PySetErr('Caractere inattendu: ' + C);
    Inc(I);
   End;
  End;
 End;
 PyTokens[PyTokCount].Kind := PtkEOF;
 PyTokens[PyTokCount].S := '';
 Inc(PyTokCount);
End;

{ Variables }

Function PyFindVar(Const N: String): Integer;
Var I: Integer;
Begin
 PyFindVar := -1;
 For I := 0 To PY_MAX_VARS - 1 Do
  If PyVars[I].InUse And (PyVars[I].Name = N) Then Begin
   PyFindVar := I; Exit;
  End;
End;

Procedure PySetVar(Const N: String; Const V: TPyValue);
Var I, F: Integer;
Begin
 I := PyFindVar(N);
 If I >= 0 Then Begin PyVars[I].Value := V; Exit; End;
 F := -1;
 For I := 0 To PY_MAX_VARS - 1 Do
  If Not PyVars[I].InUse Then Begin F := I; Break; End;
 If F < 0 Then Begin PySetErr('Trop de variables'); Exit; End;
 PyVars[F].InUse := True;
 PyVars[F].Name := Copy(N, 1, 24);
 PyVars[F].Value := V;
End;

Function PyGetVar(Const N: String; Var V: TPyValue): Boolean;
Var I: Integer;
Begin
 I := PyFindVar(N);
 If I >= 0 Then Begin V := PyVars[I].Value; PyGetVar := True; End
 Else PyGetVar := False;
End;

{ Operateurs arithmetiques }

Function PyAddOp(Const A, B: TPyValue): TPyValue;
Var NL, I: Integer;
Begin
 If (A.Kind = PvStr) And (B.Kind = PvStr) Then
  PyAddOp := PyMakeStr(A.S + B.S)
 Else If (A.Kind = PvList) And (B.Kind = PvList) Then Begin
  NL := PyNewList;
  If NL < 0 Then Begin PyAddOp := PyMakeNone; Exit; End;
  If A.L >= 0 Then
   For I := 0 To PyLists[A.L].Count - 1 Do PyListAppend(NL, PyLists[A.L].Items[I]);
  If B.L >= 0 Then
   For I := 0 To PyLists[B.L].Count - 1 Do PyListAppend(NL, PyLists[B.L].Items[I]);
  PyAddOp := PyMakeList(NL);
 End
 Else If PyIsNumeric(A) And PyIsNumeric(B) Then Begin
  If (A.Kind = PvFloat) Or (B.Kind = PvFloat) Then
   PyAddOp := PyMakeFloat(PyToFloat(A) + PyToFloat(B))
  Else
   PyAddOp := PyMakeInt(PyToInt(A) + PyToInt(B));
 End
 Else Begin
  PySetErr('+ non supporte: ' + PyTypeName(A) + ' / ' + PyTypeName(B));
  PyAddOp := PyMakeNone;
 End;
End;

Function PySubOp(Const A, B: TPyValue): TPyValue;
Begin
 If PyIsNumeric(A) And PyIsNumeric(B) Then Begin
  If (A.Kind = PvFloat) Or (B.Kind = PvFloat) Then
   PySubOp := PyMakeFloat(PyToFloat(A) - PyToFloat(B))
  Else
   PySubOp := PyMakeInt(PyToInt(A) - PyToInt(B));
 End
 Else Begin
  PySetErr('- non supporte'); PySubOp := PyMakeNone;
 End;
End;

Function PyMulOp(Const A, B: TPyValue): TPyValue;
Var I, K, N, NL: Integer; SS: String;
Begin
 If (A.Kind = PvStr) And PyIsNumeric(B) Then Begin
  N := PyToInt(B); SS := '';
  For I := 1 To N Do SS := SS + A.S;
  PyMulOp := PyMakeStr(SS);
 End
 Else If (A.Kind = PvList) And PyIsNumeric(B) Then Begin
  NL := PyNewList;
  If NL < 0 Then Begin PyMulOp := PyMakeNone; Exit; End;
  N := PyToInt(B);
  For I := 1 To N Do
   If A.L >= 0 Then
    For K := 0 To PyLists[A.L].Count - 1 Do
     PyListAppend(NL, PyLists[A.L].Items[K]);
  PyMulOp := PyMakeList(NL);
 End
 Else If PyIsNumeric(A) And PyIsNumeric(B) Then Begin
  If (A.Kind = PvFloat) Or (B.Kind = PvFloat) Then
   PyMulOp := PyMakeFloat(PyToFloat(A) * PyToFloat(B))
  Else
   PyMulOp := PyMakeInt(PyToInt(A) * PyToInt(B));
 End
 Else Begin
  PySetErr('* non supporte'); PyMulOp := PyMakeNone;
 End;
End;

Function PyDivOp(Const A, B: TPyValue): TPyValue;
Begin
 If Not (PyIsNumeric(A) And PyIsNumeric(B)) Then Begin
  PySetErr('/ non supporte'); PyDivOp := PyMakeNone; Exit;
 End;
 If PyToFloat(B) = 0.0 Then Begin
  PySetErr('Division par zero'); PyDivOp := PyMakeNone; Exit;
 End;
 PyDivOp := PyMakeFloat(PyToFloat(A) / PyToFloat(B));
End;

Function PyFloorDivOp(Const A, B: TPyValue): TPyValue;
Var F: Real; Adj: Integer;
Begin
 If Not (PyIsNumeric(A) And PyIsNumeric(B)) Then Begin
  PySetErr('// non supporte'); PyFloorDivOp := PyMakeNone; Exit;
 End;
 If PyToFloat(B) = 0.0 Then Begin
  PySetErr('Division par zero'); PyFloorDivOp := PyMakeNone; Exit;
 End;
 F := PyToFloat(A) / PyToFloat(B);
 If (F < 0) And (Frac(F) <> 0) Then Adj := 1 Else Adj := 0;
 If (A.Kind = PvFloat) Or (B.Kind = PvFloat) Then
  PyFloorDivOp := PyMakeFloat(Trunc(F) - Adj)
 Else
  PyFloorDivOp := PyMakeInt(Trunc(F) - Adj);
End;

Function PyModOp(Const A, B: TPyValue): TPyValue;
Var Ai, Bi: LongInt; Af, Bf, R: Real;
Begin
 If Not (PyIsNumeric(A) And PyIsNumeric(B)) Then Begin
  PySetErr('% non supporte'); PyModOp := PyMakeNone; Exit;
 End;
 If (A.Kind = PvFloat) Or (B.Kind = PvFloat) Then Begin
  Af := PyToFloat(A); Bf := PyToFloat(B);
  If Bf = 0 Then Begin PySetErr('Modulo par zero'); PyModOp := PyMakeNone; Exit; End;
  R := Af - Int(Af / Bf) * Bf;
  PyModOp := PyMakeFloat(R);
 End Else Begin
  Ai := PyToInt(A); Bi := PyToInt(B);
  If Bi = 0 Then Begin PySetErr('Modulo par zero'); PyModOp := PyMakeNone; Exit; End;
  PyModOp := PyMakeInt(Ai - (Ai Div Bi) * Bi);
 End;
End;

Function PyPowOp(Const A, B: TPyValue): TPyValue;
Var Af, R: Real; N, I: LongInt;
Begin
 If Not (PyIsNumeric(A) And PyIsNumeric(B)) Then Begin
  PySetErr('** non supporte'); PyPowOp := PyMakeNone; Exit;
 End;
 If B.Kind <> PvFloat Then Begin
  Af := PyToFloat(A); N := PyToInt(B); R := 1.0;
  If N >= 0 Then For I := 1 To N Do R := R * Af
  Else Begin
   If Af = 0 Then Begin PySetErr('0**negatif'); PyPowOp := PyMakeNone; Exit; End;
   For I := 1 To -N Do R := R / Af;
  End;
  If (A.Kind <> PvFloat) And (N >= 0) Then PyPowOp := PyMakeInt(Trunc(R))
  Else PyPowOp := PyMakeFloat(R);
 End Else Begin
  Af := PyToFloat(A);
  If Af <= 0 Then Begin
   PySetErr('base <= 0 avec exposant flottant'); PyPowOp := PyMakeNone; Exit;
  End;
  PyPowOp := PyMakeFloat(Exp(Ln(Af) * PyToFloat(B)));
 End;
End;

Function PyInList(Const Item: TPyValue; Const Container: TPyValue): Boolean;
Var I: Integer;
Begin
 PyInList := False;
 If Container.Kind = PvList Then Begin
  If Container.L < 0 Then Exit;
  For I := 0 To PyLists[Container.L].Count - 1 Do
   If PyEqVal(Item, PyLists[Container.L].Items[I]) Then Begin
    PyInList := True; Exit;
   End;
 End
 Else If (Container.Kind = PvStr) And (Item.Kind = PvStr) Then
  PyInList := Pos(Item.S, Container.S) > 0;
End;

{ Parser (descente recursive) }

Procedure PyAdvance;
Begin
 If PyTokPos < PyTokCount - 1 Then Inc(PyTokPos);
End;

Function PyPeek: TPyTokKind;
Begin
 If PyTokPos < PyTokCount Then PyPeek := PyTokens[PyTokPos].Kind
 Else PyPeek := PtkEOF;
End;

Function PyMatch(K: TPyTokKind): Boolean;
Begin
 If PyPeek = K Then Begin PyAdvance; PyMatch := True; End
 Else PyMatch := False;
End;

Function PyParseExpr: TPyValue;   Forward;
Function PyParseCall: TPyValue;   Forward;
Function PyParseUnary: TPyValue;  Forward;
Function PyCallBuiltin(Const N: String; Const Args: Array Of TPyValue;
                       NArgs: Integer): TPyValue; Forward;

Function PyParsePrimary: TPyValue;
Var V: TPyValue; NL: Integer;
Begin
 V := PyMakeNone;
 Case PyPeek Of
  PtkNumber: Begin
   If PyTokens[PyTokPos].IsFlt Then V := PyMakeFloat(PyTokens[PyTokPos].F)
   Else V := PyMakeInt(PyTokens[PyTokPos].N);
   PyAdvance;
  End;
  PtkString: Begin V := PyMakeStr(PyTokens[PyTokPos].S); PyAdvance; End;
  PtkTrue  : Begin V := PyMakeBool(True);  PyAdvance; End;
  PtkFalse : Begin V := PyMakeBool(False); PyAdvance; End;
  PtkNone  : Begin V := PyMakeNone;        PyAdvance; End;
  PtkLPar: Begin
   PyAdvance;
   V := PyParseExpr;
   If Not PyMatch(PtkRPar) Then PySetErr(') attendu');
  End;
  PtkLBrk: Begin
   PyAdvance;
   NL := PyNewList;
   V := PyMakeList(NL);
   If PyPeek <> PtkRBrk Then Begin
    PyListAppend(NL, PyParseExpr);
    While PyMatch(PtkComma) Do PyListAppend(NL, PyParseExpr);
   End;
   If Not PyMatch(PtkRBrk) Then PySetErr('] attendu');
  End;
 Else
  PySetErr('Expression attendue');
  V := PyMakeNone;
 End;
 PyParsePrimary := V;
End;

Function PyParseCall: TPyValue;
Var
 V, Idx : TPyValue;
 Name   : String;
 IsCall : Boolean;
 Args   : Array[0..PY_MAX_ARGS - 1] Of TPyValue;
 NArgs, K : Integer;
Begin
 IsCall := False;
 V := PyMakeNone;
 If PyPeek = PtkIdent Then Begin
  Name := PyTokens[PyTokPos].S;
  PyAdvance;
  While PyPeek = PtkDot Do Begin
   PyAdvance;
   If PyPeek <> PtkIdent Then Begin
    PySetErr('Identificateur attendu apres .'); Break;
   End;
   Name := Name + '.' + PyTokens[PyTokPos].S;
   PyAdvance;
  End;
  If PyPeek = PtkLPar Then Begin
   PyAdvance;
   NArgs := 0;
   If PyPeek <> PtkRPar Then Begin
    Args[NArgs] := PyParseExpr; Inc(NArgs);
    While PyMatch(PtkComma) And (NArgs < PY_MAX_ARGS) Do Begin
     Args[NArgs] := PyParseExpr; Inc(NArgs);
    End;
   End;
   If Not PyMatch(PtkRPar) Then PySetErr(') attendu');
   V := PyCallBuiltin(Name, Args, NArgs);
   IsCall := True;
  End Else Begin
   { Nom: variable ou constante (eventuellement qualifiee) }
   If Pos('.', Name) = 0 Then Begin
    If PyGetVar(Name, V) Then
     IsCall := True
    Else If Name = 'pi' Then Begin V := PyMakeFloat(3.141592653589793); IsCall := True; End
    Else If Name = 'e'  Then Begin V := PyMakeFloat(2.718281828459045); IsCall := True; End
    Else Begin PySetErr('Nom inconnu: ' + Name); V := PyMakeNone; IsCall := True; End;
   End Else Begin
    If      Name = 'math.pi'  Then V := PyMakeFloat(3.141592653589793)
    Else If Name = 'math.e'   Then V := PyMakeFloat(2.718281828459045)
    Else If Name = 'math.inf' Then V := PyMakeFloat(1.0e30)
    Else If Name = 'math.tau' Then V := PyMakeFloat(6.283185307179586)
    Else Begin PySetErr('Attribut inconnu: ' + Name); V := PyMakeNone; End;
    IsCall := True;
   End;
  End;
 End;
 If Not IsCall Then V := PyParsePrimary;
 While PyPeek = PtkLBrk Do Begin
  PyAdvance;
  Idx := PyParseExpr;
  If Not PyMatch(PtkRBrk) Then PySetErr('] attendu');
  If (V.Kind = PvList) And (V.L >= 0) And PyIsNumeric(Idx) Then Begin
   K := PyToInt(Idx);
   If K < 0 Then K := PyLists[V.L].Count + K;
   If (K < 0) Or (K >= PyLists[V.L].Count) Then Begin
    PySetErr('Index hors limites'); V := PyMakeNone;
   End Else V := PyLists[V.L].Items[K];
  End
  Else If (V.Kind = PvStr) And PyIsNumeric(Idx) Then Begin
   K := PyToInt(Idx);
   If K < 0 Then K := Length(V.S) + K;
   If (K < 0) Or (K >= Length(V.S)) Then Begin
    PySetErr('Index chaine hors limites'); V := PyMakeNone;
   End Else V := PyMakeStr(V.S[K+1]);
  End
  Else Begin
   PySetErr('Indexation non supportee'); V := PyMakeNone;
  End;
 End;
 PyParseCall := V;
End;

Function PyParsePower: TPyValue;
Var Base, Expn: TPyValue;
Begin
 Base := PyParseCall;
 If PyMatch(PtkDStar) Then Begin
  Expn := PyParseUnary;
  Base := PyPowOp(Base, Expn);
 End;
 PyParsePower := Base;
End;

Function PyParseUnary: TPyValue;
Var V: TPyValue;
Begin
 If PyMatch(PtkMinus) Then Begin
  V := PyParseUnary;
  If PyIsNumeric(V) Then Begin
   If V.Kind = PvFloat Then V := PyMakeFloat(-V.F)
   Else V := PyMakeInt(-PyToInt(V));
  End Else Begin PySetErr('- unaire non numerique'); V := PyMakeNone; End;
 End
 Else If PyMatch(PtkPlus) Then V := PyParseUnary
 Else V := PyParsePower;
 PyParseUnary := V;
End;

Function PyParseTerm: TPyValue;
Var Left, Right: TPyValue; Op: TPyTokKind;
Begin
 Left := PyParseUnary;
 While (PyPeek = PtkStar) Or (PyPeek = PtkSlash) Or
       (PyPeek = PtkDSlash) Or (PyPeek = PtkPercent) Do Begin
  Op := PyPeek; PyAdvance;
  Right := PyParseUnary;
  Case Op Of
   PtkStar    : Left := PyMulOp(Left, Right);
   PtkSlash   : Left := PyDivOp(Left, Right);
   PtkDSlash  : Left := PyFloorDivOp(Left, Right);
   PtkPercent : Left := PyModOp(Left, Right);
  End;
 End;
 PyParseTerm := Left;
End;

Function PyParseArith: TPyValue;
Var Left, Right: TPyValue; Op: TPyTokKind;
Begin
 Left := PyParseTerm;
 While (PyPeek = PtkPlus) Or (PyPeek = PtkMinus) Do Begin
  Op := PyPeek; PyAdvance;
  Right := PyParseTerm;
  If Op = PtkPlus Then Left := PyAddOp(Left, Right)
  Else Left := PySubOp(Left, Right);
 End;
 PyParseArith := Left;
End;

Function PyParseComp: TPyValue;
Var Left, Right: TPyValue; Op: TPyTokKind;
Begin
 Left := PyParseArith;
 While (PyPeek = PtkEqEq) Or (PyPeek = PtkNotEq) Or
       (PyPeek = PtkLT) Or (PyPeek = PtkLE) Or
       (PyPeek = PtkGT) Or (PyPeek = PtkGE) Or (PyPeek = PtkIn) Do Begin
  Op := PyPeek; PyAdvance;
  Right := PyParseArith;
  Case Op Of
   PtkEqEq : Left := PyMakeBool(PyEqVal(Left, Right));
   PtkNotEq: Left := PyMakeBool(Not PyEqVal(Left, Right));
   PtkLT   : Left := PyMakeBool(PyLtVal(Left, Right));
   PtkLE   : Left := PyMakeBool(PyLtVal(Left, Right) Or PyEqVal(Left, Right));
   PtkGT   : Left := PyMakeBool(PyLtVal(Right, Left));
   PtkGE   : Left := PyMakeBool(PyLtVal(Right, Left) Or PyEqVal(Left, Right));
   PtkIn   : Left := PyMakeBool(PyInList(Left, Right));
  End;
 End;
 PyParseComp := Left;
End;

Function PyParseNot: TPyValue;
Var V: TPyValue;
Begin
 If PyMatch(PtkNot) Then Begin
  V := PyParseNot;
  PyParseNot := PyMakeBool(Not PyTruth(V));
 End
 Else PyParseNot := PyParseComp;
End;

Function PyParseAnd: TPyValue;
Var Left, Right: TPyValue;
Begin
 Left := PyParseNot;
 While PyMatch(PtkAnd) Do Begin
  Right := PyParseNot;
  If PyTruth(Left) Then Left := Right;
 End;
 PyParseAnd := Left;
End;

Function PyParseOr: TPyValue;
Var Left, Right: TPyValue;
Begin
 Left := PyParseAnd;
 While PyMatch(PtkOr) Do Begin
  Right := PyParseAnd;
  If Not PyTruth(Left) Then Left := Right;
 End;
 PyParseOr := Left;
End;

Function PyParseExpr: TPyValue;
Var Cond, ThenV, ElseV: TPyValue;
Begin
 Cond := PyParseOr;
 If PyMatch(PtkIf) Then Begin
  ThenV := Cond;
  Cond := PyParseOr;
  If PyMatch(PtkElse) Then ElseV := PyParseExpr
  Else Begin PySetErr('else attendu'); ElseV := PyMakeNone; End;
  If PyTruth(Cond) Then PyParseExpr := ThenV
  Else PyParseExpr := ElseV;
 End Else PyParseExpr := Cond;
End;

{ Fonctions natives }

Function PyCallBuiltin(Const N: String; Const Args: Array Of TPyValue;
                       NArgs: Integer): TPyValue;
Var
 I, K, NL : Integer;
 Mn, Mx, Tmp : TPyValue;
 SS, S2 : String;
 F, Sm, Mult : Real;
 IV : LongInt;
 Code : Integer;
Begin
 PyCallBuiltin := PyMakeNone;
 { ----- Modules : math.* / time.* / os.* ----- }
 If (Length(N) > 5) And (Copy(N, 1, 5) = 'math.') Then Begin
  PyCallBuiltin := PyCallBuiltin(Copy(N, 6, 255), Args, NArgs); Exit;
 End;
 If N = 'time.time' Then Begin
  PyCallBuiltin := PyMakeFloat(PyTimeNow); Exit;
 End;
 If N = 'time.clock' Then Begin
  PyCallBuiltin := PyMakeFloat(PyTimeNow); Exit;
 End;
 If N = 'time.sleep' Then Begin
  If NArgs <> 1 Then Begin PySetErr('time.sleep() 1 arg'); Exit; End;
  Delay(Round(PyToFloat(Args[0]) * 1000));
  Exit;
 End;
 If N = 'os.getcwd' Then Begin
  GetDir(0, SS);
  PyCallBuiltin := PyMakeStr(SS); Exit;
 End;
 If N = 'floor' Then Begin
  If NArgs <> 1 Then Begin PySetErr('floor() 1 arg'); Exit; End;
  F := PyToFloat(Args[0]);
  If F >= 0 Then PyCallBuiltin := PyMakeInt(Trunc(F))
  Else If Frac(F) = 0 Then PyCallBuiltin := PyMakeInt(Trunc(F))
  Else PyCallBuiltin := PyMakeInt(Trunc(F) - 1);
  Exit;
 End;
 If N = 'ceil' Then Begin
  If NArgs <> 1 Then Begin PySetErr('ceil() 1 arg'); Exit; End;
  F := PyToFloat(Args[0]);
  If F <= 0 Then PyCallBuiltin := PyMakeInt(Trunc(F))
  Else If Frac(F) = 0 Then PyCallBuiltin := PyMakeInt(Trunc(F))
  Else PyCallBuiltin := PyMakeInt(Trunc(F) + 1);
  Exit;
 End;
 If N = 'fabs' Then Begin
  If NArgs <> 1 Then Begin PySetErr('fabs() 1 arg'); Exit; End;
  PyCallBuiltin := PyMakeFloat(Abs(PyToFloat(Args[0]))); Exit;
 End;
 If N = 'print' Then Begin
  SS := '';
  For I := 0 To NArgs - 1 Do Begin
   If I > 0 Then SS := SS + ' ';
   SS := SS + PyStrVal(Args[I]);
  End;
  ConsoleSay(SS); Exit;
 End;
 If N = 'len' Then Begin
  If NArgs <> 1 Then Begin PySetErr('len() 1 arg'); Exit; End;
  Case Args[0].Kind Of
   PvStr : PyCallBuiltin := PyMakeInt(Length(Args[0].S));
   PvList: If Args[0].L >= 0 Then
             PyCallBuiltin := PyMakeInt(PyLists[Args[0].L].Count)
           Else PyCallBuiltin := PyMakeInt(0);
  Else PySetErr('len() type non supporte');
  End; Exit;
 End;
 If N = 'str' Then Begin
  If NArgs <> 1 Then Begin PySetErr('str() 1 arg'); Exit; End;
  PyCallBuiltin := PyMakeStr(PyStrVal(Args[0])); Exit;
 End;
 If N = 'repr' Then Begin
  If NArgs <> 1 Then Begin PySetErr('repr() 1 arg'); Exit; End;
  PyCallBuiltin := PyMakeStr(PyReprVal(Args[0])); Exit;
 End;
 If N = 'int' Then Begin
  If NArgs <> 1 Then Begin PySetErr('int() 1 arg'); Exit; End;
  Case Args[0].Kind Of
   PvBool, PvInt, PvFloat: PyCallBuiltin := PyMakeInt(PyToInt(Args[0]));
   PvStr: Begin
    Val(Args[0].S, IV, Code);
    If Code <> 0 Then PySetErr('int(): conversion impossible')
    Else PyCallBuiltin := PyMakeInt(IV);
   End;
  Else PySetErr('int(): type invalide');
  End; Exit;
 End;
 If N = 'float' Then Begin
  If NArgs <> 1 Then Begin PySetErr('float() 1 arg'); Exit; End;
  Case Args[0].Kind Of
   PvBool, PvInt, PvFloat: PyCallBuiltin := PyMakeFloat(PyToFloat(Args[0]));
   PvStr: Begin
    Val(Args[0].S, F, Code);
    If Code <> 0 Then PySetErr('float(): conversion impossible')
    Else PyCallBuiltin := PyMakeFloat(F);
   End;
  Else PySetErr('float(): type invalide');
  End; Exit;
 End;
 If N = 'bool' Then Begin
  If NArgs <> 1 Then Begin PySetErr('bool() 1 arg'); Exit; End;
  PyCallBuiltin := PyMakeBool(PyTruth(Args[0])); Exit;
 End;
 If N = 'abs' Then Begin
  If NArgs <> 1 Then Begin PySetErr('abs() 1 arg'); Exit; End;
  If Args[0].Kind = PvFloat Then PyCallBuiltin := PyMakeFloat(Abs(Args[0].F))
  Else If PyIsNumeric(Args[0]) Then
   PyCallBuiltin := PyMakeInt(Abs(PyToInt(Args[0])))
  Else PySetErr('abs(): type invalide');
  Exit;
 End;
 If N = 'min' Then Begin
  If NArgs < 1 Then Begin PySetErr('min() >=1 arg'); Exit; End;
  If (NArgs = 1) And (Args[0].Kind = PvList) And (Args[0].L >= 0) Then Begin
   If PyLists[Args[0].L].Count = 0 Then Begin
    PySetErr('min() liste vide'); Exit;
   End;
   Mn := PyLists[Args[0].L].Items[0];
   For I := 1 To PyLists[Args[0].L].Count - 1 Do
    If PyLtVal(PyLists[Args[0].L].Items[I], Mn) Then
     Mn := PyLists[Args[0].L].Items[I];
   PyCallBuiltin := Mn;
  End Else Begin
   Mn := Args[0];
   For I := 1 To NArgs - 1 Do
    If PyLtVal(Args[I], Mn) Then Mn := Args[I];
   PyCallBuiltin := Mn;
  End;
  Exit;
 End;
 If N = 'max' Then Begin
  If NArgs < 1 Then Begin PySetErr('max() >=1 arg'); Exit; End;
  If (NArgs = 1) And (Args[0].Kind = PvList) And (Args[0].L >= 0) Then Begin
   If PyLists[Args[0].L].Count = 0 Then Begin
    PySetErr('max() liste vide'); Exit;
   End;
   Mx := PyLists[Args[0].L].Items[0];
   For I := 1 To PyLists[Args[0].L].Count - 1 Do
    If PyLtVal(Mx, PyLists[Args[0].L].Items[I]) Then
     Mx := PyLists[Args[0].L].Items[I];
   PyCallBuiltin := Mx;
  End Else Begin
   Mx := Args[0];
   For I := 1 To NArgs - 1 Do
    If PyLtVal(Mx, Args[I]) Then Mx := Args[I];
   PyCallBuiltin := Mx;
  End;
  Exit;
 End;
 If N = 'sum' Then Begin
  If NArgs <> 1 Then Begin PySetErr('sum() 1 arg list'); Exit; End;
  If (Args[0].Kind <> PvList) Or (Args[0].L < 0) Then Begin
   PySetErr('sum(): liste attendue'); Exit;
  End;
  Sm := 0;
  For I := 0 To PyLists[Args[0].L].Count - 1 Do
   Sm := Sm + PyToFloat(PyLists[Args[0].L].Items[I]);
  PyCallBuiltin := PyMakeFloat(Sm); Exit;
 End;
 If N = 'range' Then Begin
  NL := PyNewList;
  If NL < 0 Then Exit;
  Case NArgs Of
   1: Begin
    K := 0; I := PyToInt(Args[0]);
    While K < I Do Begin PyListAppend(NL, PyMakeInt(K)); Inc(K); End;
   End;
   2: Begin
    K := PyToInt(Args[0]); I := PyToInt(Args[1]);
    While K < I Do Begin PyListAppend(NL, PyMakeInt(K)); Inc(K); End;
   End;
   3: Begin
    K := PyToInt(Args[0]); I := PyToInt(Args[1]);
    Code := PyToInt(Args[2]);
    If Code = 0 Then Begin PySetErr('range(): step=0'); Exit; End;
    If Code > 0 Then While K < I Do Begin
     PyListAppend(NL, PyMakeInt(K)); K := K + Code; End
    Else While K > I Do Begin
     PyListAppend(NL, PyMakeInt(K)); K := K + Code; End;
   End;
  Else PySetErr('range() 1..3 args');
  End;
  PyCallBuiltin := PyMakeList(NL); Exit;
 End;
 If N = 'list' Then Begin
  NL := PyNewList; If NL < 0 Then Exit;
  If NArgs = 1 Then Begin
   If Args[0].Kind = PvList Then Begin
    If Args[0].L >= 0 Then
     For I := 0 To PyLists[Args[0].L].Count - 1 Do
      PyListAppend(NL, PyLists[Args[0].L].Items[I]);
   End
   Else If Args[0].Kind = PvStr Then
    For I := 1 To Length(Args[0].S) Do
     PyListAppend(NL, PyMakeStr(Args[0].S[I]));
  End;
  PyCallBuiltin := PyMakeList(NL); Exit;
 End;
 If N = 'sqrt' Then Begin
  If NArgs <> 1 Then Begin PySetErr('sqrt() 1 arg'); Exit; End;
  F := PyToFloat(Args[0]);
  If F < 0 Then Begin PySetErr('sqrt(<0)'); Exit; End;
  PyCallBuiltin := PyMakeFloat(Sqrt(F)); Exit;
 End;
 If (N = 'sin') And (NArgs = 1) Then Begin
  PyCallBuiltin := PyMakeFloat(Sin(PyToFloat(Args[0]))); Exit;
 End;
 If (N = 'cos') And (NArgs = 1) Then Begin
  PyCallBuiltin := PyMakeFloat(Cos(PyToFloat(Args[0]))); Exit;
 End;
 If (N = 'tan') And (NArgs = 1) Then Begin
  F := PyToFloat(Args[0]);
  PyCallBuiltin := PyMakeFloat(Sin(F) / Cos(F)); Exit;
 End;
 If (N = 'exp') And (NArgs = 1) Then Begin
  PyCallBuiltin := PyMakeFloat(Exp(PyToFloat(Args[0]))); Exit;
 End;
 If (N = 'log') Or (N = 'ln') Then Begin
  If NArgs <> 1 Then Begin PySetErr('log() 1 arg'); Exit; End;
  F := PyToFloat(Args[0]);
  If F <= 0 Then Begin PySetErr('log(<=0)'); Exit; End;
  PyCallBuiltin := PyMakeFloat(Ln(F)); Exit;
 End;
 If N = 'log10' Then Begin
  If NArgs <> 1 Then Begin PySetErr('log10() 1 arg'); Exit; End;
  F := PyToFloat(Args[0]);
  If F <= 0 Then Begin PySetErr('log10(<=0)'); Exit; End;
  PyCallBuiltin := PyMakeFloat(Ln(F) / Ln(10)); Exit;
 End;
 If N = 'pow' Then Begin
  If NArgs <> 2 Then Begin PySetErr('pow() 2 args'); Exit; End;
  PyCallBuiltin := PyPowOp(Args[0], Args[1]); Exit;
 End;
 If N = 'round' Then Begin
  If NArgs < 1 Then Begin PySetErr('round() >=1 arg'); Exit; End;
  F := PyToFloat(Args[0]);
  K := 0;
  If NArgs >= 2 Then K := PyToInt(Args[1]);
  Mult := 1;
  For I := 1 To Abs(K) Do Mult := Mult * 10;
  If K >= 0 Then F := Round(F * Mult) / Mult
  Else F := Round(F / Mult) * Mult;
  If K <= 0 Then PyCallBuiltin := PyMakeInt(Round(F))
  Else PyCallBuiltin := PyMakeFloat(F);
  Exit;
 End;
 If N = 'type' Then Begin
  If NArgs <> 1 Then Begin PySetErr('type() 1 arg'); Exit; End;
  PyCallBuiltin := PyMakeStr('<class ''' + PyTypeName(Args[0]) + '''>'); Exit;
 End;
 If N = 'hex' Then Begin
  If NArgs <> 1 Then Begin PySetErr('hex() 1 arg'); Exit; End;
  IV := PyToInt(Args[0]); SS := ''; K := Abs(IV);
  If K = 0 Then SS := '0';
  While K > 0 Do Begin
   If (K And 15) < 10 Then SS := Chr(Ord('0') + (K And 15)) + SS
   Else SS := Chr(Ord('a') + (K And 15) - 10) + SS;
   K := K Shr 4;
  End;
  If IV < 0 Then SS := '-0x' + SS Else SS := '0x' + SS;
  PyCallBuiltin := PyMakeStr(SS); Exit;
 End;
 If N = 'oct' Then Begin
  If NArgs <> 1 Then Begin PySetErr('oct() 1 arg'); Exit; End;
  IV := PyToInt(Args[0]); SS := ''; K := Abs(IV);
  If K = 0 Then SS := '0';
  While K > 0 Do Begin SS := Chr(Ord('0') + (K And 7)) + SS; K := K Shr 3; End;
  If IV < 0 Then SS := '-0o' + SS Else SS := '0o' + SS;
  PyCallBuiltin := PyMakeStr(SS); Exit;
 End;
 If N = 'bin' Then Begin
  If NArgs <> 1 Then Begin PySetErr('bin() 1 arg'); Exit; End;
  IV := PyToInt(Args[0]); SS := ''; K := Abs(IV);
  If K = 0 Then SS := '0';
  While K > 0 Do Begin SS := Chr(Ord('0') + (K And 1)) + SS; K := K Shr 1; End;
  If IV < 0 Then SS := '-0b' + SS Else SS := '0b' + SS;
  PyCallBuiltin := PyMakeStr(SS); Exit;
 End;
 If N = 'chr' Then Begin
  If NArgs <> 1 Then Begin PySetErr('chr() 1 arg'); Exit; End;
  PyCallBuiltin := PyMakeStr(Chr(PyToInt(Args[0]) And $FF)); Exit;
 End;
 If N = 'ord' Then Begin
  If (NArgs <> 1) Or (Args[0].Kind <> PvStr) Or (Length(Args[0].S) <> 1) Then
   Begin PySetErr('ord(): 1 caractere attendu'); Exit; End;
  PyCallBuiltin := PyMakeInt(Ord(Args[0].S[1])); Exit;
 End;
 If N = 'upper' Then Begin
  If (NArgs <> 1) Or (Args[0].Kind <> PvStr) Then
   Begin PySetErr('upper(): chaine attendue'); Exit; End;
  PyCallBuiltin := PyMakeStr(UpStr(Args[0].S)); Exit;
 End;
 If N = 'lower' Then Begin
  If (NArgs <> 1) Or (Args[0].Kind <> PvStr) Then
   Begin PySetErr('lower(): chaine attendue'); Exit; End;
  PyCallBuiltin := PyMakeStr(LoStr(Args[0].S)); Exit;
 End;
 If N = 'reversed' Then Begin
  If NArgs <> 1 Then Begin PySetErr('reversed() 1 arg'); Exit; End;
  If Args[0].Kind = PvList Then Begin
   NL := PyNewList;
   If (NL >= 0) And (Args[0].L >= 0) Then
    For I := PyLists[Args[0].L].Count - 1 DownTo 0 Do
     PyListAppend(NL, PyLists[Args[0].L].Items[I]);
   PyCallBuiltin := PyMakeList(NL);
  End
  Else If Args[0].Kind = PvStr Then Begin
   S2 := '';
   For I := Length(Args[0].S) DownTo 1 Do S2 := S2 + Args[0].S[I];
   PyCallBuiltin := PyMakeStr(S2);
  End
  Else PySetErr('reversed(): type invalide');
  Exit;
 End;
 If N = 'sorted' Then Begin
  If NArgs <> 1 Then Begin PySetErr('sorted() 1 arg'); Exit; End;
  If (Args[0].Kind <> PvList) Or (Args[0].L < 0) Then Begin
   PySetErr('sorted(): liste attendue'); Exit;
  End;
  NL := PyNewList;
  If NL >= 0 Then Begin
   For I := 0 To PyLists[Args[0].L].Count - 1 Do
    PyListAppend(NL, PyLists[Args[0].L].Items[I]);
   For I := 1 To PyLists[NL].Count - 1 Do Begin
    Tmp := PyLists[NL].Items[I]; K := I - 1;
    While (K >= 0) And PyLtVal(Tmp, PyLists[NL].Items[K]) Do Begin
     PyLists[NL].Items[K+1] := PyLists[NL].Items[K]; Dec(K);
    End;
    PyLists[NL].Items[K+1] := Tmp;
   End;
  End;
  PyCallBuiltin := PyMakeList(NL); Exit;
 End;
 If N = 'help' Then Begin
  ConsoleSay('Python builtins: print, len, str, repr, int, float, bool, abs,');
  ConsoleSay('  min, max, sum, range, list, sqrt, sin, cos, tan, exp, log,');
  ConsoleSay('  log10, pow, round, type, hex, oct, bin, chr, ord, upper,');
  ConsoleSay('  lower, reversed, sorted. Constantes: pi, e, True, False, None.');
  Exit;
 End;
 PySetErr('Fonction inconnue: ' + N);
End;

{ Execution }

Procedure PyExecStatement; Forward;

Procedure PyExecBlock;
Begin
 While (Not PyError) And (Not PyBreakSignal) And (Not PyContinueSignal) And
       (PyPeek <> PtkEOF) Do Begin
  PyExecStatement;
  If PyError Or PyBreakSignal Or PyContinueSignal Then Break;
  If Not PyMatch(PtkSemi) Then Break;
 End;
End;

Procedure PyExecFor;
Var
 LoopVar : String;
 IterV   : TPyValue;
 I, BodyStart, ItemsCount, Cnt : Integer;
Begin
 Inc(PyLoopDepth);
 If PyPeek <> PtkIdent Then Begin
  PySetErr('Identificateur attendu apres for'); Dec(PyLoopDepth); Exit;
 End;
 LoopVar := PyTokens[PyTokPos].S;
 PyAdvance;
 If Not PyMatch(PtkIn) Then Begin
  PySetErr('in attendu dans for'); Dec(PyLoopDepth); Exit;
 End;
 IterV := PyParseExpr;
 If PyError Then Begin Dec(PyLoopDepth); Exit; End;
 If Not PyMatch(PtkColon) Then Begin
  PySetErr(': attendu dans for'); Dec(PyLoopDepth); Exit;
 End;
 While PyMatch(PtkSemi) Do ;
 BodyStart := PyTokPos;
 If IterV.Kind = PvList Then Begin
  If IterV.L >= 0 Then Begin
   ItemsCount := PyLists[IterV.L].Count;
   For I := 0 To ItemsCount - 1 Do Begin
    PySetVar(LoopVar, PyLists[IterV.L].Items[I]);
    PyTokPos := BodyStart;
    PyContinueSignal := False;
    PyExecBlock;
    If PyError Or PyBreakSignal Then Break;
   End;
  End;
 End
 Else If IterV.Kind = PvStr Then Begin
  Cnt := Length(IterV.S);
  For I := 1 To Cnt Do Begin
   PySetVar(LoopVar, PyMakeStr(IterV.S[I]));
   PyTokPos := BodyStart;
   PyContinueSignal := False;
   PyExecBlock;
   If PyError Or PyBreakSignal Then Break;
  End;
 End
 Else PySetErr('for: type non iterable');
 PyBreakSignal := False;
 PyContinueSignal := False;
 PyTokPos := PyTokCount;  { saute le reste = fin du bloc }
 Dec(PyLoopDepth);
End;

Procedure PyExecWhile;
Var
 CondStart, BodyStart : Integer;
 Cond  : TPyValue;
 Guard : LongInt;
Begin
 Inc(PyLoopDepth);
 CondStart := PyTokPos;
 BodyStart := -1;
 Guard := 0;
 While True Do Begin
  Inc(Guard);
  If Guard > 200000 Then Begin
   PySetErr('while: limite d''iterations atteinte (200000)'); Break;
  End;
  PyTokPos := CondStart;
  Cond := PyParseExpr;
  If PyError Then Break;
  If Not PyMatch(PtkColon) Then Begin
   PySetErr(': attendu dans while'); Break;
  End;
  While PyMatch(PtkSemi) Do ;
  If BodyStart < 0 Then BodyStart := PyTokPos;
  If Not PyTruth(Cond) Then Break;
  PyTokPos := BodyStart;
  PyContinueSignal := False;
  PyExecBlock;
  If PyError Or PyBreakSignal Then Break;
 End;
 PyBreakSignal := False;
 PyContinueSignal := False;
 PyTokPos := PyTokCount;
 Dec(PyLoopDepth);
End;

Procedure PyExecStatement;
Var
 V    : TPyValue;
 Name : String;
 S    : String;
Begin
 If PyError Then Exit;
 If (PyPeek = PtkEOF) Or (PyPeek = PtkSemi) Then Exit;
 { pass : ne fait rien }
 If PyMatch(PtkPass) Then Exit;
 { break / continue }
 If PyMatch(PtkBreak) Then Begin
  If PyLoopDepth = 0 Then PySetErr('break hors boucle')
  Else PyBreakSignal := True;
  Exit;
 End;
 If PyMatch(PtkContinue) Then Begin
  If PyLoopDepth = 0 Then PySetErr('continue hors boucle')
  Else PyContinueSignal := True;
  Exit;
 End;
 { import MOD [, MOD]* : accepte silencieusement }
 If PyMatch(PtkImport) Then Begin
  If PyPeek <> PtkIdent Then Begin PySetErr('module attendu'); Exit; End;
  PyAdvance;
  While PyMatch(PtkComma) Do Begin
   If PyPeek <> PtkIdent Then Begin PySetErr('module attendu'); Exit; End;
   PyAdvance;
  End;
  Exit;
 End;
 { for IDENT in EXPR : BODY }
 If PyMatch(PtkFor) Then Begin PyExecFor; Exit; End;
 { while EXPR : BODY }
 If PyMatch(PtkWhile) Then Begin PyExecWhile; Exit; End;
 { print expr [, expr]*   (forme statement sans parentheses) }
 If PyPeek = PtkPrint Then Begin
  PyAdvance;
  If (PyPeek = PtkEOF) Or (PyPeek = PtkSemi) Then Begin ConsoleSay(''); Exit; End;
  S := '';
  V := PyParseExpr; S := PyStrVal(V);
  While PyMatch(PtkComma) Do Begin
   V := PyParseExpr; S := S + ' ' + PyStrVal(V);
  End;
  If Not PyError Then ConsoleSay(S);
  Exit;
 End;
 { Assignation: IDENT = expr }
 If (PyPeek = PtkIdent) And (PyTokPos + 1 < PyTokCount) And
    (PyTokens[PyTokPos + 1].Kind = PtkEq) Then Begin
  Name := PyTokens[PyTokPos].S;
  PyAdvance; PyAdvance;
  V := PyParseExpr;
  If Not PyError Then PySetVar(Name, V);
  Exit;
 End;
 { Expression simple : evalue et affiche son repr si non-None }
 V := PyParseExpr;
 If (Not PyError) And (V.Kind <> PvNone) Then
  ConsoleSay(PyReprVal(V));
End;

Function PyEval(Const Src: String): Boolean;
Begin
 PyError := False;
 PyErrorMsg := '';
 PyTokenize(Src);
 If Not PyError Then Begin
  PyTokPos := 0;
  Repeat
   PyExecStatement;
   If PyError Then Break;
   If Not PyMatch(PtkSemi) Then Break;
  Until PyPeek = PtkEOF;
 End;
 If PyError Then Begin
  ConsoleErr('Erreur: ' + PyErrorMsg);
  PyEval := False;
 End
 Else PyEval := True;
End;

{ --------------------------------------------------------------------- }
{ Interpreteur de commandes PyMOL                                       }
{ --------------------------------------------------------------------- }

{ Reconnait on/off/1/0/true/false/yes/no.
  Retourne 1 = on, 0 = off, -1 = invalide. }
Function ParseOnOff(Const S: String): Integer;
Var L: String;
Begin
 L := LoStr(TrimStr(S));
 If (L = 'on') Or (L = '1') Or (L = 'true') Or (L = 'yes') Then ParseOnOff := 1
 Else If (L = 'off') Or (L = '0') Or (L = 'false') Or (L = 'no') Then ParseOnOff := 0
 Else ParseOnOff := -1;
End;

Function OnOffStr(B: Boolean): String;
Begin
 If B Then OnOffStr := 'on' Else OnOffStr := 'off';
End;

{ Applique `set NAME, VAL`. Retourne True si la cle est connue. }
Function ApplySetting(Const Name, ValStr: String): Boolean;
Var
 N : String;
 B : Integer;
 F : Real;
 IV : LongInt;
 Code : Integer;
 Col : Byte;
Begin
 N := LoStr(Name);
 ApplySetting := True;
 If N = 'cartoon_discrete_colors' Then Begin
  B := ParseOnOff(ValStr);
  If B < 0 Then Begin ConsoleErr('set ' + N + ' attend on/off'); Exit; End;
  SetCartoonDiscreteColors := B = 1;
  ConsoleSay(' set ' + N + ' = ' + OnOffStr(SetCartoonDiscreteColors));
  NeedRedraw := True;
 End
 Else If N = 'cartoon_helix_color' Then Begin
  Col := ParseColorName(ValStr); SetCartoonHelixColor := Col;
  NeedRedraw := True; ConsoleSay(' set ' + N + ' = ' + ValStr);
 End
 Else If (N = 'cartoon_sheet_color') Or (N = 'cartoon_strand_color') Then Begin
  Col := ParseColorName(ValStr); SetCartoonSheetColor := Col;
  NeedRedraw := True; ConsoleSay(' set ' + N + ' = ' + ValStr);
 End
 Else If (N = 'cartoon_loop_color') Or (N = 'cartoon_coil_color') Then Begin
  Col := ParseColorName(ValStr); SetCartoonLoopColor := Col;
  NeedRedraw := True; ConsoleSay(' set ' + N + ' = ' + ValStr);
 End
 Else If N = 'ray_shadows' Then Begin
  B := ParseOnOff(ValStr);
  If B < 0 Then Begin ConsoleErr('set ' + N + ' attend on/off'); Exit; End;
  SetRayShadows := B = 1;
  ConsoleSay(' set ' + N + ' = ' + OnOffStr(SetRayShadows));
 End
 Else If N = 'orthoscopic' Then Begin
  B := ParseOnOff(ValStr);
  If B < 0 Then Begin ConsoleErr('set ' + N + ' attend on/off'); Exit; End;
  SetOrthoscopic := B = 1;
  ConsoleSay(' set ' + N + ' = ' + OnOffStr(SetOrthoscopic));
 End
 Else If N = 'sphere_scale' Then Begin
  If Not StrToReal(ValStr, F) Then Begin ConsoleErr('set ' + N + ' attend reel'); Exit; End;
  SetSphereScale := F; NeedRedraw := True;
  ConsoleSay(' set ' + N + ' = ' + RealToStr(F, 4));
 End
 Else If N = 'stick_radius' Then Begin
  If Not StrToReal(ValStr, F) Then Begin ConsoleErr('set ' + N + ' attend reel'); Exit; End;
  SetStickRadius := F; NeedRedraw := True;
  ConsoleSay(' set ' + N + ' = ' + RealToStr(F, 4));
 End
 Else If N = 'cartoon_transparency' Then Begin
  If Not StrToReal(ValStr, F) Then Begin ConsoleErr('set ' + N + ' attend reel'); Exit; End;
  SetCartoonTransparency := F; NeedRedraw := True;
  ConsoleSay(' set ' + N + ' = ' + RealToStr(F, 4));
 End
 Else If N = 'transparency' Then Begin
  If Not StrToReal(ValStr, F) Then Begin ConsoleErr('set ' + N + ' attend reel'); Exit; End;
  SetTransparency := F; NeedRedraw := True;
  ConsoleSay(' set ' + N + ' = ' + RealToStr(F, 4));
 End
 Else If N = 'ambient' Then Begin
  If Not StrToReal(ValStr, F) Then Begin ConsoleErr('set ' + N + ' attend reel'); Exit; End;
  SetAmbient := F; ConsoleSay(' set ' + N + ' = ' + RealToStr(F, 4));
 End
 Else If N = 'specular' Then Begin
  If Not StrToReal(ValStr, F) Then Begin ConsoleErr('set ' + N + ' attend reel'); Exit; End;
  SetSpecular := F; ConsoleSay(' set ' + N + ' = ' + RealToStr(F, 4));
 End
 Else If N = 'antialias' Then Begin
  Val(ValStr, IV, Code); If Code <> 0 Then Begin ConsoleErr('set ' + N + ' attend entier'); Exit; End;
  SetAntialias := IV; ConsoleSay(' set ' + N + ' = ' + IntToStr(IV));
 End
 Else If N = 'ray_trace_mode' Then Begin
  Val(ValStr, IV, Code); If Code <> 0 Then Begin ConsoleErr('set ' + N + ' attend entier'); Exit; End;
  SetRayTraceMode := IV; ConsoleSay(' set ' + N + ' = ' + IntToStr(IV));
 End
 Else If N = 'dash_width' Then Begin
  Val(ValStr, IV, Code); If Code <> 0 Then Begin ConsoleErr('set ' + N + ' attend entier'); Exit; End;
  SetDashWidth := IV; ConsoleSay(' set ' + N + ' = ' + IntToStr(IV));
 End
 Else If N = 'label_size' Then Begin
  Val(ValStr, IV, Code); If Code <> 0 Then Begin ConsoleErr('set ' + N + ' attend entier'); Exit; End;
  SetLabelSize := IV; ConsoleSay(' set ' + N + ' = ' + IntToStr(IV));
 End
 Else If N = 'hash_max' Then Begin
  Val(ValStr, IV, Code); If Code <> 0 Then Begin ConsoleErr('set ' + N + ' attend entier'); Exit; End;
  SetHashLimit := IV; ConsoleSay(' set ' + N + ' = ' + IntToStr(IV));
 End
 Else ApplySetting := False;
End;

Function QuerySetting(Const Name: String): Boolean;
Var N: String;
Begin
 N := LoStr(Name);
 QuerySetting := True;
 If      N = 'cartoon_discrete_colors' Then ConsoleSay(N + ' = ' + OnOffStr(SetCartoonDiscreteColors))
 Else If N = 'cartoon_helix_color'  Then ConsoleSay(N + ' = ' + IntToStr(SetCartoonHelixColor))
 Else If N = 'cartoon_sheet_color'  Then ConsoleSay(N + ' = ' + IntToStr(SetCartoonSheetColor))
 Else If N = 'cartoon_loop_color'   Then ConsoleSay(N + ' = ' + IntToStr(SetCartoonLoopColor))
 Else If N = 'ray_shadows'          Then ConsoleSay(N + ' = ' + OnOffStr(SetRayShadows))
 Else If N = 'orthoscopic'          Then ConsoleSay(N + ' = ' + OnOffStr(SetOrthoscopic))
 Else If N = 'sphere_scale'         Then ConsoleSay(N + ' = ' + RealToStr(SetSphereScale, 4))
 Else If N = 'stick_radius'         Then ConsoleSay(N + ' = ' + RealToStr(SetStickRadius, 4))
 Else If N = 'cartoon_transparency' Then ConsoleSay(N + ' = ' + RealToStr(SetCartoonTransparency, 4))
 Else If N = 'transparency'         Then ConsoleSay(N + ' = ' + RealToStr(SetTransparency, 4))
 Else If N = 'ambient'              Then ConsoleSay(N + ' = ' + RealToStr(SetAmbient, 4))
 Else If N = 'specular'             Then ConsoleSay(N + ' = ' + RealToStr(SetSpecular, 4))
 Else If N = 'antialias'            Then ConsoleSay(N + ' = ' + IntToStr(SetAntialias))
 Else If N = 'ray_trace_mode'       Then ConsoleSay(N + ' = ' + IntToStr(SetRayTraceMode))
 Else If N = 'dash_width'           Then ConsoleSay(N + ' = ' + IntToStr(SetDashWidth))
 Else If N = 'label_size'           Then ConsoleSay(N + ' = ' + IntToStr(SetLabelSize))
 Else If N = 'hash_max'             Then ConsoleSay(N + ' = ' + IntToStr(SetHashLimit))
 Else QuerySetting := False;
End;

Procedure ExecPyMOLCommand(Const Cmd: String);
Var
 W : Array[0..7] Of String;
 NW : Integer;
 C, Rest, Arg1, Arg2 : String;
 Ang : Real;
 I : Integer;
 Col : Byte;
Begin
 Tokenize(Cmd, W, NW, 8);
 If NW = 0 Then Exit;
 C := LoStr(W[0]);
 { commentaire }
 If (C = '') Or (C[1] = '#') Then Exit;
 If (C = 'quit') Or (C = 'exit') Or (C = 'q') Then Begin
  Running := False; Exit;
 End;
 If C = 'reset' Then Begin ZoomToFit; Exit; End;
 If C = 'reinitialize' Then Begin DoReinit; Exit; End;
 If C = 'help' Then Begin
  ConsoleSay('Commandes: load, save, select, color, show, hide,');
  ConsoleSay('  bg_color, rotate, zoom, reset, spin, ray, run, quit.');
  Exit;
 End;
 If (C = 'load') Or (C = 'open') Then Begin
  If NW < 2 Then Begin ConsoleErr('load FICHIER'); Exit; End;
  Arg1 := W[1];
  If LoadMolecule(Arg1, Mol) Then Begin
   ZoomToFit;
   ConsoleSay('Charge ' + Arg1 + ' (' + IntToStr(Mol.NAtoms) + ' atomes)');
  End
  Else ConsoleErr('Echec load: ' + Arg1);
  Exit;
 End;
 If C = 'save' Then Begin
  If NW < 2 Then Begin ConsoleErr('save FICHIER'); Exit; End;
  Arg1 := W[1];
  If (PosI('.MOL', UpStr(Arg1)) > 0) And (PosI('.MOL2', UpStr(Arg1)) = 0) Then Begin
   If SaveMOL(Arg1, Mol) Then ConsoleSay('save: ' + Arg1)
   Else ConsoleErr('Echec save.');
  End
  Else If SavePDB(Arg1, Mol) Then ConsoleSay('save: ' + Arg1)
  Else ConsoleErr('Echec save.');
  Exit;
 End;
 If C = 'run' Then Begin
  If NW < 2 Then Begin ConsoleErr('run FICHIER.PML'); Exit; End;
  If LoadPML(W[1]) Then ConsoleSay('Termine: ' + W[1])
  Else ConsoleErr('Echec script.');
  Exit;
 End;
 If C = '@' Then Begin
  If NW < 2 Then Exit;
  LoadPML(W[1]); Exit;
 End;
 If C = 'select' Then Begin
  Rest := '';
  If NW >= 3 Then Begin
   SelName := W[1];
   For I := 2 To NW - 1 Do Begin
    If Rest <> '' Then Rest := Rest + ' ';
    Rest := Rest + W[I];
   End;
  End
  Else If NW = 2 Then Begin
   Rest := W[1];
  End;
  I := SelectAtoms(Mol, Rest);
  ConsoleSay('Selection ' + SelName + ': ' + IntToStr(I) + ' atomes');
  Exit;
 End;
 If C = 'color' Then Begin
  If NW < 2 Then Begin ConsoleErr('color COULEUR [SELECTION]'); Exit; End;
  Col := ParseColorName(W[1]);
  Rest := '';
  For I := 2 To NW - 1 Do Begin
   If Rest <> '' Then Rest := Rest + ' ';
   Rest := Rest + W[I];
  End;
  If Rest = '' Then Rest := 'all';
  SelectAtoms(Mol, Rest);
  ColorSelected(Mol, Col);
  ConsoleSay('color ' + W[1] + ' ' + Rest);
  Exit;
 End;
 If (C = 'show') Or (C = 'hide') Then Begin
  If NW < 2 Then Begin ConsoleErr('show MODE [SELECTION]'); Exit; End;
  Arg1 := LoStr(W[1]);
  If Arg1 = 'lines' Then RenderMode := RENDER_LINES
  Else If Arg1 = 'sticks' Then RenderMode := RENDER_STICKS
  Else If Arg1 = 'spheres' Then RenderMode := RENDER_SPHERES
  Else If (Arg1 = 'ball') Or (Arg1 = 'ballstick') Or
          (Arg1 = 'ball-and-stick') Then RenderMode := RENDER_BALLSTICK
  Else If Arg1 = 'cartoon' Then RenderMode := RENDER_RIBBON
  Else If Arg1 = 'ribbon' Then RenderMode := RENDER_RIBBON
  Else If Arg1 = 'dots' Then RenderMode := RENDER_DOTS
  Else If Arg1 = 'hydrogens' Then Begin
   If C = 'show' Then ViewFlags := ViewFlags Or SHOW_HYDROGENS
   Else ViewFlags := ViewFlags And (Not SHOW_HYDROGENS);
  End
  Else If Arg1 = 'labels' Then Begin
   If C = 'show' Then ViewFlags := ViewFlags Or SHOW_LABELS
   Else ViewFlags := ViewFlags And (Not SHOW_LABELS);
  End
  Else ConsoleErr('Mode inconnu: ' + Arg1);
  Exit;
 End;
 If C = 'bg_color' Then Begin
  If NW < 2 Then Begin ConsoleErr('bg_color COULEUR'); Exit; End;
  BgColor := ParseColorName(W[1]);
  Exit;
 End;
 If C = 'rotate' Then Begin
  If NW < 3 Then Begin ConsoleErr('rotate AXE,DEG'); Exit; End;
  If Not StrToReal(W[2], Ang) Then Ang := 0;
  Arg1 := LoStr(W[1]);
  If Arg1 = 'x' Then RotX := RotX + Ang
  Else If Arg1 = 'y' Then RotY := RotY + Ang
  Else If Arg1 = 'z' Then RotZ := RotZ + Ang
  Else ConsoleErr('Axe inconnu');
  Exit;
 End;
 If C = 'zoom' Then Begin
  ZoomToFit; Exit;
 End;
 If C = 'spin' Then Begin
  If NW >= 2 Then Begin
   Arg2 := LoStr(W[1]);
   SimSpinning := (Arg2 = 'on') Or (Arg2 = '1');
  End
  Else SimSpinning := Not SimSpinning;
  Exit;
 End;
 If C = 'ray' Then Begin
  ConsoleSay('ray: rendu rapide (sans raytracing) effectue.');
  Exit;
 End;
 If C = 'cls' Then Begin
  ConsoleCount := 0; Exit;
 End;
 If C = 'echo' Then Begin
  Rest := '';
  For I := 1 To NW - 1 Do Begin
   If Rest <> '' Then Rest := Rest + ' ';
   Rest := Rest + W[I];
  End;
  ConsoleSay(Rest); Exit;
 End;
 { set NAME, VAL  /  set NAME = VAL  /  set NAME VAL  /  unset NAME }
 If (C = 'set') Or (C = 'unset') Then Begin
  If NW < 2 Then Begin
   ConsoleErr('Usage: set NAME, VAL    (ou unset NAME)'); Exit;
  End;
  Arg1 := W[1];
  Arg2 := '';
  If NW >= 3 Then Begin
   If W[2] = '=' Then Begin
    If NW >= 4 Then Arg2 := W[3];
   End Else Arg2 := W[2];
  End;
  If C = 'unset' Then Arg2 := '0';
  If Not ApplySetting(Arg1, Arg2) Then
   ConsoleErr('set: parametre inconnu: ' + Arg1);
  Exit;
 End;
 If C = 'get' Then Begin
  If NW < 2 Then Begin ConsoleErr('Usage: get NAME'); Exit; End;
  If Not QuerySetting(W[1]) Then
   ConsoleErr('get: parametre inconnu: ' + W[1]);
  Exit;
 End;
 { Fallback: tente l'interpreteur Python embarque }
 If Not PyEval(Cmd) Then ;
End;

{ --------------------------------------------------------------------- }
{ Saisie console                                                        }
{ --------------------------------------------------------------------- }

Procedure ConsoleHistoryAdd(Const S: String);
Var I: Integer;
Begin
 If S = '' Then Exit;
 If CmdHistCount >= MAX_HISTORY Then Begin
  For I := 0 To MAX_HISTORY - 2 Do CmdHist[I] := CmdHist[I+1];
  CmdHist[MAX_HISTORY - 1] := S;
 End
 Else Begin
  CmdHist[CmdHistCount] := S;
  Inc(CmdHistCount);
 End;
 CmdHistPos := CmdHistCount;
End;

Procedure ConsoleSubmit;
Var S: String; L: Integer; Echoed: String;
Begin
 S := ConsoleInput;
 ConsoleInput := '';
 ConsoleCursor := 0;
 If (S = '') And (Not PyContPending) Then Exit;
 If PyContPending Then Echoed := PROG_PROMPT_CONT + S
 Else Echoed := PROG_PROMPT + S;
 ConsoleEcho(Echoed);
 If S <> '' Then ConsoleHistoryAdd(S);
 { Detection de continuation : ligne se terminant par \ }
 L := Length(S);
 While (L > 0) And ((S[L] = ' ') Or (S[L] = #9)) Do Dec(L);
 If (L > 0) And (S[L] = '\') Then Begin
  S := Copy(S, 1, L - 1);
  If PyContBuf = '' Then PyContBuf := S
  Else PyContBuf := PyContBuf + #10 + S;
  PyContPending := True;
  Exit;
 End;
 If PyContPending Then Begin
  S := PyContBuf + #10 + S;
  PyContBuf := '';
  PyContPending := False;
 End;
 ExecPyMOLCommand(S);
End;

Procedure ConsoleKey(K: Char);
Begin
 If K = #13 Then ConsoleSubmit
 Else If K = #8 Then Begin
  If Length(ConsoleInput) > 0 Then
   Delete(ConsoleInput, Length(ConsoleInput), 1);
 End
 Else If K = #27 Then ConsoleInput := ''
 Else If (Ord(K) >= 32) And (Length(ConsoleInput) < MAX_CMDLINE_LEN) Then
  ConsoleInput := ConsoleInput + K;
End;

Procedure ConsoleSpecialKey(K: Char);
Begin
 Case K Of
  #72: Begin  { Up - historique precedent }
   If CmdHistPos > 0 Then Dec(CmdHistPos);
   If (CmdHistPos >= 0) And (CmdHistPos < CmdHistCount) Then
    ConsoleInput := CmdHist[CmdHistPos];
  End;
  #80: Begin  { Down - historique suivant }
   If CmdHistPos < CmdHistCount Then Inc(CmdHistPos);
   If CmdHistPos < CmdHistCount Then
    ConsoleInput := CmdHist[CmdHistPos]
   Else ConsoleInput := '';
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Boucle de redessin                                                    }
{ --------------------------------------------------------------------- }

Procedure RedrawAll;
Begin
 HideMouseCursor;
 DrawCanvas;
 DrawConsoleArea;
 DrawStatusBar;
 DrawMenuBar;
 If MenuActive Then DrawMenuPopup(ActiveMenu);
 ShowMouseCursor;
 NeedRedraw := False;
End;

{ --------------------------------------------------------------------- }
{ Gestion clavier                                                       }
{ --------------------------------------------------------------------- }

Procedure HandleKey;
Var K: Char;
Begin
 K := ReadKey;
 If K = #0 Then Begin
  K := ReadKey;
  { touches speciales }
  If MenuActive Then Begin
   Case K Of
    #75: Begin  { Left }
          Dec(ActiveMenu); If ActiveMenu < 0 Then ActiveMenu := NUM_MENUS - 1;
          ActiveMenuItem := 0; MenuNextItem(0); If ActiveMenuItem < 0 Then ActiveMenuItem := 0;
         End;
    #77: Begin  { Right }
          Inc(ActiveMenu); If ActiveMenu >= NUM_MENUS Then ActiveMenu := 0;
          ActiveMenuItem := 0;
         End;
    #72: MenuNextItem(-1);
    #80: MenuNextItem(1);
   End;
   NeedRedraw := True;
   Exit;
  End;
  { en mode console actif, fleches haut/bas = historique }
  If ConsoleVisible Then Begin
   ConsoleSpecialKey(K);
   NeedRedraw := True;
  End;
  Case K Of
   #59: ShowMessage('Aide', 'F1=Aide  F3=Recharger  Tab=Console  Echap=Menu');
   #61: DoReload;          { F3 }
   #62: SimSpinning := Not SimSpinning;        { F4 }
   #63: DispatchMenu(MENU_DISPLAY, 0);   { F5 }
   #68: Running := False;                      { F10 }
   #107: Running := False;                     { Alt+F4 }
   #71: ZoomToFit;                             { Home }
   #75: RotY := RotY - 5;
   #77: RotY := RotY + 5;
   #72: RotX := RotX - 5;
   #80: RotX := RotX + 5;
   #73: ViewScale := ViewScale * 1.1;          { PgUp }
   #81: ViewScale := ViewScale / 1.1;          { PgDn }
  End;
  NeedRedraw := True;
  Exit;
 End;
 { caractere normal }
 If MenuActive Then Begin
  Case K Of
   #27: Begin MenuActive := False; ActiveMenu := -1; End;
   #13: Begin
    If (ActiveMenuItem >= 0) And
       Not Menus[ActiveMenu].Items[ActiveMenuItem].Separator Then
     DispatchMenu(ActiveMenu, ActiveMenuItem);
    MenuActive := False;
   End;
  End;
  NeedRedraw := True;
  Exit;
 End;
 If K = #27 Then Begin
  MenuActive := True; ActiveMenu := 0; ActiveMenuItem := 0;
  NeedRedraw := True; Exit;
 End;
 If K = #9 Then Begin
  DispatchMenu(MENU_DISPLAY, 13);
  Exit;
 End;
 If ConsoleVisible Then Begin
  ConsoleKey(K);
  NeedRedraw := True;
  Exit;
 End;
 Case K Of
  'l','L': DoSetRenderMode(RENDER_LINES);
  's','S': DoSetRenderMode(RENDER_STICKS);
  'c','C': DoSetRenderMode(RENDER_SPHERES);
  'b','B': DoSetRenderMode(RENDER_BALLSTICK);
  'r','R': DoSetRenderMode(RENDER_RIBBON);
  'd','D': DoSetRenderMode(RENDER_DOTS);
  'a','A': ToggleFlag(SHOW_AXES, 'axes');
  'g','G': ToggleFlag(SHOW_GRID, 'grille');
  'h','H': ToggleFlag(SHOW_HYDROGENS, 'hydrogenes');
  'q','Q': Running := False;
  '+','=': ViewScale := ViewScale * 1.1;
  '-','_': ViewScale := ViewScale / 1.1;
  '0':     ZoomToFit;
  ' ':     SimSpinning := Not SimSpinning;
 End;
 NeedRedraw := True;
End;

{ --------------------------------------------------------------------- }
{ Gestion souris                                                        }
{ --------------------------------------------------------------------- }

Procedure HandleMouse;
Var
 Hit, ItemHit : Integer;
 Dx, Dy : Integer;
 LeftDown, RightDown, MidDown, OldLeft: Boolean;
Begin
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 LeftDown  := (MouseButtons And 1) <> 0;
 RightDown := (MouseButtons And 2) <> 0;
 MidDown   := (MouseButtons And 4) <> 0;
 OldLeft   := (OldMouseButtons And 1) <> 0;
 { Survol : menu ouvert -> changement automatique de menu / item }
 If MenuActive And ((MouseX <> OldMouseX) Or (MouseY <> OldMouseY)) Then Begin
  Hit := MenuBarHitTest(MouseX, MouseY);
  If (Hit >= 0) And (Hit <> ActiveMenu) Then Begin
   ActiveMenu := Hit; ActiveMenuItem := -1;
   NeedRedraw := True;
  End
  Else Begin
   ItemHit := MenuItemHitTest(MouseX, MouseY, ActiveMenu);
   If (ItemHit >= 0) And (ItemHit <> ActiveMenuItem) And
      (Not Menus[ActiveMenu].Items[ItemHit].Separator) Then Begin
    ActiveMenuItem := ItemHit;
    NeedRedraw := True;
   End;
  End;
 End;
 { Front bouton gauche }
 If LeftDown And (Not OldLeft) Then Begin
  Hit := MenuBarHitTest(MouseX, MouseY);
  If Hit >= 0 Then Begin
   If MenuActive And (ActiveMenu = Hit) Then Begin
    MenuActive := False; ActiveMenu := -1;
   End
   Else Begin
    MenuActive := True; ActiveMenu := Hit; ActiveMenuItem := -1;
   End;
   NeedRedraw := True;
  End
  Else If MenuActive Then Begin
   ItemHit := MenuItemHitTest(MouseX, MouseY, ActiveMenu);
   If ItemHit >= 0 Then Begin
    If Not Menus[ActiveMenu].Items[ItemHit].Separator Then
     DispatchMenu(ActiveMenu, ItemHit);
    MenuActive := False; ActiveMenu := -1;
   End
   Else Begin
    MenuActive := False; ActiveMenu := -1;
   End;
   NeedRedraw := True;
  End
  Else If (MouseY > MenuBarHeight) And (MouseY < CanvasBottom) Then Begin
   DragMode := 1;
   DragStartX := MouseX; DragStartY := MouseY;
   DragStartRotX := RotX; DragStartRotY := RotY;
  End;
 End;
 { Front bouton droit }
 If RightDown And ((OldMouseButtons And 2) = 0) Then Begin
  If (MouseY > MenuBarHeight) And (MouseY < CanvasBottom) Then Begin
   DragMode := 2;
   DragStartX := MouseX; DragStartY := MouseY;
   DragStartPanX := ViewPanX; DragStartPanY := ViewPanY;
  End;
 End;
 { Front bouton central }
 If MidDown And ((OldMouseButtons And 4) = 0) Then Begin
  DragMode := 3;
  DragStartX := MouseX; DragStartY := MouseY;
  DragStartScale := ViewScale;
 End;
 { Drag continu }
 If DragMode <> 0 Then Begin
  Dx := MouseX - DragStartX;
  Dy := MouseY - DragStartY;
  Case DragMode Of
   1: Begin
    RotY := DragStartRotY + Dx * 0.5;
    RotX := DragStartRotX + Dy * 0.5;
    NeedRedraw := True;
   End;
   2: Begin
    ViewPanX := DragStartPanX + Dx;
    ViewPanY := DragStartPanY + Dy;
    NeedRedraw := True;
   End;
   3: Begin
    If DragStartScale > 0 Then Begin
     ViewScale := DragStartScale * (1.0 + Dy * 0.01);
     If ViewScale < 0.5 Then ViewScale := 0.5;
     NeedRedraw := True;
    End;
   End;
  End;
 End;
 { Relachement }
 If (MouseButtons = 0) And (OldMouseButtons <> 0) Then DragMode := 0;
 OldMouseButtons := MouseButtons;
 OldMouseX := MouseX;
 OldMouseY := MouseY;
End;

{ --------------------------------------------------------------------- }
{ Initialisation graphique et application                               }
{ --------------------------------------------------------------------- }

Procedure InitGraphicsMode;
Var
 Driver, Mode : Integer;
 Err : Integer;
Begin
 {$IFDEF FPC}
  Driver := VGA;
  Mode := VGAHi;
 {$ELSE}
  Driver := Detect;
 {$ENDIF}
 InitGraph(Driver, Mode, '');
 If GraphResult<>grOk Then Begin
  WriteLn('Erreur d''initialisation graphique!');
  Halt(1);
 End;
 ScreenWidth := GetMaxX + 1;
 ScreenHeight := GetMaxY + 1;
End;

Procedure InitializeAll;
Begin
 ClearMolecule(Mol);
 RotX := 0; RotY := 0; RotZ := 0;
 ViewScale := 30.0;
 ViewPanX := 0; ViewPanY := 0;
 ViewDistance := 30.0;
 BgColor := Black;
 RenderMode := RENDER_BALLSTICK;
 ViewFlags := SHOW_AXES Or SHOW_INFO;
 SimSpinning := False;
 SpinSpeed := 1.5;
 SelName := 'sele';

 { Parametres PyMOL par defaut }
 SetCartoonDiscreteColors := False;
 SetCartoonHelixColor     := LightRed;
 SetCartoonSheetColor     := Yellow;
 SetCartoonLoopColor      := LightCyan;
 SetRayShadows            := True;
 SetOrthoscopic           := False;
 SetSphereScale           := 1.0;
 SetStickRadius           := 0.25;
 SetCartoonTransparency   := 0.0;
 SetTransparency          := 0.0;
 SetAmbient               := 0.10;
 SetSpecular              := 1.0;
 SetAntialias             := 1;
 SetRayTraceMode          := 0;
 SetDashWidth             := 2;
 SetLabelSize             := 14;
 SetHashLimit             := 100000;

 ConsoleCount := 0;
 ConsoleView := 0;
 ConsoleInput := '';
 ConsoleCursor := 0;
 ConsoleVisible := True;
 CmdHistCount := 0;
 CmdHistPos := 0;

 MenuActive := False;
 ActiveMenu := -1;
 ActiveMenuItem := -1;

 OldMouseButtons := 0;
 DragMode := 0;

 InitMenus;
 LayoutMenuBar;
 PyInit;

 CanvasTop := MenuBarHeight;
 If ConsoleVisible Then
  CanvasBottom := ScreenHeight - STATUSBAR_HEIGHT - CONSOLE_HEIGHT - 1
 Else
  CanvasBottom := ScreenHeight - STATUSBAR_HEIGHT - 1;

 InitMouseDriver;

 ConsoleSay(PROG_TITLE);
 ConsoleSay('Tapez "help" ou utilisez les menus.');
 ConsoleSay('Echap = barre de menus.  Tab = console on/off.');

 Running := True;
 NeedRedraw := True;
End;

Procedure CleanupAll;
Begin
 CloseGraph;
End;

Procedure MainLoop;
Var
 LastTick : LongInt;
 Now      : LongInt;
 SecH, SecL, SecS, SecMs : Word;
Begin
 GetTime(SecH, SecL, SecS, SecMs);
 LastTick := LongInt(SecH) * 3600000 + LongInt(SecL) * 60000 +
             LongInt(SecS) * 1000 + LongInt(SecMs) * 10;
 If LastTick < 0 Then LastTick := 0;
 RedrawAll;
 While Running Do Begin
  If KeyPressed Then Begin
   HandleKey;
  End;
  HandleMouse;
  GetTime(SecH, SecL, SecS, SecMs);
  Now := LongInt(SecH) * 3600000 + LongInt(SecL) * 60000 +
         LongInt(SecS) * 1000 + LongInt(SecMs) * 10;
  If SimSpinning Then Begin
   If (Now - LastTick) > 40 Then Begin
    RotY := RotY + SpinSpeed;
    LastTick := Now;
    NeedRedraw := True;
   End;
  End
  Else LastTick := Now;
  If NeedRedraw Then RedrawAll;
  Delay(15);
 End;
End;

{ ===================================================================== }
{ Programme principal                                                   }
{ ===================================================================== }

Var
 InitFile : String;

Begin
 InitGraphicsMode;
 InitializeAll;
 If ParamCount >= 1 Then Begin
  InitFile := ParamStr(1);
  If (PosI('.PML', UpStr(InitFile)) > 0) Then Begin
   If Not LoadPML(InitFile) Then
    ConsoleErr('Script introuvable: ' + InitFile);
  End
  Else Begin
   If LoadMolecule(InitFile, Mol) Then Begin
    ZoomToFit;
    ConsoleSay('Charge: ' + InitFile);
   End
   Else ConsoleErr('Echec chargement: ' + InitFile);
  End;
 End;
 MainLoop;
 CleanupAll;
End.
