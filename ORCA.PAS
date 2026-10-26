{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2026
  @version: 1.00
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal 7, Free Pascal)
  @description: Clone de l'application ORCA - logiciel de chimie quantique
                (calculs electroniques sur des molecules) avec editeur de
                fichier d'entree .INP, moteur de calcul integre, fichier de
                sortie .OUT et visualisation graphique des resultats.

  Ligne de commande :
    ORCA [fichier.INP|.XYZ]
    ORCA /B fichier.INP [sortie.OUT]   (mode lot : calcul sans interface)

  Fonctionnalites :
    - Fichier d'entree au format ORCA : lignes "! mots-cles", blocs
      "%scf ... end", "%geom ... end", "%tddft ... end", "%freq ... end",
      coordonnees "* xyz charge multiplicite ... *", "* int" (matrice Z),
      "* xyzfile charge multiplicite fichier.xyz"
    - Methodes : HF / RHF / UHF, MP2 (coeur gele), CIS (etats excites),
      DFT (LDA/SVWN, BLYP, B3LYP, PBE, PBE0) sur grille de Becke
    - Bases gaussiennes : STO-3G, 3-21G, 6-31G (H a Ne selon la base)
    - Taches : energie (SP), optimisation (Opt, TightOpt), frequences
      (Freq / NumFreq / OptFreq), thermochimie, balayage relaxe (%geom Scan),
      contraintes (%geom Constraints)
    - Analyses : charges de Mulliken et de Loewdin, ordres de liaison de
      Mayer, moment dipolaire, energies orbitalaires (HOMO / LUMO),
      spectre IR, spectre UV (CIS), coupes d'orbitales moleculaires
    - Visualisation 3D de la molecule (fils, batonnets, boules et batons,
      spheres), graphes de convergence SCF / optimisation / balayage
    - Menus deroulants adaptes de CLONEDRAW\PAINT.PAS

  Aucune dependance externe : tout le code est dans ce fichier.
}

Program ORCA;

{$IFDEF FPC}
 {$MODE TP}
 {$Q-}
 {$R-}
{$ELSE}
 {$M 60000,0,655360}
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
 PROG_TITLE      = 'ORCA Clone 1.00 - QChimie';
 PROG_VERSION    = '1.00';

 MENUBAR_H       = 20;
 TOOLBAR_H       = 24;
 STATUS_H        = 18;

 NUM_MENUS       = 8;
 MENU_FILE       = 0;
 MENU_EDIT       = 1;
 MENU_CALC       = 2;
 MENU_RES        = 3;
 MENU_VIEW       = 4;
 MENU_EXAMPLE    = 5;
 MENU_TOOLS      = 6;
 MENU_HELP       = 7;

 {$IFDEF FPC}
 MAX_BF          = 40;           { fonctions de base }
 MAX_AT          = 32;           { atomes }
 MAX_CHUNK       = 84;           { blocs d'integrales bielectroniques }
 MAX_OUT         = 20000;         { lignes de sortie conservees }
 MAX_INP         = 400;          { lignes du fichier d'entree }
 MAX_MENU_ITEMS  = 22;
 MAX_HIST        = 300;
 CIS_MAX         = 100;          { dimension maximale de la matrice CIS }
 {$ELSE}
 MAX_BF          = 14;
 MAX_AT          = 8;
 MAX_CHUNK       = 3;
 MAX_OUT         = 400;
 MAX_INP         = 80;
 MAX_MENU_ITEMS  = 18;
 MAX_HIST        = 150;
 CIS_MAX         = 30;
 {$ENDIF}

 MAX_PAIRS       = MAX_BF * (MAX_BF + 1) Div 2;
 CHUNK_SIZE      = 4096;
 MAX_PRIM        = 6;
 MAX_SHELL       = 8;
 MAX_LINE        = 100;          { longueur d'une ligne de sortie }
 MAX_CONS        = 3 * MAX_AT;   { coordonnees cartesiennes }

 ANG2BOHR        = 1.8897261246;
 HARTREE_EV      = 27.211386245;
 HARTREE_KCAL    = 627.5094740;
 HARTREE_KJ      = 2625.499638;
 HARTREE_CM      = 219474.6313;
 BOHR_DEBYE      = 2.541746473;   { u.a. -> Debye }
 FREQ_CONV       = 5140.4871;     { sqrt(Eh/(bohr^2 amu)) -> cm-1 }
 IR_CONV         = 974.8802;      { (u.a./(bohr sqrt(amu)))^2 -> km/mol }
 KB_HARTREE      = 3.166811563E-6; { constante de Boltzmann en Eh/K }

 { Methodes }
 M_HF            = 0;
 M_MP2           = 1;
 M_DFT           = 2;

 { Fonctionnelles }
 F_NONE          = 0;
 F_LDA           = 1;
 F_BLYP          = 2;
 F_B3LYP         = 3;
 F_PBE           = 4;
 F_PBE0          = 5;

 { Bases }
 B_STO3G         = 0;
 B_321G          = 1;
 B_631G          = 2;

 { Taches }
 J_SP            = 0;
 J_OPT           = 1;
 J_FREQ          = 2;
 J_OPTFREQ       = 3;
 J_SCAN          = 4;

{ --------------------------------------------------------------------- }
{ Types                                                                 }
{ --------------------------------------------------------------------- }

Type
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

 TVec3 = Array[1..3] Of Real;

 TMat  = Array[1..MAX_BF, 1..MAX_BF] Of Real;
 PMat  = ^TMat;
 TVecN = Array[1..MAX_BF] Of Real;

 TPrimPair = Record
  P          : Real;
  PX, PY, PZ : Real;
  Cm         : Real;
  E          : Array[0..26] Of Real;
 End;
 TPPArr = Array[1..40] Of TPrimPair;
 PPPArr = ^TPPArr;

 TChunk = Array[0..CHUNK_SIZE-1] Of Real;
 PChunk = ^TChunk;

 TLine  = String[MAX_LINE];
 PLine  = ^TLine;

 TCartVec = Array[1..MAX_CONS] Of Real;
 TCartMat = Array[1..MAX_CONS, 1..MAX_CONS] Of Real;
 PCartMat = ^TCartMat;

{ --------------------------------------------------------------------- }
{ Variables globales                                                    }
{ --------------------------------------------------------------------- }

Var
 { Divers }
 StatusMsg   : String;
 NeedRedraw  : Boolean;
 Running     : Boolean;
 BatchMode   : Boolean;
 UserAbort   : Boolean;
 GuiActive   : Boolean;
 ProgPct     : Integer;
 LastShowT   : LongInt;

 { Sortie }
 OutLines    : Array[0..MAX_OUT-1] Of PLine;
 NOut        : LongInt;
 OutOverflow : Boolean;
 OutFile     : Text;
 OutFileOpen : Boolean;

{ Avancement et sortie (definis avec l'interface) }
Procedure ShowProgress(Const Msg: String; Pct: Integer); Forward;
Procedure ConsoleUpdate; Forward;

{ --------------------------------------------------------------------- }
{ Utilitaires                                                          }
{ --------------------------------------------------------------------- }

Function UpStr(Const S: String): String;
Var I: Integer; R: String;
Begin
 R := S;
 For I := 1 To Length(R) Do
  If (R[I] >= 'a') And (R[I] <= 'z') Then R[I] := Chr(Ord(R[I]) - 32);
 UpStr := R;
End;

Function LoStr(Const S: String): String;
Var I: Integer; R: String;
Begin
 R := S;
 For I := 1 To Length(R) Do
  If (R[I] >= 'A') And (R[I] <= 'Z') Then R[I] := Chr(Ord(R[I]) + 32);
 LoStr := R;
End;

Function TrimStr(Const S: String): String;
Var A, B: Integer;
Begin
 A := 1; B := Length(S);
 While (A <= B) And (S[A] <= ' ') Do Inc(A);
 While (B >= A) And (S[B] <= ' ') Do Dec(B);
 TrimStr := Copy(S, A, B - A + 1);
End;

Function IStr(V: LongInt): String;
Var S: String;
Begin
 Str(V, S);
 IStr := S;
End;

Function RStr(V: Real; Dec: Integer): String;
Var S: String;
Begin
 If Abs(V) >= 1E9 Then Begin Str(V:0:0, S); RStr := S; Exit; End;
 Str(V:0:Dec, S);
 RStr := S;
End;

{ Nombre reel sans zeros inutiles }
Function RTrim(V: Real; Dec: Integer): String;
Var S: String;
Begin
 S := RStr(V, Dec);
 If Pos('.', S) > 0 Then Begin
  While (Length(S) > 1) And (S[Length(S)] = '0') Do Delete(S, Length(S), 1);
  If (Length(S) > 1) And (S[Length(S)] = '.') Then Delete(S, Length(S), 1);
 End;
 RTrim := S;
End;

Function PadR(Const S: String; W: Integer): String;
Var R: String;
Begin
 R := S;
 While Length(R) < W Do R := R + ' ';
 PadR := Copy(R, 1, W);
End;

Function PadL(Const S: String; W: Integer): String;
Var R: String;
Begin
 R := S;
 While Length(R) < W Do R := ' ' + R;
 PadL := R;
End;

{ Reel formate en largeur W avec Dec decimales }
Function FmtR(V: Real; W, Dec: Integer): String;
Begin
 FmtR := PadL(RStr(V, Dec), W);
End;

{ Notation scientifique "1.234E-05" }
Function SciStr(V: Real; Dec: Integer): String;
Var S: String;
Begin
 Str(V:Dec + 7, S);
 SciStr := TrimStr(S);
End;

Function ParseReal(Const S: String; Var V: Real): Boolean;
Var Code: Integer; R: Real; T: String; I: Integer;
Begin
 T := TrimStr(S);
 For I := 1 To Length(T) Do If T[I] = ',' Then T[I] := '.';
 If T = '' Then Begin ParseReal := False; Exit; End;
 Val(T, R, Code);
 If Code = 0 Then Begin V := R; ParseReal := True; End
 Else ParseReal := False;
End;

Function ParseInt(Const S: String; Var V: Integer): Boolean;
Var Code: Integer; R: LongInt; T: String;
Begin
 T := TrimStr(S);
 If T = '' Then Begin ParseInt := False; Exit; End;
 Val(T, R, Code);
 If (Code = 0) And (R >= -32000) And (R <= 32000) Then Begin V := R; ParseInt := True; End
 Else ParseInt := False;
End;

Function ExtractExt(Const FN: String): String;
Var I: Integer;
Begin
 I := Length(FN);
 While (I > 0) And (FN[I] <> '.') And (FN[I] <> '\') And (FN[I] <> '/') Do Dec(I);
 If (I > 0) And (FN[I] = '.') Then ExtractExt := UpStr(Copy(FN, I, 255))
 Else ExtractExt := '';
End;

Function ChangeExt(Const FN, Ext: String): String;
Var I: Integer;
Begin
 I := Length(FN);
 While (I > 0) And (FN[I] <> '.') And (FN[I] <> '\') And (FN[I] <> '/') Do Dec(I);
 If (I > 0) And (FN[I] = '.') Then ChangeExt := Copy(FN, 1, I - 1) + Ext
 Else ChangeExt := FN + Ext;
End;

Function BaseName(Const FN: String): String;
Var I: Integer;
Begin
 I := Length(FN);
 While (I > 0) And (FN[I] <> '\') And (FN[I] <> '/') Do Dec(I);
 BaseName := Copy(FN, I + 1, 255);
End;

Function FileOK(Const FN: String): Boolean;
Var F: File;
Begin
 {$I-}
 Assign(F, FN);
 Reset(F, 1);
 {$I+}
 If IOResult = 0 Then Begin
  Close(F);
  FileOK := True;
 End
 Else FileOK := False;
End;

Function Pw(X, Y: Real): Real;
Begin
 If X <= 0 Then Pw := 0 Else Pw := Exp(Y * Ln(X));
End;

Function RMin(A, B: Real): Real;
Begin
 If A < B Then RMin := A Else RMin := B;
End;

Function RMax(A, B: Real): Real;
Begin
 If A > B Then RMax := A Else RMax := B;
End;

Function IMin(A, B: LongInt): LongInt;
Begin
 If A < B Then IMin := A Else IMin := B;
End;

Function IMax(A, B: LongInt): LongInt;
Begin
 If A > B Then IMax := A Else IMax := B;
End;

Function ExpS(X: Real): Real;
Begin
 If X < -80 Then ExpS := 0 Else ExpS := Exp(X);
End;

Function ArcSinh(X: Real): Real;
Begin
 If X >= 0 Then ArcSinh := Ln(X + Sqrt(X * X + 1))
 Else ArcSinh := -Ln(-X + Sqrt(X * X + 1));
End;

Function ArcCos(X: Real): Real;
Begin
 If X >= 1 Then ArcCos := 0
 Else If X <= -1 Then ArcCos := Pi
 Else If X = 0 Then ArcCos := Pi / 2
 Else If X > 0 Then ArcCos := ArcTan(Sqrt(1 - X * X) / X)
 Else ArcCos := Pi + ArcTan(Sqrt(1 - X * X) / X);
End;

Function ISqr(I: LongInt): LongInt;
Begin
 ISqr := I * I;
End;

{ Horloge en millisecondes depuis minuit }
Function NowMs: LongInt;
Var H, M, S, S100: Word;
Begin
 GetTime(H, M, S, S100);
 NowMs := LongInt(H) * 3600000 + LongInt(M) * 60000 +
          LongInt(S) * 1000 + LongInt(S100) * 10;
End;

{ Duree ecoulee depuis T0 (en ms) en tenant compte du passage a minuit }
Function ElapsedMs(T0: LongInt): LongInt;
Var D: LongInt;
Begin
 D := NowMs - T0;
 If D < 0 Then D := D + 86400000;
 ElapsedMs := D;
End;

{ --------------------------------------------------------------------- }
{ Gestion de la souris                                                  }
{ --------------------------------------------------------------------- }

Var
 MouseX, MouseY, MouseButtons : Integer;
 OldMouseX, OldMouseY, OldMouseButtons : Integer;
 MouseAvailable : Boolean;

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
{ Tampon de sortie (fichier .OUT)                                       }
{ --------------------------------------------------------------------- }

Procedure OutClear;
Var I: LongInt;
Begin
 For I := 0 To NOut - 1 Do
  If OutLines[I] <> Nil Then Begin Dispose(OutLines[I]); OutLines[I] := Nil; End;
 NOut := 0;
 OutOverflow := False;
End;

Procedure Emit(Const S: String);
Var P: PLine;
Begin
 If OutFileOpen Then WriteLn(OutFile, S);
 If NOut < MAX_OUT Then Begin
  New(P);
  P^ := Copy(S, 1, MAX_LINE);
  OutLines[NOut] := P;
  Inc(NOut);
 End
 Else OutOverflow := True;
End;

Procedure EmitBlank;
Begin
 Emit('');
End;

Procedure EmitHeader(Const T: String);
Var L: String; I: Integer;
Begin
 EmitBlank;
 Emit(T);
 L := '';
 For I := 1 To Length(T) Do L := L + '-';
 Emit(L);
End;

Procedure EmitBar(Const T: String);
Var L: String; I: Integer;
Begin
 L := '';
 For I := 1 To 72 Do L := L + '-';
 Emit(L);
 Emit(T);
 Emit(L);
End;

{ --------------------------------------------------------------------- }
{ Elements chimiques                                                    }
{ --------------------------------------------------------------------- }

Const
 NELEM = 18;
 ElSym : Array[1..NELEM] Of String[2] = (
  'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
  'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar');
 ElName : Array[1..NELEM] Of String[12] = (
  'Hydrogene', 'Helium', 'Lithium', 'Beryllium', 'Bore', 'Carbone',
  'Azote', 'Oxygene', 'Fluor', 'Neon', 'Sodium', 'Magnesium',
  'Aluminium', 'Silicium', 'Phosphore', 'Soufre', 'Chlore', 'Argon');
 { masses isotopiques (isotope le plus abondant, comme ORCA) }
 ElMass : Array[1..NELEM] Of Real = (
  1.007825, 4.002603, 7.016004, 9.012182, 11.009305, 12.000000,
  14.003074, 15.994915, 18.998403, 19.992440, 22.989770, 23.985042,
  26.981538, 27.976927, 30.973762, 31.972071, 34.968853, 39.962383);
 { rayons covalents (angstroms) }
 ElRad : Array[1..NELEM] Of Real = (
  0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58,
  1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06);
 { couleurs d'affichage (palette VGA) }
 ElColor : Array[1..NELEM] Of Byte = (
  15, 11, 13, 2, 6, 8, 9, 12, 10, 11, 13, 10, 7, 6, 6, 14, 2, 11);
 { zeta de Slater des fonctions STO-3G (couche 1s puis couche 2sp) }
 Sto3Z1 : Array[1..10] Of Real = (
  1.24, 1.69, 2.69, 3.68, 4.68, 5.67, 6.67, 7.66, 8.65, 9.64);
 Sto3Z2 : Array[3..10] Of Real = (
  0.80, 1.15, 1.50, 1.72, 1.95, 2.25, 2.55, 2.88);

 BasisNames : Array[0..2] Of String[8] = ('STO-3G', '3-21G', '6-31G');

Function ElemZ(Const Sym: String): Integer;
Var I: Integer; S: String;
Begin
 S := UpStr(TrimStr(Sym));
 For I := 1 To NELEM Do
  If UpStr(ElSym[I]) = S Then Begin ElemZ := I; Exit; End;
 ElemZ := 0;
End;

Function ElemStr(Z: Integer): String;
Begin
 If (Z >= 1) And (Z <= NELEM) Then ElemStr := ElSym[Z] Else ElemStr := 'X';
End;

{ --------------------------------------------------------------------- }
{ Molecule                                                              }
{ --------------------------------------------------------------------- }

Var
 NAt     : Integer;
 AtZ     : Array[1..MAX_AT] Of Integer;
 AtX     : Array[1..MAX_AT, 1..3] Of Real;     { bohr }
 Charge  : Integer;
 Mult    : Integer;
 NElec   : Integer;
 MolTitle : String;

Function Dist(A, B: Integer): Real;
Begin
 Dist := Sqrt(Sqr(AtX[A,1] - AtX[B,1]) + Sqr(AtX[A,2] - AtX[B,2]) +
              Sqr(AtX[A,3] - AtX[B,3]));
End;

Function NucRepulsion: Real;
Var A, B: Integer; E: Real;
Begin
 E := 0;
 For A := 1 To NAt - 1 Do
  For B := A + 1 To NAt Do E := E + AtZ[A] * AtZ[B] / Dist(A, B);
 NucRepulsion := E;
End;

Function Bonded(A, B: Integer): Boolean;
Var R: Real;
Begin
 R := (ElRad[AtZ[A]] + ElRad[AtZ[B]]) * ANG2BOHR;
 Bonded := Dist(A, B) < 1.25 * R;
End;

{ --------------------------------------------------------------------- }
{ Bases gaussiennes                                                      }
{ --------------------------------------------------------------------- }

Var
 BasisSel : Integer;
 { coquilles temporaires de l'element en cours }
 NTS : Integer;
 TSL : Array[1..MAX_SHELL] Of Integer;
 TSN : Array[1..MAX_SHELL] Of Integer;
 TSA, TSC : Array[1..MAX_SHELL, 1..MAX_PRIM] Of Real;

 { fonctions de base }
 NBF     : Integer;
 BFAtom  : Array[1..MAX_BF] Of Integer;
 BFL     : Array[1..MAX_BF, 1..3] Of Integer;
 BFNP    : Array[1..MAX_BF] Of Integer;
 BFA, BFC : Array[1..MAX_BF, 1..MAX_PRIM] Of Real;
 BFLab   : Array[1..MAX_BF] Of String[8];

Procedure NewShell(L: Integer);
Begin
 Inc(NTS);
 TSL[NTS] := L;
 TSN[NTS] := 0;
End;

Procedure AddPrim(A, C: Real);
Begin
 Inc(TSN[NTS]);
 TSA[NTS, TSN[NTS]] := A;
 TSC[NTS, TSN[NTS]] := C;
End;

Procedure Sh1(L: Integer; A1, C1: Real);
Begin
 NewShell(L); AddPrim(A1, C1);
End;

Procedure Sh2(L: Integer; A1, C1, A2, C2: Real);
Begin
 NewShell(L); AddPrim(A1, C1); AddPrim(A2, C2);
End;

Procedure Sh3(L: Integer; A1, C1, A2, C2, A3, C3: Real);
Begin
 NewShell(L); AddPrim(A1, C1); AddPrim(A2, C2); AddPrim(A3, C3);
End;

Procedure Sh6(L: Integer; A1, C1, A2, C2, A3, C3, A4, C4, A5, C5, A6, C6: Real);
Begin
 NewShell(L);
 AddPrim(A1, C1); AddPrim(A2, C2); AddPrim(A3, C3);
 AddPrim(A4, C4); AddPrim(A5, C5); AddPrim(A6, C6);
End;

{ Coquille SP a trois primitives : une coquille s et une coquille p }
Procedure SP3(A1, S1, P1, A2, S2, P2, A3, S3, P3: Real);
Begin
 Sh3(0, A1, S1, A2, S2, A3, S3);
 Sh3(1, A1, P1, A2, P2, A3, P3);
End;

{ Coquille SP a deux primitives }
Procedure SP2(A1, S1, P1, A2, S2, P2: Real);
Begin
 Sh2(0, A1, S1, A2, S2);
 Sh2(1, A1, P1, A2, P2);
End;

Procedure SP1(A1: Real);
Begin
 Sh1(0, A1, 1.0);
 Sh1(1, A1, 1.0);
End;

Procedure Sto3Shells(Z: Integer);
Var Z1, Z2, F1, F2: Real;
Begin
 Z1 := Sto3Z1[Z];
 F1 := Z1 * Z1;
 Sh3(0, 2.227660584 * F1, 0.1543289673, 0.4057711562 * F1, 0.5353281423,
        0.1098175104 * F1, 0.4446345422);
 If Z >= 3 Then Begin
  Z2 := Sto3Z2[Z];
  F2 := Z2 * Z2;
  SP3(0.9942027296 * F2, -0.09996722919, 0.1559162750,
      0.2310313333 * F2, 0.3995128261, 0.6076837186,
      0.0751386610 * F2, 0.7001154689, 0.3919573931);
 End;
End;

Procedure Shells321(Z: Integer);
Begin
 Case Z Of
  1: Begin
   Sh2(0, 5.4471780, 0.1562850, 0.8245472, 0.9046910);
   Sh1(0, 0.1831916, 1.0);
  End;
  2: Begin
   Sh2(0, 13.6267000, 0.1752300, 1.9993500, 0.8934830);
   Sh1(0, 0.3829930, 1.0);
  End;
  6: Begin
   Sh3(0, 172.2560000, 0.0617669, 25.9109000, 0.3587940, 5.5333500, 0.7007130);
   SP2(3.6649800, -0.3958970, 0.2364600, 0.7705450, 1.2158400, 0.8606190);
   SP1(0.1958570);
  End;
  7: Begin
   Sh3(0, 242.7660000, 0.0598657, 36.4851000, 0.3529550, 7.8144900, 0.7065130);
   SP2(5.4252200, -0.4133010, 0.2379720, 1.1491500, 1.2244200, 0.8589530);
   SP1(0.2832050);
  End;
  8: Begin
   Sh3(0, 322.0370000, 0.0592394, 48.4308000, 0.3515000, 10.4206000, 0.7076580);
   SP2(7.4029400, -0.4044530, 0.2445860, 1.5762000, 1.2215600, 0.8539550);
   SP1(0.3736840);
  End;
  9: Begin
   Sh3(0, 413.8010000, 0.0585483, 62.2446000, 0.3493080, 13.4340000, 0.7096320);
   SP2(9.7775900, -0.4073270, 0.2466800, 2.0861700, 1.2231400, 0.8523210);
   SP1(0.4823830);
  End;
 End;
End;

Procedure Shells631(Z: Integer);
Begin
 Case Z Of
  1: Begin
   Sh3(0, 18.7311370, 0.03349460, 2.8253937, 0.23472695, 0.6401217, 0.81375733);
   Sh1(0, 0.1612778, 1.0);
  End;
  2: Begin
   Sh3(0, 38.4216340, 0.0237660, 5.7780300, 0.1546790, 1.2417740, 0.4696300);
   Sh1(0, 0.2979640, 1.0);
  End;
  6: Begin
   Sh6(0, 3047.5249000, 0.0018347, 457.3695100, 0.0140373, 103.9486900, 0.0688426,
          29.2101550, 0.2321844, 9.2866630, 0.4679413, 3.1639270, 0.3623120);
   SP3(7.8682724, -0.1193324, 0.0689991, 1.8812885, -0.1608542, 0.3164240,
       0.5442493, 1.1434564, 0.7443083);
   SP1(0.1687144);
  End;
  7: Begin
   Sh6(0, 4173.5110000, 0.0018348, 627.4579000, 0.0139950, 142.9021000, 0.0685870,
          40.2343300, 0.2322410, 12.8202100, 0.4690700, 4.3904370, 0.3604550);
   SP3(11.6263580, -0.1149610, 0.0675800, 2.7162800, -0.1691180, 0.3239070,
       0.7722180, 1.1458520, 0.7408950);
   SP1(0.2120313);
  End;
  8: Begin
   Sh6(0, 5484.6717000, 0.0018311, 825.2349500, 0.0139501, 188.0469600, 0.0684451,
          52.9645000, 0.2327143, 16.8975700, 0.4701930, 5.7996353, 0.3585209);
   SP3(15.5396160, -0.1107775, 0.0708743, 3.5999336, -0.1480263, 0.3397528,
       1.0137618, 1.1307670, 0.7271586);
   SP1(0.2700058);
  End;
  9: Begin
   Sh6(0, 7001.7130000, 0.0018196169, 1051.3660000, 0.01391608, 239.2857000, 0.06840532,
          67.3974400, 0.2331858, 21.5199600, 0.4712674, 7.4031010, 0.3566185);
   SP3(20.8479500, -0.1085070, 0.0716287, 4.8083080, -0.1464517, 0.3459121,
       1.3440700, 1.1286886, 0.7224700);
   SP1(0.3581514);
  End;
 End;
End;

{ Remplit les coquilles de l'element Z pour la base B ; Faux si absent }
Function GetShells(B, Z: Integer): Boolean;
Begin
 NTS := 0;
 GetShells := False;
 If (Z < 1) Or (Z > NELEM) Then Exit;
 Case B Of
  B_STO3G: If Z <= 10 Then Begin Sto3Shells(Z); GetShells := True; End;
  B_321G: If (Z = 1) Or (Z = 2) Or ((Z >= 6) And (Z <= 9)) Then Begin
   Shells321(Z); GetShells := True;
  End;
  B_631G: If (Z = 1) Or (Z = 2) Or ((Z >= 6) And (Z <= 9)) Then Begin
   Shells631(Z); GetShells := True;
  End;
 End;
 If NTS = 0 Then GetShells := False;
End;

Function BasisSupports(B, Z: Integer): Boolean;
Begin
 BasisSupports := GetShells(B, Z);
End;

Function PrimNorm(A: Real; L: Integer): Real;
Begin
 PrimNorm := Exp(0.75 * Ln(2 * A / Pi)) * Exp(L / 2 * Ln(4 * A));
End;

{ --------------------------------------------------------------------- }
{ Construction de la base et integrales (schema de McMurchie-Davidson)  }
{ --------------------------------------------------------------------- }

Const
 TWO_PI_52 = 34.986836655249725;     { 2 * pi^(5/2) }

Type
 TBoys  = Array[0..12] Of Real;
 TRTab  = Array[0..4, 0..4, 0..4, 0..4] Of Real;
 PPrim  = ^TPrimPair;

Var
 RT      : TRTab;
 PairP   : Array[1..MAX_PAIRS] Of PPPArr;
 PairNP  : Array[1..MAX_PAIRS] Of Integer;
 PairL   : Array[1..MAX_PAIRS] Of Integer;
 PairI, PairJ : Array[1..MAX_PAIRS] Of Integer;
 PairQ   : Array[1..MAX_PAIRS] Of Real;
 NPairs  : Integer;
 PairsOK : Boolean;
 ChunkP  : Array[0..MAX_CHUNK-1] Of PChunk;
 NChunk  : Integer;
 ERIReady : Boolean;
 ERIThresh : Real;

Function PIdx(I, J: Integer): Integer;
Begin
 If I >= J Then PIdx := (I * (I - 1)) Div 2 + J
 Else PIdx := (J * (J - 1)) Div 2 + I;
End;

Function PIdx2(A, B: Integer): LongInt;
Begin
 If A >= B Then PIdx2 := LongInt(A) * (A - 1) Div 2 + B
 Else PIdx2 := LongInt(B) * (B - 1) Div 2 + A;
End;

{ Ajoute une fonction de base (coquille S de TS) sur l'atome A }
Function AddFn(A, S, LX, LY, LZ: Integer; Const Lab: String): Boolean;
Var K, L: Integer;
Begin
 AddFn := False;
 If NBF >= MAX_BF Then Exit;
 Inc(NBF);
 BFAtom[NBF] := A;
 BFL[NBF, 1] := LX; BFL[NBF, 2] := LY; BFL[NBF, 3] := LZ;
 BFNP[NBF] := TSN[S];
 L := LX + LY + LZ;
 For K := 1 To TSN[S] Do Begin
  BFA[NBF, K] := TSA[S, K];
  BFC[NBF, K] := TSC[S, K] * PrimNorm(TSA[S, K], L);
 End;
 BFLab[NBF] := Lab;
 AddFn := True;
End;

Function EHerm(I, J, T: Integer; Qx, A, B: Real): Real;
Var P, Mu: Real;
Begin
 If (T < 0) Or (T > I + J) Then Begin EHerm := 0; Exit; End;
 P := A + B;
 Mu := A * B / P;
 If (I = 0) And (J = 0) And (T = 0) Then EHerm := ExpS(-Mu * Qx * Qx)
 Else If J = 0 Then
  EHerm := EHerm(I - 1, J, T - 1, Qx, A, B) / (2 * P)
         - (Mu * Qx / A) * EHerm(I - 1, J, T, Qx, A, B)
         + (T + 1) * EHerm(I - 1, J, T + 1, Qx, A, B)
 Else
  EHerm := EHerm(I, J - 1, T - 1, Qx, A, B) / (2 * P)
         + (Mu * Qx / B) * EHerm(I, J - 1, T, Qx, A, B)
         + (T + 1) * EHerm(I, J - 1, T + 1, Qx, A, B);
End;

{ Fonction de Boys par serie et recurrence descendante (precise, lente) }
Procedure BoysSlow(M: Integer; T: Real; Var F: TBoys);
Var K: Integer; Term, Sum, ET: Real;
Begin
 If T < 1E-10 Then Begin
  For K := 0 To M Do F[K] := 1 / (2 * K + 1);
 End
 Else If T > 35 Then Begin
  ET := ExpS(-T);
  F[0] := 0.5 * Sqrt(Pi / T);
  For K := 0 To M - 1 Do F[K + 1] := ((2 * K + 1) * F[K] - ET) / (2 * T);
 End
 Else Begin
  ET := Exp(-T);
  Term := 1 / (2 * M + 1);
  Sum := Term;
  K := 0;
  Repeat
   Inc(K);
   Term := Term * 2 * T / (2 * M + 2 * K + 1);
   Sum := Sum + Term;
  Until (Term < 1E-16 * Sum) Or (K > 250);
  F[M] := ET * Sum;
  For K := M DownTo 1 Do F[K - 1] := (2 * T * F[K] + ET) / (2 * K - 1);
 End;
End;

{ Table de Boys : F_m(T) pour T = i * BOYS_DT, m = 0..BOYS_NM ; interpolation de Taylor }
Const
 BOYS_NM = 10;
 {$IFDEF FPC}
 BOYS_N = 700;
 BOYS_DT = 0.05;
 {$ELSE}
 BOYS_N = 350;
 BOYS_DT = 0.1;
 {$ENDIF}

Type
 TBoysTab = Array[0..BOYS_N, 0..BOYS_NM] Of Real;
 PBoysTab = ^TBoysTab;

Var
 BoysTab : PBoysTab;

Procedure InitBoysTable;
Var I, K: Integer; F: TBoys;
Begin
 If BoysTab <> Nil Then Exit;
 New(BoysTab);
 For I := 0 To BOYS_N Do Begin
  BoysSlow(BOYS_NM, I * BOYS_DT, F);
  For K := 0 To BOYS_NM Do BoysTab^[I, K] := F[K];
 End;
End;

Procedure Boys(M: Integer; T: Real; Var F: TBoys);
Var K, I: Integer; X, D: Real;
Begin
 If T >= BOYS_N * BOYS_DT - 0.5 Then Begin
  F[0] := 0.5 * Sqrt(Pi / T);
  For K := 0 To M - 1 Do F[K + 1] := (2 * K + 1) * F[K] / (2 * T);
  Exit;
 End;
 I := Round(T / BOYS_DT);
 X := I * BOYS_DT - T;
 For K := 0 To M Do
  F[K] := BoysTab^[I, K] + X * (BoysTab^[I, K + 1] + X / 2 * (BoysTab^[I, K + 2] +
          X / 3 * (BoysTab^[I, K + 3] + X / 4 * (BoysTab^[I, K + 4] + X / 5 * (BoysTab^[I, K + 5] +
          X / 6 * BoysTab^[I, K + 6])))));
End;

{ Tenseur R_tuv^0 pour t+u+v <= L }
Procedure BuildR(L: Integer; Alpha, X, Y, Z: Real);
Var N, T, U, V, LL: Integer; F: TBoys; Fac, T2, R: Real;
Begin
 T2 := Alpha * (X * X + Y * Y + Z * Z);
 Boys(L, T2, F);
 Fac := 1;
 For N := 0 To L Do Begin
  RT[N, 0, 0, 0] := Fac * F[N];
  Fac := Fac * (-2 * Alpha);
 End;
 For LL := 1 To L Do
  For N := 0 To L - LL Do
   For T := 0 To LL Do
    For U := 0 To LL - T Do Begin
     V := LL - T - U;
     If T > 0 Then Begin
      R := X * RT[N + 1, T - 1, U, V];
      If T > 1 Then R := R + (T - 1) * RT[N + 1, T - 2, U, V];
     End
     Else If U > 0 Then Begin
      R := Y * RT[N + 1, T, U - 1, V];
      If U > 1 Then R := R + (U - 1) * RT[N + 1, T, U - 2, V];
     End
     Else Begin
      R := Z * RT[N + 1, T, U, V - 1];
      If V > 1 Then R := R + (V - 1) * RT[N + 1, T, U, V - 2];
     End;
     RT[N, T, U, V] := R;
    End;
End;

Function Ovl1D(I, J: Integer; Q, A, B: Real): Real;
Begin
 Ovl1D := EHerm(I, J, 0, Q, A, B) * Sqrt(Pi / (A + B));
End;

Function Kin1D(I, J: Integer; Q, A, B: Real): Real;
Var R: Real;
Begin
 R := -2 * B * B * Ovl1D(I, J + 2, Q, A, B) + B * (2 * J + 1) * Ovl1D(I, J, Q, A, B);
 If J >= 2 Then R := R - 0.5 * J * (J - 1) * Ovl1D(I, J - 2, Q, A, B);
 Kin1D := R;
End;

{ Recouvrement brut (sans renormalisation) entre deux fonctions de base }
Function OvlRaw(M, N: Integer): Real;
Var K, L: Integer; S, A, B: Real;
Begin
 S := 0;
 For K := 1 To BFNP[M] Do
  For L := 1 To BFNP[N] Do Begin
   A := BFA[M, K]; B := BFA[N, L];
   S := S + BFC[M, K] * BFC[N, L] *
    Ovl1D(BFL[M,1], BFL[N,1], AtX[BFAtom[M],1] - AtX[BFAtom[N],1], A, B) *
    Ovl1D(BFL[M,2], BFL[N,2], AtX[BFAtom[M],2] - AtX[BFAtom[N],2], A, B) *
    Ovl1D(BFL[M,3], BFL[N,3], AtX[BFAtom[M],3] - AtX[BFAtom[N],3], A, B);
  End;
 OvlRaw := S;
End;

Function KinRaw(M, N: Integer): Real;
Var K, L: Integer; S, A, B, SX, SY, SZ, TX, TY, TZ, QX, QY, QZ: Real;
Begin
 S := 0;
 QX := AtX[BFAtom[M],1] - AtX[BFAtom[N],1];
 QY := AtX[BFAtom[M],2] - AtX[BFAtom[N],2];
 QZ := AtX[BFAtom[M],3] - AtX[BFAtom[N],3];
 For K := 1 To BFNP[M] Do
  For L := 1 To BFNP[N] Do Begin
   A := BFA[M, K]; B := BFA[N, L];
   SX := Ovl1D(BFL[M,1], BFL[N,1], QX, A, B);
   SY := Ovl1D(BFL[M,2], BFL[N,2], QY, A, B);
   SZ := Ovl1D(BFL[M,3], BFL[N,3], QZ, A, B);
   TX := Kin1D(BFL[M,1], BFL[N,1], QX, A, B);
   TY := Kin1D(BFL[M,2], BFL[N,2], QY, A, B);
   TZ := Kin1D(BFL[M,3], BFL[N,3], QZ, A, B);
   S := S + BFC[M, K] * BFC[N, L] * (TX * SY * SZ + SX * TY * SZ + SX * SY * TZ);
  End;
 KinRaw := S;
End;

{ Construit la base pour la molecule courante ; renvoie Faux si impossible }
Function BuildBasis(Var Err: String): Boolean;
Var A, S, Ns, Np, K: Integer; Sn, Nrm: Real; OK: Boolean;
Begin
 Err := '';
 NBF := 0;
 BuildBasis := False;
 For A := 1 To NAt Do Begin
  If Not GetShells(BasisSel, AtZ[A]) Then Begin
   Err := 'La base ' + BasisNames[BasisSel] + ' n''est pas definie pour l''element ' +
          ElemStr(AtZ[A]);
   Exit;
  End;
  Ns := 0; Np := 1;
  For S := 1 To NTS Do Begin
   If TSL[S] = 0 Then Begin
    Inc(Ns);
    OK := AddFn(A, S, 0, 0, 0, IStr(Ns) + 's');
   End
   Else Begin
    Inc(Np);
    OK := AddFn(A, S, 1, 0, 0, IStr(Np) + 'px');
    If OK Then OK := AddFn(A, S, 0, 1, 0, IStr(Np) + 'py');
    If OK Then OK := AddFn(A, S, 0, 0, 1, IStr(Np) + 'pz');
   End;
   If Not OK Then Begin
    Err := 'Trop de fonctions de base (maximum ' + IStr(MAX_BF) + ')';
    Exit;
   End;
  End;
 End;
 { normalisation des fonctions contractees }
 For S := 1 To NBF Do Begin
  Sn := OvlRaw(S, S);
  Nrm := 1 / Sqrt(Sn);
  For K := 1 To BFNP[S] Do BFC[S, K] := BFC[S, K] * Nrm;
 End;
 BuildBasis := True;
End;

{ --------------------------------------------------------------------- }
{ Donnees de paires de fonctions de base                                 }
{ --------------------------------------------------------------------- }

Procedure FreePairs;
Var I: Integer;
Begin
 For I := 1 To MAX_PAIRS Do
  If PairP[I] <> Nil Then Begin
   FreeMem(PairP[I], PairNP[I] * SizeOf(TPrimPair));
   PairP[I] := Nil;
  End;
 PairsOK := False;
End;

Procedure BuildPairs;
Var
 M, N, K, L, Idx, NP, T, U, V, I1: Integer;
 A, B, P, C: Real;
 QX, QY, QZ: Real;
 EX, EY, EZ: Array[0..2] Of Real;
 PP: PPrim;
Begin
 FreePairs;
 NPairs := 0;
 For M := 1 To NBF Do
  For N := 1 To M Do Begin
   Idx := PIdx(M, N);
   If Idx > NPairs Then NPairs := Idx;
   PairI[Idx] := M; PairJ[Idx] := N;
   NP := BFNP[M] * BFNP[N];
   PairNP[Idx] := NP;
   PairL[Idx] := BFL[M,1] + BFL[M,2] + BFL[M,3] + BFL[N,1] + BFL[N,2] + BFL[N,3];
   GetMem(PairP[Idx], NP * SizeOf(TPrimPair));
   QX := AtX[BFAtom[M],1] - AtX[BFAtom[N],1];
   QY := AtX[BFAtom[M],2] - AtX[BFAtom[N],2];
   QZ := AtX[BFAtom[M],3] - AtX[BFAtom[N],3];
   I1 := 0;
   For K := 1 To BFNP[M] Do
    For L := 1 To BFNP[N] Do Begin
     Inc(I1);
     PP := PPrim(@PairP[Idx]^[I1]);
     A := BFA[M, K]; B := BFA[N, L];
     P := A + B;
     C := BFC[M, K] * BFC[N, L];
     PP^.P := P;
     PP^.PX := (A * AtX[BFAtom[M],1] + B * AtX[BFAtom[N],1]) / P;
     PP^.PY := (A * AtX[BFAtom[M],2] + B * AtX[BFAtom[N],2]) / P;
     PP^.PZ := (A * AtX[BFAtom[M],3] + B * AtX[BFAtom[N],3]) / P;
     For T := 0 To 2 Do Begin
      EX[T] := EHerm(BFL[M,1], BFL[N,1], T, QX, A, B);
      EY[T] := EHerm(BFL[M,2], BFL[N,2], T, QY, A, B);
      EZ[T] := EHerm(BFL[M,3], BFL[N,3], T, QZ, A, B);
     End;
     PP^.Cm := 0;
     For T := 0 To 2 Do
      For U := 0 To 2 Do
       For V := 0 To 2 Do Begin
        PP^.E[T * 9 + U * 3 + V] := C * EX[T] * EY[U] * EZ[V];
        If Abs(PP^.E[T * 9 + U * 3 + V]) > PP^.Cm Then PP^.Cm := Abs(PP^.E[T * 9 + U * 3 + V]);
       End;
    End;
  End;
 PairsOK := True;
End;

{ --------------------------------------------------------------------- }
{ Integrales a deux electrons                                            }
{ --------------------------------------------------------------------- }

Function ERIPairs(PA, PB: Integer): Real;
Var
 I1, I2, T, U, V, T2, U2, V2, L1, L2, LT: Integer;
 PP1, PP2: PPrim;
 P, Q, Al, Sm, E1, E2, Tot: Real;
Begin
 L1 := PairL[PA]; L2 := PairL[PB]; LT := L1 + L2;
 Tot := 0;
 For I1 := 1 To PairNP[PA] Do Begin
  PP1 := PPrim(@PairP[PA]^[I1]);
  P := PP1^.P;
  For I2 := 1 To PairNP[PB] Do Begin
   PP2 := PPrim(@PairP[PB]^[I2]);
   Q := PP2^.P;
   If PP1^.Cm * PP2^.Cm * TWO_PI_52 / (P * Q * Sqrt(P + Q)) < 1E-13 Then Continue;
   Al := P * Q / (P + Q);
   BuildR(LT, Al, PP1^.PX - PP2^.PX, PP1^.PY - PP2^.PY, PP1^.PZ - PP2^.PZ);
   Sm := 0;
   For T := 0 To L1 Do
    For U := 0 To L1 - T Do
     For V := 0 To L1 - T - U Do Begin
      E1 := PP1^.E[T * 9 + U * 3 + V];
      If E1 = 0 Then Continue;
      For T2 := 0 To L2 Do
       For U2 := 0 To L2 - T2 Do
        For V2 := 0 To L2 - T2 - U2 Do Begin
         E2 := PP2^.E[T2 * 9 + U2 * 3 + V2];
         If E2 = 0 Then Continue;
         If Odd(T2 + U2 + V2) Then
          Sm := Sm - E1 * E2 * RT[0, T + T2, U + U2, V + V2]
         Else
          Sm := Sm + E1 * E2 * RT[0, T + T2, U + U2, V + V2];
        End;
     End;
   Tot := Tot + Sm * TWO_PI_52 / (P * Q * Sqrt(P + Q));
  End;
 End;
 ERIPairs := Tot;
End;

Procedure FreeERI;
Var I: Integer;
Begin
 For I := 0 To MAX_CHUNK - 1 Do
  If ChunkP[I] <> Nil Then Begin
   Dispose(ChunkP[I]);
   ChunkP[I] := Nil;
  End;
 ERIReady := False;
End;

{ Calcule et stocke toutes les integrales (ij|kl) ; renvoie Faux si memoire }
Function ComputeERI: Boolean;
Var
 IJ, KL, C: Integer; Idx: LongInt; V: Real; NQ: LongInt;
Begin
 ComputeERI := False;
 FreeERI;
 NQ := LongInt(NPairs) * (NPairs + 1) Div 2;
 NChunk := Integer((NQ + CHUNK_SIZE - 1) Div CHUNK_SIZE);
 If NChunk > MAX_CHUNK Then Exit;
 {$IFNDEF FPC}
 If MaxAvail < LongInt(NChunk) * SizeOf(TChunk) + 40000 Then Exit;
 {$ENDIF}
 For C := 0 To NChunk - 1 Do Begin
  New(ChunkP[C]);
  FillChar(ChunkP[C]^, SizeOf(TChunk), 0);
 End;
 For IJ := 1 To NPairs Do
  PairQ[IJ] := Sqrt(Abs(ERIPairs(IJ, IJ)));
 For IJ := 1 To NPairs Do
  For KL := 1 To IJ Do Begin
   If PairQ[IJ] * PairQ[KL] < ERIThresh Then Continue;
   V := ERIPairs(IJ, KL);
   Idx := LongInt(IJ) * (IJ - 1) Div 2 + KL - 1;
   ChunkP[Idx Div CHUNK_SIZE]^[Idx Mod CHUNK_SIZE] := V;
  End;
 ERIReady := True;
 ComputeERI := True;
End;

Function ERIAt(IJ, KL: Integer): Real;
Var Idx: LongInt;
Begin
 Idx := PIdx2(IJ, KL) - 1;
 ERIAt := ChunkP[Idx Div CHUNK_SIZE]^[Idx Mod CHUNK_SIZE];
End;

Function ERI4(I, J, K, L: Integer): Real;
Begin
 ERI4 := ERIAt(PIdx(I, J), PIdx(K, L));
End;

{ --------------------------------------------------------------------- }
{ Algebre lineaire                                                      }
{ --------------------------------------------------------------------- }

Type
 {$IFDEF FPC}
 TRealArr = Array[0..16383] Of Real;
 {$ELSE}
 TRealArr = Array[0..4095] Of Real;
 {$ENDIF}
 PRealArr = ^TRealArr;
 TBig     = Array[1..CIS_MAX, 1..CIS_MAX] Of Real;
 PBig     = ^TBig;
 TCube    = Array[1..MAX_BF * MAX_BF * MAX_BF] Of Real;
 PCube    = ^TCube;

{ Diagonalisation de Jacobi d'une matrice symetrique N x N de dimension
  principale LD (A est detruite). W recoit les valeurs propres triees par
  ordre croissant (indices 0..N-1), V les vecteurs propres en colonnes. }
Procedure Jacobi(N, LD: Integer; A, V, W: PRealArr);
Var
 I, J, P, Q, K, Sweep: Integer;
 Off, Theta, T, C, S, AKP, AKQ, Tmp: Real;
 LDL: LongInt;
Begin
 LDL := LD;
 For I := 1 To N Do
  For J := 1 To N Do
   If I = J Then V^[(I - 1) * LDL + J - 1] := 1 Else V^[(I - 1) * LDL + J - 1] := 0;
 For Sweep := 1 To 60 Do Begin
  Off := 0;
  For P := 1 To N - 1 Do
   For Q := P + 1 To N Do Off := Off + Sqr(A^[(P - 1) * LDL + Q - 1]);
  If Off < 1E-26 Then Break;
  For P := 1 To N - 1 Do
   For Q := P + 1 To N Do Begin
    If Abs(A^[(P - 1) * LDL + Q - 1]) < 1E-20 Then Begin
     A^[(P - 1) * LDL + Q - 1] := 0; A^[(Q - 1) * LDL + P - 1] := 0;
     Continue;
    End;
    Theta := (A^[(Q - 1) * LDL + Q - 1] - A^[(P - 1) * LDL + P - 1]) /
             (2 * A^[(P - 1) * LDL + Q - 1]);
    If Theta >= 0 Then T := 1 / (Theta + Sqrt(1 + Theta * Theta))
    Else T := -1 / (-Theta + Sqrt(1 + Theta * Theta));
    C := 1 / Sqrt(1 + T * T);
    S := T * C;
    For K := 1 To N Do Begin
     AKP := A^[(K - 1) * LDL + P - 1];
     AKQ := A^[(K - 1) * LDL + Q - 1];
     A^[(K - 1) * LDL + P - 1] := C * AKP - S * AKQ;
     A^[(K - 1) * LDL + Q - 1] := S * AKP + C * AKQ;
    End;
    For K := 1 To N Do Begin
     AKP := A^[(P - 1) * LDL + K - 1];
     AKQ := A^[(Q - 1) * LDL + K - 1];
     A^[(P - 1) * LDL + K - 1] := C * AKP - S * AKQ;
     A^[(Q - 1) * LDL + K - 1] := S * AKP + C * AKQ;
    End;
    For K := 1 To N Do Begin
     AKP := V^[(K - 1) * LDL + P - 1];
     AKQ := V^[(K - 1) * LDL + Q - 1];
     V^[(K - 1) * LDL + P - 1] := C * AKP - S * AKQ;
     V^[(K - 1) * LDL + Q - 1] := S * AKP + C * AKQ;
    End;
   End;
 End;
 For I := 1 To N Do W^[I - 1] := A^[(I - 1) * LDL + I - 1];
 { tri croissant }
 For I := 1 To N - 1 Do Begin
  K := I;
  For J := I + 1 To N Do If W^[J - 1] < W^[K - 1] Then K := J;
  If K <> I Then Begin
   Tmp := W^[I - 1]; W^[I - 1] := W^[K - 1]; W^[K - 1] := Tmp;
   For J := 1 To N Do Begin
    Tmp := V^[(J - 1) * LDL + I - 1];
    V^[(J - 1) * LDL + I - 1] := V^[(J - 1) * LDL + K - 1];
    V^[(J - 1) * LDL + K - 1] := Tmp;
   End;
  End;
 End;
End;

{ Resolution d'un systeme lineaire N x N par elimination de Gauss
  (A : [1..N,1..N] stockee dans un tableau 10 x 10, B : second membre) }
Type
 TSmall = Array[1..12, 1..12] Of Real;
 TSmallV = Array[1..12] Of Real;

Function SolveSmall(N: Integer; Var A: TSmall; Var B: TSmallV): Boolean;
Var I, J, K, P: Integer; M, T: Real;
Begin
 SolveSmall := False;
 For I := 1 To N Do Begin
  P := I;
  For J := I + 1 To N Do If Abs(A[J, I]) > Abs(A[P, I]) Then P := J;
  If Abs(A[P, I]) < 1E-30 Then Exit;
  If P <> I Then Begin
   For K := 1 To N Do Begin T := A[I, K]; A[I, K] := A[P, K]; A[P, K] := T; End;
   T := B[I]; B[I] := B[P]; B[P] := T;
  End;
  For J := I + 1 To N Do Begin
   M := A[J, I] / A[I, I];
   For K := I To N Do A[J, K] := A[J, K] - M * A[I, K];
   B[J] := B[J] - M * B[I];
  End;
 End;
 For I := N DownTo 1 Do Begin
  T := B[I];
  For K := I + 1 To N Do T := T - A[I, K] * B[K];
  B[I] := T / A[I, I];
 End;
 SolveSmall := True;
End;

{ --------------------------------------------------------------------- }
{ Matrices de travail                                                    }
{ --------------------------------------------------------------------- }

Const
 DIIS_MAX = 8;

Var
 MS, MH, MT, MVn, MX, MXh, MFa, MFb, MDa, MDb, MCa, MCb : PMat;
 MJ, MKa, MKb, MXa, MXb, MT1, MT2, MT3, MDpx, MDpy, MDpz : PMat;
 DiisF, DiisE : Array[1..DIIS_MAX, 1..2] Of PMat;
 EpsA, EpsB : TVecN;
 MatsReady : Boolean;

Procedure AllocMats;
Var I, J: Integer;
 Procedure Mk(Var P: PMat);
 Begin
  New(P);
  FillChar(P^, SizeOf(TMat), 0);
 End;
Begin
 If MatsReady Then Exit;
 Mk(MS); Mk(MH); Mk(MT); Mk(MVn); Mk(MX); Mk(MXh); Mk(MFa); Mk(MFb);
 Mk(MDa); Mk(MDb); Mk(MCa); Mk(MCb); Mk(MJ); Mk(MKa); Mk(MKb);
 Mk(MXa); Mk(MXb); Mk(MT1); Mk(MT2); Mk(MT3); Mk(MDpx); Mk(MDpy); Mk(MDpz);
 For I := 1 To DIIS_MAX Do
  For J := 1 To 2 Do Begin
   Mk(DiisF[I, J]); Mk(DiisE[I, J]);
  End;
 MatsReady := True;
End;

Procedure MatZero(M: PMat);
Begin
 FillChar(M^, SizeOf(TMat), 0);
End;

Procedure MatCopy(Src, Dst: PMat);
Begin
 Move(Src^, Dst^, SizeOf(TMat));
End;

{ C := A * B }
Procedure MatMul(A, B, C: PMat);
Var I, J, K: Integer; S: Real;
Begin
 For I := 1 To NBF Do
  For J := 1 To NBF Do Begin
   S := 0;
   For K := 1 To NBF Do S := S + A^[I, K] * B^[K, J];
   C^[I, J] := S;
  End;
End;

{ C := A^T * B }
Procedure MatTMul(A, B, C: PMat);
Var I, J, K: Integer; S: Real;
Begin
 For I := 1 To NBF Do
  For J := 1 To NBF Do Begin
   S := 0;
   For K := 1 To NBF Do S := S + A^[K, I] * B^[K, J];
   C^[I, J] := S;
  End;
End;

{ Produit scalaire de deux matrices (trace de A^T B) }
Function MatDot(A, B: PMat): Real;
Var I, J: Integer; S: Real;
Begin
 S := 0;
 For I := 1 To NBF Do
  For J := 1 To NBF Do S := S + A^[I, J] * B^[I, J];
 MatDot := S;
End;

{ Diagonalise une matrice de Fock : F C = S C e ; C recoit les orbitales }
Procedure DiagFock(F, C: PMat; Var Eps: TVecN);
Var W: TVecN; I: Integer;
Begin
 MatMul(F, MX, MT1);
 MatTMul(MX, MT1, MT2);
 Jacobi(NBF, MAX_BF, PRealArr(MT2), PRealArr(MT3), PRealArr(@W));
 MatMul(MX, MT3, C);
 For I := 1 To NBF Do Eps[I] := W[I];
End;

{ --------------------------------------------------------------------- }
{ Parametres et etat du calcul                                           }
{ --------------------------------------------------------------------- }

Var
 Meth      : Integer;         { methode : HF, MP2, DFT }
 Func      : Integer;         { fonctionnelle }
 HFExch    : Real;            { fraction d'echange Hartree-Fock }
 UnrestrictedFlag : Boolean;
 JobType   : Integer;
 NOccA, NOccB : Integer;
 Enuc      : Real;
 EScf, EOne, ECoul, EExch, EXCdft : Real;
 EMP2, EMP2ss, EMP2os : Real;
 ETot      : Real;
 ScfConv   : Boolean;
 ScfIter   : Integer;
 ScfTolE, ScfTolD, ScfTolErr : Real;
 ScfMaxIter : Integer;
 ScfTightName : String[12];
 OrbsValid : Boolean;
 OrbN, OrbOA, OrbOB : Integer;
 HistE     : Array[0..MAX_HIST-1] Of Real;
 HistN     : Integer;
 S2Exp     : Real;
 QuietSCF  : Boolean;
 FrozenCore : Boolean;
 DoMP2     : Boolean;

{ Initialise les parametres de calcul par defaut }
Procedure DefaultParams;
Begin
 Meth := M_HF; Func := F_NONE; HFExch := 1.0;
 UnrestrictedFlag := False;
 JobType := J_SP;
 BasisSel := B_STO3G;
 ScfTolE := 1E-6; ScfTolD := 1E-6; ScfTolErr := 1E-5;
 ScfMaxIter := 125;
 ScfTightName := 'NormalSCF';
 FrozenCore := True;
 Charge := 0; Mult := 1;
End;

{ Nombre d'electrons et occupations ; renvoie Faux si incoherent }
Function SetupElectrons(Var Err: String): Boolean;
Var A, N: Integer;
Begin
 SetupElectrons := False;
 N := 0;
 For A := 1 To NAt Do N := N + AtZ[A];
 N := N - Charge;
 NElec := N;
 If N < 0 Then Begin Err := 'Nombre d''electrons negatif'; Exit; End;
 If Mult < 1 Then Begin Err := 'Multiplicite invalide'; Exit; End;
 If (N + Mult - 1) Mod 2 <> 0 Then Begin
  Err := 'Charge ' + IStr(Charge) + ' et multiplicite ' + IStr(Mult) +
         ' incompatibles avec ' + IStr(N) + ' electrons';
  Exit;
 End;
 NOccA := (N + Mult - 1) Div 2;
 NOccB := N - NOccA;
 If NOccB < 0 Then Begin Err := 'Multiplicite trop elevee'; Exit; End;
 If (NOccA > NBF) Or (NOccB > NBF) Then Begin
  Err := 'Pas assez de fonctions de base pour ' + IStr(N) + ' electrons';
  Exit;
 End;
 SetupElectrons := True;
End;

{ --------------------------------------------------------------------- }
{ Integrales a un electron                                               }
{ --------------------------------------------------------------------- }

Procedure BuildOneElectron;
Var
 IJ, I, J, K, C, T, U, V, L: Integer;
 SO, VN, DX, DY, DZ, PF, Sm, Al: Real;
 PP: PPrim;
Begin
 For IJ := 1 To NPairs Do Begin
  I := PairI[IJ]; J := PairJ[IJ];
  SO := 0; DX := 0; DY := 0; DZ := 0; VN := 0;
  For K := 1 To PairNP[IJ] Do Begin
   PP := PPrim(@PairP[IJ]^[K]);
   PF := Exp(1.5 * Ln(Pi / PP^.P));
   SO := SO + PP^.E[0] * PF;
   DX := DX + (PP^.E[9] + PP^.PX * PP^.E[0]) * PF;
   DY := DY + (PP^.E[3] + PP^.PY * PP^.E[0]) * PF;
   DZ := DZ + (PP^.E[1] + PP^.PZ * PP^.E[0]) * PF;
   L := PairL[IJ];
   For C := 1 To NAt Do Begin
    BuildR(L, PP^.P, PP^.PX - AtX[C, 1], PP^.PY - AtX[C, 2], PP^.PZ - AtX[C, 3]);
    Sm := 0;
    For T := 0 To L Do
     For U := 0 To L - T Do
      For V := 0 To L - T - U Do
       Sm := Sm + PP^.E[T * 9 + U * 3 + V] * RT[0, T, U, V];
    VN := VN - AtZ[C] * 2 * Pi / PP^.P * Sm;
   End;
  End;
  MS^[I, J] := SO; MS^[J, I] := SO;
  MVn^[I, J] := VN; MVn^[J, I] := VN;
  MDpx^[I, J] := DX; MDpx^[J, I] := DX;
  MDpy^[I, J] := DY; MDpy^[J, I] := DY;
  MDpz^[I, J] := DZ; MDpz^[J, I] := DZ;
  SO := KinRaw(I, J);
  MT^[I, J] := SO; MT^[J, I] := SO;
 End;
 For I := 1 To NBF Do
  For J := 1 To NBF Do MH^[I, J] := MT^[I, J] + MVn^[I, J];
End;

{ Orthogonalisation symetrique : X = S^(-1/2), Xh = S^(1/2) }
Function BuildOrtho: Boolean;
Var W: TVecN; I, J, K: Integer; Sm, Sp: Real;
Begin
 BuildOrtho := False;
 MatCopy(MS, MT1);
 Jacobi(NBF, MAX_BF, PRealArr(MT1), PRealArr(MT2), PRealArr(@W));
 If W[1] < 1E-7 Then Exit;
 For I := 1 To NBF Do
  For J := 1 To NBF Do Begin
   Sm := 0; Sp := 0;
   For K := 1 To NBF Do Begin
    Sm := Sm + MT2^[I, K] * MT2^[J, K] / Sqrt(W[K + 0]);
    Sp := Sp + MT2^[I, K] * MT2^[J, K] * Sqrt(W[K + 0]);
   End;
   MX^[I, J] := Sm;
   MXh^[I, J] := Sp;
  End;
 BuildOrtho := True;
End;

{ --------------------------------------------------------------------- }
{ Matrices J et K                                                        }
{ --------------------------------------------------------------------- }

Procedure BuildJK(Closed, NeedK: Boolean);
Var
 IJ, KL, I, J, K, L: Integer;
 Idx0, Idx: LongInt;
 V, Dtt, Dtt2: Real;
Begin
 MatZero(MJ); MatZero(MKa); MatZero(MKb);
 For IJ := 1 To NPairs Do Begin
  I := PairI[IJ]; J := PairJ[IJ];
  Idx0 := LongInt(IJ) * (IJ - 1) Div 2 - 1;
  For KL := 1 To IJ Do Begin
   Idx := Idx0 + KL;
   V := ChunkP[Idx Div CHUNK_SIZE]^[Idx Mod CHUNK_SIZE];
   If Abs(V) < 1E-14 Then Continue;
   K := PairI[KL]; L := PairJ[KL];
   If I = J Then V := V * 0.5;
   If K = L Then V := V * 0.5;
   If IJ = KL Then V := V * 0.5;
   { Coulomb : densite totale = Da + Db }
   Dtt := (MDa^[K, L] + MDb^[K, L]) * 2 * V;
   MJ^[I, J] := MJ^[I, J] + Dtt;
   MJ^[J, I] := MJ^[J, I] + Dtt;
   Dtt2 := (MDa^[I, J] + MDb^[I, J]) * 2 * V;
   MJ^[K, L] := MJ^[K, L] + Dtt2;
   MJ^[L, K] := MJ^[L, K] + Dtt2;
   If NeedK Then Begin
    MKa^[I, K] := MKa^[I, K] + MDa^[J, L] * V;
    MKa^[J, K] := MKa^[J, K] + MDa^[I, L] * V;
    MKa^[I, L] := MKa^[I, L] + MDa^[J, K] * V;
    MKa^[J, L] := MKa^[J, L] + MDa^[I, K] * V;
    MKa^[K, I] := MKa^[K, I] + MDa^[L, J] * V;
    MKa^[L, I] := MKa^[L, I] + MDa^[K, J] * V;
    MKa^[K, J] := MKa^[K, J] + MDa^[L, I] * V;
    MKa^[L, J] := MKa^[L, J] + MDa^[K, I] * V;
    If Not Closed Then Begin
     MKb^[I, K] := MKb^[I, K] + MDb^[J, L] * V;
     MKb^[J, K] := MKb^[J, K] + MDb^[I, L] * V;
     MKb^[I, L] := MKb^[I, L] + MDb^[J, K] * V;
     MKb^[J, L] := MKb^[J, L] + MDb^[I, K] * V;
     MKb^[K, I] := MKb^[K, I] + MDb^[L, J] * V;
     MKb^[L, I] := MKb^[L, I] + MDb^[K, J] * V;
     MKb^[K, J] := MKb^[K, J] + MDb^[L, I] * V;
     MKb^[L, J] := MKb^[L, J] + MDb^[K, I] * V;
    End;
   End;
  End;
 End;
 If Closed Then MatCopy(MKa, MKb);
End;

{ --------------------------------------------------------------------- }
{ Boucle SCF (RHF / UHF / Kohn-Sham)                                    }
{ --------------------------------------------------------------------- }

Procedure BuildXC(Closed: Boolean; Var Exc: Real); Forward;

{ Densites a partir des orbitales occupees }
Procedure BuildDensity;
Var I, J, K: Integer; S: Real;
Begin
 For I := 1 To NBF Do
  For J := 1 To I Do Begin
   S := 0;
   For K := 1 To NOccA Do S := S + MCa^[I, K] * MCa^[J, K];
   MDa^[I, J] := S; MDa^[J, I] := S;
   S := 0;
   For K := 1 To NOccB Do S := S + MCb^[I, K] * MCb^[J, K];
   MDb^[I, J] := S; MDb^[J, I] := S;
  End;
End;

{ Matrices de Fock et energie pour les densites MDa / MDb }
Procedure BuildFock(Closed: Boolean; Var E: Real);
Var I, J: Integer; Exc, Ec1, Ec2, Ek: Real; NeedK: Boolean;
Begin
 NeedK := Abs(HFExch) > 1E-12;
 BuildJK(Closed, NeedK);
 Exc := 0;
 MatZero(MXa); MatZero(MXb);
 If Meth = M_DFT Then Begin
  BuildXC(Closed, Exc);
  If Closed Then MatCopy(MXa, MXb);
 End;
 For I := 1 To NBF Do
  For J := 1 To NBF Do Begin
   MFa^[I, J] := MH^[I, J] + MJ^[I, J] - HFExch * MKa^[I, J] + MXa^[I, J];
   MFb^[I, J] := MH^[I, J] + MJ^[I, J] - HFExch * MKb^[I, J] + MXb^[I, J];
  End;
 { energies }
 EOne := 0; Ec1 := 0; Ek := 0;
 For I := 1 To NBF Do
  For J := 1 To NBF Do Begin
   EOne := EOne + (MDa^[I, J] + MDb^[I, J]) * MH^[I, J];
   Ec1 := Ec1 + (MDa^[I, J] + MDb^[I, J]) * MJ^[I, J];
   Ek := Ek + MDa^[I, J] * MKa^[I, J] + MDb^[I, J] * MKb^[I, J];
  End;
 ECoul := 0.5 * Ec1;
 EExch := -0.5 * HFExch * Ek;
 EXCdft := Exc;
 E := Enuc + EOne + ECoul + EExch + Exc;
End;

{ Orbitales de depart : orbitales precedentes orthonormalisees }
Function GuessFromOrbs: Boolean;
Var I, J, K, M, NOc: Integer; Sp: Integer; Ov: Real; CM: PMat;
Begin
 GuessFromOrbs := False;
 For Sp := 1 To 2 Do Begin
  If Sp = 1 Then Begin CM := MCa; NOc := NOccA; End
  Else Begin CM := MCb; NOc := NOccB; End;
  For I := 1 To NOc Do Begin
   { Gram-Schmidt dans la metrique S }
   For J := 1 To I - 1 Do Begin
    Ov := 0;
    For K := 1 To NBF Do
     For M := 1 To NBF Do Ov := Ov + CM^[K, J] * MS^[K, M] * CM^[M, I];
    For K := 1 To NBF Do CM^[K, I] := CM^[K, I] - Ov * CM^[K, J];
   End;
   Ov := 0;
   For K := 1 To NBF Do
    For M := 1 To NBF Do Ov := Ov + CM^[K, I] * MS^[K, M] * CM^[M, I];
   If Ov < 1E-10 Then Exit;
   Ov := 1 / Sqrt(Ov);
   For K := 1 To NBF Do CM^[K, I] := CM^[K, I] * Ov;
  End;
 End;
 GuessFromOrbs := True;
End;

{ Guess de Wolfsberg-Helmholz generalise }
Procedure GuessGWH;
Var I, J: Integer;
Begin
 For I := 1 To NBF Do
  For J := 1 To NBF Do
   If I = J Then MFa^[I, J] := MH^[I, I]
   Else MFa^[I, J] := 0.875 * MS^[I, J] * (MH^[I, I] + MH^[J, J]);
 DiagFock(MFa, MCa, EpsA);
 MatCopy(MCa, MCb);
 EpsB := EpsA;
End;

{ DIIS : erreur X^T (FDS - SDF) X }
Procedure DiisError(F, D, E: PMat);
Var I, J: Integer;
Begin
 MatMul(F, D, MT1);
 MatMul(MT1, MS, MT2);          { FDS }
 MatMul(MS, D, MT1);
 MatMul(MT1, F, MT3);           { SDF }
 For I := 1 To NBF Do
  For J := 1 To NBF Do MT2^[I, J] := MT2^[I, J] - MT3^[I, J];
 MatTMul(MX, MT2, MT1);
 MatMul(MT1, MX, E);
End;

Function MaxAbsMat(M: PMat): Real;
Var I, J: Integer; R: Real;
Begin
 R := 0;
 For I := 1 To NBF Do
  For J := 1 To NBF Do If Abs(M^[I, J]) > R Then R := Abs(M^[I, J]);
 MaxAbsMat := R;
End;

Function ScfRun(Verbose, UseGuess: Boolean): Boolean;
Var
 Iter, I, J, K, Sp, Slot, DiisN: Integer;
 Closed, Done, Guessed: Boolean;
 E, EOld, DE, RmsD, MaxD, ErrMax, Sum, Tmp: Real;
 Bm: TSmall; Rhs: TSmallV;
 ErrA, ErrB: PMat;
Begin
 ScfRun := False;
 ScfConv := False;
 Closed := (NOccA = NOccB);
 HistN := 0;
 Guessed := False;
 If UseGuess And OrbsValid And (OrbN = NBF) And (OrbOA = NOccA) And (OrbOB = NOccB) Then
  Guessed := GuessFromOrbs;
 If Not Guessed Then Begin
  GuessGWH;
 End;
 BuildDensity;
 New(ErrA); New(ErrB);
 DiisN := 0;
 EOld := 0;
 If Verbose Then Begin
  EmitBlank;
  Emit('ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]');
  Emit('               ***  Starting incremental Fock matrix formation  ***');
 End;
 Done := False;
 Iter := 0;
 While (Iter < ScfMaxIter) And Not Done Do Begin
  Inc(Iter);
  BuildFock(Closed, E);
  { erreurs DIIS }
  DiisError(MFa, MDa, ErrA);
  ErrMax := MaxAbsMat(ErrA);
  If Not Closed Then Begin
   DiisError(MFb, MDb, ErrB);
   If MaxAbsMat(ErrB) > ErrMax Then ErrMax := MaxAbsMat(ErrB);
  End;
  { pile DIIS }
  If DiisN < DIIS_MAX Then Begin Inc(DiisN); Slot := DiisN; End
  Else Begin
   { decalage }
   For K := 1 To DIIS_MAX - 1 Do
    For Sp := 1 To 2 Do Begin
     Move(DiisF[K + 1, Sp]^, DiisF[K, Sp]^, SizeOf(TMat));
     Move(DiisE[K + 1, Sp]^, DiisE[K, Sp]^, SizeOf(TMat));
    End;
   Slot := DIIS_MAX;
  End;
  Move(MFa^, DiisF[Slot, 1]^, SizeOf(TMat));
  Move(ErrA^, DiisE[Slot, 1]^, SizeOf(TMat));
  If Not Closed Then Begin
   Move(MFb^, DiisF[Slot, 2]^, SizeOf(TMat));
   Move(ErrB^, DiisE[Slot, 2]^, SizeOf(TMat));
  End;
  If DiisN >= 2 Then Begin
   For I := 1 To DiisN Do
    For J := 1 To I Do Begin
     Sum := MatDot(DiisE[I, 1], DiisE[J, 1]);
     If Not Closed Then Sum := Sum + MatDot(DiisE[I, 2], DiisE[J, 2]);
     Bm[I, J] := Sum; Bm[J, I] := Sum;
    End;
   For I := 1 To DiisN Do Begin Bm[I, DiisN + 1] := -1; Bm[DiisN + 1, I] := -1; Rhs[I] := 0; End;
   Bm[DiisN + 1, DiisN + 1] := 0; Rhs[DiisN + 1] := -1;
   If SolveSmall(DiisN + 1, Bm, Rhs) Then Begin
    MatZero(MFa);
    If Not Closed Then MatZero(MFb);
    For K := 1 To DiisN Do
     For I := 1 To NBF Do
      For J := 1 To NBF Do Begin
       MFa^[I, J] := MFa^[I, J] + Rhs[K] * DiisF[K, 1]^[I, J];
       If Not Closed Then MFb^[I, J] := MFb^[I, J] + Rhs[K] * DiisF[K, 2]^[I, J];
      End;
   End
   Else DiisN := 0;
  End;
  { nouvelles orbitales }
  DiagFock(MFa, MCa, EpsA);
  If Closed Then Begin MatCopy(MCa, MCb); EpsB := EpsA; End
  Else DiagFock(MFb, MCb, EpsB);
  MatCopy(MDa, MT3);
  BuildDensity;
  { convergence }
  MaxD := 0; RmsD := 0;
  For I := 1 To NBF Do
   For J := 1 To NBF Do Begin
    Tmp := MDa^[I, J] - MT3^[I, J];
    If Abs(Tmp) > MaxD Then MaxD := Abs(Tmp);
    RmsD := RmsD + Tmp * Tmp;
   End;
  RmsD := Sqrt(RmsD) / NBF;
  DE := E - EOld;
  If Iter = 1 Then DE := 0;
  If Iter <= MAX_HIST Then Begin HistE[Iter - 1] := E; HistN := Iter; End;
  If Verbose Then
   Emit(' ' + FmtR(Iter - 1, 3, 0) + '  ' + FmtR(E, 16, 10) + '  ' + FmtR(DE, 14, 10) + '  ' +
        SciStr(MaxD, 4) + '  ' + SciStr(RmsD, 4) + '  ' + SciStr(ErrMax, 4));
  If (Iter > 1) And (Abs(DE) < ScfTolE) And (RmsD < ScfTolD) And (ErrMax < ScfTolErr) Then Done := True;
  EOld := E;
  ShowProgress('SCF iteration ' + IStr(Iter), 0);
  If UserAbort Then Begin Done := True; End;
 End;
 Dispose(ErrA); Dispose(ErrB);
 ScfIter := Iter;
 If Done And Not UserAbort Then Begin
  { energie finale avec la densite convergee }
  BuildFock(Closed, E);
  EScf := E;
  ScfConv := True;
  ScfRun := True;
  OrbsValid := True; OrbN := NBF; OrbOA := NOccA; OrbOB := NOccB;
  If Verbose Then Begin
   EmitBlank;
   Emit('               *****************************************************');
   Emit('               *                      SUCCESS                      *');
   Emit('               *           SCF CONVERGED AFTER ' + FmtR(Iter, 3, 0) + ' CYCLES          *');
   Emit('               *****************************************************');
  End;
 End
 Else Begin
  If Verbose Then Begin
   EmitBlank;
   Emit('SCF NOT CONVERGED AFTER ' + IStr(Iter) + ' CYCLES');
  End;
  OrbsValid := False;
 End;
End;

{ --------------------------------------------------------------------- }
{ Theorie de la fonctionnelle de la densite (grille de Becke)            }
{ --------------------------------------------------------------------- }

Const
 GBLK_SIZE   = 512;
 {$IFDEF FPC}
 MAX_GBLK    = 1100;
 GRID_NR     = 36;      { points radiaux par atome }
 GRID_NTH    = 9;       { points en cos(theta) ; 2*NTH points en phi }
 {$ELSE}
 MAX_GBLK    = 12;
 GRID_NR     = 16;
 GRID_NTH    = 5;
 {$ENDIF}
 CX_LDA      = 0.9305257363491;         { (3/2)(3/(4 pi))^(1/3) }
 B88_BETA    = 0.0042;

Type
 TGridPt = Record
  X, Y, Z, W : Real;
 End;
 TGBlk  = Array[0..GBLK_SIZE-1] Of TGridPt;
 PGBlk  = ^TGBlk;
 PGridPt = ^TGridPt;

Var
 GridBlk  : Array[0..MAX_GBLK-1] Of PGBlk;
 GridN    : LongInt;
 GridOK   : Boolean;
 BFCut2   : Array[1..MAX_BF] Of Real;
 DftNAlpha, DftNBeta : Real;
 EXdft, ECdft : Real;
 GridScale : Integer;       { 1 = grille par defaut, 2 = fine }
 GridFull  : Boolean;       { capacite de la grille depassee }

Procedure FreeGrid;
Var I: Integer;
Begin
 For I := 0 To MAX_GBLK - 1 Do
  If GridBlk[I] <> Nil Then Begin Dispose(GridBlk[I]); GridBlk[I] := Nil; End;
 GridN := 0;
 GridOK := False;
End;

Procedure AddGridPt(X, Y, Z, W: Real);
Var B: Integer; O: Integer;
Begin
 B := Integer(GridN Div GBLK_SIZE);
 O := Integer(GridN Mod GBLK_SIZE);
 If B >= MAX_GBLK Then Begin GridFull := True; Exit; End;
 If GridBlk[B] = Nil Then New(GridBlk[B]);
 GridBlk[B]^[O].X := X; GridBlk[B]^[O].Y := Y; GridBlk[B]^[O].Z := Z;
 GridBlk[B]^[O].W := W;
 Inc(GridN);
End;

Function BeckeWeight(A: Integer; PX, PY, PZ: Real): Real;
Var
 B, C: Integer; Den, Pa, MU, Chi, U, AA, NU, S: Real;
 PA_: Array[1..MAX_AT] Of Real;
 RR: Array[1..MAX_AT] Of Real;
Begin
 For B := 1 To NAt Do
  RR[B] := Sqrt(Sqr(PX - AtX[B,1]) + Sqr(PY - AtX[B,2]) + Sqr(PZ - AtX[B,3]));
 Den := 0;
 For B := 1 To NAt Do Begin
  Pa := 1;
  For C := 1 To NAt Do
   If C <> B Then Begin
    MU := (RR[B] - RR[C]) / Dist(B, C);
    Chi := ElRad[AtZ[B]] / ElRad[AtZ[C]];
    U := (Chi - 1) / (Chi + 1);
    AA := U / (U * U - 1);
    If AA > 0.5 Then AA := 0.5;
    If AA < -0.5 Then AA := -0.5;
    NU := MU + AA * (1 - MU * MU);
    S := NU;
    S := 1.5 * S - 0.5 * S * S * S;
    S := 1.5 * S - 0.5 * S * S * S;
    S := 1.5 * S - 0.5 * S * S * S;
    Pa := Pa * 0.5 * (1 - S);
    If Pa < 1E-30 Then Begin Pa := 0; Break; End;
   End;
  PA_[B] := Pa;
  Den := Den + Pa;
 End;
 If Den < 1E-30 Then BeckeWeight := 0
 Else BeckeWeight := PA_[A] / Den;
End;

{ rayon de coupure des fonctions de base }
Procedure ComputeBFCut;
Var I, K: Integer; R: Real;
Begin
 For I := 1 To NBF Do Begin
  BFCut2[I] := 0;
  For K := 1 To BFNP[I] Do Begin
   R := 0.2;
   While (R < 40) And (Abs(BFC[I, K]) * (1 + R) * ExpS(-BFA[I, K] * R * R) > 1E-11) Do R := R + 0.2;
   If R * R > BFCut2[I] Then BFCut2[I] := R * R;
  End;
 End;
End;

Procedure BuildGrid;
Var
 A, I, J, K, NR, NTH, NPH: Integer;
 XL: Array[1..GRID_NTH * 2] Of Real;
 WL: Array[1..GRID_NTH * 2] Of Real;
 Xi, X, R, WR, ST, CT, Ph, PX, PY, PZ, WA, W, Z0, Z1, PP, P1, P2, P3: Real;
Begin
 FreeGrid;
 GridFull := False;
 NR := GRID_NR * GridScale;
 NTH := GRID_NTH * GridScale;
 If NTH > GRID_NTH * 2 Then NTH := GRID_NTH * 2;
 NPH := 2 * NTH;
 { Gauss-Legendre sur cos(theta) }
 For I := 1 To NTH Do Begin
  Z0 := Cos(Pi * (I - 0.25) / (NTH + 0.5));
  Repeat
   P1 := 1; P2 := 0;
   For J := 1 To NTH Do Begin
    P3 := P2; P2 := P1;
    P1 := ((2 * J - 1) * Z0 * P2 - (J - 1) * P3) / J;
   End;
   PP := NTH * (Z0 * P1 - P2) / (Z0 * Z0 - 1);
   Z1 := Z0;
   Z0 := Z1 - P1 / PP;
  Until Abs(Z0 - Z1) < 1E-13;
  XL[I] := Z0;
  WL[I] := 2 / ((1 - Z0 * Z0) * PP * PP);
 End;
 For A := 1 To NAt Do Begin
  Xi := 0.5 * (ElRad[AtZ[A]] + 0.15) * ANG2BOHR * 1.0;
  If AtZ[A] = 1 Then Xi := 0.55;
  For I := 1 To NR Do Begin
   X := Cos(Pi * I / (NR + 1));
   R := Xi * (1 + X) / (1 - X);
   WR := Pi / (NR + 1) * Sin(Pi * I / (NR + 1)) * 2 * Xi / Sqr(1 - X) * R * R;
   For J := 1 To NTH Do Begin
    CT := XL[J];
    ST := Sqrt(1 - CT * CT);
    For K := 0 To NPH - 1 Do Begin
     Ph := 2 * Pi * (K + 0.5) / NPH;
     PX := AtX[A,1] + R * ST * Cos(Ph);
     PY := AtX[A,2] + R * ST * Sin(Ph);
     PZ := AtX[A,3] + R * CT;
     WA := WR * WL[J] * 2 * Pi / NPH;
     W := WA * BeckeWeight(A, PX, PY, PZ);
     If W > 1E-15 Then AddGridPt(PX, PY, PZ, W);
    End;
   End;
  End;
 End;
  ComputeBFCut;
 GridOK := True;
End;

{ Fonctions de base et gradients au point (X,Y,Z) }
Type
 TAOVals = Array[1..MAX_BF] Of Real;
 TIdxArr = Array[0..MAX_BF-1] Of Integer;

Procedure EvalAOs(X, Y, Z: Real; Var NL: Integer; Var Idx: TIdxArr;
                  Var Phi, GX, GY, GZ: TAOVals);
Var
 M, K: Integer; DX, DY, DZ, R2, E, PX, PY, PZ, QX, QY, QZ, A, C: Real;
 SP, SX, SY, SZ: Real;
Begin
 NL := 0;
 For M := 1 To NBF Do Begin
  DX := X - AtX[BFAtom[M], 1];
  DY := Y - AtX[BFAtom[M], 2];
  DZ := Z - AtX[BFAtom[M], 3];
  R2 := DX * DX + DY * DY + DZ * DZ;
  If R2 > BFCut2[M] Then Continue;
  If BFL[M, 1] = 1 Then PX := DX Else PX := 1;
  If BFL[M, 2] = 1 Then PY := DY Else PY := 1;
  If BFL[M, 3] = 1 Then PZ := DZ Else PZ := 1;
  SP := 0; SX := 0; SY := 0; SZ := 0;
  For K := 1 To BFNP[M] Do Begin
   A := BFA[M, K];
   E := BFC[M, K] * ExpS(-A * R2);
   SP := SP + E;
   If BFL[M, 1] = 1 Then QX := 1 - 2 * A * DX * DX Else QX := -2 * A * DX;
   If BFL[M, 2] = 1 Then QY := 1 - 2 * A * DY * DY Else QY := -2 * A * DY;
   If BFL[M, 3] = 1 Then QZ := 1 - 2 * A * DZ * DZ Else QZ := -2 * A * DZ;
   SX := SX + E * QX; SY := SY + E * QY; SZ := SZ + E * QZ;
  End;
  Inc(NL);
  Idx[NL - 1] := M;
  Phi[NL] := SP * PX * PY * PZ;
  GX[NL] := SX * PY * PZ;
  GY[NL] := SY * PX * PZ;
  GZ[NL] := SZ * PX * PY;
 End;
End;

{ --------------------------------------------------------------------- }
{ Fonctionnelles d'echange et de correlation (densite d'energie)         }
{ --------------------------------------------------------------------- }

Function PW92G(Rs, A, A1, B1, B2, B3, B4: Real): Real;
Var Den, SR: Real;
Begin
 SR := Sqrt(Rs);
 Den := 2 * A * (B1 * SR + B2 * Rs + B3 * Rs * SR + B4 * Rs * Rs);
 PW92G := -2 * A * (1 + A1 * Rs) * Ln(1 + 1 / Den);
End;

Function FZeta(Z: Real): Real;
Begin
 FZeta := (Pw(1 + Z, 4 / 3) + Pw(1 - Z, 4 / 3) - 2) / (Pw(2, 4 / 3) - 2);
End;

{ Energie de correlation PW92 par electron }
Function EpsPW92(Rs, Zt: Real): Real;
Var E0, E1, AC, F, Z4: Real;
Begin
 E0 := PW92G(Rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
 E1 := PW92G(Rs, 0.01554535, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517);
 AC := -PW92G(Rs, 0.0168869, 0.11125, 10.357, 3.6231, 0.88026, 0.49671);
 F := FZeta(Zt);
 Z4 := Sqr(Sqr(Zt));
 EpsPW92 := E0 + AC * F / 1.709921 * (1 - Z4) + (E1 - E0) * F * Z4;
End;

Function VWNF(X, A, X0, B, C: Real): Real;
Var XX, X0X, Q: Real;
Begin
 XX := X * X + B * X + C;
 X0X := X0 * X0 + B * X0 + C;
 Q := Sqrt(4 * C - B * B);
 VWNF := A * (Ln(X * X / XX) + 2 * B / Q * ArcTan(Q / (2 * X + B))
          - B * X0 / X0X * (Ln(Sqr(X - X0) / XX) + 2 * (B + 2 * X0) / Q * ArcTan(Q / (2 * X + B))));
End;

{ Energie de correlation VWN5 par electron }
Function EpsVWN5(Rs, Zt: Real): Real;
Var X, EP, EF, AC, F, Z4: Real;
Begin
 X := Sqrt(Rs);
 EP := VWNF(X, 0.0310907, -0.10498, 3.72744, 12.9352);
 EF := VWNF(X, 0.01554535, -0.32500, 7.06042, 18.0578);
 AC := VWNF(X, -0.01688686394, -0.0047584, 1.13107, 13.0045);
 F := FZeta(Zt);
 Z4 := Sqr(Sqr(Zt));
 EpsVWN5 := EP + AC * F / 1.709921 * (1 - Z4) + (EF - EP) * F * Z4;
End;

Function ExSlater(Ra, Rb: Real): Real;
Begin
 ExSlater := -CX_LDA * (Pw(Ra, 4 / 3) + Pw(Rb, 4 / 3));
End;

{ correction de gradient de Becke 88 pour un spin }
Function ExB88Spin(R, S: Real): Real;
Var R43, X: Real;
Begin
 If R < 1E-20 Then Begin ExB88Spin := 0; Exit; End;
 R43 := Pw(R, 4 / 3);
 X := Sqrt(Abs(S)) / R43;
 ExB88Spin := -B88_BETA * R43 * X * X / (1 + 6 * B88_BETA * X * ArcSinh(X));
End;

Function ExPBESpin(R, S: Real): Real;
Var N, GN2, KF, SS, FX, EL: Real;
Begin
 If R < 1E-20 Then Begin ExPBESpin := 0; Exit; End;
 N := 2 * R;
 GN2 := 4 * Abs(S);
 KF := Pw(3 * Pi * Pi * N, 1 / 3);
 SS := Sqrt(GN2) / (2 * KF * N);
 FX := 1 + 0.804 - 0.804 / (1 + 0.2195149727645171 * SS * SS / 0.804);
 EL := -0.75 * Pw(3 / Pi, 1 / 3) * Pw(N, 4 / 3);
 ExPBESpin := 0.5 * EL * FX;
End;

Function EcPBE(Ra, Rb, Saa, Sab, Sbb: Real): Real;
Var
 R, Zt, Rs, Eps, Phi, KF, KS, GR2, T, A, T2, H, Beta, Gam: Real;
Begin
 R := Ra + Rb;
 If R < 1E-20 Then Begin EcPBE := 0; Exit; End;
 Beta := 0.06672455060314922;
 Gam := (1 - Ln(2.0)) / (Pi * Pi);
 Zt := (Ra - Rb) / R;
 If Zt > 1 - 1E-12 Then Zt := 1 - 1E-12;
 If Zt < -1 + 1E-12 Then Zt := -1 + 1E-12;
 Rs := Pw(3 / (4 * Pi * R), 1 / 3);
 Eps := EpsPW92(Rs, Zt);
 Phi := 0.5 * (Pw(1 + Zt, 2 / 3) + Pw(1 - Zt, 2 / 3));
 KF := Pw(3 * Pi * Pi * R, 1 / 3);
 KS := Sqrt(4 * KF / Pi);
 GR2 := Saa + 2 * Sab + Sbb;
 If GR2 < 0 Then GR2 := 0;
 T := Sqrt(GR2) / (2 * Phi * KS * R);
 T2 := T * T;
 A := Beta / Gam / (Exp(-Eps / (Gam * Phi * Phi * Phi)) - 1);
 H := Gam * Phi * Phi * Phi * Ln(1 + Beta / Gam * T2 * (1 + A * T2) / (1 + A * T2 + A * A * T2 * T2));
 EcPBE := R * (Eps + H);
End;

Function EcLYP(Ra, Rb, Saa, Sab, Sbb: Real): Real;
Const LA = 0.04918; LB = 0.132; LC = 0.2533; LD = 0.349;
Var
 R, RM13, Den, Om, Dl, CF, ST, T1, T2: Real;
Begin
 R := Ra + Rb;
 If R < 1E-20 Then Begin EcLYP := 0; Exit; End;
 RM13 := Pw(R, -1 / 3);
 Den := 1 + LD * RM13;
 Om := ExpS(-LC * RM13) / Den * Pw(R, -11 / 3);
 Dl := LC * RM13 + LD * RM13 / Den;
 CF := 0.3 * Pw(3 * Pi * Pi, 2 / 3);
 ST := Saa + 2 * Sab + Sbb;
 If ST < 0 Then ST := 0;
 T1 := -LA * 4 / Den * Ra * Rb / R;
 T2 := -LA * LB * Om * (Ra * Rb * (Pw(2, 11 / 3) * CF * (Pw(Ra, 8 / 3) + Pw(Rb, 8 / 3))
        + (47 / 18 - 7 * Dl / 18) * ST - (2.5 - Dl / 18) * (Saa + Sbb)
        - (Dl - 11) / 9 * (Ra / R * Saa + Rb / R * Sbb))
        - 2 / 3 * R * R * ST + (2 / 3 * R * R - Ra * Ra) * Sbb
        + (2 / 3 * R * R - Rb * Rb) * Saa);
 EcLYP := T1 + T2;
End;

Function EcVWN(Ra, Rb: Real): Real;
Var R, Zt, Rs: Real;
Begin
 R := Ra + Rb;
 If R < 1E-20 Then Begin EcVWN := 0; Exit; End;
 Zt := (Ra - Rb) / R;
 If Zt > 1 - 1E-12 Then Zt := 1 - 1E-12;
 If Zt < -1 + 1E-12 Then Zt := -1 + 1E-12;
 Rs := Pw(3 / (4 * Pi * R), 1 / 3);
 EcVWN := R * EpsVWN5(Rs, Zt);
End;

Function ExcDens(Ra, Rb, Saa, Sab, Sbb: Real): Real;
Begin
 Case Func Of
  F_LDA:   ExcDens := ExSlater(Ra, Rb) + EcVWN(Ra, Rb);
  F_BLYP:  ExcDens := ExSlater(Ra, Rb) + ExB88Spin(Ra, Saa) + ExB88Spin(Rb, Sbb)
                      + EcLYP(Ra, Rb, Saa, Sab, Sbb);
  F_B3LYP: ExcDens := 0.8 * ExSlater(Ra, Rb) + 0.72 * (ExB88Spin(Ra, Saa) + ExB88Spin(Rb, Sbb))
                      + 0.19 * EcVWN(Ra, Rb) + 0.81 * EcLYP(Ra, Rb, Saa, Sab, Sbb);
  F_PBE:   ExcDens := ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb) + EcPBE(Ra, Rb, Saa, Sab, Sbb);
  F_PBE0:  ExcDens := 0.75 * (ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb)) + EcPBE(Ra, Rb, Saa, Sab, Sbb);
 Else ExcDens := 0;
 End;
End;

{ Parties echange / correlation separees (pour l'impression) }
Function ExPart(Ra, Rb, Saa, Sbb: Real): Real;
Begin
 Case Func Of
  F_LDA:   ExPart := ExSlater(Ra, Rb);
  F_BLYP:  ExPart := ExSlater(Ra, Rb) + ExB88Spin(Ra, Saa) + ExB88Spin(Rb, Sbb);
  F_B3LYP: ExPart := 0.8 * ExSlater(Ra, Rb) + 0.72 * (ExB88Spin(Ra, Saa) + ExB88Spin(Rb, Sbb));
  F_PBE:   ExPart := ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb);
  F_PBE0:  ExPart := 0.75 * (ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb));
 Else ExPart := 0;
 End;
End;

{ Derivee numerique d'une variable par differences centrees }
Procedure XCDeriv(Ra, Rb, Saa, Sab, Sbb: Real; Var E, VRa, VRb, VSaa, VSab, VSbb: Real);
Var H, Ep, Em: Real;
Begin
 E := ExcDens(Ra, Rb, Saa, Sab, Sbb);
 { rho alpha }
 H := 1E-4 * Ra + 1E-14;
 Ep := ExcDens(Ra + H, Rb, Saa, Sab, Sbb);
 If Ra > 1E-12 Then Em := ExcDens(Ra - H, Rb, Saa, Sab, Sbb) Else Em := E;
 If Ra > 1E-12 Then VRa := (Ep - Em) / (2 * H) Else VRa := (Ep - Em) / H;
 { rho beta }
 H := 1E-4 * Rb + 1E-14;
 Ep := ExcDens(Ra, Rb + H, Saa, Sab, Sbb);
 If Rb > 1E-12 Then Em := ExcDens(Ra, Rb - H, Saa, Sab, Sbb) Else Em := E;
 If Rb > 1E-12 Then VRb := (Ep - Em) / (2 * H) Else VRb := (Ep - Em) / H;
 { sigma aa }
 H := 1E-4 * Saa + 1E-14 * (Ra + Rb + 1E-10);
 Ep := ExcDens(Ra, Rb, Saa + H, Sab, Sbb);
 If Saa - H >= 0 Then Begin
  Em := ExcDens(Ra, Rb, Saa - H, Sab, Sbb);
  VSaa := (Ep - Em) / (2 * H);
 End
 Else VSaa := (Ep - E) / H;
 { sigma bb }
 H := 1E-4 * Sbb + 1E-14 * (Ra + Rb + 1E-10);
 Ep := ExcDens(Ra, Rb, Saa, Sab, Sbb + H);
 If Sbb - H >= 0 Then Begin
  Em := ExcDens(Ra, Rb, Saa, Sab, Sbb - H);
  VSbb := (Ep - Em) / (2 * H);
 End
 Else VSbb := (Ep - E) / H;
 { sigma ab }
 H := 1E-4 * (Saa + Sbb + 1E-8);
 Ep := ExcDens(Ra, Rb, Saa, Sab + H, Sbb);
 Em := ExcDens(Ra, Rb, Saa, Sab - H, Sbb);
 VSab := (Ep - Em) / (2 * H);
End;



{ --- Versions a couche fermee (rho alpha = rho beta) : evaluation allegee --- }

Const
 LN2     = 0.6931471805599453;
 RS_C    = 0.6203504908994001;       { (3/(4 pi))^(1/3) }
 CF_LYP  = 2.871234000188191;        { (3/10) (3 pi^2)^(2/3) }
 P2_113  = 12.699208415745595;       { 2^(11/3) }
 KF_C    = 3.0936677262801355;       { (3 pi^2)^(1/3) }
 PBE_BETA = 0.06672455060314922;
 PBE_GAM  = 0.031090690869654895;    { (1 - ln 2) / pi^2 }

{ epsilon_c VWN5 paramagnetique pour rs donne }
Function EpsVWNp(Rs: Real): Real;
Const A = 0.0310907; X0 = -0.10498; B = 3.72744; C = 12.9352;
Var X, XX, X0X, Q, At, LXX: Real;
Begin
 X := Sqrt(Rs);
 XX := X * X + B * X + C;
 X0X := X0 * X0 + B * X0 + C;
 Q := Sqrt(4 * C - B * B);
 At := ArcTan(Q / (2 * X + B));
 LXX := Ln(XX);
 EpsVWNp := A * (2 * Ln(X) - LXX + 2 * B / Q * At
            - B * X0 / X0X * (2 * Ln(Abs(X - X0)) - LXX + 2 * (B + 2 * X0) / Q * At));
End;

{ Correction de gradient B88 pour un spin (x43 = rho_s^(4/3)) }
Function B88c(X43, S: Real): Real;
Var X: Real;
Begin
 X := Sqrt(S) / X43;
 B88c := -B88_BETA * X43 * X * X / (1 + 6 * B88_BETA * X * ArcSinh(X));
End;

{ Correlation LYP a couche fermee (rho = 2 Ra, Rm13 = rho^(-1/3)) }
Function LYPc(Ra, X43, Rm13, Saa, Sab: Real): Real;
Const LA = 0.04918; LB = 0.132; LC = 0.2533; LD = 0.349;
Var R, Den, Om, Dl, ST, T1, T2, R2, R4, R8, R11, Ra2: Real;
Begin
 R := 2 * Ra;
 Den := 1 + LD * Rm13;
 R2 := Rm13 * Rm13; R4 := R2 * R2; R8 := R4 * R4; R11 := R8 * R2 * Rm13;
 Om := ExpS(-LC * Rm13) / Den * R11;
 Dl := LC * Rm13 + LD * Rm13 / Den;
 ST := 2 * Saa + 2 * Sab;
 If ST < 0 Then ST := 0;
 Ra2 := Ra * Ra;
 T1 := -2 * LA * Ra / Den;
 T2 := -LA * LB * Om * (Ra2 * (P2_113 * CF_LYP * 2 * X43 * X43
        + (47 / 18 - 7 * Dl / 18) * ST - (2.5 - Dl / 18) * 2 * Saa
        - (Dl - 11) / 9 * Saa)
        - 2 / 3 * R * R * ST + 2 * (2 / 3 * R * R - Ra2) * Saa);
 LYPc := T1 + T2;
End;

{ Echange PBE pour un spin (Ra = rho_s, X43 = rho_s^(4/3), Kf1 = (3 pi^2 2 rho_s)^(1/3)) }
Function PBExc(Ra, X43, Kf1, S: Real): Real;
Var N43, SS, FX, EL: Real;
Begin
 N43 := 2.5198420997897464 * X43;       { (2 rho_s)^(4/3) }
 SS := 2 * Sqrt(S) / (2 * Kf1 * 2 * Ra);
 FX := 1.804 - 0.804 / (1 + 0.2195149727645171 * SS * SS / 0.804);
 EL := -0.9305257363491 * N43 / 1.2599210498948732;
 PBExc := 0.5 * EL * FX;
End;

{ Correlation PBE a couche fermee (Rs, Kf = (3 pi^2 rho)^(1/3), ST = |grad rho|^2) }
Function PBEcc(R, Rs, Kf, ST: Real): Real;
Var Eps, KS, T2, A, H: Real;
Begin
 Eps := PW92G(Rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
 KS := Sqrt(4 * Kf / Pi);
 If ST < 0 Then ST := 0;
 T2 := ST / (4 * KS * KS * R * R);
 A := PBE_BETA / PBE_GAM / (Exp(-Eps / PBE_GAM) - 1);
 H := PBE_GAM * Ln(1 + PBE_BETA / PBE_GAM * T2 * (1 + A * T2) / (1 + A * T2 + A * A * T2 * T2));
 PBEcc := R * (Eps + H);
End;

Function ExcDensC(Ra, Saa, Sab: Real): Real;
Var Lr, X43, Rm13, Rs, Kf, E: Real;
Begin
 If Ra < 1E-20 Then Begin ExcDensC := 0; Exit; End;
 Lr := Ln(Ra);
 X43 := Exp(1.3333333333333333 * Lr);
 Rm13 := Exp(-0.3333333333333333 * (Lr + LN2));
 Rs := RS_C * Rm13;
 Case Func Of
  F_LDA:   E := -2 * CX_LDA * X43 + 2 * Ra * EpsVWNp(Rs);
  F_BLYP:  E := -2 * CX_LDA * X43 + 2 * B88c(X43, Saa) + LYPc(Ra, X43, Rm13, Saa, Sab);
  F_B3LYP: E := -1.6 * CX_LDA * X43 + 1.44 * B88c(X43, Saa)
                + 0.19 * 2 * Ra * EpsVWNp(Rs) + 0.81 * LYPc(Ra, X43, Rm13, Saa, Sab);
  F_PBE, F_PBE0: Begin
   Kf := KF_C / Rm13;                           { (3 pi^2 rho)^(1/3) }
   E := PBEcc(2 * Ra, Rs, Kf, 2 * Saa + 2 * Sab);
   If Func = F_PBE Then E := E + 2 * PBExc(Ra, X43, Kf, Saa)
   Else E := E + 1.5 * PBExc(Ra, X43, Kf, Saa);
  End;
 Else E := 0;
 End;
 ExcDensC := E;
End;

Function ExPartC(Ra, Saa: Real): Real;
Var X43, Kf, Rm13: Real;
Begin
 If Ra < 1E-20 Then Begin ExPartC := 0; Exit; End;
 X43 := Pw(Ra, 4 / 3);
 Case Func Of
  F_LDA:   ExPartC := -2 * CX_LDA * X43;
  F_BLYP:  ExPartC := -2 * CX_LDA * X43 + 2 * B88c(X43, Saa);
  F_B3LYP: ExPartC := -1.6 * CX_LDA * X43 + 1.44 * B88c(X43, Saa);
  F_PBE, F_PBE0: Begin
   Kf := Pw(3 * Pi * Pi * 2 * Ra, 1 / 3);
   If Func = F_PBE Then ExPartC := 2 * PBExc(Ra, X43, Kf, Saa) Else ExPartC := 1.5 * PBExc(Ra, X43, Kf, Saa);
  End;
 Else ExPartC := 0;
 End;
End;

{ Version a couche fermee : derivees par symetrie }
Procedure XCDerivClosed(R, S: Real; Var E, VRa, VSaa, VSab: Real);
Var Ra, H, Ep, Em: Real;
Begin
 Ra := 0.5 * R;
 E := ExcDensC(Ra, S, S);
 H := 1E-4 * Ra;
 Ep := ExcDensC(Ra + H, S, S);
 Em := ExcDensC(Ra - H, S, S);
 VRa := (Ep - Em) / (4 * H);
 If Func = F_LDA Then Begin VSaa := 0; VSab := 0; Exit; End;
 H := 1E-4 * S + 1E-14 * (R + 1E-10);
 Ep := ExcDensC(Ra, S + H, S);
 If S - H >= 0 Then Begin
  Em := ExcDensC(Ra, S - H, S);
  VSaa := (Ep - Em) / (4 * H);
 End
 Else VSaa := (Ep - E) / (2 * H);
 H := 1E-4 * (S + 1E-8);
 Ep := ExcDensC(Ra, S, S + H);
 Em := ExcDensC(Ra, S, S - H);
 VSab := (Ep - Em) / (2 * H);
End;

{ Construit les potentiels d'echange-correlation Vxc (MXa, MXb) et l'energie }
Procedure BuildXC(Closed: Boolean; Var Exc: Real);
Var
 B, O, NL, I, J: Integer;
 P: LongInt;
 Idx: TIdxArr;
 Phi, GX, GY, GZ, AA, AB: TAOVals;
 Pt: PGridPt;
 Ra, Rb, Gax, Gay, Gaz, Gbx, Gby, Gbz, Saa, Sab, Sbb: Real;
 E, VRa, VRb, VSaa, VSab, VSbb, W, TA, TB, F: Real;
 QAx, QAy, QAz, QBx, QBy, QBz: Real;
Begin
 Exc := 0; DftNAlpha := 0; DftNBeta := 0; EXdft := 0; ECdft := 0;
 For P := 0 To GridN - 1 Do Begin
  B := Integer(P Div GBLK_SIZE); O := Integer(P Mod GBLK_SIZE);
  Pt := PGridPt(@GridBlk[B]^[O]);
  W := Pt^.W;
  EvalAOs(Pt^.X, Pt^.Y, Pt^.Z, NL, Idx, Phi, GX, GY, GZ);
  If NL = 0 Then Continue;
  Ra := 0; Gax := 0; Gay := 0; Gaz := 0;
  Rb := 0; Gbx := 0; Gby := 0; Gbz := 0;
  For I := 1 To NL Do Begin
   TA := 0; TB := 0;
   For J := 1 To NL Do Begin
    TA := TA + MDa^[Idx[I - 1], Idx[J - 1]] * Phi[J];
    If Not Closed Then TB := TB + MDb^[Idx[I - 1], Idx[J - 1]] * Phi[J];
   End;
   Ra := Ra + TA * Phi[I];
   Gax := Gax + 2 * TA * GX[I]; Gay := Gay + 2 * TA * GY[I]; Gaz := Gaz + 2 * TA * GZ[I];
   If Not Closed Then Begin
    Rb := Rb + TB * Phi[I];
    Gbx := Gbx + 2 * TB * GX[I]; Gby := Gby + 2 * TB * GY[I]; Gbz := Gbz + 2 * TB * GZ[I];
   End;
  End;
  If Closed Then Begin Rb := Ra; Gbx := Gax; Gby := Gay; Gbz := Gaz; End;
  DftNAlpha := DftNAlpha + W * Ra;
  DftNBeta := DftNBeta + W * Rb;
  If Ra + Rb < 1E-10 Then Continue;
  If Ra < 0 Then Ra := 0;
  If Rb < 0 Then Rb := 0;
  Saa := Gax * Gax + Gay * Gay + Gaz * Gaz;
  Sbb := Gbx * Gbx + Gby * Gby + Gbz * Gbz;
  Sab := Gax * Gbx + Gay * Gby + Gaz * Gbz;
  If Closed Then Begin
   XCDerivClosed(Ra + Rb, Saa, E, VRa, VSaa, VSab);
   VRb := VRa; VSbb := VSaa;
  End
  Else XCDeriv(Ra, Rb, Saa, Sab, Sbb, E, VRa, VRb, VSaa, VSab, VSbb);
  Exc := Exc + W * E;
  If Closed Then EXdft := EXdft + W * ExPartC(Ra, Saa) Else EXdft := EXdft + W * ExPart(Ra, Rb, Saa, Sbb);
  { vecteurs de gradient effectifs }
  QAx := 2 * VSaa * Gax + VSab * Gbx; QAy := 2 * VSaa * Gay + VSab * Gby; QAz := 2 * VSaa * Gaz + VSab * Gbz;
  QBx := 2 * VSbb * Gbx + VSab * Gax; QBy := 2 * VSbb * Gby + VSab * Gay;
  QBz := 2 * VSbb * Gbz + VSab * Gaz;
  For I := 1 To NL Do Begin
   AA[I] := W * (0.5 * VRa * Phi[I] + QAx * GX[I] + QAy * GY[I] + QAz * GZ[I]);
   If Not Closed Then AB[I] := W * (0.5 * VRb * Phi[I] + QBx * GX[I] + QBy * GY[I] + QBz * GZ[I]);
  End;
  For I := 1 To NL Do
   For J := 1 To I Do Begin
    F := Phi[I] * AA[J] + Phi[J] * AA[I];
    MXa^[Idx[I - 1], Idx[J - 1]] := MXa^[Idx[I - 1], Idx[J - 1]] + F;
    If Not Closed Then Begin
     F := Phi[I] * AB[J] + Phi[J] * AB[I];
     MXb^[Idx[I - 1], Idx[J - 1]] := MXb^[Idx[I - 1], Idx[J - 1]] + F;
    End;
   End;
 End;
 For I := 1 To NBF Do
  For J := 1 To I - 1 Do Begin
   MXa^[J, I] := MXa^[I, J];
   If Not Closed Then MXb^[J, I] := MXb^[I, J];
  End;
 If Closed Then Begin
  DftNBeta := DftNAlpha;
  EXdft := EXdft;
 End;
 ECdft := Exc - EXdft;
End;

{ --------------------------------------------------------------------- }
{ Analyses de population, moment dipolaire, <S^2>                       }
{ --------------------------------------------------------------------- }

Var
 QMull, QLoew, SpinMull : Array[1..MAX_AT] Of Real;
 BondOrd : Array[1..MAX_AT, 1..MAX_AT] Of Real;
 Dip     : TVec3;
 DipTot  : Real;
 OrbMullPct : Array[1..MAX_BF, 1..MAX_AT] Of Real;   { contribution des atomes aux MOs }

Procedure PopAnalysis;
Var A, B, I, J, K: Integer; Sm, Sl, Ss: Real;
Begin
 { Mulliken : (D S) }
 For I := 1 To NBF Do
  For J := 1 To NBF Do MT1^[I, J] := MDa^[I, J] + MDb^[I, J];
 MatMul(MT1, MS, MT2);
 For A := 1 To NAt Do Begin QMull[A] := AtZ[A]; SpinMull[A] := 0; End;
 For I := 1 To NBF Do QMull[BFAtom[I]] := QMull[BFAtom[I]] - MT2^[I, I];
 For I := 1 To NBF Do
  For J := 1 To NBF Do MT1^[I, J] := MDa^[I, J] - MDb^[I, J];
 MatMul(MT1, MS, MT3);
 For I := 1 To NBF Do SpinMull[BFAtom[I]] := SpinMull[BFAtom[I]] + MT3^[I, I];
 { Loewdin : S^1/2 D S^1/2 }
 For I := 1 To NBF Do
  For J := 1 To NBF Do MT1^[I, J] := MDa^[I, J] + MDb^[I, J];
 MatMul(MXh, MT1, MT2);
 MatMul(MT2, MXh, MT3);
 For A := 1 To NAt Do QLoew[A] := AtZ[A];
 For I := 1 To NBF Do QLoew[BFAtom[I]] := QLoew[BFAtom[I]] - MT3^[I, I];
 { Mayer : P^a S, P^b S }
 MatMul(MDa, MS, MT1);
 MatMul(MDb, MS, MT2);
 For A := 1 To NAt Do
  For B := 1 To NAt Do BondOrd[A, B] := 0;
 For I := 1 To NBF Do
  For J := 1 To NBF Do
   If BFAtom[I] <> BFAtom[J] Then
    BondOrd[BFAtom[I], BFAtom[J]] := BondOrd[BFAtom[I], BFAtom[J]] +
      2 * (MT1^[I, J] * MT1^[J, I] + MT2^[I, J] * MT2^[J, I]);
 { contributions Mulliken des atomes a chaque orbitale (alpha) }
 MatMul(MS, MCa, MT1);
 For K := 1 To NBF Do
  For A := 1 To NAt Do OrbMullPct[K, A] := 0;
 For K := 1 To NBF Do
  For I := 1 To NBF Do
   OrbMullPct[K, BFAtom[I]] := OrbMullPct[K, BFAtom[I]] + MCa^[I, K] * MT1^[I, K];
End;

Procedure ComputeDipole;
Var I, J, D: Integer; El: Real; M: PMat;
Begin
 For D := 1 To 3 Do Begin
  Case D Of 1: M := MDpx; 2: M := MDpy; Else M := MDpz; End;
  El := 0;
  For I := 1 To NBF Do
   For J := 1 To NBF Do El := El + (MDa^[I, J] + MDb^[I, J]) * M^[I, J];
  Dip[D] := -El;
  For I := 1 To NAt Do Dip[D] := Dip[D] + AtZ[I] * AtX[I, D];
 End;
 DipTot := Sqrt(Dip[1] * Dip[1] + Dip[2] * Dip[2] + Dip[3] * Dip[3]);
End;

{ <S^2> pour un determinant UHF }
Function ComputeS2: Real;
Var I, J, M, N: Integer; Ov, Sum, Sz: Real;
Begin
 Sz := (NOccA - NOccB) / 2;
 Sum := 0;
 For I := 1 To NOccA Do
  For J := 1 To NOccB Do Begin
   Ov := 0;
   For M := 1 To NBF Do
    For N := 1 To NBF Do Ov := Ov + MCa^[M, I] * MS^[M, N] * MCb^[N, J];
   Sum := Sum + Ov * Ov;
  End;
 ComputeS2 := Sz * (Sz + 1) + NOccB - Sum;
End;

{ --------------------------------------------------------------------- }
{ Transformation des integrales (MP2, CIS)                               }
{ --------------------------------------------------------------------- }

Var
 CubeA, CubeB, CubeC, CubeD : PCube;
 CubesReady : Boolean;

Procedure AllocCubes;
Begin
 If CubesReady Then Exit;
 New(CubeA); New(CubeB); New(CubeC); New(CubeD);
 CubesReady := True;
End;

Function CIx(A, B, C: Integer): LongInt;
Begin
 CIx := (LongInt(A - 1) * NBF + (B - 1)) * NBF + C;
End;

{ D[p,q,r] = (i p | q r) avec i de Ci, p de C1, q et r de C2 }
Procedure Mo3(Ci: PMat; I: Integer; C1, C2: PMat);
Var NU, LA, SI, MU, P, Q, R, N: Integer; S: Real;
Begin
 N := NBF;
 For NU := 1 To N Do
  For LA := 1 To N Do
   For SI := 1 To N Do Begin
    S := 0;
    For MU := 1 To N Do S := S + Ci^[MU, I] * ERI4(MU, NU, LA, SI);
    CubeA^[CIx(NU, LA, SI)] := S;
   End;
 For P := 1 To N Do
  For LA := 1 To N Do
   For SI := 1 To N Do Begin
    S := 0;
    For NU := 1 To N Do S := S + C1^[NU, P] * CubeA^[CIx(NU, LA, SI)];
    CubeB^[CIx(P, LA, SI)] := S;
   End;
 For P := 1 To N Do
  For Q := 1 To N Do
   For SI := 1 To N Do Begin
    S := 0;
    For LA := 1 To N Do S := S + C2^[LA, Q] * CubeB^[CIx(P, LA, SI)];
    CubeC^[CIx(P, Q, SI)] := S;
   End;
 For P := 1 To N Do
  For Q := 1 To N Do
   For R := 1 To N Do Begin
    S := 0;
    For SI := 1 To N Do S := S + C2^[SI, R] * CubeC^[CIx(P, Q, SI)];
    CubeD^[CIx(P, Q, R)] := S;
   End;
End;

Function NFrozen: Integer;
Var A, N: Integer;
Begin
 N := 0;
 If FrozenCore Then
  For A := 1 To NAt Do Begin
   If (AtZ[A] >= 3) And (AtZ[A] <= 10) Then Inc(N)
   Else If AtZ[A] >= 11 Then N := N + 5;
  End;
 NFrozen := N;
End;

{ Energie de correlation MP2 (RHF ou UHF) }
Procedure RunMP2;
Var
 I, J, A, B, NF, NA, NB, NV, NVB: Integer;
 Closed: Boolean;
 V1, V2, De, Eos, Ess, Eaa, Ebb, Eab: Real;
Begin
 AllocCubes;
 NF := NFrozen;
 Closed := (NOccA = NOccB);
 Eos := 0; Ess := 0;
 NA := NOccA; NB := NOccB;
 If Closed Then Begin
  For I := NF + 1 To NA Do Begin
   Mo3(MCa, I, MCa, MCa);
   For J := NF + 1 To NA Do
    For A := NA + 1 To NBF Do
     For B := NA + 1 To NBF Do Begin
      V1 := CubeD^[CIx(A, J, B)];
      V2 := CubeD^[CIx(B, J, A)];
      De := EpsA[I] + EpsA[J] - EpsA[A] - EpsA[B];
      Eos := Eos + V1 * V1 / De;
      Ess := Ess + V1 * (V1 - V2) / De;
     End;
   ShowProgress('MP2 : orbitale ' + IStr(I), 0);
  End;
  EMP2os := Eos; EMP2ss := Ess;
 End
 Else Begin
  Eaa := 0; Ebb := 0; Eab := 0;
  For I := NF + 1 To NA Do Begin
   Mo3(MCa, I, MCa, MCa);
   For J := NF + 1 To NA Do
    For A := NA + 1 To NBF Do
     For B := NA + 1 To NBF Do Begin
      V1 := CubeD^[CIx(A, J, B)];
      V2 := CubeD^[CIx(B, J, A)];
      De := EpsA[I] + EpsA[J] - EpsA[A] - EpsA[B];
      Eaa := Eaa + 0.25 * Sqr(V1 - V2) / De;
     End;
   Mo3(MCa, I, MCa, MCb);
   For J := NF + 1 To NB Do
    For A := NA + 1 To NBF Do
     For B := NB + 1 To NBF Do Begin
      V1 := CubeD^[CIx(A, J, B)];
      De := EpsA[I] + EpsB[J] - EpsA[A] - EpsB[B];
      Eab := Eab + V1 * V1 / De;
     End;
   ShowProgress('UMP2 : orbitale alpha ' + IStr(I), 0);
  End;
  For I := NF + 1 To NB Do Begin
   Mo3(MCb, I, MCb, MCb);
   For J := NF + 1 To NB Do
    For A := NB + 1 To NBF Do
     For B := NB + 1 To NBF Do Begin
      V1 := CubeD^[CIx(A, J, B)];
      V2 := CubeD^[CIx(B, J, A)];
      De := EpsB[I] + EpsB[J] - EpsB[A] - EpsB[B];
      Ebb := Ebb + 0.25 * Sqr(V1 - V2) / De;
     End;
   ShowProgress('UMP2 : orbitale beta ' + IStr(I), 0);
  End;
  EMP2os := Eab; EMP2ss := Eaa + Ebb;
 End;
 EMP2 := EMP2os + EMP2ss;
End;

{ --------------------------------------------------------------------- }
{ CIS : etats excites monoelectroniques (RHF)                            }
{ --------------------------------------------------------------------- }

Const
 MAX_STATES = 24;

Var
 NRoots    : Integer;
 CISTriplet : Boolean;
 CISN      : Integer;
 CISE      : Array[1..MAX_STATES] Of Real;      { energies d'excitation (Eh) }
 CISF      : Array[1..MAX_STATES] Of Real;      { forces d'oscillateur }
 CISDom    : Array[1..MAX_STATES] Of String[40];
 CISOK     : Boolean;
 CISMsg    : String;

Procedure RunCIS;
Var
 O, V, Dim, I, J, A, B, IA, JB, N, K, D, M1, M2, Cnt: Integer;
 AM, XM: PBig; W: Array[1..CIS_MAX] Of Real;
 DipMO: Array[1..3, 1..CIS_MAX] Of Real;
 TD: TVec3; Mo, X2: Real; PM: PMat;
Begin
 CISOK := False; CISN := 0; CISMsg := '';
 If NOccA <> NOccB Then Begin CISMsg := 'CIS : reference a couche fermee requise'; Exit; End;
 O := NOccA; V := NBF - O;
 Dim := O * V;
 If V < 1 Then Begin CISMsg := 'CIS : aucune orbitale virtuelle'; Exit; End;
 If Dim > CIS_MAX Then Begin
  CISMsg := 'CIS : espace d''excitation trop grand (' + IStr(Dim) + ' > ' + IStr(CIS_MAX) + ')';
  Exit;
 End;
 AllocCubes;
 New(AM); New(XM);
 FillChar(AM^, SizeOf(TBig), 0);
 For I := 1 To O Do Begin
  Mo3(MCa, I, MCa, MCa);
  For J := 1 To O Do
   For A := O + 1 To NBF Do
    For B := O + 1 To NBF Do Begin
     IA := (I - 1) * V + (A - O);
     JB := (J - 1) * V + (B - O);
     If CISTriplet Then
      AM^[IA, JB] := -CubeD^[CIx(J, A, B)]
     Else
      AM^[IA, JB] := 2 * CubeD^[CIx(A, J, B)] - CubeD^[CIx(J, A, B)];
    End;
  ShowProgress('CIS : orbitale ' + IStr(I), 0);
 End;
 For I := 1 To O Do
  For A := O + 1 To NBF Do Begin
   IA := (I - 1) * V + (A - O);
   AM^[IA, IA] := AM^[IA, IA] + EpsA[A] - EpsA[I];
  End;
 For I := 1 To Dim Do
  For J := 1 To I - 1 Do Begin
   Mo := 0.5 * (AM^[I, J] + AM^[J, I]);
   AM^[I, J] := Mo; AM^[J, I] := Mo;
  End;
 { integrales dipolaires entre orbitales occupees et virtuelles }
 For D := 1 To 3 Do Begin
  Case D Of 1: PM := MDpx; 2: PM := MDpy; Else PM := MDpz; End;
  MatMul(PM, MCa, MT1);
  For I := 1 To O Do
   For A := O + 1 To NBF Do Begin
    Mo := 0;
    For M1 := 1 To NBF Do Mo := Mo + MCa^[M1, I] * MT1^[M1, A];
    DipMO[D, (I - 1) * V + (A - O)] := Mo;
   End;
 End;
 Jacobi(Dim, CIS_MAX, PRealArr(AM), PRealArr(XM), PRealArr(@W));
 N := NRoots;
 If N > Dim Then N := Dim;
 If N > MAX_STATES Then N := MAX_STATES;
 For K := 1 To N Do Begin
  CISE[K] := W[K];
  TD[1] := 0; TD[2] := 0; TD[3] := 0;
  If Not CISTriplet Then
   For IA := 1 To Dim Do
    For D := 1 To 3 Do TD[D] := TD[D] + 1.4142135624 * XM^[IA, K] * DipMO[D, IA];
  CISF[K] := 2 / 3 * W[K] * (TD[1] * TD[1] + TD[2] * TD[2] + TD[3] * TD[3]);
  { transitions dominantes }
  CISDom[K] := '';
  Cnt := 0;
  For IA := 1 To Dim Do Begin
   X2 := XM^[IA, K] * XM^[IA, K];
   If (X2 > 0.15) And (Cnt < 2) Then Begin
    I := (IA - 1) Div V + 1;
    A := (IA - 1) Mod V + O + 1;
    If Cnt > 0 Then CISDom[K] := CISDom[K] + ', ';
    CISDom[K] := CISDom[K] + IStr(I - 1) + ' -> ' + IStr(A - 1) + ' (' + RStr(X2, 2) + ')';
    Inc(Cnt);
   End;
  End;
 End;
 Dispose(AM); Dispose(XM);
 CISN := N;
 CISOK := True;
End;

{ --------------------------------------------------------------------- }
{ Energie, gradient numerique, optimisation, frequences                  }
{ --------------------------------------------------------------------- }

Const
 MAX_CONSTR = 8;
 MAX_SCAN   = 80;
 MAX_OPT    = 120;
 THETA_ROT  = 86.61391709;       { K . amu . bohr^2 }
 LAMBDA2    = 1088.4225838;      { bohr^2 . amu . K }
 VOL_PER_K  = 919.5244285;       { bohr^3 / K a 1 atm }
 HCK        = 1.43877688;        { K par cm-1 }

Type
 TConstr = Record
  Kind  : Integer;                 { 1 liaison, 2 angle, 3 dihedre }
  A     : Array[1..4] Of Integer;
  Val   : Real;                    { valeur cible (bohr ou radian) }
  Fix   : Boolean;                 { figer la valeur de depart }
 End;

Var
 NConstr   : Integer;
 Constr    : Array[1..MAX_CONSTR] Of TConstr;
 { balayage }
 ScanOn    : Boolean;
 ScanKind  : Integer;
 ScanAt    : Array[1..4] Of Integer;
 ScanLo, ScanHi : Real;
 ScanSteps : Integer;
 ScanX, ScanE : Array[1..MAX_SCAN] Of Real;
 ScanN     : Integer;
 { optimisation }
 OptTolE, OptTolRG, OptTolMG, OptTolRS, OptTolMS : Real;
 OptMaxCyc : Integer;
 OptE      : Array[0..MAX_OPT] Of Real;
 OptN      : Integer;
 OptConv   : Boolean;
 OptName   : String[12];
 { frequences }
 NFreq     : Integer;
 FreqV     : Array[1..MAX_CONS] Of Real;
 FreqI     : Array[1..MAX_CONS] Of Real;
 FreqVec   : PCartMat;
 FreqOK    : Boolean;
 ThermoOK  : Boolean;
 ThTemp    : Real;
 ZPE, EVibT, ERotT, ETransT, EThermal, HTot, GTot, STransT, SRotT, SVibT, STot : Real;
 SymNum    : Integer;
 ThMass    : Real;
 RotConst  : TVec3;
 GradNorm  : Real;
 LastGrad  : TCartVec;
 NumHStep  : Real;
 SpaceMsg  : String;
 WorkDir   : String;

{ --------------------------------------------------------------------- }
{ Evaluation de l'energie pour la geometrie courante                     }
{ --------------------------------------------------------------------- }

Function EvalEnergy(Verbose, UseGuess: Boolean; Var E: Real; Var Err: String): Boolean;
Begin
 EvalEnergy := False;
 Err := '';
 InitBoysTable;
 BuildPairs;
 BuildOneElectron;
 If Not ComputeERI Then Begin Err := 'Memoire insuffisante pour les integrales bielectroniques'; Exit; End;
 If Not BuildOrtho Then Begin Err := 'Dependance lineaire dans la base'; Exit; End;
 Enuc := NucRepulsion;
 If Meth = M_DFT Then Begin
  BuildGrid;
  If GridFull Then Begin Err := 'Grille DFT trop grande pour cette molecule (utilisez la grille normale)'; Exit; End;
 End;
 If Not ScfRun(Verbose, UseGuess) Then Begin
  If UserAbort Then Err := 'Calcul interrompu' Else Err := 'Le SCF n''a pas converge';
  Exit;
 End;
 ETot := EScf;
 If Meth = M_MP2 Then Begin
  RunMP2;
  ETot := EScf + EMP2;
 End;
 E := ETot;
 EvalEnergy := True;
End;

{ --------------------------------------------------------------------- }
{ Coordonnees internes (contraintes)                                     }
{ --------------------------------------------------------------------- }

Function ConstrValue(Const C: TConstr): Real;
Var
 I: Integer; U, V, W, N1, N2: TVec3; Lu, Lv, Cs, X, Y: Real;
 Function Dot(Const A, B: TVec3): Real;
 Begin Dot := A[1] * B[1] + A[2] * B[2] + A[3] * B[3]; End;
 Procedure Cross(Const A, B: TVec3; Var R: TVec3);
 Begin
  R[1] := A[2] * B[3] - A[3] * B[2];
  R[2] := A[3] * B[1] - A[1] * B[3];
  R[3] := A[1] * B[2] - A[2] * B[1];
 End;
Begin
 Case C.Kind Of
  1: ConstrValue := Dist(C.A[1], C.A[2]);
  2: Begin
   For I := 1 To 3 Do Begin
    U[I] := AtX[C.A[1], I] - AtX[C.A[2], I];
    V[I] := AtX[C.A[3], I] - AtX[C.A[2], I];
   End;
   Lu := Sqrt(Dot(U, U)); Lv := Sqrt(Dot(V, V));
   Cs := Dot(U, V) / (Lu * Lv);
   ConstrValue := ArcCos(Cs);
  End;
 Else Begin
   { dihedre A1-A2-A3-A4 }
   For I := 1 To 3 Do Begin
    U[I] := AtX[C.A[2], I] - AtX[C.A[1], I];
    V[I] := AtX[C.A[3], I] - AtX[C.A[2], I];
    W[I] := AtX[C.A[4], I] - AtX[C.A[3], I];
   End;
   Cross(U, V, N1); Cross(V, W, N2);
   Lv := Sqrt(Dot(V, V));
   X := Dot(N1, N2);
   Cross(N1, N2, U);
   Y := Dot(U, V) / Lv;
   If (X = 0) And (Y = 0) Then ConstrValue := 0
   Else Begin
    If X > 0 Then ConstrValue := ArcTan(Y / X)
    Else If X < 0 Then Begin
     If Y >= 0 Then ConstrValue := ArcTan(Y / X) + Pi Else ConstrValue := ArcTan(Y / X) - Pi;
    End
    Else Begin
     If Y > 0 Then ConstrValue := Pi / 2 Else ConstrValue := -Pi / 2;
    End;
   End;
  End;
 End;
End;

Function ConstrNAt(K: Integer): Integer;
Begin
 Case K Of 1: ConstrNAt := 2; 2: ConstrNAt := 3; Else ConstrNAt := 4; End;
End;

{ difference d'angle ramenee dans [-pi, pi] pour les dihedres }
Function ConstrDiff(Const C: TConstr; Target, Cur: Real): Real;
Var D: Real;
Begin
 D := Target - Cur;
 If C.Kind = 3 Then Begin
  While D > Pi Do D := D - 2 * Pi;
  While D < -Pi Do D := D + 2 * Pi;
 End;
 ConstrDiff := D;
End;

{ gradient cartesien (numerique) de la coordonnee interne }
Procedure ConstrGrad(Const C: TConstr; Var G: TCartVec);
Var I, D, K, NA: Integer; X0, H, Vp, Vm, Df: Real;
Begin
 FillChar(G, SizeOf(G), 0);
 H := 1E-5;
 NA := ConstrNAt(C.Kind);
 For I := 1 To NA Do
  For D := 1 To 3 Do Begin
   K := (C.A[I] - 1) * 3 + D;
   X0 := AtX[C.A[I], D];
   AtX[C.A[I], D] := X0 + H; Vp := ConstrValue(C);
   AtX[C.A[I], D] := X0 - H; Vm := ConstrValue(C);
   AtX[C.A[I], D] := X0;
   Df := Vp - Vm;
   If C.Kind = 3 Then Begin
    If Df > Pi Then Df := Df - 2 * Pi;
    If Df < -Pi Then Df := Df + 2 * Pi;
   End;
   G[K] := Df / (2 * H);
  End;
End;

{ ramene la geometrie sur les valeurs cibles des contraintes }
Procedure EnforceConstraints;
Var
 It, K, I, N3: Integer; G: TCartVec; Cn, Df: Real; Done: Boolean;
Begin
 N3 := 3 * NAt;
 For It := 1 To 30 Do Begin
  Done := True;
  For K := 1 To NConstr Do Begin
   Df := ConstrDiff(Constr[K], Constr[K].Val, ConstrValue(Constr[K]));
   If Abs(Df) > 1E-9 Then Begin
    Done := False;
    ConstrGrad(Constr[K], G);
    Cn := 0;
    For I := 1 To N3 Do Cn := Cn + G[I] * G[I];
    If Cn > 1E-12 Then
     For I := 1 To N3 Do AtX[(I - 1) Div 3 + 1, (I - 1) Mod 3 + 1] :=
      AtX[(I - 1) Div 3 + 1, (I - 1) Mod 3 + 1] + G[I] * Df / Cn;
   End;
  End;
  If Done Then Break;
 End;
End;

{ --------------------------------------------------------------------- }
{ Gradient numerique                                                     }
{ --------------------------------------------------------------------- }

Var
 UB    : PCartMat;                     { base orthonormee des deplacements cartesiens }
 DuTR  : Array[1..6, 1..3] Of Real;    { derivee analytique du dipole le long des vecteurs de translation / rotation }
 NTR   : Integer;                      { nombre de vecteurs de translation / rotation }

{ Construit une base orthonormee de l'espace des deplacements cartesiens.
  Les NTR premiers vecteurs engendrent les translations et les rotations
  (invariances de l'energie), les suivants sont les deplacements internes.
  Si WithDip, DuTR recoit les derivees du dipole Dip[] (translation :
  charge ; rotation : e x mu) le long de ces vecteurs. }
Procedure BuildDispBasis(WithDip: Boolean);
Var
 N3, K, I, J, A, M, Best: Integer;
 V: TCartVec; Dv: TVec3; Nn, N0, Tmp, BestR, R2: Real;
 Used: Array[1..MAX_CONS] Of Boolean;
 Cnt: Integer;
Begin
 N3 := 3 * NAt;
 If UB = Nil Then New(UB);
 NTR := 0;
 For K := 1 To 6 Do Begin
  FillChar(V, SizeOf(V), 0);
  Dv[1] := 0; Dv[2] := 0; Dv[3] := 0;
  If K <= 3 Then Begin
   For A := 1 To NAt Do V[(A - 1) * 3 + K] := 1;
   Dv[K] := Charge;
  End
  Else Begin
   For A := 1 To NAt Do
    Case K Of
     4: Begin V[(A-1)*3+2] := -AtX[A,3]; V[(A-1)*3+3] := AtX[A,2]; End;
     5: Begin V[(A-1)*3+1] := AtX[A,3]; V[(A-1)*3+3] := -AtX[A,1]; End;
    Else Begin V[(A-1)*3+1] := -AtX[A,2]; V[(A-1)*3+2] := AtX[A,1]; End;
    End;
   Case K Of
    4: Begin Dv[2] := -Dip[3]; Dv[3] := Dip[2]; End;
    5: Begin Dv[1] := Dip[3]; Dv[3] := -Dip[1]; End;
   Else Begin Dv[1] := -Dip[2]; Dv[2] := Dip[1]; End;
   End;
  End;
  N0 := 0;
  For I := 1 To N3 Do N0 := N0 + V[I] * V[I];
  N0 := Sqrt(N0);
  If N0 < 1E-8 Then Continue;
  For M := 1 To 2 Do
   For J := 1 To NTR Do Begin
    Tmp := 0;
    For I := 1 To N3 Do Tmp := Tmp + V[I] * UB^[J, I];
    For I := 1 To N3 Do V[I] := V[I] - Tmp * UB^[J, I];
    For I := 1 To 3 Do Dv[I] := Dv[I] - Tmp * DuTR[J, I];
   End;
  Nn := 0;
  For I := 1 To N3 Do Nn := Nn + V[I] * V[I];
  Nn := Sqrt(Nn);
  If Nn < 1E-5 * N0 Then Continue;
  Inc(NTR);
  For I := 1 To N3 Do UB^[NTR, I] := V[I] / Nn;
  For I := 1 To 3 Do DuTR[NTR, I] := Dv[I] / Nn;
 End;
 { completion par des vecteurs unitaires cartesiens (residu maximal) }
 For I := 1 To N3 Do Used[I] := False;
 Cnt := NTR;
 While Cnt < N3 Do Begin
  Best := 0; BestR := -1;
  For I := 1 To N3 Do
   If Not Used[I] Then Begin
    R2 := 1;
    For J := 1 To Cnt Do R2 := R2 - Sqr(UB^[J, I]);
    If R2 > BestR Then Begin BestR := R2; Best := I; End;
   End;
  If Best = 0 Then Break;
  Used[Best] := True;
  FillChar(V, SizeOf(V), 0);
  V[Best] := 1;
  For M := 1 To 2 Do
   For J := 1 To Cnt Do Begin
    Tmp := 0;
    For I := 1 To N3 Do Tmp := Tmp + V[I] * UB^[J, I];
    For I := 1 To N3 Do V[I] := V[I] - Tmp * UB^[J, I];
   End;
  Nn := 0;
  For I := 1 To N3 Do Nn := Nn + V[I] * V[I];
  Nn := Sqrt(Nn);
  If Nn < 1E-6 Then Continue;
  Inc(Cnt);
  For I := 1 To N3 Do UB^[Cnt, I] := V[I] / Nn;
 End;
End;

Procedure ShiftAlong(K: Integer; H: Real);
Var I: Integer;
Begin
 For I := 1 To 3 * NAt Do
  AtX[(I - 1) Div 3 + 1, (I - 1) Mod 3 + 1] := AtX[(I - 1) Div 3 + 1, (I - 1) Mod 3 + 1] + H * UB^[K, I];
End;

{ Gradient numerique par differences centrees le long des NV deplacements
  internes (3N-6 evaluations doubles au lieu de 3N) }
Function NumGradient(Var G: TCartVec; Var Err: String): Boolean;
Var K, I, N3: Integer; Ep, Em, H, Gk: Real; X0: Array[1..MAX_AT, 1..3] Of Real; A, D: Integer;
Begin
 NumGradient := False;
 N3 := 3 * NAt;
 H := NumHStep;
 FillChar(G, SizeOf(G), 0);
 BuildDispBasis(False);
 For A := 1 To NAt Do For D := 1 To 3 Do X0[A, D] := AtX[A, D];
 For K := NTR + 1 To N3 Do Begin
  ShiftAlong(K, H);
  If Not EvalEnergy(False, True, Ep, Err) Then Begin
   For A := 1 To NAt Do For D := 1 To 3 Do AtX[A, D] := X0[A, D];
   Exit;
  End;
  For A := 1 To NAt Do For D := 1 To 3 Do AtX[A, D] := X0[A, D];
  ShiftAlong(K, -H);
  If Not EvalEnergy(False, True, Em, Err) Then Begin
   For A := 1 To NAt Do For D := 1 To 3 Do AtX[A, D] := X0[A, D];
   Exit;
  End;
  For A := 1 To NAt Do For D := 1 To 3 Do AtX[A, D] := X0[A, D];
  Gk := (Ep - Em) / (2 * H);
  For I := 1 To N3 Do G[I] := G[I] + Gk * UB^[K, I];
  ShowProgress('Gradient numerique ' + IStr(K - NTR) + '/' + IStr(N3 - NTR), 0);
  If UserAbort Then Begin Err := 'Calcul interrompu'; Exit; End;
 End;
 NumGradient := True;
End;

{ --------------------------------------------------------------------- }
{ Optimisation de geometrie                                              }
{ --------------------------------------------------------------------- }

Procedure ModelHessian(H: PCartMat);
Var
 N3, I, J, K, A, B, C, D1, D2: Integer;
 Rij, Rref, Rho, Rho2, R2: Real;
 V: TCartVec; U, W: TVec3; Lu, Lv, Cs, Sn: Real;
 Procedure AddOuter(Wt: Real);
 Var P, Q: Integer;
 Begin
  For P := 1 To N3 Do
   If V[P] <> 0 Then
    For Q := 1 To N3 Do
     If V[Q] <> 0 Then H^[P, Q] := H^[P, Q] + Wt * V[P] * V[Q];
 End;
Begin
 N3 := 3 * NAt;
 FillChar(H^, SizeOf(TCartMat), 0);
 For I := 1 To N3 Do H^[I, I] := 0.02;
 { liaisons (toutes les paires, ponderees par la distance) }
 For A := 1 To NAt - 1 Do
  For B := A + 1 To NAt Do Begin
   Rij := Dist(A, B);
   Rref := (ElRad[AtZ[A]] + ElRad[AtZ[B]]) * ANG2BOHR;
   Rho := Exp(0.3 * (Rref * Rref - Rij * Rij));
   If Rho < 1E-4 Then Continue;
   FillChar(V, SizeOf(V), 0);
   For D1 := 1 To 3 Do Begin
    V[(A - 1) * 3 + D1] := (AtX[A, D1] - AtX[B, D1]) / Rij;
    V[(B - 1) * 3 + D1] := -(AtX[A, D1] - AtX[B, D1]) / Rij;
   End;
   AddOuter(0.45 * Rho);
  End;
 { angles }
 For B := 1 To NAt Do
  For A := 1 To NAt Do
   For C := A + 1 To NAt Do Begin
    If (A = B) Or (C = B) Then Continue;
    Rho := Exp(0.3 * (Sqr((ElRad[AtZ[A]] + ElRad[AtZ[B]]) * ANG2BOHR) - Sqr(Dist(A, B))));
    Rho2 := Exp(0.3 * (Sqr((ElRad[AtZ[C]] + ElRad[AtZ[B]]) * ANG2BOHR) - Sqr(Dist(C, B))));
    If Rho * Rho2 < 1E-3 Then Continue;
    For D1 := 1 To 3 Do Begin
     U[D1] := AtX[A, D1] - AtX[B, D1];
     W[D1] := AtX[C, D1] - AtX[B, D1];
    End;
    Lu := Sqrt(U[1] * U[1] + U[2] * U[2] + U[3] * U[3]);
    Lv := Sqrt(W[1] * W[1] + W[2] * W[2] + W[3] * W[3]);
    Cs := (U[1] * W[1] + U[2] * W[2] + U[3] * W[3]) / (Lu * Lv);
    If Abs(Cs) > 0.9998 Then Continue;
    Sn := Sqrt(1 - Cs * Cs);
    FillChar(V, SizeOf(V), 0);
    For D1 := 1 To 3 Do Begin
     V[(A - 1) * 3 + D1] := (Cs * U[D1] / Lu - W[D1] / Lv) / (Lu * Sn);
     V[(C - 1) * 3 + D1] := (Cs * W[D1] / Lv - U[D1] / Lu) / (Lv * Sn);
     V[(B - 1) * 3 + D1] := -V[(A - 1) * 3 + D1] - V[(C - 1) * 3 + D1];
    End;
    AddOuter(0.15 * Rho * Rho2);
   End;
End;

Function ConvText(OK: Boolean): String;
Begin
 If OK Then ConvText := 'YES' Else ConvText := 'NO';
End;

Procedure PrintGeomAng(Const Title: String);
Var A: Integer;
Begin
 EmitBlank;
 Emit(Title);
 Emit('---------------------------------');
 For A := 1 To NAt Do
  Emit(PadR(ElemStr(AtZ[A]), 4) + FmtR(AtX[A, 1] / ANG2BOHR, 16, 9) +
       FmtR(AtX[A, 2] / ANG2BOHR, 16, 9) + FmtR(AtX[A, 3] / ANG2BOHR, 16, 9));
End;

Procedure SetOptTol(Const Name: String);
Var N: String;
Begin
 N := UpStr(Name);
 OptName := Name;
 If N = 'TIGHTOPT' Then Begin
  OptTolE := 1E-6; OptTolRG := 3E-5; OptTolMG := 1E-4; OptTolRS := 6E-4; OptTolMS := 1E-3;
 End
 Else If N = 'VERYTIGHTOPT' Then Begin
  OptTolE := 2E-7; OptTolRG := 2E-6; OptTolMG := 1E-5; OptTolRS := 1E-4; OptTolMS := 2E-4;
 End
 Else If N = 'LOOSEOPT' Then Begin
  OptTolE := 3E-5; OptTolRG := 5E-4; OptTolMG := 2.5E-3; OptTolRS := 1E-2; OptTolMS := 1.5E-2;
 End
 Else Begin
  OptTolE := 5E-6; OptTolRG := 1E-4; OptTolMG := 3E-4; OptTolRS := 2E-3; OptTolMS := 4E-3;
  OptName := 'NormalOpt';
 End;
End;

{ Optimisation de la geometrie (avec contraintes eventuelles).
  Verbose : impression du detail ; renvoie Faux en cas d'echec. }
Function RunOpt(Verbose: Boolean; Var Err: String): Boolean;
Var
 N3, Cyc, I, J, K, Tries: Integer;
 H, HV, Wk: PCartMat;
 X, XOld, G, GOld, S, Y, HS: TCartVec; Lam: Array[1..MAX_CONS] Of Real;
 CG: Array[1..MAX_CONSTR] Of TCartVec;
 E, EOld, EN, Trust, SNorm, SY, SHS, Tmp, RmsG, MaxG, RmsS, MaxS, EChg: Real;
 Conv, Accepted, Fail: Boolean;
 Procedure LoadX;
 Var P: Integer;
 Begin
  For P := 1 To N3 Do AtX[(P - 1) Div 3 + 1, (P - 1) Mod 3 + 1] := X[P];
 End;
 Procedure StoreX;
 Var P: Integer;
 Begin
  For P := 1 To N3 Do X[P] := AtX[(P - 1) Div 3 + 1, (P - 1) Mod 3 + 1];
 End;
 Procedure ProjectOut(Var V: TCartVec);
 Var P, Q: Integer; D: Real;
 Begin
  For P := 1 To NConstr Do Begin
   D := 0;
   For Q := 1 To N3 Do D := D + V[Q] * CG[P][Q];
   For Q := 1 To N3 Do V[Q] := V[Q] - D * CG[P][Q];
  End;
 End;
 Procedure BuildConstrVectors;
 Var P, Q, R: Integer; D, Nn: Real; Tmp2: TCartVec;
 Begin
  For P := 1 To NConstr Do Begin
   ConstrGrad(Constr[P], Tmp2);
   For R := 1 To P - 1 Do Begin
    D := 0;
    For Q := 1 To N3 Do D := D + Tmp2[Q] * CG[R][Q];
    For Q := 1 To N3 Do Tmp2[Q] := Tmp2[Q] - D * CG[R][Q];
   End;
   Nn := 0;
   For Q := 1 To N3 Do Nn := Nn + Tmp2[Q] * Tmp2[Q];
   Nn := Sqrt(Nn);
   If Nn > 1E-10 Then For Q := 1 To N3 Do CG[P][Q] := Tmp2[Q] / Nn
   Else For Q := 1 To N3 Do CG[P][Q] := 0;
  End;
 End;
 Procedure Release;
 Begin
  Dispose(H); Dispose(HV); Dispose(Wk);
 End;
Begin
 RunOpt := False;
 Err := '';
 N3 := 3 * NAt;
 New(H); New(HV); New(Wk);
 ModelHessian(H);
 OptN := 0; OptConv := False;
 Trust := 0.30;
 StoreX;
 If NConstr > 0 Then Begin EnforceConstraints; StoreX; End;
 If Not EvalEnergy(Verbose, True, E, Err) Then Begin Release; Exit; End;
 OptE[0] := E; OptN := 1;
 EOld := E;
 Conv := False; Fail := False;
 Cyc := 0;
 While Not Conv Do Begin
  Inc(Cyc);
  If Verbose Then Begin
   EmitBlank;
   Emit('*                GEOMETRY OPTIMIZATION CYCLE ' + FmtR(Cyc, 3, 0) + '                *');
   PrintGeomAng('CARTESIAN COORDINATES (ANGSTROEM)');
  End;
  If Not NumGradient(G, Err) Then Begin LoadX; Release; Exit; End;
  If NConstr > 0 Then Begin BuildConstrVectors; ProjectOut(G); End;
  RmsG := 0; MaxG := 0;
  For I := 1 To N3 Do Begin
   RmsG := RmsG + G[I] * G[I];
   If Abs(G[I]) > MaxG Then MaxG := Abs(G[I]);
  End;
  RmsG := Sqrt(RmsG / N3);
  LastGrad := G;
  If Cyc > 1 Then Begin
   For I := 1 To N3 Do Begin S[I] := X[I] - XOld[I]; Y[I] := G[I] - GOld[I]; End;
   SY := 0;
   For I := 1 To N3 Do SY := SY + S[I] * Y[I];
   For I := 1 To N3 Do Begin
    Tmp := 0;
    For J := 1 To N3 Do Tmp := Tmp + H^[I, J] * S[J];
    HS[I] := Tmp;
   End;
   SHS := 0;
   For I := 1 To N3 Do SHS := SHS + S[I] * HS[I];
   If (SY > 1E-8) And (SHS > 1E-12) Then
    For I := 1 To N3 Do
     For J := 1 To N3 Do
      H^[I, J] := H^[I, J] + Y[I] * Y[J] / SY - HS[I] * HS[J] / SHS;
  End;
  { pas de Newton dans la base propre de H }
  For I := 1 To N3 Do
   For J := 1 To N3 Do Wk^[I, J] := 0.5 * (H^[I, J] + H^[J, I]);
  Jacobi(N3, MAX_CONS, PRealArr(Wk), PRealArr(HV), PRealArr(@Lam));
  FillChar(S, SizeOf(S), 0);
  For K := 1 To N3 Do Begin
   If Abs(Lam[K]) < 1E-6 Then Continue;
   Tmp := 0;
   For I := 1 To N3 Do Tmp := Tmp + HV^[I, K] * G[I];
   Tmp := -Tmp / RMax(Abs(Lam[K]), 0.02);
   For I := 1 To N3 Do S[I] := S[I] + Tmp * HV^[I, K];
  End;
  If NConstr > 0 Then ProjectOut(S);
  SNorm := 0;
  For I := 1 To N3 Do SNorm := SNorm + S[I] * S[I];
  SNorm := Sqrt(SNorm);
  If SNorm > Trust Then Begin
   For I := 1 To N3 Do S[I] := S[I] * Trust / SNorm;
   SNorm := Trust;
  End;
  RmsS := 0; MaxS := 0;
  For I := 1 To N3 Do Begin
   RmsS := RmsS + S[I] * S[I];
   If Abs(S[I]) > MaxS Then MaxS := Abs(S[I]);
  End;
  RmsS := Sqrt(RmsS / N3);
  EChg := E - EOld;
  Conv := (Abs(EChg) < OptTolE) And (RmsG < OptTolRG) And (MaxG < OptTolMG) And
          (RmsS < OptTolRS) And (MaxS < OptTolMS);
  If Verbose Then Begin
   EmitBlank;
   Emit('FINAL SINGLE POINT ENERGY' + FmtR(E, 24, 12));
   EmitBlank;
   Emit('                     .--------------------.');
   Emit('---------------------|Geometry convergence|-------------------------');
   Emit('Item                value                   Tolerance       Converged');
   Emit('---------------------------------------------------------------------');
   Emit('Energy change      ' + FmtR(EChg, 14, 8) + '      ' + FmtR(OptTolE, 14, 8) + '        ' + ConvText(Abs(EChg) < OptTolE));
   Emit('RMS gradient       ' + FmtR(RmsG, 14, 8) + '      ' + FmtR(OptTolRG, 14, 8) + '        ' + ConvText(RmsG < OptTolRG));
   Emit('MAX gradient       ' + FmtR(MaxG, 14, 8) + '      ' + FmtR(OptTolMG, 14, 8) + '        ' + ConvText(MaxG < OptTolMG));
   Emit('RMS step           ' + FmtR(RmsS, 14, 8) + '      ' + FmtR(OptTolRS, 14, 8) + '        ' + ConvText(RmsS < OptTolRS));
   Emit('MAX step           ' + FmtR(MaxS, 14, 8) + '      ' + FmtR(OptTolMS, 14, 8) + '        ' + ConvText(MaxS < OptTolMS));
   Emit('---------------------------------------------------------------------');
  End;
  ShowProgress('Optimisation : cycle ' + IStr(Cyc) + '  E = ' + RStr(E, 8), 0);
  If UserAbort Then Begin Err := 'Calcul interrompu'; Fail := True; Break; End;
  If Conv Then Break;
  If Cyc >= OptMaxCyc Then Begin Fail := True; Break; End;
  { deplacement avec reduction du rayon de confiance si l'energie monte }
  XOld := X; GOld := G; EOld := E;
  Accepted := False; Tries := 0;
  While Not Accepted Do Begin
   Inc(Tries);
   For I := 1 To N3 Do X[I] := XOld[I] + S[I];
   LoadX;
   If NConstr > 0 Then Begin EnforceConstraints; StoreX; End;
   If Not EvalEnergy(False, True, EN, Err) Then Begin
    If UserAbort Then Begin LoadX; Release; Exit; End;
    EN := EOld + 1;
   End;
   If (EN <= EOld + 1E-7) Or (Tries >= 6) Then Accepted := True
   Else Begin
    Trust := Trust * 0.5;
    For I := 1 To N3 Do S[I] := S[I] * 0.5;
    X := XOld; LoadX;
   End;
  End;
  If (Tries = 1) And (SNorm >= 0.9 * Trust) Then Trust := RMin(Trust * 1.5, 0.5);
  E := EN;
  If OptN <= MAX_OPT Then Begin OptE[OptN] := E; Inc(OptN); End;
 End;
 { energie a la geometrie finale (met a jour les orbitales) }
 LoadX;
 If Not EvalEnergy(False, True, E, Err) Then Begin Release; Exit; End;
 If OptN <= MAX_OPT Then OptE[OptN - 1] := E;
 Release;
 OptConv := Conv;
 If Fail Then Begin
  If Err = '' Then Err := 'Nombre maximal de cycles d''optimisation atteint';
  RunOpt := False;
 End
 Else RunOpt := True;
End;

{ --------------------------------------------------------------------- }
{ Frequences harmoniques et thermochimie                                 }
{ --------------------------------------------------------------------- }

Var
 DipDer : Array[1..MAX_CONS, 1..3] Of Real;      { d mu / d x (u.a.) }
 IRTx, IRTy, IRTz : Array[1..MAX_CONS] Of Real;   { d mu / dQ par mode }
 NVib   : Integer;
 NProj  : Integer;
 Linear : Boolean;
 FreqE0 : Real;

{ Nombre de symetrie rotationnelle par recherche des rotations propres }
Function DetectSymNum: Integer;
Type
 TRot = Array[1..3, 1..3] Of Real;
Var
 Cm: TVec3; Mt, M: Real;
 R: Array[1..MAX_AT] Of TVec3;
 Ax: Array[1..400] Of TVec3; NAx: Integer;
 Rots: Array[1..120] Of TRot; NR: Integer;
 A, B, I, J, N, K: Integer; V: TVec3; L: Real;

 Procedure AddAxis(X, Y, Z: Real);
 Var Q, T: Integer; Ln2, D: Real;
 Begin
  Ln2 := Sqrt(X * X + Y * Y + Z * Z);
  If Ln2 < 1E-4 Then Exit;
  X := X / Ln2; Y := Y / Ln2; Z := Z / Ln2;
  For Q := 1 To NAx Do Begin
   D := Abs(Ax[Q][1] * X + Ax[Q][2] * Y + Ax[Q][3] * Z);
   If D > 0.9999 Then Exit;
  End;
  If NAx < 400 Then Begin
   Inc(NAx);
   Ax[NAx][1] := X; Ax[NAx][2] := Y; Ax[NAx][3] := Z;
  End;
 End;

 Function MapsOnto(Const Rm: TRot): Boolean;
 Var P, Q, T: Integer; X, Y, Z: Real; Found: Boolean;
 Begin
  MapsOnto := False;
  For P := 1 To NAt Do Begin
   X := Rm[1,1] * R[P][1] + Rm[1,2] * R[P][2] + Rm[1,3] * R[P][3];
   Y := Rm[2,1] * R[P][1] + Rm[2,2] * R[P][2] + Rm[2,3] * R[P][3];
   Z := Rm[3,1] * R[P][1] + Rm[3,2] * R[P][2] + Rm[3,3] * R[P][3];
   Found := False;
   For Q := 1 To NAt Do
    If (AtZ[Q] = AtZ[P]) And
       (Sqr(X - R[Q][1]) + Sqr(Y - R[Q][2]) + Sqr(Z - R[Q][3]) < 4E-4) Then Found := True;
   If Not Found Then Exit;
  End;
  MapsOnto := True;
 End;

 Procedure MakeRot(Const U: TVec3; Ang: Real; Var Rm: TRot);
 Var C, S, T: Real;
 Begin
  C := Cos(Ang); S := Sin(Ang); T := 1 - C;
  Rm[1,1] := C + U[1] * U[1] * T;       Rm[1,2] := U[1] * U[2] * T - U[3] * S; Rm[1,3] := U[1] * U[3] * T + U[2] * S;
  Rm[2,1] := U[2] * U[1] * T + U[3] * S; Rm[2,2] := C + U[2] * U[2] * T;       Rm[2,3] := U[2] * U[3] * T - U[1] * S;
  Rm[3,1] := U[3] * U[1] * T - U[2] * S; Rm[3,2] := U[3] * U[2] * T + U[1] * S; Rm[3,3] := C + U[3] * U[3] * T;
 End;

 Function KnownRot(Const Rm: TRot): Boolean;
 Var Q, P, S: Integer; Same: Boolean;
 Begin
  KnownRot := False;
  For Q := 1 To NR Do Begin
   Same := True;
   For P := 1 To 3 Do
    For S := 1 To 3 Do If Abs(Rots[Q][P, S] - Rm[P, S]) > 1E-3 Then Same := False;
   If Same Then Begin KnownRot := True; Exit; End;
  End;
 End;

 Procedure TryAxis(Const U: TVec3);
 Var Nn: Integer; Rm: TRot; Kk: Integer;
 Begin
  For Nn := 2 To 6 Do Begin
   MakeRot(U, 2 * Pi / Nn, Rm);
   If MapsOnto(Rm) Then
    For Kk := 1 To Nn - 1 Do Begin
     MakeRot(U, 2 * Pi * Kk / Nn, Rm);
     If Not KnownRot(Rm) And (NR < 120) Then Begin
      Inc(NR); Rots[NR] := Rm;
     End;
    End;
  End;
 End;

Begin
 If NAt < 2 Then Begin DetectSymNum := 1; Exit; End;
 Mt := 0; Cm[1] := 0; Cm[2] := 0; Cm[3] := 0;
 For A := 1 To NAt Do Begin
  M := ElMass[AtZ[A]]; Mt := Mt + M;
  For I := 1 To 3 Do Cm[I] := Cm[I] + M * AtX[A, I];
 End;
 For I := 1 To 3 Do Cm[I] := Cm[I] / Mt;
 For A := 1 To NAt Do For I := 1 To 3 Do R[A][I] := AtX[A, I] - Cm[I];
 NAx := 0; NR := 0;
 For A := 1 To NAt Do AddAxis(R[A][1], R[A][2], R[A][3]);
 For A := 1 To NAt - 1 Do
  For B := A + 1 To NAt Do Begin
   If AtZ[A] = AtZ[B] Then AddAxis(R[A][1] + R[B][1], R[A][2] + R[B][2], R[A][3] + R[B][3]);
   AddAxis(R[A][2] * R[B][3] - R[A][3] * R[B][2], R[A][3] * R[B][1] - R[A][1] * R[B][3],
           R[A][1] * R[B][2] - R[A][2] * R[B][1]);
  End;
 For I := 1 To NAx Do TryAxis(Ax[I]);
 DetectSymNum := NR + 1;
End;

Function RunFreq(Verbose: Boolean; Var Err: String): Boolean;
Var
 N3, NV, I, J, K, L, A, D: Integer;
 H, HW, HV, PM, HC: PCartMat;
 Ep, Em: Array[1..MAX_CONS] Of Real;
 E0, E1, E2, Hs, Step, Tmp: Real;
 DK: Array[1..MAX_CONS, 1..3] Of Real;
 Lam: Array[1..MAX_CONS] Of Real;
 TV: Array[1..6] Of TCartVec; NT: Integer;
 Cm: TVec3; Mt, Mi: Real; V: TCartVec; Nn: Real;
 Mass: Array[1..MAX_AT] Of Real;
 X0: Array[1..MAX_AT, 1..3] Of Real;
 Procedure Release;
 Begin
  Dispose(H); Dispose(HW); Dispose(HV); Dispose(PM); Dispose(HC);
 End;
 Procedure Reset0;
 Var P, Q: Integer;
 Begin
  For P := 1 To NAt Do For Q := 1 To 3 Do AtX[P, Q] := X0[P, Q];
 End;
Begin
 RunFreq := False;
 FreqOK := False; Err := '';
 N3 := 3 * NAt;
 Step := 0.02;
 New(H); New(HW); New(HV); New(PM); New(HC);
 FillChar(H^, SizeOf(TCartMat), 0);
 For A := 1 To NAt Do Mass[A] := ElMass[AtZ[A]];
 If Not EvalEnergy(Verbose, True, E0, Err) Then Begin Release; Exit; End;
 FreqE0 := E0;
 ComputeDipole;
 BuildDispBasis(True);
 NV := N3 - NTR;
 For A := 1 To NAt Do For D := 1 To 3 Do X0[A, D] := AtX[A, D];
 { energies et dipoles le long des deplacements internes }
 For K := 1 To NV Do Begin
  ShiftAlong(NTR + K, Step);
  If Not EvalEnergy(False, True, Ep[K], Err) Then Begin Reset0; Release; Exit; End;
  ComputeDipole;
  For L := 1 To 3 Do DK[K, L] := Dip[L];
  Reset0;
  ShiftAlong(NTR + K, -Step);
  If Not EvalEnergy(False, True, Em[K], Err) Then Begin Reset0; Release; Exit; End;
  ComputeDipole;
  For L := 1 To 3 Do DK[K, L] := (DK[K, L] - Dip[L]) / (2 * Step);
  Reset0;
  H^[K, K] := (Ep[K] + Em[K] - 2 * E0) / (Step * Step);
  ShowProgress('Frequences : deplacement ' + IStr(K) + '/' + IStr(NV), 0);
  If UserAbort Then Begin Err := 'Calcul interrompu'; Release; Exit; End;
 End;
 For K := 1 To NV - 1 Do
  For L := K + 1 To NV Do Begin
   ShiftAlong(NTR + K, Step); ShiftAlong(NTR + L, Step);
   If Not EvalEnergy(False, True, E1, Err) Then Begin Reset0; Release; Exit; End;
   Reset0;
   ShiftAlong(NTR + K, -Step); ShiftAlong(NTR + L, -Step);
   If Not EvalEnergy(False, True, E2, Err) Then Begin Reset0; Release; Exit; End;
   Reset0;
   Hs := (E1 + E2 - Ep[K] - Em[K] - Ep[L] - Em[L] + 2 * E0) / (2 * Step * Step);
   H^[K, L] := Hs; H^[L, K] := Hs;
   ShowProgress('Frequences : couple ' + IStr(K) + ',' + IStr(L), 0);
   If UserAbort Then Begin Err := 'Calcul interrompu'; Release; Exit; End;
  End;
 { retour a la geometrie de reference (orbitales a jour) }
 If Not EvalEnergy(False, True, E0, Err) Then Begin Release; Exit; End;
 { hessien cartesien : sum V_k H_kl V_l^T ; derivees du dipole cartesiennes }
 For K := 1 To NV Do
  For J := 1 To N3 Do Begin
   Tmp := 0;
   For L := 1 To NV Do Tmp := Tmp + H^[K, L] * UB^[NTR + L, J];
   HC^[K, J] := Tmp;
  End;
 For I := 1 To N3 Do
  For J := 1 To N3 Do Begin
   Tmp := 0;
   For K := 1 To NV Do Tmp := Tmp + UB^[NTR + K, I] * HC^[K, J];
   H^[I, J] := Tmp;
  End;
 For I := 1 To N3 Do
  For L := 1 To 3 Do Begin
   Tmp := 0;
   For K := 1 To NV Do Tmp := Tmp + UB^[NTR + K, I] * DK[K, L];
   For K := 1 To NTR Do Tmp := Tmp + UB^[K, I] * DuTR[K, L];
   DipDer[I, L] := Tmp;
  End;
 { ponderation par les masses }
 For I := 1 To N3 Do
  For J := 1 To N3 Do
   HW^[I, J] := H^[I, J] / Sqrt(Mass[(I - 1) Div 3 + 1] * Mass[(J - 1) Div 3 + 1]);
 { vecteurs de translation et de rotation }
 Mt := 0; Cm[1] := 0; Cm[2] := 0; Cm[3] := 0;
 For A := 1 To NAt Do Begin
  Mt := Mt + Mass[A];
  For D := 1 To 3 Do Cm[D] := Cm[D] + Mass[A] * AtX[A, D];
 End;
 For D := 1 To 3 Do Cm[D] := Cm[D] / Mt;
 NT := 0;
 For K := 1 To 6 Do Begin
  FillChar(V, SizeOf(V), 0);
  For A := 1 To NAt Do Begin
   Mi := Sqrt(Mass[A]);
   If K <= 3 Then V[(A - 1) * 3 + K] := Mi
   Else Begin
    Case K Of
     4: Begin V[(A-1)*3+2] := -Mi * (AtX[A,3] - Cm[3]); V[(A-1)*3+3] := Mi * (AtX[A,2] - Cm[2]); End;
     5: Begin V[(A-1)*3+3] := -Mi * (AtX[A,1] - Cm[1]); V[(A-1)*3+1] := Mi * (AtX[A,3] - Cm[3]); End;
    Else Begin V[(A-1)*3+1] := -Mi * (AtX[A,2] - Cm[2]); V[(A-1)*3+2] := Mi * (AtX[A,1] - Cm[1]); End;
    End;
   End;
  End;
  { orthogonalisation }
  For J := 1 To NT Do Begin
   Tmp := 0;
   For I := 1 To N3 Do Tmp := Tmp + V[I] * TV[J][I];
   For I := 1 To N3 Do V[I] := V[I] - Tmp * TV[J][I];
  End;
  Nn := 0;
  For I := 1 To N3 Do Nn := Nn + V[I] * V[I];
  Nn := Sqrt(Nn);
  If (Nn > 1E-3) And (NT < 6) Then Begin
   Inc(NT);
   For I := 1 To N3 Do TV[NT][I] := V[I] / Nn;
  End;
 End;
 NProj := NT;
 Linear := (NT = 5);
 { projection : P H P + 100 * (1 - P) }
 For I := 1 To N3 Do
  For J := 1 To N3 Do Begin
   PM^[I, J] := 0;
   If I = J Then PM^[I, J] := 1;
   For K := 1 To NT Do PM^[I, J] := PM^[I, J] - TV[K][I] * TV[K][J];
  End;
 For I := 1 To N3 Do
  For J := 1 To N3 Do Begin
   Tmp := 0;
   For K := 1 To N3 Do Tmp := Tmp + PM^[I, K] * HW^[K, J];
   H^[I, J] := Tmp;
  End;
 For I := 1 To N3 Do
  For J := 1 To N3 Do Begin
   Tmp := 0;
   For K := 1 To N3 Do Tmp := Tmp + H^[I, K] * PM^[K, J];
   HW^[I, J] := Tmp;
  End;
 For I := 1 To N3 Do
  For J := I To N3 Do Begin
   Tmp := 0;
   For K := 1 To NT Do Tmp := Tmp + 100 * TV[K][I] * TV[K][J];
   Tmp := 0.5 * (HW^[I, J] + HW^[J, I]) + Tmp;
   HW^[I, J] := Tmp; HW^[J, I] := Tmp;
  End;
 If FreqVec = Nil Then New(FreqVec);
 Jacobi(N3, MAX_CONS, PRealArr(HW), PRealArr(HV), PRealArr(@Lam));
 NVib := N3 - NT;
 NFreq := N3;
 For K := 1 To N3 Do Begin
  If K <= NVib Then Begin
   If Lam[K] >= 0 Then FreqV[K + NT] := Sqrt(Lam[K]) * FREQ_CONV
   Else FreqV[K + NT] := -Sqrt(-Lam[K]) * FREQ_CONV;
   For I := 1 To N3 Do FreqVec^[I, K + NT] := HV^[I, K];
  End;
 End;
 For K := 1 To NT Do Begin FreqV[K] := 0; For I := 1 To N3 Do FreqVec^[I, K] := TV[K][I]; End;
 { intensites IR }
 For K := 1 To N3 Do Begin
  IRTx[K] := 0; IRTy[K] := 0; IRTz[K] := 0;
  If K > NT Then
   For I := 1 To N3 Do Begin
    Tmp := FreqVec^[I, K] / Sqrt(Mass[(I - 1) Div 3 + 1]);
    IRTx[K] := IRTx[K] + DipDer[I, 1] * Tmp;
    IRTy[K] := IRTy[K] + DipDer[I, 2] * Tmp;
    IRTz[K] := IRTz[K] + DipDer[I, 3] * Tmp;
   End;
  FreqI[K] := IR_CONV * (IRTx[K] * IRTx[K] + IRTy[K] * IRTy[K] + IRTz[K] * IRTz[K]);
  If K <= NT Then FreqI[K] := 0;
 End;
 Release;
 FreqOK := True;
 RunFreq := True;
End;

{ Thermochimie (oscillateur harmonique, rotateur rigide, gaz parfait) }
Procedure RunThermo(T: Real);
Var
 K, A, I, J: Integer; Mt, Mi, X, Th: Real;
 Cm: TVec3; Inert: TMat; IW: Array[1..3] Of Real; IV: TMat; I3: TMat;
 R: TVec3; Lam2, QT, QR, KT: Real; NI: Integer; Idg: Array[1..3] Of Real;
 NR: Integer;
Begin
 ThTemp := T;
 KT := KB_HARTREE * T;
 Mt := 0; Cm[1] := 0; Cm[2] := 0; Cm[3] := 0;
 For A := 1 To NAt Do Begin
  Mi := ElMass[AtZ[A]]; Mt := Mt + Mi;
  For I := 1 To 3 Do Cm[I] := Cm[I] + Mi * AtX[A, I];
 End;
 For I := 1 To 3 Do Cm[I] := Cm[I] / Mt;
 ThMass := Mt;
 FillChar(Inert, SizeOf(Inert), 0);
 For A := 1 To NAt Do Begin
  Mi := ElMass[AtZ[A]];
  For I := 1 To 3 Do R[I] := AtX[A, I] - Cm[I];
  For I := 1 To 3 Do
   For J := 1 To 3 Do Begin
    If I = J Then
     Inert[I, I] := Inert[I, I] + Mi * (R[1] * R[1] + R[2] * R[2] + R[3] * R[3] - R[I] * R[I])
    Else Inert[I, J] := Inert[I, J] - Mi * R[I] * R[J];
   End;
 End;
 Jacobi(3, MAX_BF, PRealArr(@Inert), PRealArr(@IV), PRealArr(@IW));
 Idg[1] := IW[1]; Idg[2] := IW[2]; Idg[3] := IW[3];
 { constantes rotationnelles (cm-1) }
 For I := 1 To 3 Do
  If Idg[I] > 1E-6 Then RotConst[I] := THETA_ROT / (HCK * Idg[I]) Else RotConst[I] := 0;
 SymNum := DetectSymNum;
 ZPE := 0; EVibT := 0; SVibT := 0;
 For K := NProj + 1 To NFreq Do
  If FreqV[K] > 10 Then Begin
   Th := HCK * FreqV[K];
   X := Th / T;
   ZPE := ZPE + 0.5 * KB_HARTREE * Th;
   EVibT := EVibT + KB_HARTREE * Th / (Exp(X) - 1);
   SVibT := SVibT + KB_HARTREE * (X / (Exp(X) - 1) - Ln(1 - Exp(-X)));
  End;
 { translation }
 Lam2 := LAMBDA2 / (Mt * T);
 QT := VOL_PER_K * T / Exp(1.5 * Ln(Lam2));
 STransT := KB_HARTREE * (Ln(QT) + 2.5);
 ETransT := 1.5 * KT;
 { rotation }
 If NAt = 1 Then Begin
  SRotT := 0; ERotT := 0;
 End
 Else If Linear Then Begin
  Th := THETA_ROT / Idg[3];
  If Idg[3] < Idg[2] Then Th := THETA_ROT / Idg[2];
  Th := THETA_ROT / RMax(Idg[1], RMax(Idg[2], Idg[3]));
  QR := T / (SymNum * Th);
  SRotT := KB_HARTREE * (Ln(QR) + 1);
  ERotT := KT;
 End
 Else Begin
  QR := Sqrt(Pi) / SymNum * Exp(1.5 * Ln(T)) /
        Sqrt(THETA_ROT / Idg[1] * THETA_ROT / Idg[2] * THETA_ROT / Idg[3]);
  SRotT := KB_HARTREE * (Ln(QR) + 1.5);
  ERotT := 1.5 * KT;
 End;
 STot := STransT + SRotT + SVibT + KB_HARTREE * Ln(Mult);
 EThermal := ZPE + (EVibT) + ERotT + ETransT;
 HTot := ETot + EThermal + KT;
 GTot := HTot - T * STot;
 ThermoOK := True;
End;

{ --------------------------------------------------------------------- }
{ Fichier d'entree ORCA : lecture et interpretation                      }
{ --------------------------------------------------------------------- }

Type
 TInpLine = String[100];

Var
 InpLines    : Array[0..MAX_INP-1] Of TInpLine;
 NInp        : Integer;
 InpName     : String;
 InpModified : Boolean;
 WarnList    : Array[1..12] Of String[78];
 NWarn       : Integer;
 DoCIS       : Boolean;
 DoPrintBasis, DoPrintMOs : Boolean;
 NumTight    : Boolean;
 CalcTemp    : Real;
 CalcErr     : String;
 CalcDone    : Boolean;       { un calcul termine existe }
 CalcOK      : Boolean;
 CalcMs      : LongInt;
 KeyLine     : String;
 GeomFromFile : String;
 ZMat        : Boolean;
 GeomNeedsHdr : Boolean;

Procedure AddWarn(Const S: String);
Begin
 If NWarn < 12 Then Begin Inc(NWarn); WarnList[NWarn] := Copy(S, 1, 78); End;
End;

{ extrait le premier mot (separe par espace, tabulation ou virgule) }
Function NextTok(Var S: String): String;
Var I, J: Integer;
Begin
 I := 1;
 While (I <= Length(S)) And ((S[I] = ' ') Or (S[I] = #9) Or (S[I] = ',')) Do Inc(I);
 J := I;
 While (J <= Length(S)) And (S[J] <> ' ') And (S[J] <> #9) And (S[J] <> ',') Do Inc(J);
 NextTok := Copy(S, I, J - I);
 Delete(S, 1, J - 1);
End;

{ enleve les commentaires # }
Function StripComment(Const S: String): String;
Var P: Integer;
Begin
 P := Pos('#', S);
 If P > 0 Then StripComment := Copy(S, 1, P - 1) Else StripComment := S;
End;

Function ToAngle(Kind: Integer; V: Real): Real;
Begin
 If Kind = 1 Then ToAngle := V * ANG2BOHR Else ToAngle := V * Pi / 180;
End;

{ --------------------------------------------------------------------- }
{ Coordonnees : cartesiennes et matrice Z                                }
{ --------------------------------------------------------------------- }

Function ReadXYZFile(Const FN: String; Var Err: String): Boolean;
Var
 F: Text; L, T: String; N, I, Z: Integer; X, Y, V: Real; Code: Integer;
Begin
 ReadXYZFile := False;
 {$I-}
 Assign(F, FN); Reset(F);
 {$I+}
 If IOResult <> 0 Then Begin Err := 'Fichier introuvable : ' + FN; Exit; End;
 ReadLn(F, L);
 Val(TrimStr(L), N, Code);
 If (Code <> 0) Or (N < 1) Or (N > MAX_AT) Then Begin
  Close(F); Err := 'Fichier XYZ invalide ou trop grand (max ' + IStr(MAX_AT) + ' atomes)'; Exit;
 End;
 ReadLn(F, MolTitle);
 NAt := 0;
 For I := 1 To N Do Begin
  If Eof(F) Then Begin Close(F); Err := 'Fichier XYZ incomplet'; Exit; End;
  ReadLn(F, L);
  T := NextTok(L);
  Z := ElemZ(T);
  If Z = 0 Then Begin
   Val(T, Z, Code);
   If (Code <> 0) Or (Z < 1) Or (Z > NELEM) Then Begin Close(F); Err := 'Element inconnu : ' + T; Exit; End;
  End;
  Inc(NAt);
  AtZ[NAt] := Z;
  For Code := 1 To 3 Do Begin
   T := NextTok(L);
   If Not ParseReal(T, V) Then Begin Close(F); Err := 'Coordonnee invalide : ' + T; Exit; End;
   AtX[NAt, Code] := V * ANG2BOHR;
  End;
 End;
 Close(F);
 ReadXYZFile := True;
End;

Procedure SaveXYZ(Const FN, Comment: String);
Var F: Text; A: Integer;
Begin
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 WriteLn(F, NAt);
 WriteLn(F, Comment);
 For A := 1 To NAt Do
  WriteLn(F, PadR(ElemStr(AtZ[A]), 3), AtX[A, 1] / ANG2BOHR:16:8, AtX[A, 2] / ANG2BOHR:16:8,
          AtX[A, 3] / ANG2BOHR:16:8);
 Close(F);
End;

{ place l'atome A a partir de la matrice Z (reference NA, NB, NC) }
Procedure ZPlace(A, NA, NB, NC: Integer; R, Ang, Dih: Real);
Var
 U, V, Nv, Mv, Pa, Pb, Pc, Bc, Ab, Cr: TVec3;
 I: Integer; L, Ct, St, Cp, Sp, D1, D2, D3: Real;
Begin
 If NA = 0 Then Begin
  AtX[A, 1] := 0; AtX[A, 2] := 0; AtX[A, 3] := 0;
  Exit;
 End;
 If NB = 0 Then Begin
  AtX[A, 1] := AtX[NA, 1]; AtX[A, 2] := AtX[NA, 2]; AtX[A, 3] := AtX[NA, 3] + R;
  Exit;
 End;
 Ct := Cos(Ang); St := Sin(Ang);
 For I := 1 To 3 Do Begin Pa[I] := AtX[NA, I]; Pb[I] := AtX[NB, I]; End;
 If NC = 0 Then Begin
  { dans le plan xz }
  L := Sqrt(Sqr(Pb[1] - Pa[1]) + Sqr(Pb[2] - Pa[2]) + Sqr(Pb[3] - Pa[3]));
  For I := 1 To 3 Do U[I] := (Pb[I] - Pa[I]) / L;
  { vecteur perpendiculaire a U dans le plan xz }
  V[1] := U[3]; V[2] := 0; V[3] := -U[1];
  L := Sqrt(V[1] * V[1] + V[3] * V[3]);
  If L < 1E-8 Then Begin V[1] := 1; V[2] := 0; V[3] := 0; L := 1; End;
  For I := 1 To 3 Do AtX[A, I] := Pa[I] + R * (Ct * U[I] + St * V[I] / L);
  Exit;
 End;
 For I := 1 To 3 Do Pc[I] := AtX[NC, I];
 { NeRF : C = NA, B = NB, A = NC }
 For I := 1 To 3 Do Begin Bc[I] := Pa[I] - Pb[I]; Ab[I] := Pb[I] - Pc[I]; End;
 L := Sqrt(Bc[1] * Bc[1] + Bc[2] * Bc[2] + Bc[3] * Bc[3]);
 For I := 1 To 3 Do Bc[I] := Bc[I] / L;
 Nv[1] := Ab[2] * Bc[3] - Ab[3] * Bc[2];
 Nv[2] := Ab[3] * Bc[1] - Ab[1] * Bc[3];
 Nv[3] := Ab[1] * Bc[2] - Ab[2] * Bc[1];
 L := Sqrt(Nv[1] * Nv[1] + Nv[2] * Nv[2] + Nv[3] * Nv[3]);
 If L < 1E-8 Then Begin
  { atomes alignes : axe arbitraire }
  Nv[1] := 0; Nv[2] := 1; Nv[3] := 0;
  If Abs(Bc[2]) > 0.9 Then Begin Nv[1] := 1; Nv[2] := 0; End;
  L := 1;
 End;
 For I := 1 To 3 Do Nv[I] := Nv[I] / L;
 Mv[1] := Nv[2] * Bc[3] - Nv[3] * Bc[2];
 Mv[2] := Nv[3] * Bc[1] - Nv[1] * Bc[3];
 Mv[3] := Nv[1] * Bc[2] - Nv[2] * Bc[1];
 Cp := Cos(Dih); Sp := Sin(Dih);
 D1 := -R * Ct; D2 := R * St * Cp; D3 := R * St * Sp;
 For I := 1 To 3 Do AtX[A, I] := Pa[I] + D1 * Bc[I] + D2 * Mv[I] + D3 * Nv[I];
End;

{ --------------------------------------------------------------------- }
{ Interpretation des lignes du fichier                                   }
{ --------------------------------------------------------------------- }

Procedure SetMethodDefaults;
Begin
 Meth := M_HF; Func := F_NONE; HFExch := 1.0;
End;

Procedure ApplyKeyword(Const Kw: String);
Var K: String; P: Integer; Bn: String;
Begin
 K := UpStr(Kw);
 If K = '' Then Exit;
 If (K = 'HF') Or (K = 'RHF') Or (K = 'RKS-HF') Then Begin Meth := M_HF; Func := F_NONE; HFExch := 1.0; End
 Else If (K = 'UHF') Or (K = 'ROHF') Or (K = 'UKS') Then Begin
  If K = 'ROHF' Then AddWarn('ROHF n''est pas disponible : UHF utilise.');
  UnrestrictedFlag := True;
  If K = 'UKS' Then Begin If Meth <> M_DFT Then Begin Meth := M_DFT; Func := F_B3LYP; HFExch := 0.2; End; End
  Else Begin Meth := M_HF; Func := F_NONE; HFExch := 1.0; End;
 End
 Else If K = 'RKS' Then Begin
  If Meth <> M_DFT Then Begin Meth := M_DFT; Func := F_B3LYP; HFExch := 0.2; End;
 End
 Else If (K = 'MP2') Or (K = 'RI-MP2') Or (K = 'RHF-MP2') Or (K = 'UHF-MP2') Or (K = 'DLPNO-MP2') Then Begin
  Meth := M_MP2; Func := F_NONE; HFExch := 1.0;
  If (K = 'RI-MP2') Or (K = 'DLPNO-MP2') Then AddWarn(K + ' : MP2 canonique utilise (sans RI).');
 End
 Else If (K = 'LDA') Or (K = 'SVWN') Or (K = 'SVWN5') Or (K = 'LSD') Then Begin Meth := M_DFT; Func := F_LDA; HFExch := 0.0; End
 Else If K = 'BLYP' Then Begin Meth := M_DFT; Func := F_BLYP; HFExch := 0.0; End
 Else If (K = 'B3LYP') Or (K = 'B3LYP/G') Or (K = 'B3LYP5') Then Begin Meth := M_DFT; Func := F_B3LYP; HFExch := 0.2; End
 Else If K = 'PBE' Then Begin Meth := M_DFT; Func := F_PBE; HFExch := 0.0; End
 Else If K = 'PBE0' Then Begin Meth := M_DFT; Func := F_PBE0; HFExch := 0.25; End
 Else If (K = 'STO-3G') Or (K = 'STO3G') Then BasisSel := B_STO3G
 Else If (K = '3-21G') Or (K = '321G') Then BasisSel := B_321G
 Else If (K = '6-31G') Or (K = '631G') Then BasisSel := B_631G
 Else If (K = '6-31G*') Or (K = '6-31G**') Or (K = '6-31+G*') Or (K = '6-31G(D)') Or (K = '6-31G(D,P)') Or
         (K = '6-311G') Or (K = '6-311G*') Or (K = '6-311G**') Or (K = 'DEF2-SVP') Or (K = 'DEF2-SV(P)') Or
         (K = 'DEF2-TZVP') Or (K = 'DEF2-SVP') Or (K = 'SVP') Or (K = 'DEF2-TZVPP') Or (K = 'CC-PVDZ') Or
         (K = 'AUG-CC-PVDZ') Or (K = 'DEF2-QZVP') Or (K = 'CC-PVTZ') Or (K = 'MINIX') Or (K = 'SV') Then Begin
  BasisSel := B_631G;
  AddWarn('Base ' + Kw + ' non disponible : 6-31G utilisee.');
 End
 Else If K = 'SP' Then JobType := J_SP
 Else If (K = 'OPT') Or (K = 'NORMALOPT') Then Begin
  If JobType = J_FREQ Then JobType := J_OPTFREQ Else JobType := J_OPT;
  SetOptTol('NormalOpt');
 End
 Else If (K = 'TIGHTOPT') Or (K = 'VERYTIGHTOPT') Or (K = 'LOOSEOPT') Then Begin
  If JobType = J_FREQ Then JobType := J_OPTFREQ Else JobType := J_OPT;
  SetOptTol(K);
 End
 Else If (K = 'FREQ') Or (K = 'NUMFREQ') Or (K = 'ANFREQ') Then Begin
  If JobType = J_OPT Then JobType := J_OPTFREQ Else JobType := J_FREQ;
 End
 Else If K = 'OPTFREQ' Then JobType := J_OPTFREQ
 Else If K = 'TIGHTSCF' Then Begin ScfTolE := 1E-8; ScfTolD := 5E-9; ScfTolErr := 5E-7; ScfTightName := 'TightSCF'; End
 Else If K = 'VERYTIGHTSCF' Then Begin ScfTolE := 1E-9; ScfTolD := 1E-9; ScfTolErr := 1E-7; ScfTightName := 'VeryTightSCF'; End
 Else If (K = 'LOOSESCF') Then Begin ScfTolE := 1E-5; ScfTolD := 1E-4; ScfTolErr := 1E-3; ScfTightName := 'LooseSCF'; End
 Else If (K = 'NORMALSCF') Or (K = 'SLOWCONV') Or (K = 'VERYSLOWCONV') Then Begin
  ScfTolE := 1E-6; ScfTolD := 1E-6; ScfTolErr := 1E-5; ScfTightName := 'NormalSCF';
 End
 Else If K = 'PRINTBASIS' Then DoPrintBasis := True
 Else If K = 'PRINTMOS' Then DoPrintMOs := True
 Else If K = 'NOPRINTMOS' Then DoPrintMOs := False
 Else If K = 'NOFROZENCORE' Then FrozenCore := False
 Else If (K = 'FROZENCORE') Or (K = 'FC') Then FrozenCore := True
 Else If K = 'CIS' Then DoCIS := True
 Else If (K = 'GRID1') Or (K = 'GRID2') Or (K = 'DEFGRID1') Or (K = 'DEFGRID2') Then GridScale := 1
 Else If (K = 'GRID3') Or (K = 'GRID4') Or (K = 'GRID5') Or (K = 'GRID6') Or (K = 'DEFGRID3') Then GridScale := 2
 Else If (K = 'RI') Or (K = 'NORI') Or (K = 'RIJCOSX') Or (K = 'RIJK') Or (K = 'RIJDX') Or (K = 'AUTOAUX') Or
         (K = 'NOAUTOSTART') Or (K = 'MINIPRINT') Or (K = 'NORMALPRINT') Or (K = 'LARGEPRINT') Or
         (K = 'ENGRAD') Or (K = 'NUMGRAD') Or (K = 'MOREAD') Or (K = 'NOITER') Or (K = 'NOUSEQUALS') Or
         (K = 'CPCM') Or (K = 'XYZFILE') Or (K = 'PDBFILE') Or (K = 'KDIIS') Or (K = 'SOSCF') Or
         (K = 'UNO') Or (K = 'NBO') Or (K = 'BOHRS') Or (K = 'ANGS') Or (K = 'PRINTMOS') Then
  { options acceptees mais sans effet }
 Else If (Copy(K, 1, 3) = 'PAL') Or (Copy(K, 1, 2) = 'D3') Or (Copy(K, 1, 4) = 'DEF2') Or
         (Copy(K, 1, 3) = 'AUG') Or (Copy(K, 1, 3) = 'CC-') Then
  AddWarn('Mot-cle ' + Kw + ' ignore.')
 Else AddWarn('Mot-cle inconnu ignore : ' + Kw);
End;

{ lit un bloc %nom ... end ; renvoie les jetons dans la liste }
{$IFDEF FPC}
Const MAX_TOK = 200;
{$ELSE}
Const MAX_TOK = 90;
{$ENDIF}
Var
 Tok     : Array[1..MAX_TOK] Of String[30];
 NTok    : Integer;

Procedure TokenizeLine(Const Ln: String);
Var S, T: String; I: Integer;
Begin
 S := Ln;
 { isole les accolades et le signe = }
 I := 1;
 While I <= Length(S) Do Begin
  If (S[I] = '{') Or (S[I] = '}') Or (S[I] = '=') Then Begin
   Insert(' ', S, I + 1);
   Insert(' ', S, I);
   Inc(I, 2);
  End;
  Inc(I);
 End;
 Repeat
  T := NextTok(S);
  If T <> '' Then Begin
   If NTok < MAX_TOK Then Begin Inc(NTok); Tok[NTok] := Copy(T, 1, 30); End;
  End;
 Until T = '';
End;

Function TokNum(I: Integer; Var V: Real): Boolean;
Begin
 If (I < 1) Or (I > NTok) Then Begin TokNum := False; Exit; End;
 TokNum := ParseReal(Tok[I], V);
End;

Function TokU(I: Integer): String;
Begin
 If (I < 1) Or (I > NTok) Then TokU := '' Else TokU := UpStr(Tok[I]);
End;

Procedure ParseGeomBlock;
Var
 I, J, K, NA: Integer; V, V1, V2: Real; U: String; KindC: Integer; C: TConstr;
 Function KindOf(Const S: String): Integer;
 Begin
  If S = 'B' Then KindOf := 1 Else If S = 'A' Then KindOf := 2
  Else If S = 'D' Then KindOf := 3 Else KindOf := 0;
 End;
Begin
 I := 1;
 While I <= NTok Do Begin
  U := TokU(I);
  If U = 'MAXITER' Then Begin
   If TokNum(I + 1, V) Then OptMaxCyc := Round(V);
   Inc(I, 2);
  End
  Else If (U = 'CONVERGENCE') Or (U = 'TOL') Then Begin
   SetOptTol(TokU(I + 1) + '');
   If UpStr(TokU(I + 1)) = 'TIGHT' Then SetOptTol('TightOpt');
   If UpStr(TokU(I + 1)) = 'LOOSE' Then SetOptTol('LooseOpt');
   Inc(I, 2);
  End
  Else If U = 'SCAN' Then Begin
   KindC := KindOf(TokU(I + 1));
   If KindC = 0 Then Begin AddWarn('%geom Scan : type de coordonnee non gere'); Inc(I); Continue; End;
   NA := ConstrNAt(KindC);
   ScanKind := KindC;
   For J := 1 To NA Do Begin
    If TokNum(I + 1 + J, V) Then ScanAt[J] := Round(V) + 1;
   End;
   K := I + 2 + NA;
   If TokU(K) = '=' Then Begin
    If TokNum(K + 1, V1) And TokNum(K + 2, V2) And TokNum(K + 3, V) Then Begin
     ScanLo := V1; ScanHi := V2; ScanSteps := Round(V); ScanOn := True;
    End;
   End;
   While (I <= NTok) And (TokU(I) <> 'END') Do Inc(I);
   Inc(I);
  End
  Else If U = 'CONSTRAINTS' Then Begin
   Inc(I);
   While (I <= NTok) And (TokU(I) <> 'END') Do Begin
    If TokU(I) = '{' Then Begin
     KindC := KindOf(TokU(I + 1));
     If (KindC > 0) And (NConstr < MAX_CONSTR) Then Begin
      NA := ConstrNAt(KindC);
      C.Kind := KindC; C.Fix := True; C.Val := 0;
      For J := 1 To NA Do If TokNum(I + 1 + J, V) Then C.A[J] := Round(V) + 1;
      J := I + 2 + NA;
      If TokNum(J, V) Then Begin C.Fix := False; C.Val := ToAngle(KindC, V); End;
      Inc(NConstr); Constr[NConstr] := C;
     End;
    End;
    Inc(I);
   End;
   Inc(I);
  End
  Else Inc(I);
 End;
End;

Procedure ParseBlock(Const Name: String);
Var N: String; I: Integer; V: Real; U: String;
Begin
 N := UpStr(Name);
 If N = 'SCF' Then Begin
  I := 1;
  While I <= NTok Do Begin
   U := TokU(I);
   If U = 'MAXITER' Then Begin If TokNum(I + 1, V) Then ScfMaxIter := Round(V); Inc(I, 2); End
   Else If U = 'CONVERGENCE' Then Begin
    If TokU(I + 1) = 'TIGHT' Then ApplyKeyword('TightSCF')
    Else If TokU(I + 1) = 'VERYTIGHT' Then ApplyKeyword('VeryTightSCF')
    Else If TokU(I + 1) = 'LOOSE' Then ApplyKeyword('LooseSCF')
    Else If TokU(I + 1) = 'NORMAL' Then ApplyKeyword('NormalSCF');
    Inc(I, 2);
   End
   Else If U = 'TOLE' Then Begin If TokNum(I + 1, V) Then ScfTolE := V; Inc(I, 2); End
   Else If U = 'TOLRMSP' Then Begin If TokNum(I + 1, V) Then ScfTolD := V; Inc(I, 2); End
   Else If U = 'TOLERR' Then Begin If TokNum(I + 1, V) Then ScfTolErr := V; Inc(I, 2); End
   Else If U = 'HFTYP' Then Begin
    If (TokU(I + 1) = 'UHF') Or (TokU(I + 1) = 'UKS') Then UnrestrictedFlag := True;
    Inc(I, 2);
   End
   Else Inc(I);
  End;
 End
 Else If N = 'GEOM' Then ParseGeomBlock
 Else If (N = 'TDDFT') Or (N = 'CIS') Then Begin
  DoCIS := True;
  I := 1;
  While I <= NTok Do Begin
   U := TokU(I);
   If U = 'NROOTS' Then Begin If TokNum(I + 1, V) Then NRoots := Round(V); Inc(I, 2); End
   Else If U = 'TRIPLETS' Then Begin CISTriplet := (TokU(I + 1) = 'TRUE'); Inc(I, 2); End
   Else Inc(I);
  End;
 End
 Else If N = 'FREQ' Then Begin
  I := 1;
  While I <= NTok Do Begin
   U := TokU(I);
   If U = 'TEMP' Then Begin If TokNum(I + 1, V) Then CalcTemp := V; Inc(I, 2); End
   Else Inc(I);
  End;
 End
 Else If N = 'GRID' Then Begin
  I := 1;
  While I <= NTok Do Begin
   If TokU(I) = 'SCALE' Then Begin If TokNum(I + 1, V) Then GridScale := Round(V); Inc(I, 2); End
   Else Inc(I);
  End;
 End
 Else If (N = 'MAXCORE') Or (N = 'PAL') Or (N = 'OUTPUT') Or (N = 'METHOD') Or (N = 'BASIS') Or
         (N = 'COORDS') Or (N = 'ELPROP') Or (N = 'MDCI') Or (N = 'MP2') Then
  { ignores sans avertissement }
 Else AddWarn('Bloc %' + Name + ' ignore.');
End;

{ Interprete les lignes InpLines ; remplit parametres et geometrie }
Function ParseInput(Var Err: String): Boolean;
Var
 I, L, N, Z, Depth, GeomStart, Code: Integer;
 S, T, U, Bn, Rest: String;
 Mode: Integer;            { 0 normal, 1 bloc, 2 coordonnees xyz, 3 coordonnees int }
 BlkName: String;
 GotGeom: Boolean;
 V: Real;
 ZRec: Array[1..MAX_AT] Of Record NA, NB, NC: Integer; R, Ang, Dih: Real; End;
 Nested: Boolean;
 Function IsNested(Const Tk: String): Boolean;
 Begin
  IsNested := (Tk = 'SCAN') Or (Tk = 'CONSTRAINTS') Or (Tk = 'MODIFYINTERNAL') Or (Tk = 'INVOKE');
 End;
Begin
 ParseInput := False;
 Err := '';
 NWarn := 0;
 NAt := 0; GotGeom := False;
 Charge := 0; Mult := 1;
 DefaultParams;
 DoCIS := False; NRoots := 5; CISTriplet := False;
 DoPrintBasis := False; DoPrintMOs := True;
 NConstr := 0; ScanOn := False;
 OptMaxCyc := 60; SetOptTol('NormalOpt');
 CalcTemp := 298.15; GridScale := 1;
 ZMat := False;
 Mode := 0; I := 0;
 While I < NInp Do Begin
  S := StripComment(InpLines[I]);
  S := TrimStr(S);
  Inc(I);
  If S = '' Then Continue;
  If S[1] = '!' Then Begin
   Delete(S, 1, 1);
   Repeat
    T := NextTok(S);
    If T <> '' Then ApplyKeyword(T);
   Until T = '';
   Continue;
  End;
  If S[1] = '%' Then Begin
   Delete(S, 1, 1);
   BlkName := NextTok(S);
   NTok := 0;
   Depth := 1;
   Rest := S;
   If UpStr(BlkName) = 'MAXCORE' Then Continue;
   { bloc sur une seule ligne ou sur plusieurs lignes }
   Repeat
    TokenizeLine(Rest);
    { recherche de la fin du bloc : compte des imbrications }
    Depth := 1;
    For N := 1 To NTok Do Begin
     If IsNested(UpStr(Tok[N])) Then Inc(Depth)
     Else If UpStr(Tok[N]) = 'END' Then Dec(Depth);
    End;
    If Depth > 0 Then Begin
     If I >= NInp Then Begin Err := 'Bloc %' + BlkName + ' non termine (end manquant)'; Exit; End;
     Rest := TrimStr(StripComment(InpLines[I]));
     Inc(I);
    End;
   Until Depth <= 0;
   { le dernier END ferme le bloc : on le retire }
   If (NTok > 0) And (UpStr(Tok[NTok]) = 'END') Then Dec(NTok);
   ParseBlock(BlkName);
   Continue;
  End;
  If S[1] = '*' Then Begin
   Delete(S, 1, 1);
   S := TrimStr(S);
   T := UpStr(NextTok(S));
   If (T = 'XYZFILE') Then Begin
    Bn := NextTok(S);
    If Not ParseInt(Bn, Z) Then Begin Err := 'Charge invalide'; Exit; End;
    Charge := Z;
    Bn := NextTok(S);
    If Not ParseInt(Bn, Z) Then Begin Err := 'Multiplicite invalide'; Exit; End;
    Mult := Z;
    Bn := TrimStr(S);
    If Not ReadXYZFile(Bn, Err) Then Exit;
    GotGeom := True;
    Continue;
   End;
   If (T <> 'XYZ') And (T <> 'INT') And (T <> 'GZMT') Then Begin
    Err := 'Format de coordonnees inconnu : ' + T; Exit;
   End;
   Bn := NextTok(S);
   If Not ParseInt(Bn, Z) Then Begin Err := 'Charge invalide dans la ligne *'; Exit; End;
   Charge := Z;
   Bn := NextTok(S);
   If Not ParseInt(Bn, Z) Then Begin Err := 'Multiplicite invalide dans la ligne *'; Exit; End;
   Mult := Z;
   NAt := 0;
   ZMat := (T = 'INT') Or (T = 'GZMT');
   { lecture jusqu'a la ligne * ou END }
   While I < NInp Do Begin
    S := TrimStr(StripComment(InpLines[I]));
    Inc(I);
    If S = '' Then Continue;
    If (S[1] = '*') Or (UpStr(S) = 'END') Then Break;
    If NAt >= MAX_AT Then Begin Err := 'Trop d''atomes (maximum ' + IStr(MAX_AT) + ')'; Exit; End;
    T := NextTok(S);
    Z := ElemZ(T);
    If Z = 0 Then Begin
     Val(T, Z, Code);
     If (Code <> 0) Or (Z < 1) Or (Z > NELEM) Then Begin Err := 'Element inconnu : ' + T; Exit; End;
    End;
    Inc(NAt);
    AtZ[NAt] := Z;
    If Not ZMat Then Begin
     For N := 1 To 3 Do Begin
      T := NextTok(S);
      If Not ParseReal(T, V) Then Begin Err := 'Coordonnee invalide : ' + T; Exit; End;
      AtX[NAt, N] := V * ANG2BOHR;
     End;
    End
    Else Begin
     For N := 1 To 3 Do Begin
      T := NextTok(S);
      If Not ParseInt(T, Z) Then Begin Err := 'Reference de matrice Z invalide : ' + T; Exit; End;
      Case N Of 1: ZRec[NAt].NA := Z; 2: ZRec[NAt].NB := Z; Else ZRec[NAt].NC := Z; End;
     End;
     T := NextTok(S); If Not ParseReal(T, V) Then V := 0; ZRec[NAt].R := V * ANG2BOHR;
     T := NextTok(S); If Not ParseReal(T, V) Then V := 0; ZRec[NAt].Ang := V * Pi / 180;
     T := NextTok(S); If Not ParseReal(T, V) Then V := 0; ZRec[NAt].Dih := V * Pi / 180;
    End;
   End;
   If ZMat Then
    For N := 1 To NAt Do
     ZPlace(N, ZRec[N].NA, ZRec[N].NB, ZRec[N].NC, ZRec[N].R, ZRec[N].Ang, ZRec[N].Dih);
   GotGeom := True;
   Continue;
  End;
  AddWarn('Ligne ignoree : ' + Copy(S, 1, 50));
 End;
 If Not GotGeom Or (NAt = 0) Then Begin Err := 'Aucune geometrie dans le fichier (ligne "* xyz charge mult")'; Exit; End;
 { complete les contraintes figees avec leur valeur courante, le balayage }
 For N := 1 To NConstr Do
  If Constr[N].Fix Then Constr[N].Val := ConstrValue(Constr[N]);
 If ScanOn Then Begin
  If ScanSteps < 1 Then ScanSteps := 1;
  If ScanSteps > MAX_SCAN - 1 Then ScanSteps := MAX_SCAN - 1;
 End;
 { fonctionnelle : ajustement UKS/UHF }
 ParseInput := True;
End;

{ --------------------------------------------------------------------- }
{ Impression de la sortie (fichier .OUT)                                 }
{ --------------------------------------------------------------------- }

Function EhEv(E: Real): String;
Begin
 EhEv := FmtR(E * HARTREE_EV, 18, 5);
End;

Function MethodName: String;
Begin
 Case Meth Of
  M_MP2: If NOccA = NOccB Then MethodName := 'RHF-MP2' Else MethodName := 'UHF-MP2';
  M_DFT: Case Func Of
          F_LDA: MethodName := 'LDA (Slater + VWN5)';
          F_BLYP: MethodName := 'BLYP';
          F_B3LYP: MethodName := 'B3LYP';
          F_PBE: MethodName := 'PBE';
         Else MethodName := 'PBE0';
         End;
 Else If NOccA = NOccB Then MethodName := 'Hartree-Fock (RHF)' Else MethodName := 'Hartree-Fock (UHF)';
 End;
End;

Procedure PrintBanner;
Begin
 Emit('');
 Emit('                                 *****************');
 Emit('                                 * O   R   C   A *');
 Emit('                                 *****************');
 Emit('');
 Emit('                  --- Clone Pascal de ORCA - QChimie ' + PROG_VERSION + ' ---');
 Emit('                    Programme de chimie quantique ab initio');
 Emit('                  Moteur : Hartree-Fock, MP2, CIS et DFT (Becke)');
 Emit('');
 Emit('--------------------------------------------------------------------------------');
 Emit('                                  INPUT FILE');
 Emit('================================================================================');
End;

Procedure PrintInputEcho;
Var I: Integer;
Begin
 If InpName = '' Then Emit('NAME = ORCA.INP') Else Emit('NAME = ' + InpName);
 For I := 0 To NInp - 1 Do Emit('|' + FmtR(I + 1, 3, 0) + '> ' + InpLines[I]);
 Emit('                                 ****END OF INPUT****');
 Emit('================================================================================');
 If NWarn > 0 Then Begin
  EmitBlank;
  For I := 1 To NWarn Do Emit('WARNING: ' + WarnList[I]);
 End;
End;

Procedure PrintBoxed(Const T: String);
Var L, Pad: String; I: Integer;
Begin
 L := '';
 For I := 1 To Length(T) + 4 Do L := L + '*';
 Emit('                        ' + L);
 Emit('                        * ' + T + ' *');
 Emit('                        ' + L);
End;

Procedure PrintGeometry;
Var A, I, Na, Nb, Nc: Integer; Tot, Bd, Rv, Av, Dv: Real; C: TConstr;
Begin
 EmitBlank;
 Emit('---------------------------------');
 Emit('CARTESIAN COORDINATES (ANGSTROEM)');
 Emit('---------------------------------');
 For A := 1 To NAt Do
  Emit('  ' + PadR(ElemStr(AtZ[A]), 3) + FmtR(AtX[A, 1] / ANG2BOHR, 14, 6) +
       FmtR(AtX[A, 2] / ANG2BOHR, 12, 6) + FmtR(AtX[A, 3] / ANG2BOHR, 12, 6));
 EmitBlank;
 Emit('----------------------------');
 Emit('CARTESIAN COORDINATES (A.U.)');
 Emit('----------------------------');
 Emit('  NO LB      ZA    FRAG     MASS         X           Y           Z');
 For A := 1 To NAt Do
  Emit(FmtR(A - 1, 4, 0) + ' ' + PadR(ElemStr(AtZ[A]), 3) + FmtR(AtZ[A], 8, 1) + '    0' +
       FmtR(ElMass[AtZ[A]], 11, 3) + FmtR(AtX[A, 1], 12, 5) + FmtR(AtX[A, 2], 12, 5) +
       FmtR(AtX[A, 3], 12, 5));
 EmitBlank;
 Emit('--------------------------------');
 Emit('INTERNAL COORDINATES (ANGSTROEM)');
 Emit('--------------------------------');
 For A := 1 To NAt Do Begin
  { atomes de reference : les plus proches parmi les precedents }
  Na := 0; Nb := 0; Nc := 0;
  Bd := 1E30;
  For I := 1 To A - 1 Do If Dist(A, I) < Bd Then Begin Bd := Dist(A, I); Na := I; End;
  If Na > 0 Then Begin
   Bd := 1E30;
   For I := 1 To A - 1 Do If (I <> Na) And (Dist(Na, I) < Bd) Then Begin Bd := Dist(Na, I); Nb := I; End;
  End;
  If Nb > 0 Then Begin
   Bd := 1E30;
   For I := 1 To A - 1 Do If (I <> Na) And (I <> Nb) And (Dist(Nb, I) < Bd) Then Begin Bd := Dist(Nb, I); Nc := I; End;
  End;
  Rv := 0; Av := 0; Dv := 0;
  If Na > 0 Then Rv := Dist(A, Na) / ANG2BOHR;
  If Nb > 0 Then Begin
   C.Kind := 2; C.A[1] := A; C.A[2] := Na; C.A[3] := Nb;
   Av := ConstrValue(C) * 180 / Pi;
  End;
  If Nc > 0 Then Begin
   C.Kind := 3; C.A[1] := A; C.A[2] := Na; C.A[3] := Nb; C.A[4] := Nc;
   Dv := ConstrValue(C) * 180 / Pi;
  End;
  If Na = 0 Then Na := 1;
  If Nb = 0 Then Nb := 1;
  If Nc = 0 Then Nc := 1;
  If A = 1 Then Begin Na := 1; Nb := 1; Nc := 1; End
  Else If A = 2 Then Begin Nb := 1; Nc := 1; End
  Else If A = 3 Then Nc := 1;
  Emit('  ' + PadR(ElemStr(AtZ[A]), 3) + FmtR(Na - 1 + Ord(A = 1) * 0, 3, 0) + FmtR(Nb - 1, 4, 0) + FmtR(Nc - 1, 4, 0) +
       FmtR(Rv, 14, 6) + FmtR(Av, 12, 4) + FmtR(Dv, 12, 4));
 End;
End;

Procedure PrintBasisInfo;
Var A, S, K, NS, NP, PS, PP: Integer; Pat1, Pat2: String; Seen: Array[1..NELEM] Of Boolean;
Begin
 EmitBlank;
 Emit('-----------------------');
 Emit('BASIS SET INFORMATION');
 Emit('-----------------------');
 Emit('Basis set                                  ....  ' + BasisNames[BasisSel] + ' (gaussiennes cartesiennes)');
 FillChar(Seen, SizeOf(Seen), 0);
 For A := 1 To NAt Do
  If Not Seen[AtZ[A]] Then Begin
   Seen[AtZ[A]] := True;
   If GetShells(BasisSel, AtZ[A]) Then Begin
    NS := 0; NP := 0; PS := 0; PP := 0; Pat1 := ''; Pat2 := '';
    For S := 1 To NTS Do
     If TSL[S] = 0 Then Begin Inc(NS); PS := PS + TSN[S]; Pat1 := Pat1 + IStr(TSN[S]); End
     Else Begin Inc(NP); PP := PP + TSN[S]; Pat2 := Pat2 + IStr(TSN[S]); End;
    Emit('Type ' + PadR(ElemStr(AtZ[A]), 3) + ': ' + IStr(PS) + 's' +
         IStr(PP) + 'p contracted to ' + IStr(NS) + 's' + IStr(NP) + 'p pattern {' + Pat1 + '/' + Pat2 + '}');
   End;
  End;
 EmitBlank;
 Emit('Number of basis functions                   ....  ' + IStr(NBF));
 If DoPrintBasis Then Begin
  EmitBlank;
  Emit('Fonctions de base :');
  For K := 1 To NBF Do
   Emit('  ' + FmtR(K - 1, 3, 0) + '  ' + PadR(IStr(BFAtom[K] - 1) + ElemStr(AtZ[BFAtom[K]]), 5) + BFLab[K] +
        '  (' + IStr(BFNP[K]) + ' primitives)');
 End;
End;

Procedure PrintScfSettings;
Begin
 EmitBlank;
 Emit('------------');
 Emit('SCF SETTINGS');
 Emit('------------');
 Emit('Hamiltonian:');
 If Meth = M_DFT Then Begin
  Emit(' Density Functional     Method          .... DFT(GTOs)');
  Emit(' Functional              XC              .... ' + MethodName);
  Emit(' Exchange HF fraction    ScalHFX         .... ' + RStr(HFExch, 6));
 End
 Else Emit(' Ab initio Hamiltonian  Method          .... Hartree-Fock(GTOs)');
 EmitBlank;
 Emit('General Settings:');
 Emit(' Total Charge           Charge          .... ' + FmtR(Charge, 6, 0));
 Emit(' Multiplicity           Mult            .... ' + FmtR(Mult, 6, 0));
 Emit(' Number of Electrons    NEL             .... ' + FmtR(NElec, 6, 0));
 Emit(' Basis Dimension        Dim             .... ' + FmtR(NBF, 6, 0));
 Emit(' Nuclear Repulsion      ENuc            .... ' + FmtR(NucRepulsion, 18, 10) + ' Eh');
 EmitBlank;
 Emit('SCF convergence criteria (' + ScfTightName + '):');
 Emit(' Energy Change           TolE            .... ' + SciStr(ScfTolE, 2) + ' Eh');
 Emit(' 1-El. energy change     TolRMSP         .... ' + SciStr(ScfTolD, 2));
 Emit(' Orbital Gradient        TolG            .... ' + SciStr(ScfTolErr, 2));
 Emit(' Max. Iterations         MaxIter         .... ' + FmtR(ScfMaxIter, 6, 0));
 If Meth = M_DFT Then Begin
  EmitBlank;
  Emit(' Grille d''integration : Becke (radiale Gauss-Chebyshev, angulaire ' + IStr(GRID_NTH * GridScale) + 'x' +
       IStr(2 * GRID_NTH * GridScale) + ')');
 End;
End;

Procedure PrintOrbitalEnergies;
Var K: Integer; Occ: Real; Homo, Lumo: Real;
Begin
 EmitBlank;
 Emit('----------------');
 Emit('ORBITAL ENERGIES');
 Emit('----------------');
 If NOccA = NOccB Then Begin
  Emit('  NO   OCC          E(Eh)            E(eV)');
  For K := 1 To NBF Do Begin
   If K <= NOccA Then Occ := 2 Else Occ := 0;
   Emit(FmtR(K - 1, 4, 0) + FmtR(Occ, 8, 4) + FmtR(EpsA[K], 16, 6) + FmtR(EpsA[K] * HARTREE_EV, 16, 4));
  End;
  If (NOccA > 0) And (NOccA < NBF) Then Begin
   EmitBlank;
   Emit('HOMO-LUMO gap : ' + RStr((EpsA[NOccA + 1] - EpsA[NOccA]) * HARTREE_EV, 4) + ' eV');
  End;
 End
 Else Begin
  Emit('                      SPIN UP ORBITALS');
  Emit('  NO   OCC          E(Eh)            E(eV)');
  For K := 1 To NBF Do Begin
   If K <= NOccA Then Occ := 1 Else Occ := 0;
   Emit(FmtR(K - 1, 4, 0) + FmtR(Occ, 8, 4) + FmtR(EpsA[K], 16, 6) + FmtR(EpsA[K] * HARTREE_EV, 16, 4));
  End;
  EmitBlank;
  Emit('                      SPIN DOWN ORBITALS');
  Emit('  NO   OCC          E(Eh)            E(eV)');
  For K := 1 To NBF Do Begin
   If K <= NOccB Then Occ := 1 Else Occ := 0;
   Emit(FmtR(K - 1, 4, 0) + FmtR(Occ, 8, 4) + FmtR(EpsB[K], 16, 6) + FmtR(EpsB[K] * HARTREE_EV, 16, 4));
  End;
 End;
End;

Procedure PrintMOCoefficients;
Var K, J, I0, Cnt: Integer; L1, L2, L3, L4: String;
Begin
 EmitBlank;
 Emit('-----------------');
 Emit('MOLECULAR ORBITALS');
 Emit('-----------------');
 I0 := 1;
 While I0 <= NBF Do Begin
  Cnt := NBF - I0 + 1;
  If Cnt > 6 Then Cnt := 6;
  L1 := '                  '; L2 := '          '; L3 := '                  '; L4 := '';
  For J := 0 To Cnt - 1 Do Begin
   L1 := L1 + FmtR(EpsA[I0 + J], 10, 5);
   If I0 + J <= NOccA Then L2 := L2 + FmtR(2 - Ord(NOccA <> NOccB), 10, 5)
   Else L2 := L2 + FmtR(0, 10, 5);
   L3 := L3 + FmtR(I0 + J - 1, 10, 0);
  End;
  Emit(L3);
  Emit(L1);
  Emit('                 ' + '--------------------------------------------------------------');
  For K := 1 To NBF Do Begin
   L4 := FmtR(K - 1, 3, 0) + ' ' + PadR(IStr(BFAtom[K] - 1) + ElemStr(AtZ[BFAtom[K]]), 4) + PadR(BFLab[K], 6) + '  ';
   For J := 0 To Cnt - 1 Do L4 := L4 + FmtR(MCa^[K, I0 + J], 10, 5);
   Emit(L4);
  End;
  EmitBlank;
  Inc(I0, 6);
 End;
End;

Procedure PrintSCFResults;
Var Tk, Vir, Sum: Real; A, B: Integer; Line: String;
Begin
 EmitBlank;
 Emit('----------------');
 Emit('TOTAL SCF ENERGY');
 Emit('----------------');
 EmitBlank;
 Emit('Total Energy       :' + FmtR(EScf, 20, 12) + ' Eh  ' + EhEv(EScf) + ' eV');
 EmitBlank;
 Emit('Components:');
 Emit('Nuclear Repulsion  :' + FmtR(Enuc, 20, 12) + ' Eh  ' + EhEv(Enuc) + ' eV');
 Emit('Electronic Energy  :' + FmtR(EScf - Enuc, 20, 12) + ' Eh  ' + EhEv(EScf - Enuc) + ' eV');
 Emit('One Electron Energy:' + FmtR(EOne, 20, 12) + ' Eh  ' + EhEv(EOne) + ' eV');
 Emit('Two Electron Energy:' + FmtR(ECoul + EExch + EXCdft, 20, 12) + ' Eh  ' + EhEv(ECoul + EExch + EXCdft) + ' eV');
 If Meth = M_DFT Then Begin
  Emit('  Coulomb          :' + FmtR(ECoul, 20, 12) + ' Eh');
  Emit('  Exact exchange   :' + FmtR(EExch, 20, 12) + ' Eh');
  Emit('  E(XC)            :' + FmtR(EXCdft, 20, 12) + ' Eh');
  Emit('  E(X) / E(C)      :' + FmtR(EXdft, 20, 12) + ' Eh' + FmtR(ECdft, 20, 12) + ' Eh');
  Emit('N(Alpha)           :' + FmtR(DftNAlpha, 20, 10) + ' electrons');
  Emit('N(Beta)            :' + FmtR(DftNBeta, 20, 10) + ' electrons');
  Emit('N(Total)           :' + FmtR(DftNAlpha + DftNBeta, 20, 10) + ' electrons');
 End;
 Tk := 0;
 For A := 1 To NBF Do
  For B := 1 To NBF Do Tk := Tk + (MDa^[A, B] + MDb^[A, B]) * MT^[A, B];
 EmitBlank;
 Emit('Virial components:');
 Emit('Potential Energy   :' + FmtR(EScf - Tk, 20, 12) + ' Eh');
 Emit('Kinetic Energy     :' + FmtR(Tk, 20, 12) + ' Eh');
 If Tk <> 0 Then Emit('Virial Ratio       :' + FmtR(-(EScf - Tk) / Tk, 20, 12));
 If NOccA <> NOccB Then Begin
  EmitBlank;
  Emit('Expectation value of <S**2>     : ' + RStr(S2Exp, 6));
  Emit('Ideal value S*(S+1) for S=' + RStr((Mult - 1) / 2, 1) + '   : ' + RStr((Mult - 1) / 2 * ((Mult - 1) / 2 + 1), 6));
  Emit('Deviation                       : ' + RStr(S2Exp - (Mult - 1) / 2 * ((Mult - 1) / 2 + 1), 6));
 End;
 PrintOrbitalEnergies;
 If DoPrintMOs And (NBF <= 24) Then PrintMOCoefficients;
 { populations }
 EmitBlank;
 Emit('-----------------------');
 Emit('MULLIKEN ATOMIC CHARGES');
 If NOccA <> NOccB Then Begin Emit('   (et populations de spin)'); End;
 Emit('-----------------------');
 Sum := 0;
 For A := 1 To NAt Do Begin
  Line := FmtR(A - 1, 4, 0) + ' ' + PadR(ElemStr(AtZ[A]), 2) + ':' + FmtR(QMull[A], 12, 6);
  If NOccA <> NOccB Then Line := Line + FmtR(SpinMull[A], 12, 6);
  Emit(Line);
  Sum := Sum + QMull[A];
 End;
 Emit('Sum of atomic charges:' + FmtR(Sum, 12, 7));
 EmitBlank;
 Emit('----------------------');
 Emit('LOEWDIN ATOMIC CHARGES');
 Emit('----------------------');
 Sum := 0;
 For A := 1 To NAt Do Begin
  Emit(FmtR(A - 1, 4, 0) + ' ' + PadR(ElemStr(AtZ[A]), 2) + ':' + FmtR(QLoew[A], 12, 6));
  Sum := Sum + QLoew[A];
 End;
 Emit('Sum of atomic charges:' + FmtR(Sum, 12, 7));
 If NAt > 1 Then Begin
  EmitBlank;
  Emit('--------------------------');
  Emit('MAYER POPULATION ANALYSIS');
  Emit('--------------------------');
  Emit('  Mayer bond orders larger than 0.1');
  For A := 1 To NAt - 1 Do
   For B := A + 1 To NAt Do
    If BondOrd[A, B] > 0.1 Then
     Emit('B(' + FmtR(A - 1, 3, 0) + '-' + PadR(ElemStr(AtZ[A]), 2) + ' ,' + FmtR(B - 1, 3, 0) + '-' +
          PadR(ElemStr(AtZ[B]), 2) + ' ) :' + FmtR(BondOrd[A, B], 8, 4));
 End;
 EmitBlank;
 Emit('-------------');
 Emit('DIPOLE MOMENT');
 Emit('-------------');
 Emit('                                X             Y             Z');
 Emit('Total Dipole Moment    :' + FmtR(Dip[1], 14, 9) + FmtR(Dip[2], 14, 9) + FmtR(Dip[3], 14, 9));
 Emit('Magnitude (a.u.)       :' + FmtR(DipTot, 14, 9));
 Emit('Magnitude (Debye)      :' + FmtR(DipTot * BOHR_DEBYE, 14, 9));
End;

Procedure PrintMP2Results;
Var Scs: Real;
Begin
 EmitBlank;
 Emit('-----------------------------------------------');
 Emit('MP2 CORRELATION ENERGY' + PadR('', 0));
 Emit('-----------------------------------------------');
 If NFrozen > 0 Then Emit('Orbitales de coeur gelees : ' + IStr(NFrozen))
 Else Emit('Aucune orbitale de coeur gelee');
 Emit('Opposite-spin component          ' + FmtR(EMP2os, 18, 10) + ' Eh');
 Emit('Same-spin component              ' + FmtR(EMP2ss, 18, 10) + ' Eh');
 Scs := 1.2 * EMP2os + EMP2ss / 3;
 Emit('E(SCS-MP2) correlation           ' + FmtR(Scs, 18, 10) + ' Eh');
 Emit('E(MP2) correlation energy        ' + FmtR(EMP2, 18, 10) + ' Eh');
 EmitBlank;
 Emit('Initial E(tot)  (SCF)            ' + FmtR(EScf, 18, 10) + ' Eh');
 Emit('Final   E(MP2)                   ' + FmtR(EScf + EMP2, 18, 10) + ' Eh');
End;

Procedure PrintCISResults;
Var K: Integer; Ev, Cm, Nm: Real;
Begin
 EmitBlank;
 If CISTriplet Then Emit('CIS - ETATS TRIPLETS') Else Emit('CIS - ETATS SINGULETS');
 Emit('---------------------------------------------------------------------------');
 If Not CISOK Then Begin Emit(CISMsg); Exit; End;
 Emit('STATE   ENERGY(Eh)      (eV)       (cm-1)    WAVELENGTH(nm)    fosc   TRANSITIONS');
 For K := 1 To CISN Do Begin
  Ev := CISE[K] * HARTREE_EV;
  Cm := CISE[K] * HARTREE_CM;
  If CISE[K] > 1E-9 Then Nm := 1E7 / Cm Else Nm := 0;
  Emit(FmtR(K, 4, 0) + FmtR(CISE[K], 14, 6) + FmtR(Ev, 10, 4) + FmtR(Cm, 12, 1) + FmtR(Nm, 14, 1) +
       FmtR(CISF[K], 12, 6) + '   ' + CISDom[K]);
 End;
 EmitBlank;
 Emit('Energie totale de l''etat excite 1 : ' + RStr(EScf + CISE[1], 8) + ' Eh');
End;

Procedure PrintFreqResults;
Var K: Integer;
Begin
 EmitBlank;
 Emit('-----------------------');
 Emit('VIBRATIONAL FREQUENCIES');
 Emit('-----------------------');
 EmitBlank;
 Emit('Scaling factor for frequencies =  1.000000000  (already applied!)');
 EmitBlank;
 For K := 1 To NFreq Do Begin
  Emit(FmtR(K - 1, 6, 0) + ':' + FmtR(FreqV[K], 14, 2) + ' cm**-1' + PadR('', 0));
 End;
 EmitBlank;
 Emit('-----------');
 Emit('IR SPECTRUM');
 Emit('-----------');
 EmitBlank;
 Emit(' Mode   freq (cm**-1)   Int (km/mol)   T**2 (a.u.)      TX        TY        TZ');
 Emit('--------------------------------------------------------------------------------');
 For K := NProj + 1 To NFreq Do
  Emit(FmtR(K - 1, 4, 0) + ':' + FmtR(FreqV[K], 14, 2) + FmtR(FreqI[K], 14, 2) +
       FmtR(FreqI[K] / IR_CONV, 14, 6) + FmtR(IRTx[K], 10, 6) + FmtR(IRTy[K], 10, 6) + FmtR(IRTz[K], 10, 6));
 If FreqV[NProj + 1] < -20 Then Begin
  EmitBlank;
  Emit('***imaginary mode***  : la geometrie n''est pas un minimum de la surface d''energie');
 End;
End;

Procedure PrintThermoResults;
Var Rn: Real;
Begin
 EmitBlank;
 Emit('--------------------------');
 Emit('THERMOCHEMISTRY AT ' + RStr(ThTemp, 2) + 'K');
 Emit('--------------------------');
 EmitBlank;
 Emit('Temperature         ... ' + FmtR(ThTemp, 10, 2) + ' K');
 Emit('Pressure            ... ' + FmtR(1, 10, 2) + ' atm');
 Emit('Total Mass          ... ' + FmtR(ThMass, 10, 3) + ' AMU');
 Emit('Rotational symmetry number  ... ' + IStr(SymNum));
 Emit('Rotational constants in cm-1: ' + FmtR(RotConst[1], 12, 6) + FmtR(RotConst[2], 12, 6) + FmtR(RotConst[3], 12, 6));
 EmitBlank;
 Emit('-------------------');
 Emit('INNER ENERGY');
 Emit('-------------------');
 Emit('Electronic energy                ' + FmtR(ETot, 18, 8) + ' Eh');
 Emit('Zero point energy                ' + FmtR(ZPE, 18, 8) + ' Eh  ' + FmtR(ZPE * HARTREE_KCAL, 10, 2) + ' kcal/mol');
 Emit('Thermal vibrational correction   ' + FmtR(EVibT, 18, 8) + ' Eh');
 Emit('Thermal rotational correction    ' + FmtR(ERotT, 18, 8) + ' Eh');
 Emit('Thermal translational correction ' + FmtR(ETransT, 18, 8) + ' Eh');
 Emit('-----------------------------------------------------------------------');
 Emit('Total thermal energy             ' + FmtR(ETot + EThermal, 18, 8) + ' Eh');
 EmitBlank;
 Emit('-------');
 Emit('ENTHALPY');
 Emit('-------');
 Emit('Thermal Enthalpy correction      ' + FmtR(KB_HARTREE * ThTemp, 18, 8) + ' Eh');
 Emit('Total Enthalpy                   ' + FmtR(HTot, 18, 8) + ' Eh');
 EmitBlank;
 Emit('-------');
 Emit('ENTROPY');
 Emit('-------');
 Emit('Electronic entropy               ' + FmtR(KB_HARTREE * Ln(Mult) * ThTemp, 18, 8) + ' Eh');
 Emit('Vibrational entropy              ' + FmtR(SVibT * ThTemp, 18, 8) + ' Eh');
 Emit('Rotational entropy               ' + FmtR(SRotT * ThTemp, 18, 8) + ' Eh');
 Emit('Translational entropy            ' + FmtR(STransT * ThTemp, 18, 8) + ' Eh');
 Emit('-----------------------------------------------------------------------');
 Emit('Final entropy term               ' + FmtR(STot * ThTemp, 18, 8) + ' Eh');
 EmitBlank;
 Emit('-------');
 Emit('GIBBS FREE ENERGY');
 Emit('-------');
 Emit('Final Gibbs free energy          ' + FmtR(GTot, 18, 8) + ' Eh');
 Emit('G-E(el)                          ' + FmtR(GTot - ETot, 18, 8) + ' Eh  ' +
      FmtR((GTot - ETot) * HARTREE_KCAL, 10, 2) + ' kcal/mol');
End;

{ --------------------------------------------------------------------- }
{ Balayage relaxe de surface                                             }
{ --------------------------------------------------------------------- }

Function ScanValueToUser(Kind: Integer; V: Real): Real;
Begin
 If Kind = 1 Then ScanValueToUser := V / ANG2BOHR Else ScanValueToUser := V * 180 / Pi;
End;

Function RunScan(Var Err: String): Boolean;
Var
 K, I, SaveN: Integer; Tgt, V0, V1: Real; C: TConstr; OK: Boolean;
 Unit_: String;
Begin
 RunScan := False;
 ScanN := 0;
 SaveN := NConstr;
 If SaveN >= MAX_CONSTR Then Begin Err := 'Trop de contraintes'; Exit; End;
 C.Kind := ScanKind; C.Fix := False;
 For I := 1 To 4 Do C.A[I] := ScanAt[I];
 For I := 1 To ConstrNAt(ScanKind) Do
  If (C.A[I] < 1) Or (C.A[I] > NAt) Then Begin Err := 'Numero d''atome invalide dans %geom Scan'; Exit; End;
 If ScanKind = 1 Then Unit_ := 'Angstroem' Else Unit_ := 'degres';
 V0 := ToAngle(ScanKind, ScanLo); V1 := ToAngle(ScanKind, ScanHi);
 Inc(NConstr);
 For K := 0 To ScanSteps Do Begin
  If ScanSteps > 0 Then Tgt := V0 + (V1 - V0) * K / ScanSteps Else Tgt := V0;
  C.Val := Tgt;
  Constr[NConstr] := C;
  EmitBlank;
  Emit('         *************************************************************');
  Emit('         *               RELAXED SURFACE SCAN STEP ' + FmtR(K + 1, 3, 0) + ' OF ' + FmtR(ScanSteps + 1, 3, 0) + '          *');
  Emit('         *                                                           *');
  Emit('         *   Constraint = ' + FmtR(ScanValueToUser(ScanKind, Tgt), 10, 5) + ' ' + PadR(Unit_, 10) + '                       *');
  Emit('         *************************************************************');
  If Not RunOpt(False, Err) Then Begin
   If UserAbort Then Begin NConstr := SaveN; Exit; End;
   Emit('WARNING: ' + Err + ' (etape ' + IStr(K + 1) + ')');
   Err := '';
  End;
  Inc(ScanN);
  ScanX[ScanN] := ScanValueToUser(ScanKind, ConstrValue(C));
  ScanE[ScanN] := ETot;
  Emit('  Etape ' + IStr(K + 1) + ' : coordonnee = ' + RStr(ScanX[ScanN], 5) + '   E = ' + RStr(ScanE[ScanN], 8) + ' Eh');
  ShowProgress('Balayage ' + IStr(K + 1) + '/' + IStr(ScanSteps + 1), 0);
  If UserAbort Then Begin Err := 'Calcul interrompu'; NConstr := SaveN; Exit; End;
 End;
 NConstr := SaveN;
 EmitBlank;
 Emit('       *************************************************************');
 Emit('       *           RELAXED SURFACE SCAN : RESULTATS                *');
 Emit('       *************************************************************');
 EmitBlank;
 Emit('The Calculated Surface using the ''Actual Energy''');
 For K := 1 To ScanN Do
  Emit('   ' + FmtR(ScanX[K], 12, 6) + '   ' + FmtR(ScanE[K], 18, 10));
 RunScan := True;
End;

{ --------------------------------------------------------------------- }
{ Execution du travail decrit par le fichier d'entree                     }
{ --------------------------------------------------------------------- }

Function TimeText(Ms: LongInt): String;
Var S, M, H: LongInt;
Begin
 S := Ms Div 1000; H := S Div 3600; M := (S Mod 3600) Div 60;
 TimeText := '0 days ' + IStr(H) + ' hours ' + IStr(M) + ' minutes ' + IStr(S Mod 60) +
             ' seconds ' + IStr(Ms Mod 1000) + ' msec';
End;

Procedure FinishScf;
Begin
 PopAnalysis;
 ComputeDipole;
 If NOccA <> NOccB Then S2Exp := ComputeS2 Else S2Exp := 0;
End;

Procedure PrintFinalEnergy;
Begin
 EmitBlank;
 Emit('-------------------------   --------------------');
 Emit('FINAL SINGLE POINT ENERGY' + FmtR(ETot, 24, 12));
 Emit('-------------------------   --------------------');
End;

Function ExecuteJob(Var Err: String): Boolean;
Var
 T0: LongInt; E: Real; Ok: Boolean; A: Integer;
Begin
 ExecuteJob := False;
 CalcOK := False; CalcDone := False; UserAbort := False;
 ThermoOK := False; FreqOK := False; CISOK := False; ScanN := 0; OptN := 0; HistN := 0;
 OrbsValid := False; NWarn := 0;
 T0 := NowMs;
 CalcMs := 0;
 OutClear;
 AllocMats;
 ERIThresh := 1E-11;
 NumHStep := 0.005;
 If Not ParseInput(Err) Then Begin
  PrintBanner; PrintInputEcho;
  EmitBlank; Emit('ERREUR : ' + Err);
  EmitBlank; Emit('                             ****ORCA FINISHED BY ERROR****');
  CalcErr := Err;
  Exit;
 End;
 PrintBanner; PrintInputEcho;
 If ScanOn Then JobType := J_SCAN;
 { methode et spin }
 If Not BuildBasis(Err) Then Begin
  EmitBlank; Emit('ERREUR : ' + Err); Emit('                             ****ORCA FINISHED BY ERROR****');
  CalcErr := Err; Exit;
 End;
 If Not SetupElectrons(Err) Then Begin
  EmitBlank; Emit('ERREUR : ' + Err); Emit('                             ****ORCA FINISHED BY ERROR****');
  CalcErr := Err; Exit;
 End;
 If (NOccA <> NOccB) And (Not UnrestrictedFlag) And (Meth = M_HF) Then
  AddWarn('Systeme a couche ouverte : calcul UHF (non restreint).');
 If JobType <> J_SP Then Begin
  { derivees numeriques : seuils serres }
  If ScfTolE > 1E-9 Then ScfTolE := 1E-9;
  If ScfTolD > 1E-8 Then ScfTolD := 1E-8;
  If ScfTolErr > 1E-7 Then ScfTolErr := 1E-7;
  ScfTightName := 'VeryTightSCF';
 End;
 Case JobType Of
  J_SP: PrintBoxed('Single Point Calculation');
  J_OPT: PrintBoxed('Geometry Optimization Run');
  J_FREQ: PrintBoxed('Frequency Calculation');
  J_OPTFREQ: PrintBoxed('Optimization + Frequency Run');
 Else PrintBoxed('Relaxed Surface Scan');
 End;
 PrintGeometry;
 PrintBasisInfo;
 PrintScfSettings;
 If Meth = M_DFT Then Begin
  If NOccA <> NOccB Then Emit(' Reference determinant : Kohn-Sham non restreint (UKS)');
 End;
 Ok := True;
 Case JobType Of
  J_SP, J_FREQ: Begin
   Ok := EvalEnergy(True, False, E, Err);
   If Ok Then Begin
    FinishScf; PrintSCFResults;
    If Meth = M_MP2 Then PrintMP2Results;
    If DoCIS Then Begin RunCIS; PrintCISResults; End;
    PrintFinalEnergy;
    If JobType = J_FREQ Then Begin
     Ok := RunFreq(False, Err);
     If Ok Then Begin
      PrintFreqResults;
      RunThermo(CalcTemp);
      PrintThermoResults;
     End;
    End;
   End;
  End;
  J_OPT, J_OPTFREQ: Begin
   Ok := RunOpt(True, Err);
   If Ok Then Begin
    EmitBlank;
    Emit('                    *****************************************************');
    Emit('                    *** THE OPTIMIZATION HAS CONVERGED (' + IStr(OptN - 1) + ' cycles) ***');
    Emit('                    *****************************************************');
    PrintGeometry;
    EmitBlank;
    Emit('--------------------------------------------------------------------------------');
    Emit('                FINAL ENERGY EVALUATION AT THE STATIONARY POINT');
    Emit('--------------------------------------------------------------------------------');
    FinishScf; PrintSCFResults;
    If Meth = M_MP2 Then PrintMP2Results;
    If DoCIS Then Begin RunCIS; PrintCISResults; End;
    PrintFinalEnergy;
    SaveXYZ(ChangeExt(InpName, '.xyz'), 'Geometrie optimisee par ORCA Clone - E = ' + RStr(ETot, 8) + ' Eh');
    If JobType = J_OPTFREQ Then Begin
     Ok := RunFreq(False, Err);
     If Ok Then Begin
      PrintFreqResults;
      RunThermo(CalcTemp);
      PrintThermoResults;
     End;
    End;
   End;
  End;
 Else Begin
   Ok := RunScan(Err);
   If Ok Then Begin
    FinishScf;
    PrintGeometry;
   End;
  End;
 End;
 CalcMs := ElapsedMs(T0);
 EmitBlank;
 If Ok Then Begin
  Emit('Total run time: ' + TimeText(CalcMs));
  Emit('                             ****ORCA TERMINATED NORMALLY****');
  CalcOK := True;
  CalcDone := True;
 End
 Else Begin
  EmitBlank;
  Emit('ERREUR : ' + Err);
  Emit('Total run time: ' + TimeText(CalcMs));
  Emit('                             ****ORCA FINISHED BY ERROR****');
  CalcErr := Err;
  If NAt > 0 Then CalcDone := OrbsValid;
 End;
 ExecuteJob := Ok;
End;

{ --------------------------------------------------------------------- }
{ Fichiers : entree, sortie                                              }
{ --------------------------------------------------------------------- }

Function LoadInputFile(Const FN: String): Boolean;
Var F: Text; L: String; N: Integer;
Begin
 LoadInputFile := False;
 {$I-}
 Assign(F, FN); Reset(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 N := 0;
 While Not Eof(F) And (N < MAX_INP) Do Begin
  ReadLn(F, L);
  { enleve CR eventuel et tabulations }
  While (Length(L) > 0) And (L[Length(L)] < ' ') Do Delete(L, Length(L), 1);
  InpLines[N] := Copy(L, 1, 100);
  Inc(N);
 End;
 Close(F);
 NInp := N;
 InpName := FN;
 InpModified := False;
 LoadInputFile := True;
End;

Function SaveInputFile(Const FN: String): Boolean;
Var F: Text; I: Integer;
Begin
 SaveInputFile := False;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 For I := 0 To NInp - 1 Do WriteLn(F, InpLines[I]);
 Close(F);
 InpName := FN;
 InpModified := False;
 SaveInputFile := True;
End;

Function SaveOutputFile(Const FN: String): Boolean;
Var F: Text; I: LongInt;
Begin
 SaveOutputFile := False;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Exit;
 For I := 0 To NOut - 1 Do WriteLn(F, OutLines[I]^);
 If OutOverflow Then WriteLn(F, '... (sortie tronquee dans la memoire de l''interface)');
 Close(F);
 SaveOutputFile := True;
End;

{ --------------------------------------------------------------------- }
{ Interface graphique : menus, boites de dialogue, fichiers, graphiques  }
{ (adaptes de CLONEDRAW\PAINT.PAS)                                       }
{ --------------------------------------------------------------------- }

Const
 {$IFDEF FPC}
 MAXPT    = 320;
 {$ELSE}
 MAXPT    = 160;
 {$ENDIF}
 {$IFDEF FPC}
 MAX_ROWS = 24;
 {$ELSE}
 MAX_ROWS = 16;
 {$ENDIF}

Var
 ScreenWidth, ScreenHeight : Integer;
 Menus       : Array[0..NUM_MENUS-1] Of TMenu;
 MenuBarHeight : Integer;
 MenuActive  : Boolean;
 ActiveMenu, ActiveMenuItem : Integer;

Function PtInRect(MX, MY, X1, Y1, X2, Y2: Integer): Boolean;
Begin
 PtInRect := (MX >= X1) And (MX <= X2) And (MY >= Y1) And (MY <= Y2);
End;

Procedure Bevel(X1, Y1, X2, Y2: Integer; Down: Boolean);
Begin
 If Down Then SetColor(DarkGray) Else SetColor(White);
 Line(X1, Y1, X2 - 1, Y1); Line(X1, Y1, X1, Y2 - 1);
 If Down Then SetColor(White) Else SetColor(DarkGray);
 Line(X2, Y1, X2, Y2); Line(X1, Y2, X2, Y2);
End;

Procedure DrawButton(X1, Y1, X2, Y2: Integer; Const Caption: String;
                     Focused: Boolean);
Var TX, TY: Integer;
Begin
 SetFillStyle(SolidFill, LightGray);
 Bar(X1 + 1, Y1 + 1, X2 - 1, Y2 - 1);
 SetColor(White);
 Line(X1, Y1, X2, Y1);
 Line(X1, Y1, X1, Y2);
 SetColor(DarkGray);
 Line(X1 + 1, Y2 - 1, X2 - 1, Y2 - 1);
 Line(X2 - 1, Y1 + 1, X2 - 1, Y2 - 1);
 SetColor(Black);
 Line(X1, Y2, X2, Y2);
 Line(X2, Y1, X2, Y2);
 SetTextStyle(DefaultFont, HorizDir, 1);
 TX := (X1 + X2 - TextWidth(Caption)) Div 2;
 TY := (Y1 + Y2 - TextHeight(Caption)) Div 2;
 SetColor(Black);
 OutTextXY(TX, TY, Caption);
 If Focused Then Begin
  SetColor(Black);
  Rectangle(X1 + 3, Y1 + 3, X2 - 3, Y2 - 3);
 End;
End;

Procedure DrawDialogBox(X1, Y1, X2, Y2: Integer; Const Title: String);
Begin
 SetFillStyle(SolidFill, LightGray);
 Bar(X1, Y1, X2, Y2);
 SetColor(White);
 Rectangle(X1, Y1, X2, Y2);
 SetColor(DarkGray);
 Rectangle(X1 + 1, Y1 + 1, X2 - 1, Y2 - 1);
 SetFillStyle(SolidFill, Blue);
 Bar(X1 + 2, Y1 + 2, X2 - 2, Y1 + 16);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(White);
 OutTextXY(X1 + 6, Y1 + 5, Title);
 SetColor(Black);
End;

{ --------------------------------------------------------------------- }
{ Barre de menus                                                        }
{ --------------------------------------------------------------------- }

Procedure LayoutMenuBar;
Var I, X, W : Integer;
Begin
 X := 6;
 For I := 0 To NUM_MENUS - 1 Do Begin
  W := 8 * Length(Menus[I].Title);
  Menus[I].X := X;
  Menus[I].Y := 0;
  X := X + W + 14;
 End;
 MenuBarHeight := MENUBAR_H;
End;

Procedure DrawMenuBar;
Var
 I, W : Integer;
Begin
 LayoutMenuBar;
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetFillStyle(SolidFill, LightGray);
 Bar(0, 0, ScreenWidth - 1, MenuBarHeight - 1);
 SetColor(DarkGray);
 Line(0, MenuBarHeight - 1, ScreenWidth - 1, MenuBarHeight - 1);
 For I := 0 To NUM_MENUS - 1 Do Begin
  W := 8 * Length(Menus[I].Title);
  If MenuActive And (ActiveMenu = I) Then Begin
   SetFillStyle(SolidFill, Blue);
   Bar(Menus[I].X - 3, Menus[I].Y + 2, Menus[I].X + W + 3, Menus[I].Y + MENUBAR_H - 3);
   SetColor(White);
  End
  Else SetColor(Black);
  OutTextXY(Menus[I].X, Menus[I].Y + 6, Menus[I].Title);
 End;
End;

Procedure DrawMenuPopup(MenuIdx: Integer);
Var
 I, ItemHeight, X, Y, W, H, Wmax, WSc : Integer;
 Cap, Sc : String;
Begin
 If (MenuIdx < 0) Or (MenuIdx >= NUM_MENUS) Then Exit;
 SetTextStyle(DefaultFont, HorizDir, 1);
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
 X := Menus[MenuIdx].X - 3;
 Y := Menus[MenuIdx].Y + MENUBAR_H;
 W := Wmax;
 H := Menus[MenuIdx].ItemCount * ItemHeight + 4;
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
    Bar(X + 2, Y + 2 + I * ItemHeight, X + W - 2, Y + 2 + (I + 1) * ItemHeight);
    SetColor(White);
   End
   Else Begin
    If Menus[MenuIdx].Items[I].Enabled Then SetColor(Black)
    Else SetColor(DarkGray);
   End;
   If Menus[MenuIdx].Items[I].Checked Then
    OutTextXY(X + 5, Y + 4 + I * ItemHeight, #251);
   OutTextXY(X + 22, Y + 4 + I * ItemHeight, Menus[MenuIdx].Items[I].Caption);
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
  X1 := Menus[I].X - 3;
  X2 := Menus[I].X + 8 * Length(Menus[I].Title) + 3;
  If (Mx >= X1) And (Mx <= X2) Then Begin
   MenuBarHitTest := I; Exit;
  End;
 End;
End;

Function MenuItemHitTest(Mx, My, MenuIdx: Integer): Integer;
Var Idx, Top, Bot, X, W: Integer;
Begin
 MenuItemHitTest := -1;
 If (MenuIdx < 0) Or (MenuIdx >= NUM_MENUS) Then Exit;
 Top := Menus[MenuIdx].Y + MENUBAR_H;
 Bot := Top + Menus[MenuIdx].ItemCount * 16 + 4;
 X := Menus[MenuIdx].X - 3;
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
{ Boites de dialogue                                                   }
{ --------------------------------------------------------------------- }

Const
 {$IFDEF FPC}
 MAX_TL = 300;
 {$ELSE}
 MAX_TL = 110;
 {$ENDIF}

Var
 DlgOldLeft : Boolean;
 DlgHappened : Boolean;
 TLines     : Array[0..MAX_TL - 1] Of String[72];
 NTL        : Integer;

Procedure DlgBegin;
Var K: Char;
Begin
 While KeyPressed Do Begin
  K := ReadKey;
  If K = #0 Then If KeyPressed Then K := ReadKey;
 End;
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 DlgOldLeft := (MouseButtons And 1) <> 0;
End;

{ Interroge clavier et souris. Key : 0=rien, 1..255=ASCII, 256+code=touche etendue }
Function DlgPoll(Var Key, MX, MY: Integer; Var Click: Boolean): Boolean;
Var K: Char; L: Boolean;
Begin
 Key := 0; Click := False;
 If KeyPressed Then Begin
  K := ReadKey;
  If K = #0 Then Begin
   If KeyPressed Then Begin K := ReadKey; Key := 256 + Ord(K); End;
  End
  Else Key := Ord(K);
 End;
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 MX := MouseX; MY := MouseY;
 L := (MouseButtons And 1) <> 0;
 Click := L And Not DlgOldLeft;
 DlgOldLeft := L;
 DlgPoll := (Key <> 0) Or Click;
End;

Procedure DlgEnd;
Begin
 DlgHappened := True;
 OldMouseButtons := MouseButtons;
 OldMouseX := MouseX; OldMouseY := MouseY;
 NeedRedraw := True;
End;

Function WordLine(Const S: String; N: Integer): String;
{ N-ieme segment (base 1) d'une chaine separee par '|' }
Var I, K, St: Integer;
Begin
 K := 1; St := 1;
 For I := 1 To Length(S) + 1 Do Begin
  If (I > Length(S)) Or (S[I] = '|') Then Begin
   If K = N Then Begin WordLine := Copy(S, St, I - St); Exit; End;
   Inc(K); St := I + 1;
  End;
 End;
 WordLine := '';
End;

Function CountLines(Const S: String): Integer;
Var I, N: Integer;
Begin
 N := 1;
 For I := 1 To Length(S) Do If S[I] = '|' Then Inc(N);
 CountLines := N;
End;

Procedure ShowMessage(Const Title, Msg: String);
Const BW = 70; BH = 22;
Var
 X1, Y1, X2, Y2, BX1, BY1, BX2, BY2, W, H, I, NL, Key, MX, MY : Integer;
 Click, Done : Boolean;
Begin
 NL := CountLines(Msg);
 W := TextWidth(Title) + 40;
 For I := 1 To NL Do
  If TextWidth(WordLine(Msg, I)) + 40 > W Then W := TextWidth(WordLine(Msg, I)) + 40;
 If W < 240 Then W := 240;
 If W > ScreenWidth - 20 Then W := ScreenWidth - 20;
 H := 70 + NL * 12;
 X1 := (ScreenWidth - W) Div 2; Y1 := (ScreenHeight - H) Div 2;
 X2 := X1 + W; Y2 := Y1 + H;
 BX1 := (X1 + X2 - BW) Div 2; BY1 := Y2 - BH - 8;
 BX2 := BX1 + BW; BY2 := BY1 + BH;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(Black);
 For I := 1 To NL Do OutTextXY(X1 + 20, Y1 + 28 + (I - 1) * 12, WordLine(Msg, I));
 DrawButton(BX1, BY1, BX2, BY2, 'OK', True);
 ShowMouseCursor;
 DlgBegin;
 Done := False;
 While Not Done Do Begin
  If DlgPoll(Key, MX, MY, Click) Then Begin
   If (Key = 13) Or (Key = 27) Or (Key = 32) Then Done := True;
   If Click And PtInRect(MX, MY, BX1, BY1, BX2, BY2) Then Done := True;
  End;
  Delay(10);
 End;
 DlgEnd;
End;

Function ConfirmBox(Const Title, Msg: String): Boolean;
Const BW = 70; BH = 22;
Var
 X1, Y1, X2, Y2, W, H, I, NL, Key, MX, MY, BY1, BY2, B1X, B2X : Integer;
 Click, Done, Res : Boolean;
Begin
 NL := CountLines(Msg);
 W := TextWidth(Title) + 40;
 For I := 1 To NL Do
  If TextWidth(WordLine(Msg, I)) + 40 > W Then W := TextWidth(WordLine(Msg, I)) + 40;
 If W < 260 Then W := 260;
 If W > ScreenWidth - 20 Then W := ScreenWidth - 20;
 H := 70 + NL * 12;
 X1 := (ScreenWidth - W) Div 2; Y1 := (ScreenHeight - H) Div 2;
 X2 := X1 + W; Y2 := Y1 + H;
 BY1 := Y2 - BH - 8; BY2 := BY1 + BH;
 B1X := (X1 + X2) Div 2 - BW - 8; B2X := (X1 + X2) Div 2 + 8;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(Black);
 For I := 1 To NL Do OutTextXY(X1 + 20, Y1 + 28 + (I - 1) * 12, WordLine(Msg, I));
 DrawButton(B1X, BY1, B1X + BW, BY2, 'Oui', True);
 DrawButton(B2X, BY1, B2X + BW, BY2, 'Non', False);
 ShowMouseCursor;
 DlgBegin;
 Done := False; Res := False;
 While Not Done Do Begin
  If DlgPoll(Key, MX, MY, Click) Then Begin
   If (Key = 13) Or (Key = Ord('o')) Or (Key = Ord('O')) Then Begin Res := True; Done := True; End;
   If (Key = 27) Or (Key = Ord('n')) Or (Key = Ord('N')) Then Begin Res := False; Done := True; End;
   If Click And PtInRect(MX, MY, B1X, BY1, B1X + BW, BY2) Then Begin Res := True; Done := True; End;
   If Click And PtInRect(MX, MY, B2X, BY1, B2X + BW, BY2) Then Begin Res := False; Done := True; End;
  End;
  Delay(10);
 End;
 DlgEnd;
 ConfirmBox := Res;
End;

Function InputString(Const Title, Prompt: String; Var Value: String): Boolean;
Const BW = 80; BH = 22;
Var
 X1, Y1, X2, Y2, W, BBY1, BBY2, BOk1, BOk2, BCa1, BCa2, Key, MX, MY : Integer;
 Click, Done, OK, First : Boolean;
Begin
 W := 460;
 If W > ScreenWidth - 10 Then W := ScreenWidth - 10;
 X1 := (ScreenWidth - W) Div 2; Y1 := (ScreenHeight - 140) Div 2;
 X2 := X1 + W; Y2 := Y1 + 140;
 BBY1 := Y2 - BH - 10; BBY2 := BBY1 + BH;
 BOk1 := X1 + W Div 2 - BW - 10; BOk2 := BOk1 + BW;
 BCa1 := X1 + W Div 2 + 10; BCa2 := BCa1 + BW;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(Black);
 OutTextXY(X1 + 15, Y1 + 30, Copy(Prompt, 1, (W - 30) Div 8));
 DrawButton(BOk1, BBY1, BOk2, BBY2, 'OK', False);
 DrawButton(BCa1, BBY1, BCa2, BBY2, 'Annuler', False);
 ShowMouseCursor;
 DlgBegin;
 Done := False; OK := False; First := True;
 While Not Done Do Begin
  HideMouseCursor;
  SetFillStyle(SolidFill, White);
  Bar(X1 + 15, Y1 + 50, X2 - 15, Y1 + 70);
  SetColor(Black);
  Rectangle(X1 + 14, Y1 + 49, X2 - 14, Y1 + 71);
  If Length(Value) > (W - 50) Div 8 Then
   OutTextXY(X1 + 18, Y1 + 56, Copy(Value, Length(Value) - (W - 50) Div 8 + 1, 255) + '_')
  Else OutTextXY(X1 + 18, Y1 + 56, Value + '_');
  ShowMouseCursor;
  Repeat
   If DlgPoll(Key, MX, MY, Click) Then Break;
   Delay(10);
  Until False;
  If Click Then Begin
   If PtInRect(MX, MY, BOk1, BBY1, BOk2, BBY2) Then Begin Done := True; OK := True; End
   Else If PtInRect(MX, MY, BCa1, BBY1, BCa2, BBY2) Then Begin Done := True; OK := False; End;
  End;
  If Key = 13 Then Begin Done := True; OK := True; End
  Else If Key = 27 Then Begin Done := True; OK := False; End
  Else If Key = 8 Then Begin
   If Length(Value) > 0 Then Delete(Value, Length(Value), 1);
   First := False;
  End
  Else If (Key >= 32) And (Key < 256) And (Length(Value) < 250) Then Begin
   If First Then Value := '';
   First := False;
   Value := Value + Chr(Key);
  End
  Else If Key > 0 Then First := False;
 End;
 DlgEnd;
 InputString := OK;
End;

{ Liste / visionneuse de texte dans TLines[0..NTL-1] ; renvoie l'indice choisi }
Function ListDialog(Const Title: String; Selectable: Boolean): Integer;
Const BW = 80; BH = 22; LH = 12;
Var
 X1, Y1, X2, Y2, LX1, LY1, LX2, LY2, MaxLn, Top, Sel, I, Y : Integer;
 B1X, B2X, BY1, BY2, Key, MX, MY, Row : Integer;
 Click, Done, Redraw : Boolean;
 Res : Integer;
 LastT : LongInt;
 LastRow : Integer;
Begin
 X1 := 40; Y1 := 24; X2 := ScreenWidth - 40; Y2 := ScreenHeight - 24;
 LX1 := X1 + 8; LY1 := Y1 + 24; LX2 := X2 - 24; LY2 := Y2 - 40;
 MaxLn := (LY2 - LY1 - 4) Div LH;
 BY1 := Y2 - BH - 8; BY2 := BY1 + BH;
 B1X := X2 - 2 * (BW + 10); B2X := X2 - (BW + 10);
 Top := 0; Sel := 0; Res := -1; Done := False; Redraw := True;
 LastT := 0; LastRow := -2;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 If Selectable Then Begin
  DrawButton(B1X, BY1, B1X + BW, BY2, 'OK', True);
  DrawButton(B2X, BY1, B2X + BW, BY2, 'Annuler', False);
 End
 Else DrawButton(B2X, BY1, B2X + BW, BY2, 'Fermer', True);
 DrawButton(LX2 + 2, LY1, LX2 + 18, LY1 + 16, '^', False);
 DrawButton(LX2 + 2, LY2 - 16, LX2 + 18, LY2, 'v', False);
 ShowMouseCursor;
 DlgBegin;
 While Not Done Do Begin
  If Redraw Then Begin
   HideMouseCursor;
   SetFillStyle(SolidFill, White);
   Bar(LX1, LY1, LX2, LY2);
   SetColor(Black);
   Rectangle(LX1, LY1, LX2, LY2);
   SetTextStyle(DefaultFont, HorizDir, 1);
   For I := 0 To MaxLn - 1 Do Begin
    If Top + I >= NTL Then Break;
    Y := LY1 + 3 + I * LH;
    If Selectable And (Top + I = Sel) Then Begin
     SetFillStyle(SolidFill, Blue);
     Bar(LX1 + 1, Y - 1, LX2 - 1, Y + LH - 2);
     SetColor(White);
    End
    Else SetColor(Black);
    OutTextXY(LX1 + 4, Y, TLines[Top + I]);
   End;
   ShowMouseCursor;
   Redraw := False;
  End;
  If DlgPoll(Key, MX, MY, Click) Then Begin
   Case Key Of
    256 + 72: Begin If Selectable Then Begin If Sel > 0 Then Dec(Sel); If Sel < Top Then Top := Sel; End Else If Top > 0 Then Dec(Top); Redraw := True; End;
    256 + 80: Begin If Selectable Then Begin If Sel < NTL - 1 Then Inc(Sel); If Sel >= Top + MaxLn Then Top := Sel - MaxLn + 1; End Else If Top < NTL - MaxLn Then Inc(Top); Redraw := True; End;
    256 + 73: Begin Dec(Top, MaxLn); If Top < 0 Then Top := 0; If Selectable Then Sel := Top; Redraw := True; End;
    256 + 81: Begin
     Inc(Top, MaxLn); If Top > NTL - MaxLn Then Top := NTL - MaxLn; If Top < 0 Then Top := 0;
     If Selectable Then Sel := Top; Redraw := True;
    End;
    256 + 71: Begin Top := 0; Sel := 0; Redraw := True; End;
    256 + 79: Begin
     Sel := NTL - 1; Top := NTL - MaxLn; If Top < 0 Then Top := 0; Redraw := True;
    End;
    13: Begin If Selectable Then Res := Sel Else Res := 0; Done := True; End;
    27: Begin Res := -1; Done := True; End;
   End;
   If Click Then Begin
    If PtInRect(MX, MY, B2X, BY1, B2X + BW, BY2) Then Begin Res := -1; Done := True; End
    Else If Selectable And PtInRect(MX, MY, B1X, BY1, B1X + BW, BY2) Then Begin Res := Sel; Done := True; End
    Else If PtInRect(MX, MY, LX2 + 2, LY1, LX2 + 18, LY1 + 16) Then Begin If Top > 0 Then Dec(Top, 3); If Top < 0 Then Top := 0; Redraw := True; End
    Else If PtInRect(MX, MY, LX2 + 2, LY2 - 16, LX2 + 18, LY2) Then Begin
     If Top < NTL - MaxLn Then Inc(Top, 3); If Top > NTL - MaxLn Then Top := NTL - MaxLn; If Top < 0 Then Top := 0; Redraw := True;
    End
    Else If Selectable And PtInRect(MX, MY, LX1, LY1, LX2, LY2) Then Begin
     Row := Top + (MY - LY1 - 3) Div LH;
     If (Row >= 0) And (Row < NTL) Then Begin
      If (Row = LastRow) And (NowMs - LastT < 600) And (NowMs >= LastT) Then Begin
       Res := Row; Done := True;
      End;
      Sel := Row; LastRow := Row; LastT := NowMs;
      Redraw := True;
     End;
    End;
   End;
  End;
  Delay(10);
 End;
 DlgEnd;
 ListDialog := Res;
End;

{ --------------------------------------------------------------------- }
{ Boite d'ouverture de fichier                                          }
{ --------------------------------------------------------------------- }

Const
 {$IFDEF FPC}
 FD_MAX_FILES = 256;
 {$ELSE}
 FD_MAX_FILES = 128;
 {$ENDIF}

Type
 TFDEntry = Record
  Name  : String[14];
  IsDir : Boolean;
 End;

Var
 FDFiles : Array[0..FD_MAX_FILES - 1] Of TFDEntry;
 FDCount, FDTop, FDSel : Integer;
 FDMask, FDFileName : String;

Procedure FD_ScanMask(Const Mask: String);
Var Sr: SearchRec;
Begin
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
End;

Procedure FD_ScanDir(Const Mask: String);
Var
 Sr : SearchRec;
 Tmp : TFDEntry;
 I, J, P : Integer;
 M : String;
Begin
 FDCount := 0;
 FindFirst('*', Directory, Sr);
 While (DosError = 0) And (FDCount < FD_MAX_FILES) Do Begin
  If ((Sr.Attr And Directory) <> 0) And (Sr.Name <> '.') Then Begin
   FDFiles[FDCount].Name := Sr.Name;
   FDFiles[FDCount].IsDir := True;
   Inc(FDCount);
  End;
  FindNext(Sr);
 End;
 {$IFDEF FPC}FindClose(Sr);{$ENDIF}
 M := Mask;
 Repeat
  P := Pos(';', M);
  If P > 0 Then Begin FD_ScanMask(Copy(M, 1, P - 1)); Delete(M, 1, P); End
  Else Begin FD_ScanMask(M); M := ''; End;
 Until M = '';
 For I := 1 To FDCount - 1 Do Begin
  Tmp := FDFiles[I];
  J := I - 1;
  While (J >= 0) And
        ((FDFiles[J].IsDir < Tmp.IsDir) Or
         ((FDFiles[J].IsDir = Tmp.IsDir) And (UpStr(FDFiles[J].Name) > UpStr(Tmp.Name)))) Do Begin
   FDFiles[J + 1] := FDFiles[J];
   Dec(J);
  End;
  FDFiles[J + 1] := Tmp;
 End;
End;

Procedure FD_DrawList(X1, Y1, X2, Y2: Integer);
Var I, Y, MaxLn : Integer; S : String;
Begin
 SetFillStyle(SolidFill, White);
 Bar(X1, Y1, X2, Y2);
 SetColor(Black);
 Rectangle(X1, Y1, X2, Y2);
 SetTextStyle(DefaultFont, HorizDir, 1);
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
 X1, Y1, X2, Y2, LX1, LY1, LX2, LY2, BOu1, BAn1, BBY1, BBY2, MaxLn : Integer;
 Key, MX, MY, SimK : Integer;
 Click, Done, Result : Boolean;
 OldDir, CurDir : String;
 NewMask : String;
 LastT : LongInt;
 LastRow, Row : Integer;
Begin
 FileOpenDialog := False;
 GetDir(0, OldDir);
 FDMask := InitialMask;
 FD_ScanDir(FDMask);
 FDTop := 0; FDSel := 0;
 X1 := 50; Y1 := 40; X2 := ScreenWidth - 50; Y2 := ScreenHeight - 40;
 LX1 := X1 + 10; LY1 := Y1 + 50; LX2 := X2 - 10; LY2 := Y2 - 90;
 MaxLn := (LY2 - LY1 - 2) Div 10;
 BBY1 := Y2 - BH - 8; BBY2 := BBY1 + BH;
 BOu1 := X2 - 2 * (BW + 8) - 10;
 BAn1 := X2 - (BW + 8) - 10;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, Title);
 GetDir(0, CurDir);
 SetColor(Black);
 OutTextXY(X1 + 10, Y1 + 22, 'Dossier : ' + Copy(CurDir, 1, 52));
 OutTextXY(X1 + 10, Y1 + 35, 'Masque : ' + FDMask + '   (Tab = changer le masque)');
 OutTextXY(X1 + 10, Y2 - 80, 'Fichier:');
 DrawButton(BOu1, BBY1, BOu1 + BW, BBY2, 'Ouvrir', True);
 DrawButton(BAn1, BBY1, BAn1 + BW, BBY2, 'Annuler', False);
 ShowMouseCursor;
 Done := False; Result := False; FDFileName := '';
 LastT := 0; LastRow := -2;
 DlgBegin;
 While Not Done Do Begin
  HideMouseCursor;
  SetFillStyle(SolidFill, White);
  Bar(X1 + 80, Y2 - 83, X2 - 10, Y2 - 70);
  SetColor(Black);
  Rectangle(X1 + 80, Y2 - 84, X2 - 9, Y2 - 69);
  OutTextXY(X1 + 84, Y2 - 81, FDFileName + '_');
  FD_DrawList(LX1, LY1, LX2, LY2);
  ShowMouseCursor;
  SimK := 0;
  Repeat
   If DlgPoll(Key, MX, MY, Click) Then Begin
    If Click Then Begin
     If PtInRect(MX, MY, BOu1, BBY1, BOu1 + BW, BBY2) Then SimK := 13
     Else If PtInRect(MX, MY, BAn1, BBY1, BAn1 + BW, BBY2) Then SimK := 27
     Else If PtInRect(MX, MY, LX1, LY1, LX2, LY2) Then Begin
      Row := FDTop + (MY - LY1 - 2) Div 10;
      If (Row >= 0) And (Row < FDCount) Then Begin
       If (Row = LastRow) And (NowMs - LastT < 600) And (NowMs >= LastT) Then SimK := 13;
       FDSel := Row; LastRow := Row; LastT := NowMs;
       FDFileName := '';
      End;
     End;
    End;
    If SimK = 0 Then SimK := Key;
    If SimK = 0 Then SimK := -1;
    Break;
   End;
   Delay(10);
  Until False;
  Case SimK Of
   -1: ;
   256 + 72: Begin If FDSel > 0 Then Dec(FDSel); If FDSel < FDTop Then FDTop := FDSel; End;
   256 + 80: Begin
    If FDSel < FDCount - 1 Then Inc(FDSel);
    If FDSel >= FDTop + MaxLn Then FDTop := FDSel - MaxLn + 1;
   End;
   256 + 73: Begin Dec(FDSel, MaxLn); If FDSel < 0 Then FDSel := 0; FDTop := FDSel; End;
   256 + 81: Begin
    Inc(FDSel, MaxLn); If FDSel >= FDCount Then FDSel := FDCount - 1;
    FDTop := FDSel - MaxLn + 1; If FDTop < 0 Then FDTop := 0;
   End;
   256 + 71: Begin FDSel := 0; FDTop := 0; End;
   256 + 79: Begin
    FDSel := FDCount - 1; FDTop := FDSel - MaxLn + 1; If FDTop < 0 Then FDTop := 0;
   End;
   13: Begin
    If FDFileName <> '' Then Begin
     FileName := FDFileName; Result := True; Done := True;
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
      Bar(X1 + 10, Y1 + 20, X2 - 5, Y1 + 31);
      SetColor(Black);
      OutTextXY(X1 + 10, Y1 + 22, 'Dossier : ' + Copy(CurDir, 1, 52));
      ShowMouseCursor;
     End
     Else Begin
      FileName := FDFiles[FDSel].Name; Result := True; Done := True;
     End;
    End;
   End;
   27: Begin Result := False; Done := True; End;
   9: Begin
    NewMask := FDMask;
    If InputString('Masque', 'Entrer le masque (ex: *.CDT;*.MOL):', NewMask) Then Begin
     FDMask := NewMask; FD_ScanDir(FDMask); FDTop := 0; FDSel := 0;
    End;
    HideMouseCursor;
    DrawDialogBox(X1, Y1, X2, Y2, Title);
    SetColor(Black);
    OutTextXY(X1 + 10, Y1 + 22, 'Dossier : ' + Copy(CurDir, 1, 52));
    OutTextXY(X1 + 10, Y1 + 35, 'Masque : ' + FDMask + '   (Tab = changer le masque)');
    OutTextXY(X1 + 10, Y2 - 80, 'Fichier:');
    DrawButton(BOu1, BBY1, BOu1 + BW, BBY2, 'Ouvrir', True);
    DrawButton(BAn1, BBY1, BAn1 + BW, BBY2, 'Annuler', False);
    ShowMouseCursor;
   End;
   8: Begin If Length(FDFileName) > 0 Then Delete(FDFileName, Length(FDFileName), 1); End;
  Else
   If (SimK >= 32) And (SimK < 256) And (Length(FDFileName) < 60) Then
    FDFileName := FDFileName + Chr(SimK);
  End;
 End;
 If Not Result Then ChDir(OldDir);
 DlgEnd;
 FileOpenDialog := Result;
End;

Function FileSaveDialog(Const Title: String; Var FileName: String): Boolean;
Var Cur: String;
Begin
 Cur := FileName;
 If InputString(Title, 'Nom du fichier :', Cur) And (Cur <> '') Then Begin
  FileName := Cur;
  FileSaveDialog := True;
 End
 Else FileSaveDialog := False;
End;


{ --------------------------------------------------------------------- }
{ Boite de parametres generique                                         }
{ --------------------------------------------------------------------- }

Var
 PN     : Integer;
 PLabel : Array[0..MAX_ROWS-1] Of String[44];
 PVal   : Array[0..MAX_ROWS-1] Of Real;
 PMin   : Array[0..MAX_ROWS-1] Of Real;
 PMax   : Array[0..MAX_ROWS-1] Of Real;
 PKind  : Array[0..MAX_ROWS-1] Of Byte;   { 0 reel, 1 entier, 2 choix, 3 oui/non, 4 titre, 5 texte }
 PDec   : Array[0..MAX_ROWS-1] Of Byte;
 PEnum  : Array[0..MAX_ROWS-1] Of String[90];

Procedure PReset;
Begin
 PN := 0;
End;

Procedure PAdd(Const L: String; K: Byte; V, Mn, Mx: Real; Dc: Byte; Const E: String);
Begin
 If PN >= MAX_ROWS Then Exit;
 PLabel[PN] := L; PKind[PN] := K; PVal[PN] := V; PMin[PN] := Mn; PMax[PN] := Mx;
 PDec[PN] := Dc; PEnum[PN] := E;
 Inc(PN);
End;

Procedure PAddR(Const L: String; V, Mn, Mx: Real; Dc: Integer);
Begin PAdd(L, 0, V, Mn, Mx, Dc, ''); End;

Procedure PAddI(Const L: String; V, Mn, Mx: LongInt);
Begin PAdd(L, 1, V, Mn, Mx, 0, ''); End;

Procedure PAddE(Const L: String; V: Integer; Const E: String);
Begin PAdd(L, 2, V, 0, 0, 0, E); End;

Procedure PAddB(Const L: String; V: Boolean);
Begin If V Then PAdd(L, 3, 1, 0, 1, 0, '') Else PAdd(L, 3, 0, 0, 1, 0, ''); End;

Procedure PAddH(Const L: String);
Begin PAdd(L, 4, 0, 0, 0, 0, ''); End;

Procedure PAddT(Const L, Txt: String);
Begin PAdd(L, 5, 0, 0, 0, 0, Txt); End;

Function PValStr(I: Integer): String;
Begin
 Case PKind[I] Of
  0: PValStr := RTrim(PVal[I], PDec[I]);
  1: PValStr := IStr(Round(PVal[I]));
  2: PValStr := WordLine(PEnum[I], Round(PVal[I]) + 1);
  3: If PVal[I] > 0.5 Then PValStr := 'Oui' Else PValStr := 'Non';
  5: PValStr := PEnum[I];
 Else PValStr := '';
 End;
End;

Function PBool(I: Integer): Boolean;
Begin PBool := PVal[I] > 0.5; End;

Function PInt(I: Integer): LongInt;
Begin PInt := Round(PVal[I]); End;

Procedure PEdit(I: Integer);
Var S: String; V: Real;
Begin
 Case PKind[I] Of
  0, 1: Begin
   S := PValStr(I);
   If InputString(PLabel[I], 'Nouvelle valeur (' + RTrim(PMin[I], 4) + ' a ' + RTrim(PMax[I], 4) + ') :', S) Then Begin
    If ParseReal(S, V) Then Begin
     If (PMax[I] > PMin[I]) And ((V < PMin[I]) Or (V > PMax[I])) Then
      ShowMessage('Valeur refusee', 'La valeur doit etre comprise entre|' + RTrim(PMin[I], 4) + ' et ' + RTrim(PMax[I], 4))
     Else Begin
      If PKind[I] = 1 Then V := Round(V);
      PVal[I] := V;
     End;
    End
    Else ShowMessage('Valeur refusee', 'Ce n''est pas un nombre valide.');
   End;
  End;
  2: Begin
   PVal[I] := PVal[I] + 1;
   If PVal[I] > CountLines(PEnum[I]) - 1 Then PVal[I] := 0;
  End;
  3: If PVal[I] > 0.5 Then PVal[I] := 0 Else PVal[I] := 1;
 End;
End;

Function ParamDialog(Const Title, Hint: String): Boolean;
Const BW = 90; BH = 22; RH = 16;
Var
 X1, Y1, X2, Y2, RY, I, Sel, Top, Vis, Key, MX, MY, BOk1, BCa1, BBY1, BBY2, Row, LabW : Integer;
 Click, Done, Res, Redraw : Boolean;
 LastT : LongInt; LastRow : Integer;

 Procedure NextSel(Step: Integer);
 Var K: Integer;
 Begin
  K := Sel;
  Repeat
   K := K + Step;
   If (K < 0) Or (K >= PN) Then Exit;
  Until PKind[K] < 4;
  Sel := K;
  If Sel < Top Then Top := Sel;
  If Sel >= Top + Vis Then Top := Sel - Vis + 1;
 End;

Begin
 X1 := 50; Y1 := 26; X2 := ScreenWidth - 50; Y2 := ScreenHeight - 26;
 RY := Y1 + 24;
 Vis := (Y2 - Y1 - 24 - 52) Div RH;
 BBY2 := Y2 - 8; BBY1 := BBY2 - BH;
 BOk1 := X2 - 2 * (BW + 10); BCa1 := X2 - (BW + 10);
 LabW := 0;
 For I := 0 To PN - 1 Do If 8 * Length(PLabel[I]) + 28 > LabW Then LabW := 8 * Length(PLabel[I]) + 28;
 If LabW < 170 Then LabW := 170;
 If LabW > 360 Then LabW := 360;
 Sel := 0; Top := 0; Done := False; Res := False; Redraw := True;
 LastT := 0; LastRow := -1;
 If (PN > 0) And (PKind[0] >= 4) Then NextSel(1);
 DlgBegin;
 While Not Done Do Begin
  If Redraw Then Begin
   HideMouseCursor;
   DrawDialogBox(X1, Y1, X2, Y2, Title);
   SetTextStyle(DefaultFont, HorizDir, 1);
   For I := 0 To Vis - 1 Do Begin
    Row := Top + I;
    If Row >= PN Then Break;
    If PKind[Row] = 4 Then Begin
     SetColor(Blue);
     OutTextXY(X1 + 10, RY + I * RH + 4, PLabel[Row]);
     SetColor(DarkGray);
     Line(X1 + 10, RY + I * RH + RH - 2, X2 - 10, RY + I * RH + RH - 2);
    End
    Else Begin
     If Row = Sel Then Begin
      SetFillStyle(SolidFill, Blue);
      Bar(X1 + 6, RY + I * RH, X2 - 6, RY + I * RH + RH - 2);
      SetColor(White);
     End
     Else SetColor(Black);
     OutTextXY(X1 + 12, RY + I * RH + 4, PLabel[Row]);
     If PKind[Row] < 4 Then Begin
      If Row <> Sel Then Begin
       SetFillStyle(SolidFill, White);
       Bar(X1 + LabW, RY + I * RH, X2 - 12, RY + I * RH + RH - 3);
       SetColor(Black);
      End;
      OutTextXY(X1 + LabW + 4, RY + I * RH + 4, Copy(PValStr(Row), 1, (X2 - X1 - LabW - 20) Div 8));
     End
     Else Begin
      SetColor(Black);
      If Row = Sel Then SetColor(White);
      OutTextXY(X1 + LabW + 4, RY + I * RH + 4, Copy(PValStr(Row), 1, (X2 - X1 - LabW - 20) Div 8));
     End;
    End;
   End;
   SetColor(DarkGray);
   If Hint <> '' Then OutTextXY(X1 + 10, BBY1 - 14, Copy(Hint, 1, (X2 - X1 - 20) Div 8));
   If PN > Vis Then OutTextXY(X1 + 10, BBY1 + 6, 'Haut/bas : defiler (' + IStr(Sel + 1) + '/' + IStr(PN) + ')')
   Else OutTextXY(X1 + 10, BBY1 + 6, 'Entree : modifier  Tab : OK  Echap : annuler');
   DrawButton(BOk1, BBY1, BOk1 + BW, BBY2, 'OK', True);
   DrawButton(BCa1, BBY1, BCa1 + BW, BBY2, 'Annuler', False);
   ShowMouseCursor;
   Redraw := False;
  End;
  If DlgPoll(Key, MX, MY, Click) Then Begin
   Case Key Of
    256 + 72: Begin NextSel(-1); Redraw := True; End;
    256 + 80: Begin NextSel(1); Redraw := True; End;
    256 + 73: Begin For I := 1 To Vis Do NextSel(-1); Redraw := True; End;
    256 + 81: Begin For I := 1 To Vis Do NextSel(1); Redraw := True; End;
    256 + 75: If (Sel < PN) And (PKind[Sel] = 2) Then Begin
     PVal[Sel] := PVal[Sel] - 1;
     If PVal[Sel] < 0 Then PVal[Sel] := CountLines(PEnum[Sel]) - 1;
     Redraw := True;
    End;
    256 + 77: If (Sel < PN) And ((PKind[Sel] = 2) Or (PKind[Sel] = 3)) Then Begin PEdit(Sel); Redraw := True; End;
    32: If (Sel < PN) And ((PKind[Sel] = 2) Or (PKind[Sel] = 3)) Then Begin PEdit(Sel); Redraw := True; End;
    13: Begin
     If (Sel < PN) And ((PKind[Sel] = 0) Or (PKind[Sel] = 1)) Then Begin PEdit(Sel); Redraw := True; End
     Else If (Sel < PN) And ((PKind[Sel] = 2) Or (PKind[Sel] = 3)) Then Begin PEdit(Sel); Redraw := True; End
     Else Begin Res := True; Done := True; End;
    End;
    27: Begin Res := False; Done := True; End;
    9: Begin Res := True; Done := True; End;
   End;
   If Click Then Begin
    If PtInRect(MX, MY, BOk1, BBY1, BOk1 + BW, BBY2) Then Begin Res := True; Done := True; End
    Else If PtInRect(MX, MY, BCa1, BBY1, BCa1 + BW, BBY2) Then Begin Res := False; Done := True; End
    Else If PtInRect(MX, MY, X1 + 6, RY, X2 - 6, RY + Vis * RH) Then Begin
     Row := Top + (MY - RY) Div RH;
     If (Row >= 0) And (Row < PN) And (PKind[Row] < 4) Then Begin
      Sel := Row;
      If (PKind[Row] >= 2) And (PKind[Row] <= 3) And (MX >= X1 + LabW) Then PEdit(Row)
      Else If (Row = LastRow) And (NowMs - LastT < 600) And (NowMs >= LastT) Then PEdit(Row);
      LastRow := Row; LastT := NowMs;
      Redraw := True;
     End;
    End;
   End;
  End;
  Delay(10);
 End;
 DlgEnd;
 ParamDialog := Res;
End;

{ --------------------------------------------------------------------- }
{ Lignes de texte et visionneuse                                        }
{ --------------------------------------------------------------------- }

Procedure TClear;
Begin
 NTL := 0;
End;

Procedure TAdd(Const S: String);
Begin
 If NTL < MAX_TL Then Begin
  TLines[NTL] := Copy(S, 1, 72);
  Inc(NTL);
 End;
End;

Procedure TView(Const Title: String);
Var R: Integer;
Begin
 R := ListDialog(Title, False);
End;

{ --------------------------------------------------------------------- }
{ Graphiques                                                            }
{ --------------------------------------------------------------------- }

Const
 {$IFDEF FPC}
 MAXSER = 5;
 {$ELSE}
 MAXSER = 3;
 {$ENDIF}
 SerColors : Array[0..4] Of Byte = (1, 4, 2, 5, 6);

Var
 ChNS    : Integer;
 ChN     : Array[0..MAXSER-1] Of Integer;
 ChX, ChY : Array[0..MAXSER-1, 0..MAXPT-1] Of Real;
 ChName  : Array[0..MAXSER-1] Of String[32];
 ChMark  : Array[0..MAXSER-1] Of Boolean;
 ChTitle, ChXL, ChYL : String[64];
 ChLogX  : Boolean;
 ChFixY0 : Boolean;
 ChLogMin : Real;

Procedure ChartReset(Const T, XL, YL: String; LogX: Boolean);
Var I: Integer;
Begin
 ChNS := 0; ChTitle := T; ChXL := XL; ChYL := YL; ChLogX := LogX;
 ChFixY0 := True;
 If LogX Then ChLogMin := 0.1 Else ChLogMin := 0;
 For I := 0 To MAXSER - 1 Do Begin ChN[I] := 0; ChMark[I] := False; ChName[I] := ''; End;
End;

Function ChartAddSeries(Const Name: String; Mark: Boolean): Integer;
Begin
 If ChNS >= MAXSER Then Begin ChartAddSeries := MAXSER - 1; Exit; End;
 ChName[ChNS] := Name; ChMark[ChNS] := Mark; ChN[ChNS] := 0;
 ChartAddSeries := ChNS;
 Inc(ChNS);
End;

Procedure ChartPoint(S: Integer; X, Y: Real);
Begin
 If (S < 0) Or (S >= MAXSER) Or (ChN[S] >= MAXPT) Then Exit;
 If ChLogX And (X < ChLogMin) Then Exit;
 ChX[S, ChN[S]] := X; ChY[S, ChN[S]] := Y;
 Inc(ChN[S]);
End;

Procedure NiceRange(Lo, Hi: Real; Var NLo, NHi, Step: Real);
Var Raw, Mag, R: Real;
Begin
 If Hi - Lo < 1E-12 Then Begin Hi := Lo + 1; End;
 Raw := (Hi - Lo) / 6;
 Mag := 1;
 While Raw >= 10 * Mag Do Mag := Mag * 10;
 While Raw < Mag Do Mag := Mag / 10;
 R := Raw / Mag;
 If R < 1.5 Then Step := Mag
 Else If R < 3 Then Step := 2 * Mag
 Else If R < 7 Then Step := 5 * Mag
 Else Step := 10 * Mag;
 NLo := Int(Lo / Step) * Step;
 If NLo > Lo + 1E-12 Then NLo := NLo - Step;
 NHi := Int(Hi / Step) * Step;
 If NHi < Hi - 1E-12 Then NHi := NHi + Step;
End;

Procedure ChartExport;
Var F: Text; FN: String; S, I: Integer;
Begin
 FN := 'GRAPHE.CSV';
 If Not InputString('Exporter le graphique', 'Nom du fichier CSV :', FN) Then Exit;
 If FN = '' Then Exit;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire ' + FN); Exit; End;
 WriteLn(F, ChTitle);
 For S := 0 To ChNS - 1 Do Begin
  WriteLn(F, ChName[S]);
  WriteLn(F, ChXL, ';', ChYL);
  For I := 0 To ChN[S] - 1 Do WriteLn(F, ChX[S, I]:0:6, ';', ChY[S, I]:0:6);
 End;
 Close(F);
End;

Procedure ChartDialog;
Const BW = 90; BH = 22;
Var
 X1, Y1, X2, Y2, PX1, PY1, PX2, PY2, BBY1, BBY2, BOk1, BCs1 : Integer;
 S, I, Key, MX, MY, OldMX, OldMY, XP, YP, XP0, YP0 : Integer;
 XLo, XHi, YLo, YHi, XSt, YSt, V, XMin, XMax, YMin, YMax : Real;
 Click, Done, Any : Boolean;
 BestD, D : Real; BS, BI : Integer; Info : String;

 Function SXp(V: Real): Integer;
 Begin
  If ChLogX Then Begin
   If V <= 0 Then V := XLo;
   SXp := PX1 + Round((Ln(V) - Ln(XLo)) / (Ln(XHi) - Ln(XLo)) * (PX2 - PX1));
  End
  Else SXp := PX1 + Round((V - XLo) / (XHi - XLo) * (PX2 - PX1));
 End;

 Function SYp(V: Real): Integer;
 Begin
  SYp := PY2 - Round((V - YLo) / (YHi - YLo) * (PY2 - PY1));
 End;

Begin
 X1 := 24; Y1 := 16; X2 := ScreenWidth - 24; Y2 := ScreenHeight - 16;
 PX1 := X1 + 62; PX2 := X2 - 18; PY1 := Y1 + 40; PY2 := Y2 - 100;
 BBY2 := Y2 - 8; BBY1 := BBY2 - BH;
 BOk1 := X2 - BW - 12; BCs1 := BOk1 - BW - 10;
 { bornes }
 Any := False; XMin := 1E30; XMax := -1E30; YMin := 1E30; YMax := -1E30;
 For S := 0 To ChNS - 1 Do
  For I := 0 To ChN[S] - 1 Do Begin
   Any := True;
   If ChLogX And (ChX[S, I] <= 0) Then Continue;
   If ChX[S, I] < XMin Then XMin := ChX[S, I];
   If ChX[S, I] > XMax Then XMax := ChX[S, I];
   If ChY[S, I] < YMin Then YMin := ChY[S, I];
   If ChY[S, I] > YMax Then YMax := ChY[S, I];
  End;
 If Not Any Then Begin XMin := 0; XMax := 1; YMin := 0; YMax := 1; End;
 If ChFixY0 And (YMin > 0) Then YMin := 0;
 If YMax <= YMin Then YMax := YMin + 1;
 NiceRange(YMin, YMax, YLo, YHi, YSt);
 If ChLogX Then Begin
  XLo := 1; While XLo > XMin Do XLo := XLo / 10;
  While XLo * 10 <= XMin Do XLo := XLo * 10;
  XHi := XLo;
  While XHi < XMax Do XHi := XHi * 10;
  If XHi <= XLo Then XHi := XLo * 10;
  XSt := 10;
 End
 Else Begin
  If XMax <= XMin Then XMax := XMin + 1;
  NiceRange(XMin, XMax, XLo, XHi, XSt);
 End;
 HideMouseCursor;
 DrawDialogBox(X1, Y1, X2, Y2, ChTitle);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetFillStyle(SolidFill, White);
 Bar(PX1, PY1, PX2, PY2);
 { grille et graduations }
 SetColor(Black);
 V := YLo;
 While V <= YHi + YSt * 0.001 Do Begin
  YP := SYp(V);
  SetColor(LightGray); Line(PX1, YP, PX2, YP);
  SetColor(Black); Line(PX1 - 3, YP, PX1, YP);
  OutTextXY(PX1 - 6 - 8 * Length(RTrim(V, 3)), YP - 3, RTrim(V, 3));
  V := V + YSt;
 End;
 If ChLogX Then Begin
  V := XLo;
  While V <= XHi * 1.001 Do Begin
   XP := SXp(V);
   SetColor(LightGray); Line(XP, PY1, XP, PY2);
   SetColor(Black); Line(XP, PY2, XP, PY2 + 3);
   OutTextXY(XP - 4 * Length(RTrim(V, 3)), PY2 + 6, RTrim(V, 3));
   V := V * 10;
  End;
 End
 Else Begin
  V := XLo;
  While V <= XHi + XSt * 0.001 Do Begin
   XP := SXp(V);
   SetColor(LightGray); Line(XP, PY1, XP, PY2);
   SetColor(Black); Line(XP, PY2, XP, PY2 + 3);
   OutTextXY(XP - 4 * Length(RTrim(V, 3)), PY2 + 6, RTrim(V, 3));
   V := V + XSt;
  End;
 End;
 SetColor(Black);
 Rectangle(PX1, PY1, PX2, PY2);
 OutTextXY((PX1 + PX2) Div 2 - 4 * Length(ChXL), PY2 + 20, ChXL);
 OutTextXY(PX1 - 54, PY1 - 14, ChYL);
 { courbes }
 SetViewPort(PX1 + 1, PY1 + 1, PX2 - 1, PY2 - 1, True);
 For S := 0 To ChNS - 1 Do Begin
  SetColor(SerColors[S]);
  For I := 0 To ChN[S] - 1 Do Begin
   XP := SXp(ChX[S, I]) - PX1 - 1; YP := SYp(ChY[S, I]) - PY1 - 1;
   If ChMark[S] Then Begin
    Rectangle(XP - 2, YP - 2, XP + 2, YP + 2);
   End
   Else If I > 0 Then
    Line(XP0, YP0, XP, YP);
   XP0 := XP; YP0 := YP;
  End;
 End;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
 { legende }
 For S := 0 To ChNS - 1 Do Begin
  SetFillStyle(SolidFill, SerColors[S]);
  Bar(X1 + 14 + (S Mod 2) * 280, PY2 + 36 + (S Div 2) * 12, X1 + 26 + (S Mod 2) * 280, PY2 + 42 + (S Div 2) * 12);
  SetColor(Black);
  OutTextXY(X1 + 32 + (S Mod 2) * 280, PY2 + 36 + (S Div 2) * 12, Copy(ChName[S], 1, 32));
 End;
 DrawButton(BCs1, BBY1, BCs1 + BW, BBY2, 'Exporter CSV', False);
 DrawButton(BOk1, BBY1, BOk1 + BW, BBY2, 'Fermer', True);
 ShowMouseCursor;
 DlgBegin;
 Done := False; OldMX := -1; OldMY := -1; Info := '';
 While Not Done Do Begin
  If DlgPoll(Key, MX, MY, Click) Then Begin
   If (Key = 13) Or (Key = 27) Or (Key = 32) Then Done := True;
   If Click Then Begin
    If PtInRect(MX, MY, BOk1, BBY1, BOk1 + BW, BBY2) Then Done := True
    Else If PtInRect(MX, MY, BCs1, BBY1, BCs1 + BW, BBY2) Then Begin
     ChartExport;
     HideMouseCursor;
     SetFillStyle(SolidFill, LightGray);
     Bar(0, 0, ScreenWidth - 1, ScreenHeight - 1);
     ShowMouseCursor;
     Done := True;
    End;
   End;
  End;
  If (MX <> OldMX) Or (MY <> OldMY) Then Begin
   OldMX := MX; OldMY := MY;
   If PtInRect(MX, MY, PX1, PY1, PX2, PY2) Then Begin
    BestD := 1E30; BS := -1; BI := -1;
    For S := 0 To ChNS - 1 Do
     For I := 0 To ChN[S] - 1 Do Begin
      D := Sqr(SXp(ChX[S, I]) - MX) + Sqr(SYp(ChY[S, I]) - MY);
      If D < BestD Then Begin BestD := D; BS := S; BI := I; End;
     End;
    If (BS >= 0) And (BestD < 900) Then
     Info := ChName[BS] + ' : ' + RTrim(ChX[BS, BI], 4) + ' ; ' + RTrim(ChY[BS, BI], 4)
    Else Info := '';
   End
   Else Info := '';
   HideMouseCursor;
   SetFillStyle(SolidFill, LightGray);
   Bar(X1 + 14, BBY1 + 2, BCs1 - 8, BBY2);
   SetColor(Blue);
   OutTextXY(X1 + 14, BBY1 + 6, Copy(Info, 1, (BCs1 - X1 - 30) Div 8));
   ShowMouseCursor;
  End;
  Delay(10);
 End;
 DlgEnd;
End;

{ --------------------------------------------------------------------- }
{ Disposition de la fenetre principale                                   }
{ --------------------------------------------------------------------- }

Const
 ED_X1 = 0;    ED_Y1 = 44;   ED_X2 = 307;  ED_Y2 = 461;
 VW_X1 = 310;  VW_Y1 = 44;   VW_X2 = 639;  VW_Y2 = 299;
 CN_X1 = 310;  CN_Y1 = 302;  CN_X2 = 639;  CN_Y2 = 461;
 ED_CHARS = 37;
 ED_ROWS  = 39;

 TB_RUN = 1;  TB_OUT = 2;  TB_MOL = 3;  TB_ORB = 4;  TB_SCF = 5;
 TB_OPT = 6;  TB_IR = 7;   TB_UV = 8;   TB_SCAN = 9; TB_EX = 10;
 TB_HELP = 11;
 MAX_TB = 12;

 STYLE_WIRE = 0; STYLE_STICK = 1; STYLE_BALL = 2; STYLE_CPK = 3;

Var
 { editeur }
 CurLn, CurCol, TopLn, LeftCol : Integer;
 InsertMode : Boolean;
 { vue }
 ViewStyle  : Integer;
 ShowSym, ShowNum, ShowChg, ShowAxes : Boolean;
 RotM       : Array[1..3, 1..3] Of Real;
 ZoomF      : Real;
 DragOn     : Boolean;
 DragX, DragY : Integer;
 { geometries (apercu de l'entree et resultat) }
 PrevN      : Integer;
 PrevZ      : Array[1..MAX_AT] Of Integer;
 PrevX      : Array[1..MAX_AT, 1..3] Of Real;
 PrevOK     : Boolean;
 PrevErr    : String;
 ResN       : Integer;
 ResZ       : Array[1..MAX_AT] Of Integer;
 ResX       : Array[1..MAX_AT, 1..3] Of Real;
 ResCharge, ResMult : Integer;
 ResMethod  : String;
 ResMeth    : Integer;
 ResStale   : Boolean;
 HaveResult : Boolean;
 { animation d'un mode de vibration }
 VibMode    : Integer;
 VibPhase   : Real;
 VibLast    : LongInt;
 { barre d'outils }
 NTB        : Integer;
 TBId       : Array[0..MAX_TB-1] Of Integer;
 TBX1, TBX2 : Array[0..MAX_TB-1] Of Integer;
 TBCap      : Array[0..MAX_TB-1] Of String[6];
 TBHint     : Array[0..MAX_TB-1] Of String[60];
 CalcRunning : Boolean;
 BusyTick   : Integer;
 SearchText : String;
 LastConsoleN : LongInt;

{ --------------------------------------------------------------------- }
{ Editeur de texte (fichier d'entree)                                    }
{ --------------------------------------------------------------------- }

Procedure EnsureOneLine;
Begin
 If NInp < 1 Then Begin NInp := 1; InpLines[0] := ''; End;
End;

Procedure MarkModified;
Begin
 InpModified := True;
 ResStale := True;
End;

Procedure EdScroll;
Begin
 If CurLn < TopLn Then TopLn := CurLn;
 If CurLn >= TopLn + ED_ROWS Then TopLn := CurLn - ED_ROWS + 1;
 If CurCol < LeftCol Then LeftCol := CurCol;
 If CurCol >= LeftCol + ED_CHARS Then LeftCol := CurCol - ED_CHARS + 1;
 If TopLn < 0 Then TopLn := 0;
 If LeftCol < 0 Then LeftCol := 0;
End;

Procedure EdClamp;
Begin
 EnsureOneLine;
 If CurLn < 0 Then CurLn := 0;
 If CurLn >= NInp Then CurLn := NInp - 1;
 If CurCol < 0 Then CurCol := 0;
 If CurCol > Length(InpLines[CurLn]) Then CurCol := Length(InpLines[CurLn]);
End;

Procedure EdInsertLine(At: Integer; Const S: String);
Var I: Integer;
Begin
 If NInp >= MAX_INP Then Exit;
 For I := NInp DownTo At + 1 Do InpLines[I] := InpLines[I - 1];
 InpLines[At] := Copy(S, 1, 100);
 Inc(NInp);
End;

Procedure EdDeleteLine(At: Integer);
Var I: Integer;
Begin
 If (At < 0) Or (At >= NInp) Then Exit;
 For I := At To NInp - 2 Do InpLines[I] := InpLines[I + 1];
 Dec(NInp);
 EnsureOneLine;
End;

Procedure EdTypeChar(C: Char);
Var L: String;
Begin
 EdClamp;
 L := InpLines[CurLn];
 If Length(L) >= 100 Then Exit;
 If InsertMode Or (CurCol >= Length(L)) Then Insert(C, L, CurCol + 1)
 Else L[CurCol + 1] := C;
 InpLines[CurLn] := L;
 Inc(CurCol);
 MarkModified;
End;

Procedure EdEnter;
Var L, Rest: String; Ind: Integer;
Begin
 EdClamp;
 L := InpLines[CurLn];
 Rest := Copy(L, CurCol + 1, 100);
 Ind := 0;
 While (Ind < Length(L)) And (L[Ind + 1] = ' ') And (Ind < CurCol) Do Inc(Ind);
 Rest := Copy('                                        ', 1, Ind) + Rest;
 InpLines[CurLn] := Copy(L, 1, CurCol);
 EdInsertLine(CurLn + 1, Rest);
 Inc(CurLn); CurCol := Ind;
 MarkModified;
End;

Procedure EdBackspace;
Var L: String;
Begin
 EdClamp;
 If CurCol > 0 Then Begin
  L := InpLines[CurLn];
  Delete(L, CurCol, 1);
  InpLines[CurLn] := L;
  Dec(CurCol);
  MarkModified;
 End
 Else If CurLn > 0 Then Begin
  CurCol := Length(InpLines[CurLn - 1]);
  InpLines[CurLn - 1] := Copy(InpLines[CurLn - 1] + InpLines[CurLn], 1, 100);
  EdDeleteLine(CurLn);
  Dec(CurLn);
  MarkModified;
 End;
End;

Procedure EdDelete;
Var L: String;
Begin
 EdClamp;
 L := InpLines[CurLn];
 If CurCol < Length(L) Then Begin
  Delete(L, CurCol + 1, 1);
  InpLines[CurLn] := L;
  MarkModified;
 End
 Else If CurLn < NInp - 1 Then Begin
  InpLines[CurLn] := Copy(L + InpLines[CurLn + 1], 1, 100);
  EdDeleteLine(CurLn + 1);
  MarkModified;
 End;
End;

Procedure EdClear;
Begin
 NInp := 1; InpLines[0] := '';
 CurLn := 0; CurCol := 0; TopLn := 0; LeftCol := 0;
 MarkModified;
End;

Procedure EdAppend(Const S: String);
Begin
 If (NInp = 1) And (InpLines[0] = '') Then Begin InpLines[0] := Copy(S, 1, 100); End
 Else If NInp < MAX_INP Then Begin InpLines[NInp] := Copy(S, 1, 100); Inc(NInp); End;
End;

{ --------------------------------------------------------------------- }
{ Rotation de la vue                                                     }
{ --------------------------------------------------------------------- }

Procedure ResetView;
Var C1, S1, C2, S2: Real; I, J: Integer;
Begin
 { rotation de -0.55 rad autour de X puis 0.45 rad autour de Y }
 C1 := Cos(-0.55); S1 := Sin(-0.55);
 C2 := Cos(0.45); S2 := Sin(0.45);
 For I := 1 To 3 Do For J := 1 To 3 Do RotM[I, J] := 0;
 { RotY * RotX }
 RotM[1, 1] := C2;       RotM[1, 2] := S2 * S1;  RotM[1, 3] := S2 * C1;
 RotM[2, 1] := 0;        RotM[2, 2] := C1;       RotM[2, 3] := -S1;
 RotM[3, 1] := -S2;      RotM[3, 2] := C2 * S1;  RotM[3, 3] := C2 * C1;
 ZoomF := 1.0;
End;

Procedure RotateView(DX, DY: Real);
Var C, S: Real; R: Array[1..3, 1..3] Of Real; I, J, K: Integer; T: Array[1..3, 1..3] Of Real;
Begin
 { rotation autour de Y (ecran) de DX puis autour de X de DY }
 C := Cos(DX); S := Sin(DX);
 T[1,1] := C; T[1,2] := 0; T[1,3] := S;
 T[2,1] := 0; T[2,2] := 1; T[2,3] := 0;
 T[3,1] := -S; T[3,2] := 0; T[3,3] := C;
 For I := 1 To 3 Do For J := 1 To 3 Do Begin
  R[I, J] := 0;
  For K := 1 To 3 Do R[I, J] := R[I, J] + T[I, K] * RotM[K, J];
 End;
 For I := 1 To 3 Do For J := 1 To 3 Do RotM[I, J] := R[I, J];
 C := Cos(DY); S := Sin(DY);
 T[1,1] := 1; T[1,2] := 0; T[1,3] := 0;
 T[2,1] := 0; T[2,2] := C; T[2,3] := -S;
 T[3,1] := 0; T[3,2] := S; T[3,3] := C;
 For I := 1 To 3 Do For J := 1 To 3 Do Begin
  R[I, J] := 0;
  For K := 1 To 3 Do R[I, J] := R[I, J] + T[I, K] * RotM[K, J];
 End;
 For I := 1 To 3 Do For J := 1 To 3 Do RotM[I, J] := R[I, J];
End;

{ --------------------------------------------------------------------- }
{ Apercu de la geometrie du fichier d'entree                             }
{ --------------------------------------------------------------------- }

Var
 SvN, SvCharge, SvMult : Integer;
 SvZ : Array[1..MAX_AT] Of Integer;
 SvX : Array[1..MAX_AT, 1..3] Of Real;
 Sv2N : Integer;
 Sv2Z : Array[1..MAX_AT] Of Integer;
 Sv2X : Array[1..MAX_AT, 1..3] Of Real;

{ ParseInput modifie les variables globales du calcul : on sauvegarde
  la geometrie du calcul termine pour pouvoir la restaurer ensuite }
Procedure PreviewInput;
Var Err: String; A, D: Integer;
Begin
 EnsureOneLine;
 { sauvegarde }
 Sv2N := NAt;
 For A := 1 To NAt Do Begin
  Sv2Z[A] := AtZ[A];
  For D := 1 To 3 Do Sv2X[A, D] := AtX[A, D];
 End;
 SvCharge := Charge; SvMult := Mult;
 PrevOK := ParseInput(Err);
 If PrevOK Then Begin
  PrevN := NAt;
  For A := 1 To NAt Do Begin
   PrevZ[A] := AtZ[A];
   For D := 1 To 3 Do PrevX[A, D] := AtX[A, D];
  End;
  PrevErr := '';
 End
 Else Begin PrevN := 0; PrevErr := Err; End;
 { restauration }
 NAt := Sv2N;
 For A := 1 To NAt Do Begin
  AtZ[A] := Sv2Z[A];
  For D := 1 To 3 Do AtX[A, D] := Sv2X[A, D];
 End;
 Charge := SvCharge; Mult := SvMult;
End;

Procedure CaptureResults;
Var A, D: Integer;
Begin
 ResN := NAt;
 For A := 1 To NAt Do Begin
  ResZ[A] := AtZ[A];
  For D := 1 To 3 Do ResX[A, D] := AtX[A, D];
 End;
 ResCharge := Charge; ResMult := Mult;
 ResMethod := MethodName + ' / ' + BasisNames[BasisSel];
 ResMeth := Meth;
 HaveResult := True;
 ResStale := False;
End;

{ Active la geometrie du resultat dans AtX / AtZ (pour les coupes d'orbitales) }
Procedure UseResultGeom;
Var A, D: Integer;
Begin
 SvN := NAt;
 For A := 1 To NAt Do Begin
  SvZ[A] := AtZ[A];
  For D := 1 To 3 Do SvX[A, D] := AtX[A, D];
 End;
 NAt := ResN;
 For A := 1 To ResN Do Begin
  AtZ[A] := ResZ[A];
  For D := 1 To 3 Do AtX[A, D] := ResX[A, D];
 End;
End;

Procedure UsePrevGeom;
Var A, D: Integer;
Begin
 SvN := NAt;
 For A := 1 To NAt Do Begin
  SvZ[A] := AtZ[A];
  For D := 1 To 3 Do SvX[A, D] := AtX[A, D];
 End;
 NAt := PrevN;
 For A := 1 To PrevN Do Begin
  AtZ[A] := PrevZ[A];
  For D := 1 To 3 Do AtX[A, D] := PrevX[A, D];
 End;
End;

Procedure RestoreGeom;
Var A, D: Integer;
Begin
 NAt := SvN;
 For A := 1 To SvN Do Begin
  AtZ[A] := SvZ[A];
  For D := 1 To 3 Do AtX[A, D] := SvX[A, D];
 End;
End;

{ --------------------------------------------------------------------- }
{ Dessin de la molecule                                                  }
{ --------------------------------------------------------------------- }

{ Formule brute de la molecule (Src = 1 resultat, 0 apercu) }
Function MolFormula(Src: Integer): String;
Var Cnt: Array[1..NELEM] Of Integer; I, N, E: Integer; S: String;
 Procedure Add(E: Integer);
 Begin
  If Cnt[E] > 0 Then Begin
   S := S + ElSym[E];
   If Cnt[E] > 1 Then S := S + IStr(Cnt[E]);
  End;
 End;
Begin
 FillChar(Cnt, SizeOf(Cnt), 0);
 If Src = 1 Then Begin
  N := ResN;
  For I := 1 To N Do If (ResZ[I] >= 1) And (ResZ[I] <= NELEM) Then Inc(Cnt[ResZ[I]]);
 End
 Else Begin
  N := PrevN;
  For I := 1 To N Do If (PrevZ[I] >= 1) And (PrevZ[I] <= NELEM) Then Inc(Cnt[PrevZ[I]]);
 End;
 S := '';
 If Cnt[6] > 0 Then Begin Add(6); Add(1); End Else Add(1);
 For E := 2 To NELEM Do If E <> 6 Then Add(E);
 MolFormula := S;
End;

Procedure DrawSphere(CX, CY, R, Col: Integer);
Begin
 If R < 1 Then R := 1;
 SetFillStyle(SolidFill, Col);
 SetColor(Col);
 FillEllipse(CX, CY, R, R);
 If R >= 4 Then Begin
  SetColor(Black);
  Circle(CX, CY, R);
  SetFillStyle(SolidFill, White);
  SetColor(White);
  FillEllipse(CX - R Div 3, CY - R Div 3, R Div 4 + 1, R Div 4 + 1);
 End;
End;

Procedure DrawMolecule;
Const MAXB = 160;
Var
 N, I, J, K, NB, NItems, BX, BY, CXs, CYs, X1, Y1, X2, Y2, R, Rn, Col: Integer;
 Z: Array[1..MAX_AT] Of Integer;
 X: Array[1..MAX_AT, 1..3] Of Real;
 P: Array[1..MAX_AT, 1..3] Of Real;
 SX, SY: Array[1..MAX_AT] Of Integer;
 Cm: Array[1..3] Of Real; Rmax, Sc, D, Mx, DispAmp, Ph: Real;
 BA, BB: Array[1..MAXB] Of Integer;
 ItK, ItA, ItB: Array[1..MAX_AT + MAXB] Of Integer;
 ItZ: Array[1..MAX_AT + MAXB] Of Real;
 Tz: Real; Tk, Ta, Tb: Integer; Msg, Lab: String;
 UseRes: Boolean; V: Array[1..3] Of Real; DV: Array[1..MAX_AT, 1..3] Of Real; MaxD, Fc: Real;
 RadF: Real; Thick: Integer;
Begin
 SetViewPort(VW_X1, VW_Y1, VW_X2, VW_Y2, True);
 SetFillStyle(SolidFill, Black);
 Bar(0, 0, VW_X2 - VW_X1, VW_Y2 - VW_Y1);
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
 UseRes := HaveResult And Not ResStale;
 If UseRes Then Begin
  N := ResN;
  For I := 1 To N Do Begin Z[I] := ResZ[I]; For J := 1 To 3 Do X[I, J] := ResX[I, J]; End;
 End
 Else If PrevOK Then Begin
  N := PrevN;
  For I := 1 To N Do Begin Z[I] := PrevZ[I]; For J := 1 To 3 Do X[I, J] := PrevX[I, J]; End;
 End
 Else Begin
  SetColor(LightRed);
  SetTextStyle(DefaultFont, HorizDir, 1);
  OutTextXY(VW_X1 + 8, VW_Y1 + 30, 'Geometrie non valide :');
  SetColor(LightGray);
  OutTextXY(VW_X1 + 8, VW_Y1 + 44, Copy(PrevErr, 1, 40));
  OutTextXY(VW_X1 + 8, VW_Y1 + 54, Copy(PrevErr, 41, 40));
  Exit;
 End;
 If N < 1 Then Exit;
 { deplacement du mode de vibration }
 FillChar(DV, SizeOf(DV), 0);
 If UseRes And (VibMode > 0) And FreqOK And (VibMode <= NFreq) And (FreqVec <> Nil) Then Begin
  MaxD := 0;
  For I := 1 To N Do
   For J := 1 To 3 Do Begin
    DV[I, J] := FreqVec^[(I - 1) * 3 + J, VibMode] / Sqrt(ElMass[Z[I]]);
    If Abs(DV[I, J]) > MaxD Then MaxD := Abs(DV[I, J]);
   End;
  If MaxD > 1E-9 Then Fc := 0.6 * Sin(VibPhase) / MaxD Else Fc := 0;
  For I := 1 To N Do
   For J := 1 To 3 Do X[I, J] := X[I, J] + Fc * DV[I, J];
 End;
 Cm[1] := 0; Cm[2] := 0; Cm[3] := 0;
 For I := 1 To N Do For J := 1 To 3 Do Cm[J] := Cm[J] + X[I, J] / N;
 Rmax := 1;
 For I := 1 To N Do Begin
  D := Sqrt(Sqr(X[I,1] - Cm[1]) + Sqr(X[I,2] - Cm[2]) + Sqr(X[I,3] - Cm[3]));
  If D > Rmax Then Rmax := D;
 End;
 Sc := ZoomF * 0.5 * (VW_Y2 - VW_Y1 - 30) / (Rmax + 1.3);
 CXs := (VW_X1 + VW_X2) Div 2;
 CYs := (VW_Y1 + VW_Y2) Div 2 + 6;
 For I := 1 To N Do Begin
  For J := 1 To 3 Do V[J] := X[I, J] - Cm[J];
  For J := 1 To 3 Do P[I, J] := RotM[J, 1] * V[1] + RotM[J, 2] * V[2] + RotM[J, 3] * V[3];
  SX[I] := CXs + Round(P[I, 1] * Sc);
  SY[I] := CYs - Round(P[I, 2] * Sc);
 End;
 { liaisons }
 NB := 0;
 For I := 1 To N - 1 Do
  For J := I + 1 To N Do Begin
   D := Sqrt(Sqr(X[I,1] - X[J,1]) + Sqr(X[I,2] - X[J,2]) + Sqr(X[I,3] - X[J,3]));
   If D < 1.25 * (ElRad[Z[I]] + ElRad[Z[J]]) * ANG2BOHR Then
    If NB < MAXB Then Begin Inc(NB); BA[NB] := I; BB[NB] := J; End;
  End;
 { elements a dessiner, tries par profondeur }
 NItems := 0;
 For I := 1 To NB Do Begin
  Inc(NItems); ItK[NItems] := 1; ItA[NItems] := BA[I]; ItB[NItems] := BB[I];
  ItZ[NItems] := 0.5 * (P[BA[I], 3] + P[BB[I], 3]) - 0.0001;
 End;
 For I := 1 To N Do Begin
  Inc(NItems); ItK[NItems] := 2; ItA[NItems] := I; ItB[NItems] := 0; ItZ[NItems] := P[I, 3];
 End;
 For I := 2 To NItems Do Begin
  Tz := ItZ[I]; Tk := ItK[I]; Ta := ItA[I]; Tb := ItB[I];
  J := I - 1;
  While (J >= 1) And (ItZ[J] > Tz) Do Begin
   ItZ[J + 1] := ItZ[J]; ItK[J + 1] := ItK[J]; ItA[J + 1] := ItA[J]; ItB[J + 1] := ItB[J];
   Dec(J);
  End;
  ItZ[J + 1] := Tz; ItK[J + 1] := Tk; ItA[J + 1] := Ta; ItB[J + 1] := Tb;
 End;
 SetViewPort(VW_X1, VW_Y1, VW_X2, VW_Y2, True);
 For K := 1 To NItems Do Begin
  If ItK[K] = 1 Then Begin
   I := ItA[K]; J := ItB[K];
   X1 := SX[I] - VW_X1; Y1 := SY[I] - VW_Y1; X2 := SX[J] - VW_X1; Y2 := SY[J] - VW_Y1;
   BX := (X1 + X2) Div 2; BY := (Y1 + Y2) Div 2;
   Case ViewStyle Of
    STYLE_WIRE: Thick := 1;
    STYLE_STICK: Thick := 3;
    STYLE_BALL: Thick := 3;
   Else Thick := 1;
   End;
   If ViewStyle <> STYLE_CPK Then Begin
    If ViewStyle = STYLE_BALL Then Begin
     SetLineStyle(SolidLn, 0, ThickWidth);
     SetColor(LightGray);
     Line(X1, Y1, X2, Y2);
    End
    Else Begin
     If Thick = 3 Then SetLineStyle(SolidLn, 0, ThickWidth) Else SetLineStyle(SolidLn, 0, NormWidth);
     SetColor(ElColor[Z[I]]); Line(X1, Y1, BX, BY);
     SetColor(ElColor[Z[J]]); Line(BX, BY, X2, Y2);
    End;
    SetLineStyle(SolidLn, 0, NormWidth);
   End;
  End
  Else Begin
   I := ItA[K];
   X1 := SX[I] - VW_X1; Y1 := SY[I] - VW_Y1;
   Col := ElColor[Z[I]];
   Case ViewStyle Of
    STYLE_WIRE: Begin
     R := 0;
     { atomes isoles uniquement }
     Rn := 0;
     For J := 1 To NB Do If (BA[J] = I) Or (BB[J] = I) Then Inc(Rn);
     If Rn = 0 Then Begin SetColor(Col); Circle(X1, Y1, 3); End;
    End;
    STYLE_STICK: Begin
     Rn := 0;
     For J := 1 To NB Do If (BA[J] = I) Or (BB[J] = I) Then Inc(Rn);
     If Rn = 0 Then DrawSphere(X1, Y1, 4, Col);
    End;
    STYLE_BALL: Begin
     R := Round(0.36 * ElRad[Z[I]] * ANG2BOHR * Sc) + 2;
     If R < 3 Then R := 3;
     DrawSphere(X1, Y1, R, Col);
    End;
   Else Begin
     R := Round(0.95 * ElRad[Z[I]] * ANG2BOHR * Sc) + 2;
     DrawSphere(X1, Y1, R, Col);
    End;
   End;
  End;
 End;
 { etiquettes }
 SetTextStyle(DefaultFont, HorizDir, 1);
 For I := 1 To N Do Begin
  Lab := '';
  If ShowSym Then Lab := ElemStr(Z[I]);
  If ShowNum Then Lab := Lab + IStr(I - 1);
  If ShowChg And UseRes And (I <= NAt) Then Lab := Lab + ' ' + RStr(QMull[I], 2);
  If Lab <> '' Then Begin
   SetColor(Yellow);
   OutTextXY(SX[I] - VW_X1 + 5, SY[I] - VW_Y1 - 12, Lab);
  End;
 End;
 { axes }
 If ShowAxes Then Begin
  X1 := 26; Y1 := VW_Y2 - VW_Y1 - 26;
  For K := 1 To 3 Do Begin
   Case K Of 1: SetColor(LightRed); 2: SetColor(LightGreen); Else SetColor(LightBlue); End;
   Line(X1, Y1, X1 + Round(RotM[1, K] * 18), Y1 - Round(RotM[2, K] * 18));
   Case K Of 1: Lab := 'x'; 2: Lab := 'y'; Else Lab := 'z'; End;
   OutTextXY(X1 + Round(RotM[1, K] * 22) - 3, Y1 - Round(RotM[2, K] * 22) - 3, Lab);
  End;
 End;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
 { titre }
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(White);
 If UseRes Then Begin
  Msg := MolFormula(1) + '  E = ' + RStr(ETot, 6) + ' Eh';
  If (VibMode > 0) And FreqOK Then Msg := 'Mode ' + IStr(VibMode - 1) + ' : ' + RStr(FreqV[VibMode], 1) + ' cm-1';
 End
 Else Msg := MolFormula(0) + '  (apercu de l''entree)';
 OutTextXY(VW_X1 + 6, VW_Y1 + 4, Copy(Msg, 1, 40));
End;

{ --------------------------------------------------------------------- }
{ Console (fin de la sortie)                                             }
{ --------------------------------------------------------------------- }

Const LH = 10;

Function Collapse(Const S: String): String;
Var R: String; I, Run: Integer;
Begin
 R := ''; Run := 0;
 For I := 1 To Length(S) Do Begin
  If S[I] = ' ' Then Inc(Run) Else Run := 0;
  If Run <= 2 Then R := R + S[I];
 End;
 While (Length(R) > 0) And (R[Length(R)] = ' ') Do Delete(R, Length(R), 1);
 Collapse := R;
End;

Procedure DrawConsole;
Var Rows, I, First: LongInt; Y: Integer; S: String;
Begin
 SetViewPort(CN_X1, CN_Y1, CN_X2, CN_Y2, True);
 SetFillStyle(SolidFill, Black);
 Bar(0, 0, CN_X2 - CN_X1, CN_Y2 - CN_Y1);
 SetFillStyle(SolidFill, DarkGray);
 Bar(0, 0, CN_X2 - CN_X1, 12);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(White);
 OutTextXY(6, 3, 'Console ORCA  (' + IStr(NOut) + ' lignes)');
 Rows := (CN_Y2 - CN_Y1 - 16) Div LH;
 First := NOut - Rows;
 If First < 0 Then First := 0;
 SetColor(LightGray);
 Y := 16;
 For I := First To NOut - 1 Do Begin
  S := OutLines[I]^;
  If Pos('ERREUR', S) = 1 Then SetColor(LightRed)
  Else If Pos('WARNING', S) = 1 Then SetColor(Yellow)
  Else If (Pos('****ORCA TERMINATED', S) > 0) Then SetColor(LightGreen)
  Else SetColor(LightGray);
  OutTextXY(4, Y, Copy(Collapse(S), 1, (CN_X2 - CN_X1 - 8) Div 8));
  Inc(Y, LH);
 End;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
End;

Procedure ConsoleUpdate;
Begin
 If Not GuiActive Then Exit;
 HideMouseCursor;
 DrawConsole;
 ShowMouseCursor;
End;

{ --------------------------------------------------------------------- }
{ Dessin de l'editeur                                                    }
{ --------------------------------------------------------------------- }

Procedure DrawEditor;
Var I, Y, Row: Integer; L, S, T: String; C: Byte; Hdr: String;
Begin
 SetFillStyle(SolidFill, White);
 Bar(ED_X1, ED_Y1 + 14, ED_X2, ED_Y2);
 SetFillStyle(SolidFill, Blue);
 Bar(ED_X1, ED_Y1, ED_X2, ED_Y1 + 13);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(White);
 Hdr := 'Entree : ' + BaseName(InpName);
 If InpName = '' Then Hdr := 'Entree : (sans nom)';
 If InpModified Then Hdr := Hdr + ' *';
 OutTextXY(ED_X1 + 4, ED_Y1 + 3, Copy(Hdr, 1, 26));
 Hdr := 'L' + IStr(CurLn + 1) + ' C' + IStr(CurCol + 1);
 OutTextXY(ED_X2 - 8 * Length(Hdr) - 4, ED_Y1 + 3, Hdr);
 SetViewPort(ED_X1, ED_Y1 + 14, ED_X2, ED_Y2, True);
 For Row := 0 To ED_ROWS - 1 Do Begin
  I := TopLn + Row;
  If I >= NInp Then Break;
  Y := 3 + Row * 10;
  L := InpLines[I];
  S := Copy(L, LeftCol + 1, ED_CHARS);
  T := TrimStr(L);
  If Length(T) > 0 Then Begin
   Case T[1] Of
    '!': C := Blue;
    '%': C := Magenta;
    '#': C := DarkGray;
    '*': C := Green;
   Else C := Black;
   End;
  End
  Else C := Black;
  SetColor(C);
  OutTextXY(4, Y, S);
  If (I = CurLn) And ((CurCol - LeftCol) >= 0) And ((CurCol - LeftCol) < ED_CHARS) Then Begin
   SetFillStyle(SolidFill, Black);
   If InsertMode Then Bar(4 + (CurCol - LeftCol) * 8, Y - 1, 5 + (CurCol - LeftCol) * 8, Y + 8)
   Else Bar(4 + (CurCol - LeftCol) * 8, Y - 1, 11 + (CurCol - LeftCol) * 8, Y + 8);
   If (Not InsertMode) And (CurCol < Length(L)) Then Begin
    SetColor(White);
    OutTextXY(4 + (CurCol - LeftCol) * 8, Y, L[CurCol + 1]);
   End;
  End;
 End;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
End;

{ --------------------------------------------------------------------- }
{ Barre d'outils                                                         }
{ --------------------------------------------------------------------- }

Procedure AddTB(Id: Integer; Const Cap, Hint: String; Var X: Integer);
Var W: Integer;
Begin
 W := 8 * Length(Cap) + 8;
 If W < 22 Then W := 22;
 TBId[NTB] := Id; TBCap[NTB] := Cap; TBHint[NTB] := Hint;
 TBX1[NTB] := X; TBX2[NTB] := X + W - 1;
 X := X + W + 1;
 Inc(NTB);
End;

Procedure InitToolbar;
Var X: Integer;
Begin
 NTB := 0; X := 3;
 AddTB(TB_RUN, 'Exec', 'Executer le calcul (F5)', X);
 AddTB(TB_OUT, 'Sortie', 'Afficher le fichier de sortie .OUT (F3)', X);
 X := X + 6;
 AddTB(TB_MOL, 'Style', 'Changer le style d''affichage de la molecule', X);
 AddTB(TB_ORB, 'Orb.', 'Diagramme des orbitales moleculaires (F6)', X);
 X := X + 6;
 AddTB(TB_SCF, 'SCF', 'Convergence du SCF', X);
 AddTB(TB_OPT, 'Opt.', 'Courbe de l''optimisation de geometrie', X);
 AddTB(TB_IR, 'IR', 'Spectre infrarouge (F4)', X);
 AddTB(TB_UV, 'UV', 'Spectre UV-visible (CIS)', X);
 AddTB(TB_SCAN, 'Scan', 'Surface d''energie potentielle (balayage)', X);
 X := X + 6;
 AddTB(TB_EX, 'Exemp.', 'Charger un exemple', X);
 AddTB(TB_HELP, '?', 'Aide (F1)', X);
End;

Function ToolbarHit(MX, MY: Integer): Integer;
Var I: Integer;
Begin
 ToolbarHit := -1;
 If (MY < MENUBAR_H) Or (MY >= MENUBAR_H + TOOLBAR_H) Then Exit;
 For I := 0 To NTB - 1 Do
  If (MX >= TBX1[I]) And (MX <= TBX2[I]) Then Begin ToolbarHit := I; Exit; End;
End;

Function TBEnabled(Id: Integer): Boolean;
Begin
 Case Id Of
  TB_OUT: TBEnabled := NOut > 0;
  TB_ORB: TBEnabled := HaveResult And OrbsValid;
  TB_SCF: TBEnabled := HistN > 1;
  TB_OPT: TBEnabled := OptN > 1;
  TB_IR: TBEnabled := FreqOK;
  TB_UV: TBEnabled := CISOK;
  TB_SCAN: TBEnabled := ScanN > 1;
 Else TBEnabled := True;
 End;
End;

Procedure DrawToolbar;
Var I, Y1, Y2: Integer;
Begin
 Y1 := MENUBAR_H; Y2 := MENUBAR_H + TOOLBAR_H - 1;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, Y1, ScreenWidth - 1, Y2);
 SetColor(DarkGray);
 Line(0, Y2, ScreenWidth - 1, Y2);
 SetTextStyle(DefaultFont, HorizDir, 1);
 For I := 0 To NTB - 1 Do Begin
  SetFillStyle(SolidFill, LightGray);
  Bar(TBX1[I], Y1 + 2, TBX2[I], Y2 - 2);
  Bevel(TBX1[I], Y1 + 2, TBX2[I], Y2 - 2, False);
  If TBEnabled(TBId[I]) Then SetColor(Black) Else SetColor(DarkGray);
  OutTextXY(TBX1[I] + 4, Y1 + 8, TBCap[I]);
 End;
End;

{ --------------------------------------------------------------------- }
{ Barre d'etat                                                           }
{ --------------------------------------------------------------------- }

Procedure DrawStatusBar;
Var Y, W, MaxC, K: Integer; S, R: String;
Begin
 Y := ScreenHeight - STATUS_H;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, Y, ScreenWidth - 1, ScreenHeight - 1);
 SetColor(DarkGray);
 Line(0, Y, ScreenWidth - 1, Y);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(Black);
 If CalcRunning Then Begin
  W := 120;
  SetFillStyle(SolidFill, White);
  Bar(ScreenWidth - W - 8, Y + 4, ScreenWidth - 8, Y + STATUS_H - 4);
  SetFillStyle(SolidFill, Blue);
  K := (BusyTick * 12) Mod (W - 30);
  Bar(ScreenWidth - W - 8 + K, Y + 4, ScreenWidth - W - 8 + K + 30, Y + STATUS_H - 4);
  SetColor(Black);
  Rectangle(ScreenWidth - W - 8, Y + 4, ScreenWidth - 8, Y + STATUS_H - 4);
  OutTextXY(6, Y + 5, Copy(StatusMsg, 1, (ScreenWidth - W - 30) Div 8));
 End
 Else Begin
  R := '';
  If HaveResult And Not ResStale And CalcOK Then R := 'E = ' + RStr(ETot, 6) + ' Eh';
  OutTextXY(6, Y + 5, Copy(StatusMsg, 1, (ScreenWidth - 8 * Length(R) - 30) Div 8));
  SetColor(Blue);
  OutTextXY(ScreenWidth - 8 * Length(R) - 8, Y + 5, R);
 End;
End;

{ Mise a jour pendant les calculs : console, barre d'etat, touche Echap }
Procedure ShowProgress(Const Msg: String; Pct: Integer);
Var K: Char; T: LongInt;
Begin
 If BatchMode Or Not GuiActive Then Exit;
 While KeyPressed Do Begin
  K := ReadKey;
  If K = #27 Then UserAbort := True
  Else If K = #0 Then If KeyPressed Then K := ReadKey;
 End;
 T := NowMs;
 If ElapsedMs(LastShowT) < 250 Then Exit;
 LastShowT := T;
 If Msg <> '' Then StatusMsg := Msg + '   (Echap = arreter)';
 Inc(BusyTick);
 HideMouseCursor;
 If NOut <> LastConsoleN Then Begin DrawConsole; LastConsoleN := NOut; End;
 DrawStatusBar;
 ShowMouseCursor;
End;

Procedure RedrawAll;
Begin
 HideMouseCursor;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, MENUBAR_H + TOOLBAR_H, ScreenWidth - 1, ScreenHeight - 1);
 DrawEditor;
 DrawMolecule;
 DrawConsole;
 DrawToolbar;
 DrawStatusBar;
 DrawMenuBar;
 If MenuActive Then DrawMenuPopup(ActiveMenu);
 ShowMouseCursor;
 NeedRedraw := False;
End;

{ --------------------------------------------------------------------- }
{ Visionneuse du fichier de sortie                                      }
{ --------------------------------------------------------------------- }

Function ContainsText(Const S, Sub: String): Boolean;
Begin
 ContainsText := Pos(UpStr(Sub), UpStr(S)) > 0;
End;

Function OutFileName: String;
Begin
 If InpName = '' Then OutFileName := 'ORCA.OUT' Else OutFileName := BaseName(ChangeExt(InpName, '.out'));
End;

Procedure OutputViewer;
Const LH = 10;
Var
 X1, Y1, X2, Y2, LX1, LY1, LX2, LY2, MaxLn, Top, Left, I, Y, Key, MX, MY, Found, Cols: Integer;
 Click, Done, Redraw: Boolean;
 S, Q: String; FN: String; Hi: Integer;
Begin
 X1 := 0; Y1 := 0; X2 := ScreenWidth - 1; Y2 := ScreenHeight - 1;
 LX1 := 4; LY1 := 16; LX2 := X2 - 4; LY2 := Y2 - 16;
 MaxLn := (LY2 - LY1 - 4) Div LH;
 Cols := (LX2 - LX1 - 8) Div 8;
 Top := 0; Left := 0; Hi := -1; Done := False; Redraw := True;
 HideMouseCursor;
 SetFillStyle(SolidFill, LightGray);
 Bar(X1, Y1, X2, Y2);
 SetFillStyle(SolidFill, Blue);
 Bar(X1, Y1, X2, Y1 + 13);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(White);
 OutTextXY(6, 3, 'Sortie ORCA : ' + OutFileName + '   (' + IStr(NOut) + ' lignes)');
 SetColor(Black);
 OutTextXY(6, Y2 - 11, 'Fleches/PgUp/PgDn/Debut/Fin F=chercher N=suivant S=enregistrer Echap=fermer');
 ShowMouseCursor;
 DlgBegin;
 While Not Done Do Begin
  If Redraw Then Begin
   HideMouseCursor;
   SetFillStyle(SolidFill, Black);
   Bar(LX1, LY1, LX2, LY2);
   SetTextStyle(DefaultFont, HorizDir, 1);
   For I := 0 To MaxLn - 1 Do Begin
    If Top + I >= NOut Then Break;
    Y := LY1 + 3 + I * LH;
    S := OutLines[Top + I]^;
    If Top + I = Hi Then Begin
     SetFillStyle(SolidFill, Blue);
     Bar(LX1 + 1, Y - 1, LX2 - 1, Y + LH - 2);
     SetColor(White);
    End
    Else If Pos('ERREUR', S) = 1 Then SetColor(LightRed)
    Else If (Pos('FINAL SINGLE POINT', S) = 1) Or (Pos('Total Energy', S) = 1) Then SetColor(Yellow)
    Else If (Pos('---', S) = 1) Or (Pos('***', S) > 0) Or (Pos('===', S) = 1) Then SetColor(LightCyan)
    Else If Pos('WARNING', S) = 1 Then SetColor(Yellow)
    Else SetColor(LightGray);
    OutTextXY(LX1 + 4, Y, Copy(S, Left + 1, Cols));
   End;
   If NOut = 0 Then Begin SetColor(LightGray); OutTextXY(LX1 + 8, LY1 + 8, '(aucune sortie : executez un calcul avec F5)'); End;
   ShowMouseCursor;
   Redraw := False;
  End;
  If DlgPoll(Key, MX, MY, Click) Then Begin
   Case Key Of
    256 + 72: If Top > 0 Then Begin Dec(Top); Redraw := True; End;
    256 + 80: If Top < NOut - MaxLn Then Begin Inc(Top); Redraw := True; End;
    256 + 73: Begin Dec(Top, MaxLn - 1); If Top < 0 Then Top := 0; Redraw := True; End;
    256 + 81: Begin Inc(Top, MaxLn - 1); If Top > NOut - MaxLn Then Top := NOut - MaxLn; If Top < 0 Then Top := 0; Redraw := True; End;
    256 + 71: Begin Top := 0; Left := 0; Redraw := True; End;
    256 + 79: Begin Top := NOut - MaxLn; If Top < 0 Then Top := 0; Redraw := True; End;
    256 + 75: If Left > 0 Then Begin Dec(Left, 8); If Left < 0 Then Left := 0; Redraw := True; End;
    256 + 77: If Left < 40 Then Begin Inc(Left, 8); Redraw := True; End;
    Ord('f'), Ord('F'), 6, Ord('n'), Ord('N'): Begin
     If (Key = Ord('f')) Or (Key = Ord('F')) Or (Key = 6) Or (SearchText = '') Then Begin
      Q := SearchText;
      If InputString('Rechercher', 'Texte a chercher dans la sortie :', Q) Then SearchText := Q;
      HideMouseCursor;
      SetFillStyle(SolidFill, LightGray);
      Bar(X1, Y1, X2, Y2);
      SetFillStyle(SolidFill, Blue);
      Bar(X1, Y1, X2, Y1 + 13);
      SetColor(White);
      OutTextXY(6, 3, 'Sortie ORCA : ' + OutFileName + '   (' + IStr(NOut) + ' lignes)');
      SetColor(Black);
      OutTextXY(6, Y2 - 11, 'Fleches/PgUp/PgDn/Debut/Fin F=chercher N=suivant S=enregistrer Echap=fermer');
      ShowMouseCursor;
      DlgBegin;
     End;
     If SearchText <> '' Then Begin
      Found := -1;
      For I := Hi + 1 To NOut - 1 Do
       If ContainsText(OutLines[I]^, SearchText) Then Begin Found := I; Break; End;
      If (Found < 0) And (Hi >= 0) Then
       For I := 0 To Hi Do
        If ContainsText(OutLines[I]^, SearchText) Then Begin Found := I; Break; End;
      If Found >= 0 Then Begin
       Hi := Found;
       If (Hi < Top) Or (Hi >= Top + MaxLn) Then Begin Top := Hi - 3; If Top < 0 Then Top := 0; End;
      End
      Else ShowMessage('Recherche', 'Texte introuvable : ' + SearchText);
      Redraw := True;
     End;
    End;
    Ord('s'), Ord('S'): Begin
     FN := ChangeExt(InpName, '.out');
     If InpName = '' Then FN := 'ORCA.OUT';
     If FileSaveDialog('Enregistrer la sortie', FN) Then Begin
      If SaveOutputFile(FN) Then ShowMessage('Sortie', 'Fichier enregistre : ' + FN)
      Else ShowMessage('Erreur', 'Impossible d''ecrire ' + FN);
     End;
     HideMouseCursor;
     SetFillStyle(SolidFill, LightGray);
     Bar(X1, Y1, X2, Y2);
     SetFillStyle(SolidFill, Blue);
     Bar(X1, Y1, X2, Y1 + 13);
     SetColor(White);
     OutTextXY(6, 3, 'Sortie ORCA : ' + OutFileName + '   (' + IStr(NOut) + ' lignes)');
     SetColor(Black);
     OutTextXY(6, Y2 - 11, 'Fleches/PgUp/PgDn/Debut/Fin F=chercher N=suivant S=enregistrer Echap=fermer');
     ShowMouseCursor;
     DlgBegin;
     Redraw := True;
    End;
    27, 13: Done := True;
   End;
  End;
  Delay(10);
 End;
 DlgEnd;
End;

{ --------------------------------------------------------------------- }
{ Resumes textuels                                                       }
{ --------------------------------------------------------------------- }

Function NeedResult: Boolean;
Begin
 If Not HaveResult Or Not CalcDone Then Begin
  ShowMessage('Aucun resultat', 'Executez d''abord un calcul (F5).');
  NeedResult := False;
 End
 Else NeedResult := True;
End;

Procedure ShowEnergySummary;
Var HL: Real;
Begin
 If Not NeedResult Then Exit;
 TClear;
 TAdd('Molecule        : ' + MolFormula(1) + '   (' + IStr(ResN) + ' atomes)');
 TAdd('Methode / base  : ' + ResMethod);
 TAdd('Charge / mult.  : ' + IStr(ResCharge) + ' / ' + IStr(ResMult));
 TAdd('Electrons       : ' + IStr(NElec) + '   Fonctions de base : ' + IStr(NBF));
 TAdd('');
 TAdd('Repulsion nucleaire  : ' + FmtR(Enuc, 18, 10) + ' Eh');
 TAdd('Energie mono-elec.   : ' + FmtR(EOne, 18, 10) + ' Eh');
 TAdd('Coulomb              : ' + FmtR(ECoul, 18, 10) + ' Eh');
 TAdd('Echange exact        : ' + FmtR(EExch, 18, 10) + ' Eh');
 If ResMeth = M_DFT Then TAdd('Echange-correlation  : ' + FmtR(EXCdft, 18, 10) + ' Eh');
 TAdd('Energie SCF          : ' + FmtR(EScf, 18, 10) + ' Eh');
 If ResMeth = M_MP2 Then Begin
  TAdd('Correlation MP2      : ' + FmtR(EMP2, 18, 10) + ' Eh');
 End;
 TAdd('Energie totale       : ' + FmtR(ETot, 18, 10) + ' Eh');
 TAdd('                       ' + FmtR(ETot * HARTREE_EV, 18, 4) + ' eV');
 TAdd('                       ' + FmtR(ETot * HARTREE_KCAL, 18, 3) + ' kcal/mol');
 TAdd('                       ' + FmtR(ETot * HARTREE_KJ, 18, 3) + ' kJ/mol');
 TAdd('');
 TAdd('Iterations SCF       : ' + IStr(ScfIter));
 If (NOccA > 0) And (NOccA < NBF) And (NOccA = NOccB) Then Begin
  TAdd('HOMO                 : ' + FmtR(EpsA[NOccA], 12, 5) + ' Eh  (' + RStr(EpsA[NOccA] * HARTREE_EV, 3) + ' eV)');
  TAdd('LUMO                 : ' + FmtR(EpsA[NOccA + 1], 12, 5) + ' Eh  (' + RStr(EpsA[NOccA + 1] * HARTREE_EV, 3) + ' eV)');
  HL := (EpsA[NOccA + 1] - EpsA[NOccA]) * HARTREE_EV;
  TAdd('Ecart HOMO-LUMO      : ' + RStr(HL, 3) + ' eV');
 End;
 If NOccA <> NOccB Then TAdd('<S**2>               : ' + RStr(S2Exp, 5));
 TAdd('Moment dipolaire     : ' + RStr(DipTot * BOHR_DEBYE, 4) + ' Debye');
 If OptN > 1 Then TAdd('Cycles d''optimisation: ' + IStr(OptN - 1));
 If FreqOK Then TAdd('Frequences calculees : ' + IStr(NFreq - NProj) + ' modes de vibration');
 TAdd('Temps de calcul      : ' + TimeText(CalcMs));
 TView('Resume des energies');
End;

Procedure ShowPopulations;
Var A, B: Integer;
Begin
 If Not NeedResult Then Exit;
 TClear;
 TAdd('Atome      Mulliken     Loewdin     Spin(Mull.)');
 TAdd('------------------------------------------------');
 For A := 1 To ResN Do
  TAdd(FmtR(A - 1, 3, 0) + ' ' + PadR(ElemStr(ResZ[A]), 4) + FmtR(QMull[A], 11, 5) + FmtR(QLoew[A], 12, 5) +
       FmtR(SpinMull[A], 12, 5));
 TAdd('');
 TAdd('Ordres de liaison de Mayer (> 0.1) :');
 For A := 1 To ResN - 1 Do
  For B := A + 1 To ResN Do
   If BondOrd[A, B] > 0.1 Then
    TAdd('  ' + IStr(A - 1) + ElemStr(ResZ[A]) + ' - ' + IStr(B - 1) + ElemStr(ResZ[B]) + '   ' + RStr(BondOrd[A, B], 4));
 TAdd('');
 TAdd('Moment dipolaire (u.a.) : ' + RStr(Dip[1], 5) + '  ' + RStr(Dip[2], 5) + '  ' + RStr(Dip[3], 5));
 TAdd('Module : ' + RStr(DipTot, 5) + ' u.a. = ' + RStr(DipTot * BOHR_DEBYE, 4) + ' Debye');
 TView('Charges, populations et dipole');
End;

Procedure ShowFrequencies;
Var K, Sel: Integer;
Begin
 If Not NeedResult Then Exit;
 If Not FreqOK Then Begin
  ShowMessage('Frequences', 'Aucun calcul de frequences.|Ajoutez le mot-cle Freq a la ligne "!".');
  Exit;
 End;
 TClear;
 For K := NProj + 1 To NFreq Do
  TAdd('Mode ' + FmtR(K - 1, 3, 0) + ' : ' + FmtR(FreqV[K], 10, 2) + ' cm-1   IR = ' + FmtR(FreqI[K], 8, 2) + ' km/mol');
 If NProj > 0 Then TAdd('(' + IStr(NProj) + ' modes de translation/rotation exclus)');
 Sel := ListDialog('Frequences : Entree = animer le mode', True);
 If (Sel >= 0) And (Sel < NFreq - NProj) Then Begin
  VibMode := NProj + 1 + Sel; VibPhase := 0; VibLast := NowMs;
 End;
End;

Procedure ShowThermo;
Begin
 If Not NeedResult Then Exit;
 If Not ThermoOK Then Begin
  ShowMessage('Thermochimie', 'Aucun calcul de frequences.|Ajoutez le mot-cle Freq a la ligne "!".');
  Exit;
 End;
 TClear;
 TAdd('Temperature        : ' + RStr(ThTemp, 2) + ' K     Pression : 1 atm');
 TAdd('Nombre de symetrie : ' + IStr(SymNum));
 TAdd('Const. rotation(cm-1): ' + RStr(RotConst[1], 4) + '  ' + RStr(RotConst[2], 4) + '  ' + RStr(RotConst[3], 4));
 TAdd('');
 TAdd('Energie electronique : ' + FmtR(ETot, 16, 8) + ' Eh');
 TAdd('Energie point zero   : ' + FmtR(ZPE, 16, 8) + ' Eh  (' + RStr(ZPE * HARTREE_KCAL, 2) + ' kcal/mol)');
 TAdd('Correction vibration : ' + FmtR(EVibT, 16, 8) + ' Eh');
 TAdd('Correction rotation  : ' + FmtR(ERotT, 16, 8) + ' Eh');
 TAdd('Correction translat. : ' + FmtR(ETransT, 16, 8) + ' Eh');
 TAdd('Energie thermique    : ' + FmtR(ETot + EThermal, 16, 8) + ' Eh');
 TAdd('Enthalpie H          : ' + FmtR(HTot, 16, 8) + ' Eh');
 TAdd('Terme entropique T*S : ' + FmtR(STot * ThTemp, 16, 8) + ' Eh');
 TAdd('Energie libre G      : ' + FmtR(GTot, 16, 8) + ' Eh');
 TAdd('G - E(el)            : ' + RStr((GTot - ETot) * HARTREE_KCAL, 3) + ' kcal/mol');
 TView('Thermochimie (oscillateur harmonique, rotateur rigide)');
End;

Procedure ShowExcited;
Var K: Integer; Ev, Nm: Real;
Begin
 If Not NeedResult Then Exit;
 If Not CISOK Then Begin
  If CISMsg <> '' Then ShowMessage('Etats excites', CISMsg)
  Else ShowMessage('Etats excites', 'Aucun calcul CIS.|Ajoutez le mot-cle CIS a la ligne "!".');
  Exit;
 End;
 TClear;
 TAdd('Etat   E(eV)    lambda(nm)     f        transition');
 TAdd('--------------------------------------------------------');
 For K := 1 To CISN Do Begin
  Ev := CISE[K] * HARTREE_EV;
  If CISE[K] > 1E-9 Then Nm := 1239.84193 / Ev Else Nm := 0;
  TAdd(FmtR(K, 3, 0) + FmtR(Ev, 9, 3) + FmtR(Nm, 11, 1) + FmtR(CISF[K], 10, 4) + '   ' + Copy(CISDom[K], 1, 28));
 End;
 TView('Etats excites CIS');
End;

{ --------------------------------------------------------------------- }
{ Graphiques                                                             }
{ --------------------------------------------------------------------- }

Procedure ShowSCFChart;
Var S, I: Integer; Last, V: Real;
Begin
 If HistN < 2 Then Begin ShowMessage('Convergence SCF', 'Aucun historique SCF disponible.'); Exit; End;
 ChartReset('Convergence du SCF', 'Iteration', 'log10 |E - E(final)| (Eh)', False);
 ChFixY0 := False;
 S := ChartAddSeries('Erreur sur l''energie', False);
 Last := HistE[HistN - 1];
 For I := 0 To HistN - 2 Do Begin
  V := Abs(HistE[I] - Last);
  If V < 1E-13 Then V := 1E-13;
  ChartPoint(S, I, Ln(V) / Ln(10));
 End;
 ChartDialog;
End;

Procedure ShowOptChart;
Var S, I: Integer;
Begin
 If OptN < 2 Then Begin ShowMessage('Optimisation', 'Aucune optimisation de geometrie.'); Exit; End;
 ChartReset('Optimisation de geometrie', 'Cycle', 'Energie (Eh)', False);
 ChFixY0 := False;
 S := ChartAddSeries('E (Eh)', False);
 For I := 0 To OptN - 1 Do ChartPoint(S, I, OptE[I]);
 S := ChartAddSeries('points', True);
 For I := 0 To OptN - 1 Do ChartPoint(S, I, OptE[I]);
 ChartDialog;
End;

Procedure ShowScanChart;
Var S, I: Integer; Emin: Real;
Begin
 If ScanN < 2 Then Begin ShowMessage('Balayage', 'Aucun balayage de surface.|Utilisez %geom Scan ... end end.'); Exit; End;
 Emin := ScanE[1];
 For I := 2 To ScanN Do If ScanE[I] < Emin Then Emin := ScanE[I];
 If ScanKind = 1 Then
  ChartReset('Surface d''energie potentielle', 'Distance (Angstroem)', 'E relative (kcal/mol)', False)
 Else ChartReset('Surface d''energie potentielle', 'Angle (degres)', 'E relative (kcal/mol)', False);
 ChFixY0 := True;
 S := ChartAddSeries('Energie relaxee', False);
 For I := 1 To ScanN Do ChartPoint(S, ScanX[I], (ScanE[I] - Emin) * HARTREE_KCAL);
 S := ChartAddSeries('points calcules', True);
 For I := 1 To ScanN Do ChartPoint(S, ScanX[I], (ScanE[I] - Emin) * HARTREE_KCAL);
 ChartDialog;
End;

Procedure ShowIRChart;
Var S, I, K: Integer; X, Y, W: Real;
Begin
 If Not NeedResult Then Exit;
 If Not FreqOK Then Begin ShowMessage('Spectre IR', 'Aucun calcul de frequences.|Ajoutez le mot-cle Freq a la ligne "!".'); Exit; End;
 ChartReset('Spectre infrarouge (HF/DFT harmonique)', 'Nombre d''onde (cm-1)', 'Intensite (km/mol)', False);
 ChFixY0 := True;
 S := ChartAddSeries('Spectre (largeur 30 cm-1)', False);
 W := 15;
 For I := 0 To 150 Do Begin
  X := I * 30;
  Y := 0;
  For K := NProj + 1 To NFreq Do
   If FreqV[K] > 0 Then Y := Y + FreqI[K] * W * W / (Sqr(X - FreqV[K]) + W * W);
  ChartPoint(S, X, Y);
 End;
 S := ChartAddSeries('Bandes', True);
 For K := NProj + 1 To NFreq Do If FreqV[K] > 0 Then ChartPoint(S, FreqV[K], FreqI[K]);
 ChartDialog;
End;

Procedure ShowUVChart;
Var S, I, K: Integer; X, Y, Sg, Emax: Real;
Begin
 If Not NeedResult Then Exit;
 If Not CISOK Then Begin ShowMessage('Spectre UV', 'Aucun calcul CIS.|Ajoutez le mot-cle CIS a la ligne "!".'); Exit; End;
 ChartReset('Spectre UV-visible (CIS)', 'Energie (eV)', 'Force d''oscillateur (elargie)', False);
 ChFixY0 := True;
 Emax := CISE[CISN] * HARTREE_EV + 2;
 Sg := 0.25;
 S := ChartAddSeries('Spectre (gaussiennes 0.25 eV)', False);
 For I := 0 To 150 Do Begin
  X := Emax * I / 150;
  Y := 0;
  For K := 1 To CISN Do Y := Y + CISF[K] * Exp(-Sqr(X - CISE[K] * HARTREE_EV) / (2 * Sg * Sg));
  ChartPoint(S, X, Y);
 End;
 S := ChartAddSeries('Etats', True);
 For K := 1 To CISN Do ChartPoint(S, CISE[K] * HARTREE_EV, CISF[K]);
 ChartDialog;
End;

{ --------------------------------------------------------------------- }
{ Orbitales moleculaires : diagramme et coupes                           }
{ --------------------------------------------------------------------- }

Procedure OrbitalCut(Orb, Spin: Integer);
Const CELL = 5; NC = 70;
Var
 X1, Y1, X2, Y2, PX1, PY1, I, J, K, A, Plane, Key, MX, MY, Q, NL: Integer; MaxA: Real;
 Click, Done, Redraw, Dens: Boolean;
 Cm, E1, E2, E3, U, V, P: TVec3; Cov: TMat; W: TVecN; Ev: TMat;
 Ext, Offs, Psi, Rho, Su, Sv, Dd: Real;
 Idx: TIdxArr; Phi, GX, GY, GZ: TAOVals; Cm2: PMat;
 Col: Integer; Lbl: String; Mx2: Array[1..3] Of TVec3;
 Rmax, Dist2: Real; CurOrb: Integer;
 Function Level(V: Real): Integer;
 Begin
  V := Abs(V);
  If V > 0.15 Then Level := 3 Else If V > 0.06 Then Level := 2 Else If V > 0.02 Then Level := 1 Else Level := 0;
 End;
Begin
 UseResultGeom;
 ComputeBFCut;
 { axes principaux du nuage d'atomes }
 For J := 1 To 3 Do Cm[J] := 0;
 For A := 1 To NAt Do For J := 1 To 3 Do Cm[J] := Cm[J] + AtX[A, J] / NAt;
 FillChar(Cov, SizeOf(Cov), 0);
 For A := 1 To NAt Do
  For I := 1 To 3 Do
   For J := 1 To 3 Do Cov[I, J] := Cov[I, J] + (AtX[A, I] - Cm[I]) * (AtX[A, J] - Cm[J]);
 { le plan XY du molecule : vecteurs propres de variance decroissante }
 Jacobi(3, MAX_BF, PRealArr(@Cov), PRealArr(@Ev), PRealArr(@W));
 For J := 1 To 3 Do Begin
  E1[J] := Ev[J, 3]; E2[J] := Ev[J, 2]; E3[J] := Ev[J, 1];
 End;
 Rmax := 1;
 For A := 1 To NAt Do Begin
  Dd := Sqrt(Sqr(AtX[A,1] - Cm[1]) + Sqr(AtX[A,2] - Cm[2]) + Sqr(AtX[A,3] - Cm[3]));
  If Dd > Rmax Then Rmax := Dd;
 End;
 Ext := Rmax + 3.0;
 Plane := 1; Offs := 0; Dens := False; CurOrb := Orb;
 If Spin = 2 Then Cm2 := MCb Else Cm2 := MCa;
 { orbitale nulle dans le plan moleculaire (pi) : decale le plan }
 MaxA := 0;
 For I := 0 To 12 Do
  For J := 0 To 12 Do Begin
   Su := (I - 6) / 6 * Ext; Sv := (J - 6) / 6 * Ext;
   For K := 1 To 3 Do Mx2[1][K] := Cm[K] + Su * E1[K] + Sv * E2[K];
   EvalAOs(Mx2[1][1], Mx2[1][2], Mx2[1][3], NL, Idx, Phi, GX, GY, GZ);
   Psi := 0;
   For K := 1 To NL Do Psi := Psi + Cm2^[Idx[K - 1], CurOrb] * Phi[K];
   If Abs(Psi) > MaxA Then MaxA := Abs(Psi);
  End;
 If MaxA < 0.01 Then Offs := 0.8;
 X1 := 40; Y1 := 16; X2 := ScreenWidth - 40; Y2 := ScreenHeight - 16;
 PX1 := X1 + 14; PY1 := Y1 + 30;
 Done := False; Redraw := True;
 DlgBegin;
 While Not Done Do Begin
  If Redraw Then Begin
   HideMouseCursor;
   If Dens Then Lbl := 'Densite electronique totale'
   Else Lbl := 'Orbitale ' + IStr(CurOrb - 1);
   DrawDialogBox(X1, Y1, X2, Y2, Lbl + ' (' + ResMethod + ')');
   SetTextStyle(DefaultFont, HorizDir, 1);
   Case Plane Of
    1: Begin U := E1; V := E2; P := E3; End;
    2: Begin U := E1; V := E3; P := E2; End;
   Else Begin U := E2; V := E3; P := E1; End;
   End;
   For I := 0 To NC - 1 Do
    For J := 0 To NC - 1 Do Begin
     Su := (I - NC / 2 + 0.5) / NC * 2 * Ext;
     Sv := (NC / 2 - J - 0.5) / NC * 2 * Ext;
     For K := 1 To 3 Do Mx2[1][K] := Cm[K] + Su * U[K] + Sv * V[K] + Offs * P[K];
     EvalAOs(Mx2[1][1], Mx2[1][2], Mx2[1][3], NL, Idx, Phi, GX, GY, GZ);
     If Dens Then Begin
      Rho := 0;
      For Q := 1 To NOccA Do Begin
       Psi := 0;
       For K := 1 To NL Do Psi := Psi + MCa^[Idx[K - 1], Q] * Phi[K];
       Rho := Rho + Psi * Psi;
      End;
      For Q := 1 To NOccB Do Begin
       Psi := 0;
       For K := 1 To NL Do Psi := Psi + MCb^[Idx[K - 1], Q] * Phi[K];
       Rho := Rho + Psi * Psi;
      End;
      If Rho > 1.0 Then Col := Red Else If Rho > 0.3 Then Col := LightRed
      Else If Rho > 0.1 Then Col := Yellow Else If Rho > 0.03 Then Col := LightGreen
      Else If Rho > 0.008 Then Col := LightCyan Else If Rho > 0.002 Then Col := LightGray Else Col := White;
     End
     Else Begin
      Psi := 0;
      For K := 1 To NL Do Psi := Psi + Cm2^[Idx[K - 1], CurOrb] * Phi[K];
      Case Level(Psi) Of
       0: Col := White;
       1: If Psi > 0 Then Col := LightMagenta Else Col := LightCyan;
       2: If Psi > 0 Then Col := LightRed Else Col := LightBlue;
      Else If Psi > 0 Then Col := Red Else Col := Blue;
      End;
     End;
     SetFillStyle(SolidFill, Col);
     Bar(PX1 + I * CELL, PY1 + J * CELL, PX1 + I * CELL + CELL - 1, PY1 + J * CELL + CELL - 1);
    End;
   SetColor(Black);
   Rectangle(PX1 - 1, PY1 - 1, PX1 + NC * CELL, PY1 + NC * CELL);
   { atomes }
   For A := 1 To NAt Do Begin
    Su := 0; Sv := 0; Dd := 0;
    For K := 1 To 3 Do Begin
     Su := Su + (AtX[A, K] - Cm[K]) * U[K];
     Sv := Sv + (AtX[A, K] - Cm[K]) * V[K];
     Dd := Dd + (AtX[A, K] - Cm[K]) * P[K];
    End;
    If Abs(Dd - Offs) < 0.8 Then Begin
     I := PX1 + Round((Su / (2 * Ext) + 0.5) * NC * CELL);
     J := PY1 + Round((0.5 - Sv / (2 * Ext)) * NC * CELL);
     If (I > PX1) And (I < PX1 + NC * CELL) And (J > PY1) And (J < PY1 + NC * CELL) Then Begin
      SetColor(Black);
      SetFillStyle(SolidFill, ElColor[AtZ[A]]);
      FillEllipse(I, J, 4, 4);
      Circle(I, J, 4);
      OutTextXY(I + 6, J - 4, ElemStr(AtZ[A]) + IStr(A - 1));
     End;
    End;
   End;
   { informations }
   SetColor(Black);
   OutTextXY(PX1 + NC * CELL + 14, PY1, 'Plan : ' + IStr(Plane) + '/3  (P)');
   OutTextXY(PX1 + NC * CELL + 14, PY1 + 14, 'Decalage : ' + RStr(Offs, 2) + ' bohr');
   OutTextXY(PX1 + NC * CELL + 14, PY1 + 28, '(Haut/Bas)');
   If Not Dens Then Begin
    OutTextXY(PX1 + NC * CELL + 14, PY1 + 48, 'Orbitale : ' + IStr(CurOrb - 1));
    OutTextXY(PX1 + NC * CELL + 14, PY1 + 62, '(Gauche/Droite)');
    If CurOrb <= NOccA Then OutTextXY(PX1 + NC * CELL + 14, PY1 + 76, 'occupee')
    Else OutTextXY(PX1 + NC * CELL + 14, PY1 + 76, 'virtuelle');
    OutTextXY(PX1 + NC * CELL + 14, PY1 + 90, 'E = ' + RStr(EpsA[CurOrb], 4) + ' Eh');
    OutTextXY(PX1 + NC * CELL + 14, PY1 + 108, 'Rouge : positive');
    OutTextXY(PX1 + NC * CELL + 14, PY1 + 120, 'Bleu  : negative');
   End;
   OutTextXY(PX1 + NC * CELL + 14, PY1 + 146, 'D : densite');
   OutTextXY(PX1 + NC * CELL + 14, PY1 + 158, 'Echap : fermer');
   ShowMouseCursor;
   Redraw := False;
  End;
  If DlgPoll(Key, MX, MY, Click) Then Begin
   Case Key Of
    27, 13: Done := True;
    Ord('p'), Ord('P'): Begin Inc(Plane); If Plane > 3 Then Plane := 1; Redraw := True; End;
    Ord('d'), Ord('D'): Begin Dens := Not Dens; Redraw := True; End;
    256 + 72: Begin Offs := Offs + 0.25; Redraw := True; End;
    256 + 80: Begin Offs := Offs - 0.25; Redraw := True; End;
    256 + 77: If CurOrb < NBF Then Begin Inc(CurOrb); Redraw := True; End;
    256 + 75: If CurOrb > 1 Then Begin Dec(CurOrb); Redraw := True; End;
   End;
  End;
  Delay(10);
 End;
 RestoreGeom;
 DlgEnd;
End;

Procedure OrbitalDialog;
Const BW = 90; BH = 22;
Var
 X1, Y1, X2, Y2, I, K, Sel, Key, MX, MY, BCl, BCu, BBY1, BBY2, YL, LY1, LY2, Top, MaxLn, Spin: Integer;
 Click, Done, Redraw: Boolean;
 Eo: TVecN; NO, LastLab: Integer; ELo, EHi, Sc, Y: Real; S: String; A: Integer;
 Cm2: PMat;
Begin
 If Not NeedResult Then Exit;
 If Not OrbsValid Then Begin ShowMessage('Orbitales', 'Orbitales indisponibles (calcul non converge).'); Exit; End;
 X1 := 30; Y1 := 14; X2 := ScreenWidth - 30; Y2 := ScreenHeight - 14;
 BBY2 := Y2 - 8; BBY1 := BBY2 - BH;
 BCl := X2 - BW - 12; BCu := BCl - BW - 10;
 LY1 := Y1 + 40; LY2 := BBY1 - 70;
 MaxLn := (LY2 - LY1 - 4) Div 10;
 Spin := 1; Sel := NOccA; If Sel < 1 Then Sel := 1;
 Top := 0; Done := False; Redraw := True;
 DlgBegin;
 While Not Done Do Begin
  If Redraw Then Begin
   HideMouseCursor;
   If Spin = 1 Then Begin Eo := EpsA; NO := NOccA; End Else Begin Eo := EpsB; NO := NOccB; End;
   If Sel >= Top + MaxLn Then Top := Sel - MaxLn + 1;
   If Sel < Top Then Top := Sel - 1;
   If Top < 0 Then Top := 0;
   DrawDialogBox(X1, Y1, X2, Y2, 'Orbitales moleculaires - ' + ResMethod);
   SetTextStyle(DefaultFont, HorizDir, 1);
   SetColor(Black);
   If NOccA <> NOccB Then Begin
    If Spin = 1 Then OutTextXY(X1 + 10, Y1 + 22, 'Spin alpha (Tab : changer)') Else OutTextXY(X1 + 10, Y1 + 22, 'Spin beta (Tab : changer)');
   End
   Else OutTextXY(X1 + 10, Y1 + 22, 'NO OCC   E(Eh)     E(eV)');
   SetFillStyle(SolidFill, White);
   Bar(X1 + 8, LY1, X1 + 260, LY2);
   Rectangle(X1 + 8, LY1, X1 + 260, LY2);
   For I := 0 To MaxLn - 1 Do Begin
    K := Top + I + 1;
    If K > NBF Then Break;
    If K = Sel Then Begin
     SetFillStyle(SolidFill, Blue); Bar(X1 + 9, LY1 + 2 + I * 10, X1 + 259, LY1 + 11 + I * 10); SetColor(White);
    End
    Else SetColor(Black);
    If K <= NO Then S := ' 2' Else S := ' 0';
    If NOccA <> NOccB Then If K <= NO Then S := ' 1';
    S := FmtR(K - 1, 3, 0) + S + FmtR(Eo[K], 11, 4) + FmtR(Eo[K] * HARTREE_EV, 9, 2);
    If (K = NO) Then S := S + ' HOMO' Else If K = NO + 1 Then S := S + ' LUMO';
    OutTextXY(X1 + 12, LY1 + 3 + I * 10, Copy(S, 1, 31));
   End;
   { diagramme des niveaux }
   ELo := Eo[1]; EHi := Eo[NBF];
   If NO >= 1 Then ELo := Eo[NO] - 2.2;
   If NO < NBF Then EHi := Eo[NO + 1] + 2.2 Else EHi := Eo[NO] + 1;
   If Eo[1] > ELo Then ELo := Eo[1];
   If EHi <= ELo Then EHi := ELo + 1;
   Sc := (LY2 - LY1 - 20) / (EHi - ELo);
   SetFillStyle(SolidFill, White);
   Bar(X1 + 280, LY1, X2 - 12, LY2);
   SetColor(Black);
   Rectangle(X1 + 280, LY1, X2 - 12, LY2);
   LastLab := -100;
   For K := 1 To NBF Do
    If (Eo[K] >= ELo) And (Eo[K] <= EHi) Then Begin
     YL := LY2 - 10 - Round((Eo[K] - ELo) * Sc);
     If K = Sel Then SetColor(Red) Else If K <= NO Then SetColor(Blue) Else SetColor(DarkGray);
     Line(X1 + 330, YL, X1 + 400, YL);
     Line(X1 + 330, YL - 1, X1 + 400, YL - 1);
     If K <= NO Then Begin
      SetColor(Blue);
      If NOccA = NOccB Then Begin
       Line(X1 + 355, YL - 8, X1 + 355, YL - 2); Line(X1 + 352, YL - 6, X1 + 355, YL - 8); Line(X1 + 358, YL - 6, X1 + 355, YL - 8);
       Line(X1 + 372, YL - 2, X1 + 372, YL - 8); Line(X1 + 369, YL - 4, X1 + 372, YL - 2); Line(X1 + 375, YL - 4, X1 + 372, YL - 2);
      End
      Else Begin
       If Spin = 1 Then Begin Line(X1 + 365, YL - 8, X1 + 365, YL - 2); Line(X1 + 362, YL - 6, X1 + 365, YL - 8); Line(X1 + 368, YL - 6, X1 + 365, YL - 8); End
       Else Begin Line(X1 + 365, YL - 2, X1 + 365, YL - 8); Line(X1 + 362, YL - 4, X1 + 365, YL - 2); Line(X1 + 368, YL - 4, X1 + 365, YL - 2); End;
      End;
     End;
     SetColor(Black);
     If (Abs(LastLab - YL) >= 9) Or (K = Sel) Or (K = NO) Or (K = NO + 1) Then Begin
      OutTextXY(X1 + 408, YL - 4, IStr(K - 1) + '  ' + RStr(Eo[K], 3));
      LastLab := YL;
     End;
     If K = NO Then OutTextXY(X1 + 286, YL - 4, 'HOMO');
     If K = NO + 1 Then OutTextXY(X1 + 286, YL - 4, 'LUMO');
    End;
   { composition de l'orbitale selectionnee }
   SetColor(Black);
   S := 'Orbitale ' + IStr(Sel - 1) + ' : ';
   If Spin = 1 Then Begin
    For A := 1 To ResN Do
     If OrbMullPct[Sel, A] > 0.05 Then S := S + ElemStr(ResZ[A]) + IStr(A - 1) + ' ' + RStr(100 * OrbMullPct[Sel, A], 0) + '%  ';
   End
   Else S := S + '(composition : spin alpha uniquement)';
   OutTextXY(X1 + 10, LY2 + 8, Copy(S, 1, (X2 - X1 - 20) Div 8));
   OutTextXY(X1 + 10, LY2 + 22, 'Haut/Bas : choisir   Entree : coupe de l''orbitale   Echap : fermer');
   DrawButton(BCu, BBY1, BCu + BW, BBY2, 'Coupe', False);
   DrawButton(BCl, BBY1, BCl + BW, BBY2, 'Fermer', True);
   ShowMouseCursor;
   Redraw := False;
  End;
  If DlgPoll(Key, MX, MY, Click) Then Begin
   Case Key Of
    27: Done := True;
    256 + 72: If Sel > 1 Then Begin Dec(Sel); Redraw := True; End;
    256 + 80: If Sel < NBF Then Begin Inc(Sel); Redraw := True; End;
    256 + 73: Begin Dec(Sel, 5); If Sel < 1 Then Sel := 1; Redraw := True; End;
    256 + 81: Begin Inc(Sel, 5); If Sel > NBF Then Sel := NBF; Redraw := True; End;
    9: If NOccA <> NOccB Then Begin Spin := 3 - Spin; Redraw := True; End;
    13: Begin OrbitalCut(Sel, Spin); Redraw := True; End;
   End;
   If Click Then Begin
    If PtInRect(MX, MY, BCl, BBY1, BCl + BW, BBY2) Then Done := True
    Else If PtInRect(MX, MY, BCu, BBY1, BCu + BW, BBY2) Then Begin OrbitalCut(Sel, Spin); Redraw := True; End
    Else If PtInRect(MX, MY, X1 + 8, LY1, X1 + 260, LY2) Then Begin
     K := Top + (MY - LY1 - 2) Div 10 + 1;
     If (K >= 1) And (K <= NBF) Then Begin Sel := K; Redraw := True; End;
    End;
   End;
  End;
  Delay(10);
 End;
 DlgEnd;
End;

{ --------------------------------------------------------------------- }
{ Assistant de calcul et insertion de mots-cles                          }
{ --------------------------------------------------------------------- }

Function FindBangLine: Integer;
Var I: Integer; T: String;
Begin
 FindBangLine := -1;
 For I := 0 To NInp - 1 Do Begin
  T := TrimStr(InpLines[I]);
  If (Length(T) > 0) And (T[1] = '!') Then Begin FindBangLine := I; Exit; End;
 End;
End;

Function FindCoordLine: Integer;
Var I: Integer; S, T: String;
Begin
 FindCoordLine := -1;
 For I := 0 To NInp - 1 Do Begin
  S := TrimStr(InpLines[I]);
  If (Length(S) > 0) And (S[1] = '*') Then Begin
   Delete(S, 1, 1); T := UpStr(NextTok(S));
   If (T = 'XYZ') Or (T = 'INT') Or (T = 'XYZFILE') Or (T = 'GZMT') Then Begin FindCoordLine := I; Exit; End;
  End;
 End;
End;

Procedure AddKeywords(Const Kw: String);
Var B: Integer;
Begin
 B := FindBangLine;
 If B < 0 Then Begin EdInsertLine(0, '! ' + Kw); CurLn := 0; End
 Else InpLines[B] := Copy(InpLines[B] + ' ' + Kw, 1, 100);
 MarkModified;
End;

Procedure AssistantDialog;
Var
 Err: String; B, C, Mi, Ji, Si, I: Integer; Line, Rest, T, U: String; Kw: String;
Begin
 PreviewInput;
 { valeurs actuelles }
 Case Meth Of
  M_HF: Mi := 0;
  M_MP2: Mi := 1;
 Else Mi := 1 + Func;
 End;
 If (Meth = M_DFT) And (Func = F_NONE) Then Mi := 0;
 Case JobType Of
  J_OPT, J_SCAN: Ji := 1;
  J_FREQ: Ji := 2;
  J_OPTFREQ: Ji := 3;
 Else Ji := 0;
 End;
 If ScfTolE <= 1E-9 Then Si := 2 Else If ScfTolE <= 1E-8 Then Si := 1 Else Si := 0;
 PReset;
 PAddH('Methode et base');
 PAddE('Methode', Mi, 'HF|MP2|LDA|BLYP|B3LYP|PBE|PBE0');
 PAddE('Base gaussienne', BasisSel, 'STO-3G|3-21G|6-31G');
 PAddH('Type de calcul');
 PAddE('Tache', Ji, 'Energie (SP)|Optimisation (Opt)|Frequences (Freq)|Opt + Freq');
 PAddB('Etats excites CIS', DoCIS);
 PAddH('Etat electronique');
 PAddI('Charge totale', Charge, -5, 5);
 PAddI('Multiplicite (2S+1)', Mult, 1, 7);
 PAddH('Convergence');
 PAddE('Convergence SCF', Si, 'NormalSCF|TightSCF|VeryTightSCF');
 PAddE('Grille DFT', GridScale - 1, 'Normale|Fine (lente)');
 If Not ParamDialog('Assistant de calcul ORCA', 'Met a jour la ligne "!" et la ligne "* xyz" du fichier.') Then Exit;
 Kw := '! ';
 Case PInt(1) Of
  0: Kw := Kw + 'HF';
  1: Kw := Kw + 'MP2';
  2: Kw := Kw + 'LDA';
  3: Kw := Kw + 'BLYP';
  4: Kw := Kw + 'B3LYP';
  5: Kw := Kw + 'PBE';
 Else Kw := Kw + 'PBE0';
 End;
 If (PInt(1) = 0) And (PInt(8) > 1) Then Kw := '! UHF';
 Kw := Kw + ' ' + BasisNames[PInt(2)];
 Case PInt(4) Of
  1: Kw := Kw + ' Opt';
  2: Kw := Kw + ' Freq';
  3: Kw := Kw + ' Opt Freq';
 End;
 Case PInt(10) Of
  1: Kw := Kw + ' TightSCF';
  2: Kw := Kw + ' VeryTightSCF';
 End;
 If PBool(5) Then Kw := Kw + ' CIS';
 If PInt(11) = 1 Then Kw := Kw + ' Grid3';
 B := FindBangLine;
 If B >= 0 Then InpLines[B] := Copy(Kw, 1, 100)
 Else Begin EdInsertLine(0, Kw); End;
 C := FindCoordLine;
 If C >= 0 Then Begin
  Line := TrimStr(InpLines[C]); Delete(Line, 1, 1);
  U := NextTok(Line);
  If UpStr(U) = 'XYZFILE' Then Begin
   T := NextTok(Line); T := NextTok(Line);
   InpLines[C] := '* xyzfile ' + IStr(PInt(7)) + ' ' + IStr(PInt(8)) + ' ' + TrimStr(Line);
  End
  Else InpLines[C] := '* ' + LoStr(U) + ' ' + IStr(PInt(7)) + ' ' + IStr(PInt(8));
 End;
 MarkModified;
End;

Procedure InsertKeywordsDialog;
Const NK = 22;
 KwList : Array[0..NK-1] Of String[16] = (
  'HF', 'UHF', 'MP2', 'LDA', 'BLYP', 'B3LYP', 'PBE', 'PBE0', 'STO-3G', '3-21G', '6-31G',
  'Opt', 'TightOpt', 'Freq', 'CIS', 'TightSCF', 'VeryTightSCF', 'LooseSCF', 'NoFrozenCore',
  'PrintBasis', 'Grid3', 'NoPrintMOs');
 KwDesc : Array[0..NK-1] Of String[44] = (
  'Hartree-Fock (RHF)', 'Hartree-Fock non restreint', 'Perturbation Moller-Plesset du 2e ordre',
  'DFT : echange de Slater + correlation VWN5', 'DFT : Becke88 + LYP', 'DFT hybride B3LYP (20% HF)',
  'DFT : Perdew-Burke-Ernzerhof', 'DFT hybride PBE0 (25% HF)', 'Base minimale STO-3G (H a Ne)',
  'Base 3-21G (H, He, C a F)', 'Base 6-31G (H, He, C a F)', 'Optimisation de geometrie',
  'Optimisation avec seuils serres', 'Frequences harmoniques et thermochimie', 'Etats excites CIS',
  'Convergence SCF serree', 'Convergence SCF tres serree', 'Convergence SCF lache',
  'MP2 sans coeur gele', 'Imprimer les fonctions de base', 'Grille DFT fine', 'Ne pas imprimer les MOs');
Var I, Sel: Integer;
Begin
 TClear;
 For I := 0 To NK - 1 Do TAdd(PadR(KwList[I], 14) + ' ' + KwDesc[I]);
 Sel := ListDialog('Inserer un mot-cle dans la ligne "!"', True);
 If (Sel >= 0) And (Sel < NK) Then AddKeywords(KwList[Sel]);
End;

{ --------------------------------------------------------------------- }
{ Exemples tires des tutoriels                                           }
{ --------------------------------------------------------------------- }

Const
 NUM_EX = 12;
 ExTitle : Array[0..NUM_EX-1] Of String[40] = (
  'H2 : energie ponctuelle HF/STO-3G',
  'Eau : energie HF/6-31G + populations',
  'Eau : optimisation + frequences',
  'Eau : B3LYP/6-31G',
  'Eau : MP2/6-31G',
  'Radical OH : UHF/6-31G',
  'Ethylene : optimisation HF/STO-3G',
  'Eau : etats excites CIS',
  'H2 : balayage de la liaison',
  'Eau : matrice Z (* int)',
  'Ammoniac : opt + freq (STO-3G)',
  'Formaldehyde : PBE0/3-21G');

Procedure LoadExample(N: Integer);
Begin
 NInp := 0;
 Case N Of
  0: Begin
   EdAppend('# Molecule H2 - energie ponctuelle');
   EdAppend('! HF STO-3G TightSCF');
   EdAppend('* xyz 0 1');
   EdAppend('H   0.000000   0.000000   0.000000');
   EdAppend('H   0.000000   0.000000   0.740000');
   EdAppend('*');
  End;
  1: Begin
   EdAppend('# Eau : energie, orbitales, charges et dipole');
   EdAppend('! HF 6-31G TightSCF');
   EdAppend('* xyz 0 1');
   EdAppend('O   0.000000   0.000000   0.000000');
   EdAppend('H   0.757000   0.586000   0.000000');
   EdAppend('H  -0.757000   0.586000   0.000000');
   EdAppend('*');
  End;
  2: Begin
   EdAppend('# Eau : optimisation de geometrie puis frequences');
   EdAppend('! HF 6-31G Opt Freq');
   EdAppend('* xyz 0 1');
   EdAppend('O   0.000000   0.000000   0.000000');
   EdAppend('H   0.800000   0.600000   0.000000');
   EdAppend('H  -0.700000   0.650000   0.100000');
   EdAppend('*');
  End;
  3: Begin
   EdAppend('# Eau : DFT hybride B3LYP (grille de Becke)');
   EdAppend('! B3LYP 6-31G TightSCF');
   EdAppend('* xyz 0 1');
   EdAppend('O   0.000000   0.000000   0.000000');
   EdAppend('H   0.757000   0.586000   0.000000');
   EdAppend('H  -0.757000   0.586000   0.000000');
   EdAppend('*');
  End;
  4: Begin
   EdAppend('# Eau : energie de correlation MP2 (coeur gele)');
   EdAppend('! MP2 6-31G TightSCF');
   EdAppend('* xyz 0 1');
   EdAppend('O   0.000000   0.000000   0.000000');
   EdAppend('H   0.757000   0.586000   0.000000');
   EdAppend('H  -0.757000   0.586000   0.000000');
   EdAppend('*');
  End;
  5: Begin
   EdAppend('# Radical hydroxyle : couche ouverte (doublet)');
   EdAppend('! UHF 6-31G TightSCF');
   EdAppend('* xyz 0 2');
   EdAppend('O   0.000000   0.000000   0.000000');
   EdAppend('H   0.000000   0.000000   0.970000');
   EdAppend('*');
  End;
  6: Begin
   EdAppend('# Ethylene : optimisation de geometrie');
   EdAppend('! HF STO-3G Opt');
   EdAppend('* xyz 0 1');
   EdAppend('C   0.000000   0.000000   0.670000');
   EdAppend('C   0.000000   0.000000  -0.670000');
   EdAppend('H   0.000000   0.930000   1.240000');
   EdAppend('H   0.000000  -0.930000   1.240000');
   EdAppend('H   0.000000   0.930000  -1.240000');
   EdAppend('H   0.000000  -0.930000  -1.240000');
   EdAppend('*');
  End;
  7: Begin
   EdAppend('# Eau : etats excites singulets CIS et spectre UV');
   EdAppend('! HF 6-31G TightSCF CIS');
   EdAppend('%tddft');
   EdAppend('  nroots 6');
   EdAppend('end');
   EdAppend('* xyz 0 1');
   EdAppend('O   0.000000   0.000000   0.000000');
   EdAppend('H   0.757000   0.586000   0.000000');
   EdAppend('H  -0.757000   0.586000   0.000000');
   EdAppend('*');
  End;
  8: Begin
   EdAppend('# H2 : balayage relaxe de la distance H-H');
   EdAppend('! HF 6-31G Opt');
   EdAppend('%geom Scan');
   EdAppend('  B 0 1 = 0.50, 2.00, 12');
   EdAppend('  end');
   EdAppend('end');
   EdAppend('* xyz 0 1');
   EdAppend('H   0.000000   0.000000   0.000000');
   EdAppend('H   0.000000   0.000000   0.740000');
   EdAppend('*');
  End;
  9: Begin
   EdAppend('# Eau en coordonnees internes (matrice Z)');
   EdAppend('! HF 3-21G Opt');
   EdAppend('* int 0 1');
   EdAppend('O  0 0 0  0.000  0.00  0.00');
   EdAppend('H  1 0 0  0.960  0.00  0.00');
   EdAppend('H  1 2 0  0.960  104.50  0.00');
   EdAppend('*');
  End;
  10: Begin
   EdAppend('# Ammoniac : geometrie et vibrations');
   EdAppend('! HF STO-3G TightSCF Opt Freq');
   EdAppend('* xyz 0 1');
   EdAppend('N   0.000000   0.000000   0.120000');
   EdAppend('H   0.000000   0.940000  -0.280000');
   EdAppend('H   0.814000  -0.470000  -0.280000');
   EdAppend('H  -0.814000  -0.470000  -0.280000');
   EdAppend('*');
  End;
  11: Begin
   EdAppend('# Formaldehyde : DFT hybride PBE0');
   EdAppend('! PBE0 3-21G TightSCF');
   EdAppend('* xyz 0 1');
   EdAppend('C   0.000000   0.000000  -0.600000');
   EdAppend('O   0.000000   0.000000   0.600000');
   EdAppend('H   0.000000   0.940000  -1.180000');
   EdAppend('H   0.000000  -0.940000  -1.180000');
   EdAppend('*');
  End;
 End;
 EnsureOneLine;
 CurLn := 0; CurCol := 0; TopLn := 0; LeftCol := 0;
 InpName := '';
 InpModified := False;
 ResStale := True;
 PreviewInput;
 ResetView;
 NeedRedraw := True;
End;

Procedure ExamplesDialog;
Var I, Sel: Integer;
Begin
 TClear;
 For I := 0 To NUM_EX - 1 Do TAdd(ExTitle[I]);
 Sel := ListDialog('Exemples de fichiers d''entree', True);
 If (Sel >= 0) And (Sel < NUM_EX) Then LoadExample(Sel);
End;

{ --------------------------------------------------------------------- }
{ Outils                                                                 }
{ --------------------------------------------------------------------- }

Procedure ShowMolInfo;
Var A, D, N: Integer; M, Mt: Real; Cm: TVec3; Zs: Array[1..MAX_AT] Of Integer; Xs: Array[1..MAX_AT, 1..3] Of Real; Src: Integer;
Begin
 If HaveResult And Not ResStale Then Begin
  N := ResN; Src := 1;
  For A := 1 To N Do Begin Zs[A] := ResZ[A]; For D := 1 To 3 Do Xs[A, D] := ResX[A, D]; End;
 End
 Else If PrevOK Then Begin
  N := PrevN; Src := 0;
  For A := 1 To N Do Begin Zs[A] := PrevZ[A]; For D := 1 To 3 Do Xs[A, D] := PrevX[A, D]; End;
 End
 Else Begin ShowMessage('Molecule', 'La geometrie du fichier est invalide.'); Exit; End;
 Mt := 0; For D := 1 To 3 Do Cm[D] := 0;
 For A := 1 To N Do Begin
  M := ElMass[Zs[A]]; Mt := Mt + M;
  For D := 1 To 3 Do Cm[D] := Cm[D] + M * Xs[A, D];
 End;
 For D := 1 To 3 Do Cm[D] := Cm[D] / Mt;
 TClear;
 TAdd('Formule brute      : ' + MolFormula(Src));
 TAdd('Nombre d''atomes    : ' + IStr(N));
 TAdd('Masse moleculaire  : ' + RStr(Mt, 4) + ' g/mol');
 TAdd('Centre de masse (A): ' + RStr(Cm[1] / ANG2BOHR, 4) + '  ' + RStr(Cm[2] / ANG2BOHR, 4) + '  ' + RStr(Cm[3] / ANG2BOHR, 4));
 TAdd('');
 TAdd('Atomes (angstroems) :');
 For A := 1 To N Do
  TAdd(FmtR(A - 1, 3, 0) + ' ' + PadR(ElemStr(Zs[A]), 3) + FmtR(Xs[A,1] / ANG2BOHR, 11, 5) +
       FmtR(Xs[A,2] / ANG2BOHR, 11, 5) + FmtR(Xs[A,3] / ANG2BOHR, 11, 5));
 TView('Informations moleculaires');
End;

Procedure MeasureDialog;
Var S, T: String; Ids: Array[1..4] Of Integer; N, I, V: Integer; C: TConstr; R: Real; Res: Boolean;
Begin
 Res := HaveResult And Not ResStale;
 If Not Res And Not PrevOK Then Begin
  ShowMessage('Mesures', 'La geometrie du fichier est invalide.'); Exit;
 End;
 S := '0 1';
 If Not InputString('Mesures geometriques', 'Numeros d''atomes (2=distance 3=angle 4=dihedre) :', S) Then Exit;
 N := 0;
 Repeat
  T := NextTok(S);
  If T <> '' Then Begin
   If ParseInt(T, V) And (N < 4) Then Begin Inc(N); Ids[N] := V + 1; End;
  End;
 Until T = '';
 If N < 2 Then Begin ShowMessage('Mesures', 'Entrez 2, 3 ou 4 numeros d''atomes.'); Exit; End;
 If Res Then UseResultGeom Else UsePrevGeom;
 For I := 1 To N Do
  If (Ids[I] < 1) Or (Ids[I] > NAt) Then Begin
   RestoreGeom; ShowMessage('Mesures', 'Numero d''atome invalide.'); Exit;
  End;
 C.Kind := N - 1;
 For I := 1 To N Do C.A[I] := Ids[I];
 R := ConstrValue(C);
 RestoreGeom;
 If N = 2 Then ShowMessage('Distance', 'Distance = ' + RStr(R / ANG2BOHR, 4) + ' angstroems')
 Else If N = 3 Then ShowMessage('Angle', 'Angle = ' + RStr(R * 180 / Pi, 3) + ' degres')
 Else ShowMessage('Dihedre', 'Angle dihedre = ' + RStr(R * 180 / Pi, 3) + ' degres');
End;

Procedure ConvertEnergy;
Var S: String; V: Real;
Begin
 S := '-1.0';
 If Not InputString('Convertisseur', 'Energie en hartree (Eh) :', S) Then Exit;
 If Not ParseReal(S, V) Then Begin ShowMessage('Convertisseur', 'Nombre invalide.'); Exit; End;
 ShowMessage('Conversion d''energie', RStr(V, 8) + ' Eh|' + RStr(V * HARTREE_EV, 5) + ' eV|' +
             RStr(V * HARTREE_KCAL, 3) + ' kcal/mol|' + RStr(V * HARTREE_KJ, 3) + ' kJ/mol|' +
             RStr(V * HARTREE_CM, 1) + ' cm-1');
End;

Procedure ShowElements;
Var I: Integer;
Begin
 TClear;
 TAdd('Z  Sym  Nom           Masse     Rayon cov.(A)  Bases');
 TAdd('-----------------------------------------------------------');
 For I := 1 To NELEM Do Begin
  TAdd(FmtR(I, 2, 0) + '  ' + PadR(ElSym[I], 4) + PadR(ElName[I], 13) + FmtR(ElMass[I], 9, 3) + FmtR(ElRad[I], 10, 2) + '      ' +
       Copy('S', 1, Ord(GetShells(B_STO3G, I))) + Copy('3', 1, Ord(GetShells(B_321G, I))) + Copy('6', 1, Ord(GetShells(B_631G, I))));
 End;
 TAdd('');
 TAdd('Bases : S = STO-3G, 3 = 3-21G, 6 = 6-31G');
 TView('Elements chimiques supportes');
End;

{ --------------------------------------------------------------------- }
{ Aide                                                                   }
{ --------------------------------------------------------------------- }

Procedure ShowHelp;
Begin
 TClear;
 TAdd('ORCA Clone - chimie quantique ab initio (QChimie)');
 TAdd('');
 TAdd('Fichier d''entree (.INP) a gauche, molecule en haut a droite,');
 TAdd('console ORCA en bas a droite.  F5 = executer le calcul.');
 TAdd('');
 TAdd('Format du fichier :');
 TAdd('  ! HF 6-31G TightSCF Opt Freq        (mots-cles)');
 TAdd('  %scf MaxIter 150 end                (blocs)');
 TAdd('  * xyz 0 1  ...  *                   (charge et multiplicite)');
 TAdd('  * int 0 1  ...  *                   (matrice Z)');
 TAdd('  * xyzfile 0 1 fichier.xyz');
 TAdd('');
 TAdd('Methodes : HF UHF MP2 LDA BLYP B3LYP PBE PBE0, CIS');
 TAdd('Bases    : STO-3G (H-Ne)  3-21G  6-31G (H, He, C-F)');
 TAdd('Taches   : SP (defaut) Opt Freq Opt Freq, %geom Scan');
 TAdd('');
 TAdd('Touches : F1 aide   F2 assistant   F3 sortie   F4 IR');
 TAdd('          F5/Ctrl+R executer   F6 orbitales   F10 menus');
 TAdd('          Ctrl+N/O/S nouveau/ouvrir/enregistrer');
 TAdd('          Ctrl+Y supprimer la ligne   Ctrl+Q quitter');
 TAdd('Souris  : glisser dans la vue pour tourner la molecule');
 TAdd('          touches + et - pour le zoom, R pour reinitialiser');
 TAdd('Echap interrompt un calcul en cours.');
 TView('Aide - ORCA Clone');
End;

Procedure ShowKeywordHelp;
Begin
 TClear;
 TAdd('Mots-cles de la ligne "!" pris en charge :');
 TAdd('');
 TAdd('Methodes  HF RHF UHF MP2 LDA SVWN BLYP B3LYP PBE PBE0');
 TAdd('Bases     STO-3G  3-21G  6-31G');
 TAdd('          (def2-SVP, 6-31G*, cc-pVDZ : 6-31G utilisee)');
 TAdd('Taches    SP Opt TightOpt LooseOpt Freq NumFreq');
 TAdd('SCF       NormalSCF TightSCF VeryTightSCF LooseSCF');
 TAdd('Autres    CIS FrozenCore NoFrozenCore PrintBasis');
 TAdd('          PrintMOs NoPrintMOs Grid1..Grid5');
 TAdd('');
 TAdd('Blocs :');
 TAdd('  %scf MaxIter n  Convergence tight|verytight  end');
 TAdd('  %geom MaxIter n  Convergence tight  end');
 TAdd('  %geom Scan B i j = debut, fin, n  end  end');
 TAdd('  %geom Constraints { B i j C } { A i j k C } end end');
 TAdd('  %tddft NRoots n  Triplets true  end');
 TAdd('  %freq Temp 298.15 end');
 TAdd('Les numeros d''atomes des blocs commencent a 0 (comme ORCA).');
 TView('Mots-cles');
End;

Procedure ShowAbout;
Begin
 ShowMessage('A propos', 'ORCA Clone 1.00 - QChimie|Clone du logiciel de chimie quantique ORCA|' +
             'Auteur : Sylvain Maltais (support@gladir.com)|' +
             'Moteur : HF, MP2, CIS, DFT (Becke), bases STO-3G/3-21G/6-31G|https://www.gladir.com/qchimie');
End;

{ --------------------------------------------------------------------- }
{ Menus                                                                 }
{ --------------------------------------------------------------------- }

Procedure ZeroMenu(Var M: TMenu);
Var I: Integer;
Begin
 For I := 0 To MAX_MENU_ITEMS - 1 Do Begin
  M.Items[I].Caption := '';
  M.Items[I].Shortcut := '';
  M.Items[I].Enabled := True;
  M.Items[I].Checked := False;
  M.Items[I].Separator := False;
 End;
End;

Procedure MI(Var Mn: TMenu; Idx: Integer; Const Cap, Sc: String);
Begin
 Mn.Items[Idx].Caption := Cap;
 Mn.Items[Idx].Shortcut := Sc;
End;

Procedure MSep(Var Mn: TMenu; Idx: Integer);
Begin
 Mn.Items[Idx].Separator := True;
End;

Procedure InitMenus;
Var I: Integer;
Begin
 ZeroMenu(Menus[MENU_FILE]);
 Menus[MENU_FILE].Title := 'Fichier';
 Menus[MENU_FILE].ItemCount := 10;
 MI(Menus[MENU_FILE], 0, 'Nouveau', 'Ctrl+N');
 MI(Menus[MENU_FILE], 1, 'Ouvrir un fichier .inp...', 'Ctrl+O');
 MI(Menus[MENU_FILE], 2, 'Importer une geometrie XYZ...', '');
 MI(Menus[MENU_FILE], 3, 'Enregistrer', 'Ctrl+S');
 MI(Menus[MENU_FILE], 4, 'Enregistrer sous...', '');
 MSep(Menus[MENU_FILE], 5);
 MI(Menus[MENU_FILE], 6, 'Enregistrer la sortie (.out)...', '');
 MI(Menus[MENU_FILE], 7, 'Exporter la geometrie (.xyz)...', '');
 MSep(Menus[MENU_FILE], 8);
 MI(Menus[MENU_FILE], 9, 'Quitter', 'Ctrl+Q');

 ZeroMenu(Menus[MENU_EDIT]);
 Menus[MENU_EDIT].Title := 'Edition';
 Menus[MENU_EDIT].ItemCount := 9;
 MI(Menus[MENU_EDIT], 0, 'Assistant de calcul...', 'F2');
 MI(Menus[MENU_EDIT], 1, 'Inserer un mot-cle...', '');
 MI(Menus[MENU_EDIT], 2, 'Geometrie finale -> entree', '');
 MSep(Menus[MENU_EDIT], 3);
 MI(Menus[MENU_EDIT], 4, 'Supprimer la ligne', 'Ctrl+Y');
 MI(Menus[MENU_EDIT], 5, 'Tout effacer', '');
 MI(Menus[MENU_EDIT], 6, 'Aller a la ligne...', 'Ctrl+G');
 MI(Menus[MENU_EDIT], 7, 'Mode insertion / ecrasement', 'Inser');
 MI(Menus[MENU_EDIT], 8, 'Commenter / decommenter', 'Ctrl+K');

 ZeroMenu(Menus[MENU_CALC]);
 Menus[MENU_CALC].Title := 'Calcul';
 Menus[MENU_CALC].ItemCount := 6;
 MI(Menus[MENU_CALC], 0, 'Executer', 'F5');
 MSep(Menus[MENU_CALC], 1);
 MI(Menus[MENU_CALC], 2, 'Options SCF et DFT...', '');
 MI(Menus[MENU_CALC], 3, 'Options d''optimisation...', '');
 MI(Menus[MENU_CALC], 4, 'Etats excites (CIS)...', '');
 MI(Menus[MENU_CALC], 5, 'Temperature thermochimie...', '');

 ZeroMenu(Menus[MENU_RES]);
 Menus[MENU_RES].Title := 'Resultats';
 Menus[MENU_RES].ItemCount := 17;
 MI(Menus[MENU_RES], 0, 'Sortie ORCA (.out)', 'F3');
 MI(Menus[MENU_RES], 1, 'Resume des energies', '');
 MI(Menus[MENU_RES], 2, 'Orbitales moleculaires', 'F6');
 MI(Menus[MENU_RES], 3, 'Coupe d''orbitale / densite', '');
 MI(Menus[MENU_RES], 4, 'Charges, populations et dipole', '');
 MSep(Menus[MENU_RES], 5);
 MI(Menus[MENU_RES], 6, 'Frequences (animation)', '');
 MI(Menus[MENU_RES], 7, 'Spectre infrarouge', 'F4');
 MI(Menus[MENU_RES], 8, 'Thermochimie', '');
 MSep(Menus[MENU_RES], 9);
 MI(Menus[MENU_RES], 10, 'Etats excites CIS', '');
 MI(Menus[MENU_RES], 11, 'Spectre UV-visible', '');
 MSep(Menus[MENU_RES], 12);
 MI(Menus[MENU_RES], 13, 'Convergence du SCF', '');
 MI(Menus[MENU_RES], 14, 'Courbe d''optimisation', '');
 MI(Menus[MENU_RES], 15, 'Surface d''energie (scan)', '');
 MI(Menus[MENU_RES], 16, 'Resume du calcul', '');

 ZeroMenu(Menus[MENU_VIEW]);
 Menus[MENU_VIEW].Title := 'Affichage';
 Menus[MENU_VIEW].ItemCount := 16;
 MI(Menus[MENU_VIEW], 0, 'Fil de fer', '');
 MI(Menus[MENU_VIEW], 1, 'Batonnets', '');
 MI(Menus[MENU_VIEW], 2, 'Boules et batons', '');
 MI(Menus[MENU_VIEW], 3, 'Spheres (CPK)', '');
 MSep(Menus[MENU_VIEW], 4);
 MI(Menus[MENU_VIEW], 5, 'Symboles des atomes', '');
 MI(Menus[MENU_VIEW], 6, 'Numeros des atomes', '');
 MI(Menus[MENU_VIEW], 7, 'Charges de Mulliken', '');
 MI(Menus[MENU_VIEW], 8, 'Axes', '');
 MSep(Menus[MENU_VIEW], 9);
 MI(Menus[MENU_VIEW], 10, 'Reinitialiser la vue', 'F9');
 MI(Menus[MENU_VIEW], 11, 'Agrandir', 'F7');
 MI(Menus[MENU_VIEW], 12, 'Reduire', 'F8');
 MSep(Menus[MENU_VIEW], 13);
 MI(Menus[MENU_VIEW], 14, 'Animer un mode de vibration', '');
 MI(Menus[MENU_VIEW], 15, 'Arreter l''animation', '');

 ZeroMenu(Menus[MENU_EXAMPLE]);
 Menus[MENU_EXAMPLE].Title := 'Exemples';
 Menus[MENU_EXAMPLE].ItemCount := NUM_EX;
 For I := 0 To NUM_EX - 1 Do MI(Menus[MENU_EXAMPLE], I, ExTitle[I], '');

 ZeroMenu(Menus[MENU_TOOLS]);
 Menus[MENU_TOOLS].Title := 'Outils';
 Menus[MENU_TOOLS].ItemCount := 5;
 MI(Menus[MENU_TOOLS], 0, 'Informations moleculaires', '');
 MI(Menus[MENU_TOOLS], 1, 'Mesures geometriques...', '');
 MI(Menus[MENU_TOOLS], 2, 'Convertisseur d''energie...', '');
 MI(Menus[MENU_TOOLS], 3, 'Elements supportes', '');
 MI(Menus[MENU_TOOLS], 4, 'Mots-cles ORCA supportes', '');

 ZeroMenu(Menus[MENU_HELP]);
 Menus[MENU_HELP].Title := 'Aide';
 Menus[MENU_HELP].ItemCount := 3;
 MI(Menus[MENU_HELP], 0, 'Aide', 'F1');
 MI(Menus[MENU_HELP], 1, 'Mots-cles', '');
 MI(Menus[MENU_HELP], 2, 'A propos...', '');
End;

Procedure UpdateMenuChecks;
Var I: Integer;
Begin
 For I := 0 To 3 Do Menus[MENU_VIEW].Items[I].Checked := (ViewStyle = I);
 Menus[MENU_VIEW].Items[5].Checked := ShowSym;
 Menus[MENU_VIEW].Items[6].Checked := ShowNum;
 Menus[MENU_VIEW].Items[7].Checked := ShowChg;
 Menus[MENU_VIEW].Items[8].Checked := ShowAxes;
 Menus[MENU_EDIT].Items[7].Checked := InsertMode;
 Menus[MENU_RES].Items[1].Enabled := HaveResult;
 Menus[MENU_RES].Items[2].Enabled := HaveResult;
 Menus[MENU_RES].Items[3].Enabled := HaveResult;
 Menus[MENU_RES].Items[4].Enabled := HaveResult;
 Menus[MENU_RES].Items[6].Enabled := FreqOK;
 Menus[MENU_RES].Items[7].Enabled := FreqOK;
 Menus[MENU_RES].Items[8].Enabled := ThermoOK;
 Menus[MENU_RES].Items[10].Enabled := CISOK;
 Menus[MENU_RES].Items[11].Enabled := CISOK;
 Menus[MENU_RES].Items[13].Enabled := HistN > 1;
 Menus[MENU_RES].Items[14].Enabled := OptN > 1;
 Menus[MENU_RES].Items[15].Enabled := ScanN > 1;
 Menus[MENU_VIEW].Items[14].Enabled := FreqOK;
 Menus[MENU_VIEW].Items[15].Enabled := VibMode > 0;
End;

{ --------------------------------------------------------------------- }
{ Actions de fichier                                                     }
{ --------------------------------------------------------------------- }

Function ConfirmDiscard: Boolean;
Begin
 If InpModified Then
  ConfirmDiscard := ConfirmBox('Fichier modifie', 'Le fichier d''entree a ete modifie.|Abandonner les modifications ?')
 Else ConfirmDiscard := True;
End;

Procedure AfterEdit;
Begin
 PreviewInput;
 EdClamp;
 EdScroll;
End;

Procedure DoNew;
Begin
 If Not ConfirmDiscard Then Exit;
 NInp := 0;
 EdAppend('# Nouveau calcul ORCA');
 EdAppend('! HF STO-3G');
 EdAppend('* xyz 0 1');
 EdAppend('*');
 InpName := '';
 InpModified := False;
 CurLn := 1; CurCol := 0; TopLn := 0; LeftCol := 0;
 ResStale := True;
 AfterEdit;
 StatusMsg := 'Nouveau fichier d''entree';
End;

Procedure DoOpen;
Var FN: String;
Begin
 If Not ConfirmDiscard Then Exit;
 FN := '';
 If FileOpenDialog('Ouvrir un fichier d''entree ORCA', FN, '*.INP;*.TXT;*.IN') Then Begin
  If LoadInputFile(FN) Then Begin
   CurLn := 0; CurCol := 0; TopLn := 0; LeftCol := 0;
   ResStale := True;
   EnsureOneLine;
   AfterEdit;
   ResetView;
   StatusMsg := 'Fichier ouvert : ' + FN;
  End
  Else ShowMessage('Erreur', 'Impossible de lire ' + FN);
 End;
End;

Procedure DoSaveAs;
Var FN: String;
Begin
 FN := InpName;
 If FN = '' Then FN := 'ORCA.INP';
 If FileSaveDialog('Enregistrer le fichier d''entree', FN) Then Begin
  If ExtractExt(FN) = '' Then FN := FN + '.inp';
  If SaveInputFile(FN) Then StatusMsg := 'Enregistre : ' + FN
  Else ShowMessage('Erreur', 'Impossible d''ecrire ' + FN);
 End;
End;

Procedure DoSave;
Begin
 If InpName = '' Then DoSaveAs
 Else If SaveInputFile(InpName) Then StatusMsg := 'Enregistre : ' + InpName
 Else ShowMessage('Erreur', 'Impossible d''ecrire ' + InpName);
End;

Type
 TZArr = Array[0..MAX_AT-1] Of Integer;
 TXArr = Array[0..3*MAX_AT-1] Of Real;

{ Insere une ligne avant le bloc de coordonnees (ou a la fin) }
Procedure InsertBeforeCoords(Const S: String);
Var C: Integer;
Begin
 C := FindCoordLine;
 If C < 0 Then C := NInp;
 EdInsertLine(C, S);
End;

{ Remplace le bloc de coordonnees par la geometrie donnee }
Procedure ReplaceGeometry(N: Integer; Const Z: TZArr; Const X: TXArr; Ch, Mu: Integer);
Var C, I, J: Integer; Hdr: String;
Begin
 C := FindCoordLine;
 If C >= 0 Then Begin
  { supprime l'ancien bloc jusqu'a la ligne * de fin }
  Hdr := TrimStr(InpLines[C]);
  I := C + 1;
  While (I < NInp) And (Copy(TrimStr(InpLines[I]), 1, 1) <> '*') And (UpStr(TrimStr(InpLines[I])) <> 'END') Do Inc(I);
  If I < NInp Then Inc(I);
  For J := C To I - 1 Do EdDeleteLine(C);
 End;
 If C < 0 Then C := NInp;
 EdInsertLine(C, '* xyz ' + IStr(Ch) + ' ' + IStr(Mu));
 For I := 0 To N - 1 Do
  EdInsertLine(C + 1 + I, PadR(ElemStr(Z[I]), 3) + FmtR(X[I * 3] / ANG2BOHR, 14, 8) +
               FmtR(X[I * 3 + 1] / ANG2BOHR, 14, 8) + FmtR(X[I * 3 + 2] / ANG2BOHR, 14, 8));
 EdInsertLine(C + 1 + N, '*');
 MarkModified;
End;

Procedure DoImportXYZ;
Var FN, Err: String; Z: TZArr; X: TXArr; I, SvNa, A, D: Integer;
    Sz: Array[1..MAX_AT] Of Integer; Sx: Array[1..MAX_AT, 1..3] Of Real;
Begin
 FN := '';
 If Not FileOpenDialog('Importer une geometrie XYZ', FN, '*.XYZ') Then Exit;
 { sauvegarde de la geometrie courante }
 SvNa := NAt;
 For A := 1 To NAt Do Begin Sz[A] := AtZ[A]; For D := 1 To 3 Do Sx[A, D] := AtX[A, D]; End;
 If ReadXYZFile(FN, Err) Then Begin
  For I := 1 To NAt Do Begin
   Z[I - 1] := AtZ[I];
   For D := 1 To 3 Do X[(I - 1) * 3 + D - 1] := AtX[I, D];
  End;
  I := NAt;
  NAt := SvNa;
  For A := 1 To NAt Do Begin AtZ[A] := Sz[A]; For D := 1 To 3 Do AtX[A, D] := Sx[A, D]; End;
  PreviewInput;
  If PrevOK Then ReplaceGeometry(I, Z, X, Charge, Mult) Else ReplaceGeometry(I, Z, X, 0, 1);
  AfterEdit;
  StatusMsg := 'Geometrie importee : ' + FN;
 End
 Else Begin
  NAt := SvNa;
  For A := 1 To NAt Do Begin AtZ[A] := Sz[A]; For D := 1 To 3 Do AtX[A, D] := Sx[A, D]; End;
  ShowMessage('Erreur', Err);
 End;
End;

Procedure DoExportXYZ;
Var FN: String;
Begin
 If Not (HaveResult And Not ResStale) And Not PrevOK Then Begin ShowMessage('Erreur', 'Geometrie invalide.'); Exit; End;
 FN := ChangeExt(InpName, '.xyz');
 If FN = '.xyz' Then FN := 'ORCA.XYZ';
 If Not FileSaveDialog('Exporter la geometrie XYZ', FN) Then Exit;
 If HaveResult And Not ResStale Then UseResultGeom Else UsePrevGeom;
 SaveXYZ(FN, 'ORCA Clone - ' + ResMethod);
 RestoreGeom;
 StatusMsg := 'Geometrie exportee : ' + FN;
End;

Procedure DoSaveOutput;
Var FN: String;
Begin
 If NOut = 0 Then Begin ShowMessage('Sortie', 'Aucune sortie : executez un calcul.'); Exit; End;
 FN := ChangeExt(InpName, '.out');
 If FN = '.out' Then FN := 'ORCA.OUT';
 If FileSaveDialog('Enregistrer la sortie', FN) Then Begin
  If SaveOutputFile(FN) Then StatusMsg := 'Sortie enregistree : ' + FN
  Else ShowMessage('Erreur', 'Impossible d''ecrire ' + FN);
 End;
End;

Procedure DoFinalToInput;
Var Z: TZArr; X: TXArr; I, D: Integer;
Begin
 If Not NeedResult Then Exit;
 For I := 1 To ResN Do Begin
  Z[I - 1] := ResZ[I];
  For D := 1 To 3 Do X[(I - 1) * 3 + D - 1] := ResX[I, D];
 End;
 ReplaceGeometry(ResN, Z, X, ResCharge, ResMult);
 AfterEdit;
 StatusMsg := 'Geometrie finale copiee dans l''entree';
End;

{ --------------------------------------------------------------------- }
{ Execution                                                              }
{ --------------------------------------------------------------------- }

Procedure DoRun;
Var Err, FN: String; Ok: Boolean;
Begin
 If CalcRunning Then Exit;
 EnsureOneLine;
 CalcRunning := True;
 UserAbort := False;
 StatusMsg := 'Calcul en cours...';
 BusyTick := 0;
 GuiActive := True;
 NeedRedraw := True;
 { fichiers de sortie }
 FN := ChangeExt(InpName, '.out');
 If InpName = '' Then FN := 'ORCA.OUT';
 RedrawAll;
 Ok := ExecuteJob(Err);
 CalcRunning := False;
 If Ok Or OrbsValid Then CaptureResults Else HaveResult := False;
 ResStale := False;
 VibMode := 0;
 If Ok Then StatusMsg := 'Calcul termine en ' + RStr(CalcMs / 1000, 2) + ' s'
 Else StatusMsg := 'Calcul interrompu : ' + Err;
 SaveOutputFile(FN);
 If Not Ok Then ShowMessage('Calcul ORCA', Copy(Err, 1, 60) + '|Voir la sortie (F3) pour le detail.');
 NeedRedraw := True;
End;

{ --------------------------------------------------------------------- }
{ Parametres divers                                                      }
{ --------------------------------------------------------------------- }

Procedure OptionsSCF;
Var B: Integer; Line: String;
Begin
 PreviewInput;
 PReset;
 PAddH('Convergence SCF');
 PAddI('Nombre maximal d''iterations', ScfMaxIter, 10, 500);
 PAddE('Seuils de convergence', 0, 'NormalSCF|TightSCF|VeryTightSCF|LooseSCF');
 PAddH('Integration DFT');
 PAddE('Grille de Becke', GridScale - 1, 'Normale|Fine (lente)');
 If Not ParamDialog('Options SCF et DFT', 'Ajoute un bloc %scf au fichier d''entree.') Then Exit;
 Line := '%scf MaxIter ' + IStr(PInt(1));
 Case PInt(2) Of
  1: Line := Line + ' Convergence tight';
  2: Line := Line + ' Convergence verytight';
 End;
 Line := Line + ' end';
 InsertBeforeCoords(Line);
 If PInt(4) = 1 Then AddKeywords('Grid3');
 If PInt(2) = 3 Then AddKeywords('LooseSCF');
 MarkModified; AfterEdit;
End;

Procedure OptionsOpt;
Var Line: String;
Begin
 PReset;
 PAddH('Optimisation de geometrie');
 PAddI('Cycles maximum', 60, 5, 300);
 PAddE('Seuils', 0, 'NormalOpt|TightOpt|LooseOpt');
 If Not ParamDialog('Options d''optimisation', 'Ajoute un bloc %geom au fichier d''entree.') Then Exit;
 Line := '%geom MaxIter ' + IStr(PInt(1));
 Case PInt(2) Of
  1: Line := Line + ' Convergence tight';
  2: Line := Line + ' Convergence loose';
 End;
 Line := Line + ' end';
 InsertBeforeCoords(Line);
 MarkModified; AfterEdit;
End;

Procedure OptionsCIS;
Begin
 PReset;
 PAddH('Etats excites CIS (couche fermee, HF)');
 PAddI('Nombre d''etats (NRoots)', 5, 1, 20);
 PAddB('Etats triplets', False);
 If Not ParamDialog('Etats excites', 'Ajoute le mot-cle CIS et un bloc %tddft.') Then Exit;
 AddKeywords('CIS');
 If PBool(2) Then InsertBeforeCoords('%tddft nroots ' + IStr(PInt(1)) + ' triplets true end')
 Else InsertBeforeCoords('%tddft nroots ' + IStr(PInt(1)) + ' end');
 MarkModified; AfterEdit;
End;

Procedure OptionsTemp;
Var S: String; V: Real;
Begin
 S := RTrim(CalcTemp, 2);
 If InputString('Temperature', 'Temperature de la thermochimie (K) :', S) And ParseReal(S, V) And (V > 1) And (V < 5000) Then Begin
  InsertBeforeCoords('%freq Temp ' + RTrim(V, 2) + ' end');
  MarkModified; AfterEdit;
 End;
End;

Procedure CommentLine;
Var L: String;
Begin
 EdClamp;
 L := InpLines[CurLn];
 If (Length(L) > 0) And (L[1] = '#') Then Delete(L, 1, 1)
 Else L := Copy('#' + L, 1, 100);
 InpLines[CurLn] := L;
 MarkModified; AfterEdit;
End;

Procedure GotoLine;
Var S: String; V: Real;
Begin
 S := IStr(CurLn + 1);
 If InputString('Aller a la ligne', 'Numero de ligne (1 a ' + IStr(NInp) + ') :', S) And ParseReal(S, V) Then Begin
  CurLn := Round(V) - 1; CurCol := 0; EdClamp; EdScroll;
 End;
End;

Procedure AnimateDialog;
Begin
 If Not NeedResult Then Exit;
 If Not FreqOK Then Begin ShowMessage('Animation', 'Aucun calcul de frequences.'); Exit; End;
 ShowFrequencies;
End;

{ --------------------------------------------------------------------- }
{ Dispatch des menus et de la barre d'outils                             }
{ --------------------------------------------------------------------- }

Procedure CycleStyle;
Begin
 ViewStyle := (ViewStyle + 1) Mod 4;
 NeedRedraw := True;
End;

Procedure Zoom(F: Real);
Begin
 ZoomF := ZoomF * F;
 If ZoomF < 0.2 Then ZoomF := 0.2;
 If ZoomF > 6 Then ZoomF := 6;
 NeedRedraw := True;
End;

Procedure DispatchMenu(M, I: Integer);
Begin
 UpdateMenuChecks;
 RedrawAll;
 NeedRedraw := True;
 Case M Of
  MENU_FILE: Case I Of
   0: DoNew;
   1: DoOpen;
   2: DoImportXYZ;
   3: DoSave;
   4: DoSaveAs;
   6: DoSaveOutput;
   7: DoExportXYZ;
   9: Begin If ConfirmDiscard Then Running := False; End;
  End;
  MENU_EDIT: Case I Of
   0: AssistantDialog;
   1: InsertKeywordsDialog;
   2: DoFinalToInput;
   4: Begin EdClamp; EdDeleteLine(CurLn); MarkModified; AfterEdit; End;
   5: If ConfirmBox('Tout effacer', 'Effacer tout le fichier d''entree ?') Then Begin EdClear; AfterEdit; End;
   6: GotoLine;
   7: InsertMode := Not InsertMode;
   8: CommentLine;
  End;
  MENU_CALC: Case I Of
   0: DoRun;
   2: OptionsSCF;
   3: OptionsOpt;
   4: OptionsCIS;
   5: OptionsTemp;
  End;
  MENU_RES: Case I Of
   0: OutputViewer;
   1: ShowEnergySummary;
   2: OrbitalDialog;
   3: Begin If NeedResult Then If OrbsValid Then OrbitalCut(NOccA, 1); End;
   4: ShowPopulations;
   6: ShowFrequencies;
   7: ShowIRChart;
   8: ShowThermo;
   10: ShowExcited;
   11: ShowUVChart;
   13: ShowSCFChart;
   14: ShowOptChart;
   15: ShowScanChart;
   16: ShowEnergySummary;
  End;
  MENU_VIEW: Case I Of
   0..3: Begin ViewStyle := I; End;
   5: ShowSym := Not ShowSym;
   6: ShowNum := Not ShowNum;
   7: ShowChg := Not ShowChg;
   8: ShowAxes := Not ShowAxes;
   10: ResetView;
   11: Zoom(1.25);
   12: Zoom(0.8);
   14: AnimateDialog;
   15: VibMode := 0;
  End;
  MENU_EXAMPLE: Begin
   If ConfirmDiscard Then LoadExample(I);
  End;
  MENU_TOOLS: Case I Of
   0: ShowMolInfo;
   1: MeasureDialog;
   2: ConvertEnergy;
   3: ShowElements;
   4: ShowKeywordHelp;
  End;
  MENU_HELP: Case I Of
   0: ShowHelp;
   1: ShowKeywordHelp;
   2: ShowAbout;
  End;
 End;
 NeedRedraw := True;
End;

Procedure DoToolbar(Id: Integer);
Begin
 Case Id Of
  TB_RUN: DoRun;
  TB_OUT: OutputViewer;
  TB_MOL: CycleStyle;
  TB_ORB: OrbitalDialog;
  TB_SCF: ShowSCFChart;
  TB_OPT: ShowOptChart;
  TB_IR: ShowIRChart;
  TB_UV: ShowUVChart;
  TB_SCAN: ShowScanChart;
  TB_EX: Begin If ConfirmDiscard Then ExamplesDialog; End;
  TB_HELP: ShowHelp;
 End;
 NeedRedraw := True;
End;

{ --------------------------------------------------------------------- }
{ Clavier et souris                                                      }
{ --------------------------------------------------------------------- }

Procedure RedrawEditorOnly;
Begin
 HideMouseCursor;
 DrawEditor;
 DrawMolecule;
 DrawStatusBar;
 ShowMouseCursor;
End;

Procedure HandleKey;
Var K: Char; MZ, MN, P: Integer; Edited: Boolean;
Begin
 K := ReadKey;
 Edited := False;
 If K = #0 Then Begin
  K := ReadKey;
  If MenuActive Then Begin
   Case K Of
    #75: Begin
     Dec(ActiveMenu); If ActiveMenu < 0 Then ActiveMenu := NUM_MENUS - 1;
     ActiveMenuItem := 0; MenuNextItem(0); If ActiveMenuItem < 0 Then ActiveMenuItem := 0;
    End;
    #77: Begin
     Inc(ActiveMenu); If ActiveMenu >= NUM_MENUS Then ActiveMenu := 0;
     ActiveMenuItem := 0; MenuNextItem(0); If ActiveMenuItem < 0 Then ActiveMenuItem := 0;
    End;
    #72: MenuNextItem(-1);
    #80: MenuNextItem(1);
    #68, #67: Begin MenuActive := False; ActiveMenu := -1; End;
   End;
   NeedRedraw := True;
   Exit;
  End;
  Case K Of
   #59: ShowHelp;
   #60: AssistantDialog;
   #61: OutputViewer;
   #62: ShowIRChart;
   #63: DoRun;
   #64: OrbitalDialog;
   #65: Zoom(1.25);
   #66: Zoom(0.8);
   #67: ResetView;
   #68: Begin UpdateMenuChecks; MenuActive := True; ActiveMenu := 0; ActiveMenuItem := 0; MenuNextItem(0); End;
   #45: Begin If ConfirmDiscard Then Running := False; End;
   #72: Begin Dec(CurLn); EdClamp; End;
   #80: Begin Inc(CurLn); EdClamp; End;
   #75: Begin If CurCol > 0 Then Dec(CurCol) Else If CurLn > 0 Then Begin Dec(CurLn); CurCol := 100; End; EdClamp; End;
   #77: Begin
    If CurCol < Length(InpLines[CurLn]) Then Inc(CurCol)
    Else If CurLn < NInp - 1 Then Begin Inc(CurLn); CurCol := 0; End;
    EdClamp;
   End;
   #71: CurCol := 0;
   #79: CurCol := Length(InpLines[CurLn]);
   #73: Begin Dec(CurLn, ED_ROWS - 1); EdClamp; End;
   #81: Begin Inc(CurLn, ED_ROWS - 1); EdClamp; End;
   #83: Begin EdDelete; Edited := True; End;
   #82: InsertMode := Not InsertMode;
  End;
  EdClamp; EdScroll;
  If Edited Then AfterEdit;
  NeedRedraw := True;
  Exit;
 End;
 If MenuActive Then Begin
  Case K Of
   #27: Begin MenuActive := False; ActiveMenu := -1; End;
   #13: Begin
    If (ActiveMenuItem >= 0) And Not Menus[ActiveMenu].Items[ActiveMenuItem].Separator Then Begin
     MZ := ActiveMenu; MN := ActiveMenuItem;
     MenuActive := False; ActiveMenu := -1;
     DispatchMenu(MZ, MN);
    End;
   End;
  End;
  NeedRedraw := True;
  Exit;
 End;
 Case K Of
  #14: DoNew;
  #15: DoOpen;
  #19: DoSave;
  #17: Begin If ConfirmDiscard Then Running := False; End;
  #18: DoRun;
  #25: Begin EdClamp; EdDeleteLine(CurLn); MarkModified; Edited := True; End;
  #7: GotoLine;
  #11: CommentLine;
  #27: Begin If VibMode > 0 Then VibMode := 0; End;
  #13: Begin EdEnter; Edited := True; End;
  #8: Begin EdBackspace; Edited := True; End;
  #9: Begin EdTypeChar(' '); EdTypeChar(' '); Edited := True; End;
 Else
  If (K >= ' ') Then Begin EdTypeChar(K); Edited := True; End;
 End;
 EdClamp; EdScroll;
 If Edited Then Begin
  AfterEdit;
  RedrawEditorOnly;
 End
 Else NeedRedraw := True;
End;

Procedure HandleMouse;
Var
 Hit, ItemHit, T, MZ, Row, Col: Integer;
 LeftDown, OldLeft, Moved: Boolean;
 DX, DY: Integer;
Begin
 GetMouseLocal(MouseX, MouseY, MouseButtons);
 LeftDown := (MouseButtons And 1) <> 0;
 OldLeft := (OldMouseButtons And 1) <> 0;
 Moved := (MouseX <> OldMouseX) Or (MouseY <> OldMouseY);
 DlgHappened := False;
 If MenuActive And Moved Then Begin
  Hit := MenuBarHitTest(MouseX, MouseY);
  If (Hit >= 0) And (Hit <> ActiveMenu) Then Begin
   ActiveMenu := Hit; ActiveMenuItem := -1; NeedRedraw := True;
  End
  Else Begin
   ItemHit := MenuItemHitTest(MouseX, MouseY, ActiveMenu);
   If (ItemHit >= 0) And (ItemHit <> ActiveMenuItem) And
      (Not Menus[ActiveMenu].Items[ItemHit].Separator) Then Begin
    ActiveMenuItem := ItemHit; NeedRedraw := True;
   End;
  End;
 End;
 If Moved And Not MenuActive And Not DragOn Then Begin
  T := ToolbarHit(MouseX, MouseY);
  If T >= 0 Then Begin
   If StatusMsg <> TBHint[T] Then Begin StatusMsg := TBHint[T]; HideMouseCursor; DrawStatusBar; ShowMouseCursor; End;
  End;
 End;
 { rotation a la souris dans la vue }
 If DragOn Then Begin
  If LeftDown Then Begin
   If Moved Then Begin
    DX := MouseX - DragX; DY := MouseY - DragY;
    RotateView(DX * 0.02, DY * 0.02);
    DragX := MouseX; DragY := MouseY;
    HideMouseCursor; DrawMolecule; ShowMouseCursor;
   End;
  End
  Else DragOn := False;
 End;
 If LeftDown And Not OldLeft Then Begin
  Hit := MenuBarHitTest(MouseX, MouseY);
  If Hit >= 0 Then Begin
   If MenuActive And (ActiveMenu = Hit) Then Begin MenuActive := False; ActiveMenu := -1; End
   Else Begin MenuActive := True; ActiveMenu := Hit; ActiveMenuItem := -1; UpdateMenuChecks; End;
   NeedRedraw := True;
  End
  Else If MenuActive Then Begin
   ItemHit := MenuItemHitTest(MouseX, MouseY, ActiveMenu);
   If (ItemHit >= 0) And Not Menus[ActiveMenu].Items[ItemHit].Separator Then Begin
    MZ := ActiveMenu;
    MenuActive := False; ActiveMenu := -1;
    DispatchMenu(MZ, ItemHit);
   End
   Else If ItemHit < 0 Then Begin MenuActive := False; ActiveMenu := -1; End;
   NeedRedraw := True;
  End
  Else If (MouseY >= MENUBAR_H) And (MouseY < MENUBAR_H + TOOLBAR_H) Then Begin
   T := ToolbarHit(MouseX, MouseY);
   If (T >= 0) And TBEnabled(TBId[T]) Then DoToolbar(TBId[T]);
  End
  Else If PtInRect(MouseX, MouseY, ED_X1, ED_Y1 + 14, ED_X2, ED_Y2) Then Begin
   Row := (MouseY - ED_Y1 - 14 - 2) Div 10 + TopLn;
   Col := (MouseX - 4) Div 8 + LeftCol;
   If Row >= NInp Then Row := NInp - 1;
   CurLn := Row; CurCol := Col;
   EdClamp; EdScroll;
   NeedRedraw := True;
  End
  Else If PtInRect(MouseX, MouseY, VW_X1, VW_Y1, VW_X2, VW_Y2) Then Begin
   DragOn := True; DragX := MouseX; DragY := MouseY;
  End;
 End;
 If Not DlgHappened Then OldMouseButtons := MouseButtons;
 OldMouseX := MouseX; OldMouseY := MouseY;
End;

{ --------------------------------------------------------------------- }
{ Mode lot                                                               }
{ --------------------------------------------------------------------- }

Procedure BatchRun;
Var FN, OutFN, Err: String; Ok: Boolean;
Begin
 BatchMode := True;
 FN := ParamStr(2);
 If FN = '' Then Begin WriteLn('Syntaxe : ORCA /B fichier.inp [sortie.out]'); Halt(1); End;
 If Not LoadInputFile(FN) Then Begin WriteLn('Fichier introuvable : ', FN); Halt(1); End;
 If ParamCount >= 3 Then OutFN := ParamStr(3) Else OutFN := ChangeExt(FN, '.out');
 Assign(OutFile, OutFN);
 {$I-}
 Rewrite(OutFile);
 {$I+}
 If IOResult <> 0 Then Begin WriteLn('Impossible d''ecrire ', OutFN); Halt(1); End;
 OutFileOpen := True;
 WriteLn('ORCA Clone - mode lot : ', FN);
 Ok := ExecuteJob(Err);
 Close(OutFile);
 OutFileOpen := False;
 If Ok Then Begin
  WriteLn('Energie finale : ', ETot:0:10, ' Eh');
  WriteLn('Sortie : ', OutFN);
 End
 Else WriteLn('ERREUR : ', Err);
End;

{ --------------------------------------------------------------------- }
{ Initialisation et boucle principale                                    }
{ --------------------------------------------------------------------- }

Procedure InitializeGraphics;
Var
 Driver, Mode : Integer;
 ErrCode : Integer;
Begin
 {$IFDEF FPC}
  Driver := VGA;
  Mode := VGAHi;
 {$ELSE}
  Driver := Detect;
  Mode := VGAHi;
 {$ENDIF}
 InitGraph(Driver, Mode, '');
 ErrCode := GraphResult;
 If ErrCode = grOk Then Begin
  SetColor(White);
  SetLineStyle(0, 0, 1);
 End
 Else Begin
  WriteLn('Erreur graphique : ', GraphErrorMsg(ErrCode));
  Halt;
 End;
 ScreenWidth := GetMaxX + 1;
 ScreenHeight := GetMaxY + 1;
End;

Procedure InitializeAll;
Begin
 InitMenus;
 InitToolbar;
 LayoutMenuBar;
 MenuActive := False; ActiveMenu := -1; ActiveMenuItem := -1;
 OldMouseButtons := 0; OldMouseX := -1; OldMouseY := -1;
 ProgPct := -1;
 InitMouseDriver;
 Running := True;
 NeedRedraw := True;
 ViewStyle := STYLE_BALL;
 ShowSym := True; ShowNum := False; ShowChg := False; ShowAxes := True;
 InsertMode := True;
 ResetView;
 StatusMsg := 'Pret. F5 = executer le calcul, F1 = aide';
 GuiActive := True;
End;

Procedure MainLoop;
Var LastAnim: LongInt;
Begin
 RedrawAll;
 LastAnim := NowMs;
 While Running Do Begin
  If KeyPressed Then HandleKey;
  HandleMouse;
  If (VibMode > 0) And FreqOK And Not MenuActive And Not DragOn Then Begin
   If ElapsedMs(LastAnim) > 70 Then Begin
    LastAnim := NowMs;
    VibPhase := VibPhase + 0.5;
    If VibPhase > 2 * Pi Then VibPhase := VibPhase - 2 * Pi;
    HideMouseCursor; DrawMolecule; ShowMouseCursor;
   End;
  End;
  If MenuActive Then UpdateMenuChecks;
  If NeedRedraw Then Begin RedrawAll; LastAnim := NowMs; End;
  Delay(8);
 End;
End;

{ ===================================================================== }
{ Programme principal                                                   }
{ ===================================================================== }

Var Arg1: String;
BEGIN
 NOut := 0; OutFileOpen := False; NInp := 0; InpName := ''; InpModified := False;
 DefaultParams;
 GridScale := 1;
 ERIThresh := 1E-11;
 Arg1 := ParamStr(1);
 If (Arg1 = '/?') Or (Arg1 = '--help') Or (Arg1 = '-h') Then Begin
  WriteLn('ORCA - Clone du logiciel de chimie quantique ORCA (calculs ab initio et DFT).');
  WriteLn;
  WriteLn('Syntaxe : ORCA [fichier.INP | fichier.XYZ]');
  WriteLn('          ORCA /B fichier.INP [sortie.OUT]');
  WriteLn;
  WriteLn(' fichier.INP   Fichier d''entree ORCA a ouvrir dans l''editeur');
  WriteLn(' fichier.XYZ   Geometrie a importer (calcul HF/STO-3G par defaut)');
  WriteLn(' /B            Mode lot : execute le calcul sans interface graphique');
  Halt;
 End;
 If (UpStr(Arg1) = '/B') Or (UpStr(Arg1) = '-B') Then Begin
  BatchRun;
  Halt;
 End;
 InitializeGraphics;
 InitializeAll;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
 LoadExample(1);
 If ParamCount >= 1 Then Begin
  If FileOK(Arg1) Then Begin
   If ExtractExt(Arg1) = '.XYZ' Then Begin
    NInp := 0;
    EdAppend('# Geometrie importee de ' + BaseName(Arg1));
    EdAppend('! HF STO-3G TightSCF');
    EdAppend('* xyzfile 0 1 ' + Arg1);
    EdAppend('*');
    InpName := '';
    PreviewInput;
   End
   Else If LoadInputFile(Arg1) Then Begin
    EnsureOneLine;
    PreviewInput;
    StatusMsg := 'Fichier ouvert : ' + Arg1;
   End;
  End
  Else StatusMsg := 'Fichier introuvable : ' + Arg1;
 End;
 InpModified := False;
 MainLoop;
 CloseGraph;
END.
