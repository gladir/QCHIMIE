{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2026
  @version: 1.00
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal 7, Free Pascal)
  @description: Clone de la ligne de commande CP2K - logiciel de chimie
                quantique et de physique de l'etat solide pour les
                simulations atomistiques (molecules, solides, liquides,
                systemes periodiques).

  Ligne de commande :
    CP2K [options] [-i] fichier.inp [-o sortie.out]

  Options (comme CP2K) :
    -i fichier           fichier d'entree (ou premier argument sans option)
    -o fichier           fichier de sortie (sinon sortie standard)
    -c, --check          verifie seulement la syntaxe du fichier d'entree
    -e, --echo           affiche l'entree apres pretraitement et analyse
    -E, --set nom=valeur definit une variable du preprocesseur (@SET)
    -r, --run            lance le calcul (comportement par defaut)
    -v, --version        affiche la version
    --license            affiche la licence
    --xml                decrit les sections et mots-cles reconnus (XML)
    -h, --help           affiche l'aide

  Fichier d'entree au format CP2K :
    - sections &NOM ... &END NOM, mots-cles, commentaires # et !
    - preprocesseur : @SET, @IF / @ELSE / @ENDIF, @INCLUDE, @PRINT,
      $VARIABLE, unites entre crochets ([angstrom], [fs], [kcalmol]...)
    - &GLOBAL (PROJECT, RUN_TYPE, PRINT_LEVEL, SEED)
    - &FORCE_EVAL (METHOD FIST | QS) avec &SUBSYS (&CELL, &COORD, &KIND,
      &TOPOLOGY), &MM (&FORCEFIELD, &POISSON/&EWALD) et &DFT (&QS, &XC,
      &SCF)
    - &MOTION (&GEO_OPT, &CELL_OPT, &MD avec &THERMOSTAT et &BAROSTAT,
      &PRINT), &VIBRATIONAL_ANALYSIS, &DEBUG, &EXT_RESTART

  Types d'execution (RUN_TYPE) :
    ENERGY, ENERGY_FORCE, GEO_OPT, CELL_OPT, MD, VIBRATIONAL_ANALYSIS,
    DEBUG (et leurs synonymes)

  Moteurs de calcul :
    - FIST : champ de force classique (Lennard-Jones, Williams, charges
      ponctuelles avec sommation d'Ewald, liaisons, angles et torsions),
      forces et tenseur de contraintes analytiques
    - QS (Quickstep) : Hartree-Fock et DFT (LDA, BLYP, B3LYP, PBE, PBE0)
      tout-electron avec bases gaussiennes STO-3G (H a Ne), 3-21G et
      6-31G (H, He et C a F), forces par differences finies
    - Dynamique moleculaire NVE, NVT (Nose-Hoover, CSVR) et NPT_I,
      optimisation de geometrie et de cellule (BFGS), analyse
      vibrationnelle, fichiers .xyz, .ener, .restart et .cell

  Aucune dependance externe : tout le code est dans ce fichier.
}

Program CP2K;

{$IFDEF FPC}
 {$MODE TP}
 {$Q-}
 {$R-}
{$ELSE}
 {$M 60000,0,655360}
{$ENDIF}

Uses DOS;

{ --------------------------------------------------------------------- }
{ Constantes generales                                                  }
{ --------------------------------------------------------------------- }

Const
 PROG_TITLE      = 'CP2K Clone 1.00 - QChimie';
 PROG_VERSION    = '1.00';
 CP2K_VERSION    = '2024.1';

 {$IFDEF FPC}
 MAX_BF          = 40;           { fonctions de base }
 MAX_AT          = 32;           { atomes du moteur quantique }
 MAX_CHUNK       = 84;           { blocs d'integrales bielectroniques }
 MAX_HIST        = 300;
 CIS_MAX         = 40;
 MAX_ATOM        = 4000;         { atomes du systeme }
 MAX_KIND        = 32;           { sortes d'atomes }
 MAX_SEC         = 500;          { sections du fichier d'entree }
 MAX_KW          = 4000;         { mots-cles }
 MAX_LINES       = 6000;         { lignes du fichier d'entree }
 MAX_BONDS       = 8000;
 MAX_BENDS       = 12000;
 MAX_TORS        = 20000;
 MAX_EXCL        = 40;
 MAX_PTAB        = 60;           { parametres de champ de force }
 MAX_KVEC        = 6000;
 MAX_VIB         = 120;          { atomes pour l'analyse vibrationnelle }
 {$ELSE}
 MAX_BF          = 14;
 MAX_AT          = 8;
 MAX_CHUNK       = 3;
 MAX_HIST        = 150;
 CIS_MAX         = 10;
 MAX_ATOM        = 40;
 MAX_KIND        = 12;
 MAX_SEC         = 90;
 MAX_KW          = 400;
 MAX_LINES       = 500;
 MAX_BONDS       = 60;
 MAX_BENDS       = 100;
 MAX_TORS        = 120;
 MAX_EXCL        = 20;
 MAX_PTAB        = 16;
 MAX_KVEC        = 300;
 MAX_VIB         = 12;
 {$ENDIF}

 MAX_PAIRS       = MAX_BF * (MAX_BF + 1) Div 2;
 CHUNK_SIZE      = 4096;
 MAX_PRIM        = 6;
 MAX_SHELL       = 8;
 MAX_LINE        = 100;
 MAX_CONS        = 3 * MAX_AT;

 { constantes physiques (CODATA 2018, comme CP2K) }
 BOHR_ANG        = 0.529177210903;       { angstrom par bohr }
 ANG2BOHR        = 1.8897261246257;
 HARTREE_EV      = 27.211386245988;
 HARTREE_KCAL    = 627.5094740631;
 HARTREE_KJ      = 2625.4996394799;
 HARTREE_CM      = 219474.6313632;
 BOHR_DEBYE      = 2.541746473;          { u.a. -> Debye }
 KB_AU           = 3.1668115634556E-6;   { Eh par kelvin }
 AMU_AU          = 1822.888486209;       { u.m.a. en masses electroniques }
 FS_AU           = 41.341373336493;      { unites atomiques de temps par fs }
 PRESS_AU_PA     = 2.942101569713E13;    { pascal par unite atomique }
 BAR_AU          = 1.0E5 / PRESS_AU_PA;
 GPA_AU          = 1.0E9 / PRESS_AU_PA;
 FREQ_FAC        = 219474.6313632;       { Eh -> cm-1 }

 { Methodes du moteur quantique }
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

 J_SP            = 0;

 { Types d'execution }
 RT_NONE         = 0;
 RT_ENERGY       = 1;
 RT_FORCE        = 2;
 RT_GEOOPT       = 3;
 RT_CELLOPT      = 4;
 RT_MD           = 5;
 RT_VIB          = 6;
 RT_DEBUG        = 7;

 { Methodes de calcul des forces }
 FE_FIST         = 1;
 FE_QS           = 2;

 { Niveaux d'impression }
 PL_SILENT       = 0;
 PL_LOW          = 1;
 PL_MEDIUM       = 2;
 PL_HIGH         = 3;
 PL_DEBUG        = 4;

{ --------------------------------------------------------------------- }
{ Types                                                                 }
{ --------------------------------------------------------------------- }

Type
 TVec3 = Array[1..3] Of Real;
 TMat3 = Array[1..3, 1..3] Of Real;

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

 TXYZ = Array[1..MAX_ATOM, 1..3] Of Real;

 { mot-cle du fichier d'entree }
 TKw = Record
  U, O : String[48];     { nom en majuscules / nom d'origine }
  V    : String[250];    { reste de la ligne }
  Sec  : Integer;
  Line : Integer;
 End;
 PKw = ^TKw;

 PStr = ^String;

{ --------------------------------------------------------------------- }
{ Variables globales (etat general)                                     }
{ --------------------------------------------------------------------- }

Var
 UserAbort   : Boolean;
 OutF        : Text;
 OutOpen     : Boolean;
 NWarn       : Integer;
 PrintLevel  : Integer;
 T0Run       : LongInt;
 StartDate   : String;
 InputName   : String;
 OutputName  : String;
 ProjectName : String;
 CurDir      : String;
 WallTime    : Real;          { duree maximale en secondes (0 = illimitee) }

Procedure ShowProgress(Const Msg: String; Pct: Integer); Forward;

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
 If (Code = 0) And (R >= -32000) And (R <= 32000) Then Begin V := R; ParseInt :=
   True; End
 Else ParseInt := False;
End;

Function ExtractExt(Const FN: String): String;
Var I: Integer;
Begin
 I := Length(FN);
 While (I > 0) And (FN[I] <> '.') And (FN[I] <> '\') And (FN[I] <> '/') Do
   Dec(I);
 If (I > 0) And (FN[I] = '.') Then ExtractExt := UpStr(Copy(FN, I, 255))
 Else ExtractExt := '';
End;

Function ChangeExt(Const FN, Ext: String): String;
Var I: Integer;
Begin
 I := Length(FN);
 While (I > 0) And (FN[I] <> '.') And (FN[I] <> '\') And (FN[I] <> '/') Do
   Dec(I);
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

{ --------------------------------------------------------------------- }
{ Utilitaires supplementaires                                           }
{ --------------------------------------------------------------------- }

Function Rep(C: Char; N: Integer): String;
Var S: String; I: Integer;
Begin
 S := '';
 For I := 1 To N Do S := S + C;
 Rep := S;
End;

Function Pw10(E: Integer): Real;
Begin
 Pw10 := Exp(E * Ln(10.0));
End;

{ Notation scientifique a la CP2K : 1.234E-05 }
Function EStr(V: Real; Nd: Integer): String;
Var E: Integer; M: Real; S, Sg, Es: String;
Begin
 Sg := '';
 E := 0;
 M := 0;
 If V <> 0 Then Begin
  If V < 0 Then Begin Sg := '-'; V := -V; End;
  E := Trunc(Ln(V) / Ln(10.0));
  M := V / Pw10(E);
  While M >= 10 Do Begin M := M / 10; Inc(E); End;
  While M < 1 Do Begin M := M * 10; Dec(E); End;
  S := RStr(M, Nd);
  If Copy(S, 1, 2) = '10' Then Begin M := 1; Inc(E); End;
 End;
 S := RStr(M, Nd);
 If E < 0 Then Begin Es := '-'; E := -E; End Else Es := '+';
 If E < 10 Then Es := Es + '0' + IStr(E) Else Es := Es + IStr(E);
 EStr := Sg + S + 'E' + Es;
End;

{ Reel en largeur W avec Dec decimales (format Fw.d) }
Function FR(V: Real; W, Nd: Integer): String;
Var S: String;
Begin
 If Abs(V) >= 1E15 Then Begin FR := PadL(EStr(V, 6), W); Exit; End;
 S := RStr(V, Nd);
 If S = '-0' Then S := '0';
 If (Nd > 0) And (Copy(S, 1, 2) = '-0') Then Begin
  { evite -0.000 }
  If Abs(V) < 0.5 * Pw10(-Nd) Then Delete(S, 1, 1);
 End;
 FR := PadL(S, W);
End;

Function Cu(X: Real): Real;
Begin
 Cu := X * X * X;
End;

{ Partie entiere par defaut (plancher) }
Function Fl(X: Real): Real;
Begin
 If X < Int(X) Then Fl := Int(X) - 1 Else Fl := Int(X);
End;

Function ATan2(Y, X: Real): Real;
Begin
 If X > 0 Then ATan2 := ArcTan(Y / X)
 Else If X < 0 Then Begin
  If Y >= 0 Then ATan2 := ArcTan(Y / X) + Pi Else ATan2 := ArcTan(Y / X) - Pi;
 End
 Else Begin
  If Y > 0 Then ATan2 := Pi / 2
  Else If Y < 0 Then ATan2 := -Pi / 2
  Else ATan2 := 0;
 End;
End;

Function ISt(V: LongInt; W: Integer): String;
Begin
 ISt := PadL(IStr(V), W);
End;

Function BoolStr(B: Boolean): String;
Begin
 If B Then BoolStr := 'T' Else BoolStr := 'F';
End;

Function YesNo(B: Boolean): String;
Begin
 If B Then YesNo := 'YES' Else YesNo := 'NO';
End;

Function Z2(V: LongInt): String;
Begin
 If V < 10 Then Z2 := '0' + IStr(V) Else Z2 := IStr(V);
End;

Function Z3(V: LongInt): String;
Begin
 Z3 := Rep('0', 3 - Length(IStr(V))) + IStr(V);
End;

{ Date et heure : AAAA-MM-JJ HH:MM:SS.mmm }
Function DateTimeStr: String;
Var Y, M, D, DW, H, Mi, S, S100: Word;
Begin
 GetDate(Y, M, D, DW);
 GetTime(H, Mi, S, S100);
 DateTimeStr := IStr(Y) + '-' + Z2(M) + '-' + Z2(D) + ' ' + Z2(H) + ':' +
                Z2(Mi) + ':' + Z2(S) + '.' + Z3(S100 * 10);
End;

{ Temps ecoule depuis le debut, en secondes }
Function RunSec: Real;
Begin
 RunSec := ElapsedMs(T0Run) / 1000.0;
End;

{ --------------------------------------------------------------------- }
{ Sortie                                                                }
{ --------------------------------------------------------------------- }

Var
 RedirOn : Boolean;
 RedirF  : Text;

Procedure Emit(Const S: String);
Begin
 If RedirOn Then WriteLn(RedirF, S)
 Else If OutOpen Then WriteLn(OutF, S);
End;

Procedure EmitBlank;
Begin
 Emit('');
End;

Procedure ShowProgress(Const Msg: String; Pct: Integer);
Begin
End;

{ Message sur la console (erreur standard) }
Procedure ConMsg(Const S: String);
Begin
 {$IFDEF FPC}
 WriteLn(StdErr, S);
 {$ELSE}
 WriteLn(S);
 {$ENDIF}
End;

{ Ligne d'information : prefixe| libelle ................ valeur }
Procedure Info(Const Pfx, Lbl, Val: String);
Var S: String; N: Integer;
Begin
 S := ' ' + Pfx + '| ' + Lbl;
 N := 80 - Length(S) - Length(Val);
 If N < 1 Then N := 1;
 Emit(S + Rep(' ', N) + Val);
End;

Function Stars: String;
Begin
 Stars := ' ' + Rep('*', 79);
End;

{ Texte centre dans une boite *** ... *** }
Procedure BoxText(Const T: String);
Var N, L: Integer; S: String;
Begin
 L := Length(T);
 N := (73 - L) Div 2;
 If N < 0 Then N := 0;
 S := Rep(' ', N) + T;
 S := S + Rep(' ', 73 - Length(S));
 Emit(' ***' + Copy(S, 1, 73) + '***');
End;

Procedure BoxBanner(Const T1, T2: String);
Begin
 Emit(Stars);
 BoxText(T1);
 If T2 <> '' Then BoxText(T2);
 Emit(Stars);
End;

Procedure Warn(Const Msg: String);
Begin
 Inc(NWarn);
 Emit(' *** WARNING in CP2K.PAS :: ' + Msg + ' ***');
 If OutOpen Then ConMsg('AVERTISSEMENT : ' + Msg);
End;

{ Ligne de l'encadre ABORT }
Procedure AbortLine(Const T: String);
Var S: String;
Begin
 S := T + Rep(' ', 77);
 Emit(' *' + Copy(S, 1, 77) + '*');
End;

{ Arret avec l'encadre ABORT de CP2K }
Procedure Abort(Const Msg: String);
Var M, L: String; N, NL: Integer; First: Boolean;
Begin
 If OutOpen Then Begin
  EmitBlank;
  Emit(' ' + Rep('*', 79));
  AbortLine('   ___');
  AbortLine('  /   \');
  AbortLine(' [ABORT]');
  M := Msg;
  First := True;
  NL := 0;
  While M <> '' Do Begin
   N := 62;
   If Length(M) > N Then Begin
    While (N > 20) And (M[N] <> ' ') Do Dec(N);
    If M[N] <> ' ' Then N := 62;
   End;
   L := Copy(M, 1, N);
   M := Copy(M, N + 1, 255);
   While (M <> '') And (M[1] = ' ') Do Delete(M, 1, 1);
   If First Then AbortLine('  \___/   ' + L) Else AbortLine('    |     ' + L);
   First := False;
   Inc(NL);
  End;
  If NL = 1 Then AbortLine('    |');
  AbortLine('  O/|');
  AbortLine(' /| |');
  AbortLine(' / \' + Rep(' ', 62) + 'CP2K.PAS');
  Emit(' ' + Rep('*', 79));
  EmitBlank;
  EmitBlank;
  If RedirOn Then Begin Close(RedirF); RedirOn := False; End;
  Close(OutF);
  OutOpen := False;
 End;
 ConMsg('ERREUR : ' + Msg);
 Halt(1);
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
 AtBas   : Array[1..MAX_AT] Of Integer;
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
   Sh6(0, 3047.5249000, 0.0018347, 457.3695100, 0.0140373, 103.9486900,
     0.0688426,
          29.2101550, 0.2321844, 9.2866630, 0.4679413, 3.1639270, 0.3623120);
   SP3(7.8682724, -0.1193324, 0.0689991, 1.8812885, -0.1608542, 0.3164240,
       0.5442493, 1.1434564, 0.7443083);
   SP1(0.1687144);
  End;
  7: Begin
   Sh6(0, 4173.5110000, 0.0018348, 627.4579000, 0.0139950, 142.9021000,
     0.0685870,
          40.2343300, 0.2322410, 12.8202100, 0.4690700, 4.3904370, 0.3604550);
   SP3(11.6263580, -0.1149610, 0.0675800, 2.7162800, -0.1691180, 0.3239070,
       0.7722180, 1.1458520, 0.7408950);
   SP1(0.2120313);
  End;
  8: Begin
   Sh6(0, 5484.6717000, 0.0018311, 825.2349500, 0.0139501, 188.0469600,
     0.0684451,
          52.9645000, 0.2327143, 16.8975700, 0.4701930, 5.7996353, 0.3585209);
   SP3(15.5396160, -0.1107775, 0.0708743, 3.5999336, -0.1480263, 0.3397528,
       1.0137618, 1.1307670, 0.7271586);
   SP1(0.2700058);
  End;
  9: Begin
   Sh6(0, 7001.7130000, 0.0018196169, 1051.3660000, 0.01391608, 239.2857000,
     0.06840532,
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

{ Table de Boys : F_m(T) pour T = i * BOYS_DT, m = 0..BOYS_NM ; interpolation de
  Taylor }
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
Var K, I: Integer; X: Real;
Begin
 If T >= BOYS_N * BOYS_DT - 0.5 Then Begin
  F[0] := 0.5 * Sqrt(Pi / T);
  For K := 0 To M - 1 Do F[K + 1] := (2 * K + 1) * F[K] / (2 * T);
  Exit;
 End;
 I := Round(T / BOYS_DT);
 X := I * BOYS_DT - T;
 For K := 0 To M Do
  F[K] := BoysTab^[I, K] + X * (BoysTab^[I, K + 1] + X / 2 * (BoysTab^[I, K + 2]
    +
          X / 3 * (BoysTab^[I, K + 3] + X / 4 * (BoysTab^[I, K + 4] + X / 5 *
            (BoysTab^[I, K + 5] +
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
 R := -2 * B * B * Ovl1D(I, J + 2, Q, A, B) + B * (2 * J + 1) * Ovl1D(I, J, Q,
   A, B);
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
   S := S + BFC[M, K] * BFC[N, L] * (TX * SY * SZ + SX * TY * SZ + SX * SY *
     TZ);
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
  If Not GetShells(AtBas[A], AtZ[A]) Then Begin
   Err := 'La base ' + BasisNames[AtBas[A]] +
     ' n''est pas definie pour l''element ' +
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
   PairL[Idx] := BFL[M,1] + BFL[M,2] + BFL[M,3] + BFL[N,1] + BFL[N,2] +
     BFL[N,3];
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
        If Abs(PP^.E[T * 9 + U * 3 + V]) > PP^.Cm Then PP^.Cm := Abs(PP^.E[T * 9
          + U * 3 + V]);
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
   If PP1^.Cm * PP2^.Cm * TWO_PI_52 / (P * Q * Sqrt(P + Q)) < 1E-13 Then
     Continue;
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
   If I = J Then V^[(I - 1) * LDL + J - 1] := 1 Else V^[(I - 1) * LDL + J - 1]
     := 0;
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
 FrozenCore : Boolean;

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
 SO, VN, DX, DY, DZ, PF, Sm: Real;
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
    BuildR(L, PP^.P, PP^.PX - AtX[C, 1], PP^.PY - AtX[C, 2], PP^.PZ - AtX[C,
      3]);
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
Var I, J: Integer; Exc, Ec1, Ek: Real; NeedK: Boolean;
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

Var ScfTi: LongInt;

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
 If UseGuess And OrbsValid And (OrbN = NBF) And (OrbOA = NOccA) And (OrbOB =
   NOccB) Then
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
  Emit(' SCF WAVEFUNCTION OPTIMIZATION');
  EmitBlank;
  Emit('  Step     Update method      Time    Convergence         Total ' +
    'energy    Change');
  Emit('  ' + Rep('-', 77));
  ScfTi := NowMs;
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
   For I := 1 To DiisN Do Begin Bm[I, DiisN + 1] := -1; Bm[DiisN + 1, I] := -1;
     Rhs[I] := 0; End;
   Bm[DiisN + 1, DiisN + 1] := 0; Rhs[DiisN + 1] := -1;
   If SolveSmall(DiisN + 1, Bm, Rhs) Then Begin
    MatZero(MFa);
    If Not Closed Then MatZero(MFb);
    For K := 1 To DiisN Do
     For I := 1 To NBF Do
      For J := 1 To NBF Do Begin
       MFa^[I, J] := MFa^[I, J] + Rhs[K] * DiisF[K, 1]^[I, J];
       If Not Closed Then MFb^[I, J] := MFb^[I, J] + Rhs[K] * DiisF[K, 2]^[I,
         J];
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
  If Verbose Then Begin
   Emit(ISt(Iter, 6) + ' ' + PadR('DIIS/Diag.', 11) + ' ' + PadL(EStr(0.4, 2),
     8) + FR(ElapsedMs(ScfTi) / 1000.0, 7, 1) +
        FR(ErrMax, 15, 8) + FR(E, 21, 10) + PadL(EStr(DE, 2), 10));
   ScfTi := NowMs;
  End;
  If (Iter > 1) And (Abs(DE) < ScfTolE) And (RmsD < ScfTolD) And (ErrMax <
    ScfTolErr) Then Done := True;
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
   Emit(' *** SCF run converged in ' + ISt(Iter, 4) + ' steps ***');
  End;
 End
 Else Begin
  If Verbose Then Begin
   EmitBlank;
   Emit(' *** SCF run NOT converged ***');
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
   While (R < 40) And (Abs(BFC[I, K]) * (1 + R) * ExpS(-BFA[I, K] * R * R) >
     1E-11) Do R := R + 0.2;
   If R * R > BFCut2[I] Then BFCut2[I] := R * R;
  End;
 End;
End;

Function GridNRad: Integer;
Begin
 GridNRad := GRID_NR + (GridScale - 1) * (GRID_NR Div 2);
End;

Function GridNTheta: Integer;
Var N: Integer;
Begin
 N := GRID_NTH + (GridScale - 1) * 4;
 If N > GRID_NTH + 8 Then N := GRID_NTH + 8;
 GridNTheta := N;
End;

Procedure BuildGrid;
Var
 A, I, J, K, NR, NTH, NPH: Integer;
 XL: Array[1..GRID_NTH + 8] Of Real;
 WL: Array[1..GRID_NTH + 8] Of Real;
 Xi, X, R, WR, ST, CT, Ph, PX, PY, PZ, WA, W, Z0, Z1, PP, P1, P2, P3: Real;
Begin
 FreeGrid;
 GridFull := False;
 NR := GridNRad;
 NTH := GridNTheta;
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
 M, K: Integer; DX, DY, DZ, R2, E, PX, PY, PZ, QX, QY, QZ, A: Real;
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
          - B * X0 / X0X * (Ln(Sqr(X - X0) / XX) + 2 * (B + 2 * X0) / Q *
            ArcTan(Q / (2 * X + B))));
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
 H := Gam * Phi * Phi * Phi * Ln(1 + Beta / Gam * T2 * (1 + A * T2) / (1 + A *
   T2 + A * A * T2 * T2));
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
 T2 := -LA * LB * Om * (Ra * Rb * (Pw(2, 11 / 3) * CF * (Pw(Ra, 8 / 3) + Pw(Rb,
   8 / 3))
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
  F_B3LYP: ExcDens := 0.8 * ExSlater(Ra, Rb) + 0.72 * (ExB88Spin(Ra, Saa) +
    ExB88Spin(Rb, Sbb))
                      + 0.19 * EcVWN(Ra, Rb) + 0.81 * EcLYP(Ra, Rb, Saa, Sab,
                        Sbb);
  F_PBE:   ExcDens := ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb) + EcPBE(Ra, Rb,
    Saa, Sab, Sbb);
  F_PBE0:  ExcDens := 0.75 * (ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb)) +
    EcPBE(Ra, Rb, Saa, Sab, Sbb);
 Else ExcDens := 0;
 End;
End;

{ Parties echange / correlation separees (pour l'impression) }
Function ExPart(Ra, Rb, Saa, Sbb: Real): Real;
Begin
 Case Func Of
  F_LDA:   ExPart := ExSlater(Ra, Rb);
  F_BLYP:  ExPart := ExSlater(Ra, Rb) + ExB88Spin(Ra, Saa) + ExB88Spin(Rb, Sbb);
  F_B3LYP: ExPart := 0.8 * ExSlater(Ra, Rb) + 0.72 * (ExB88Spin(Ra, Saa) +
    ExB88Spin(Rb, Sbb));
  F_PBE:   ExPart := ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb);
  F_PBE0:  ExPart := 0.75 * (ExPBESpin(Ra, Saa) + ExPBESpin(Rb, Sbb));
 Else ExPart := 0;
 End;
End;

{ Derivee numerique d'une variable par differences centrees }
Procedure XCDeriv(Ra, Rb, Saa, Sab, Sbb: Real; Var E, VRa, VRb, VSaa, VSab,
  VSbb: Real);
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
            - B * X0 / X0X * (2 * Ln(Abs(X - X0)) - LXX + 2 * (B + 2 * X0) / Q *
              At));
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

{ Echange PBE pour un spin (Ra = rho_s, X43 = rho_s^(4/3), Kf1 = (3 pi^2 2
  rho_s)^(1/3)) }
Function PBExc(Ra, X43, Kf1, S: Real): Real;
Var N43, SS, FX, EL: Real;
Begin
 N43 := 2.5198420997897464 * X43;       { (2 rho_s)^(4/3) }
 SS := 2 * Sqrt(S) / (2 * Kf1 * 2 * Ra);
 FX := 1.804 - 0.804 / (1 + 0.2195149727645171 * SS * SS / 0.804);
 EL := -0.9305257363491 * N43 / 1.2599210498948732;
 PBExc := 0.5 * EL * FX;
End;

{ Correlation PBE a couche fermee (Rs, Kf = (3 pi^2 rho)^(1/3), ST = |grad
  rho|^2) }
Function PBEcc(R, Rs, Kf, ST: Real): Real;
Var Eps, KS, T2, A, H: Real;
Begin
 Eps := PW92G(Rs, 0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294);
 KS := Sqrt(4 * Kf / Pi);
 If ST < 0 Then ST := 0;
 T2 := ST / (4 * KS * KS * R * R);
 A := PBE_BETA / PBE_GAM / (Exp(-Eps / PBE_GAM) - 1);
 H := PBE_GAM * Ln(1 + PBE_BETA / PBE_GAM * T2 * (1 + A * T2) / (1 + A * T2 + A
   * A * T2 * T2));
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
  F_BLYP:  E := -2 * CX_LDA * X43 + 2 * B88c(X43, Saa) + LYPc(Ra, X43, Rm13,
    Saa, Sab);
  F_B3LYP: E := -1.6 * CX_LDA * X43 + 1.44 * B88c(X43, Saa)
                + 0.19 * 2 * Ra * EpsVWNp(Rs) + 0.81 * LYPc(Ra, X43, Rm13, Saa,
                  Sab);
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
Var X43, Kf: Real;
Begin
 If Ra < 1E-20 Then Begin ExPartC := 0; Exit; End;
 X43 := Pw(Ra, 4 / 3);
 Case Func Of
  F_LDA:   ExPartC := -2 * CX_LDA * X43;
  F_BLYP:  ExPartC := -2 * CX_LDA * X43 + 2 * B88c(X43, Saa);
  F_B3LYP: ExPartC := -1.6 * CX_LDA * X43 + 1.44 * B88c(X43, Saa);
  F_PBE, F_PBE0: Begin
   Kf := Pw(3 * Pi * Pi * 2 * Ra, 1 / 3);
   If Func = F_PBE Then ExPartC := 2 * PBExc(Ra, X43, Kf, Saa) Else ExPartC :=
     1.5 * PBExc(Ra, X43, Kf, Saa);
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
   Gax := Gax + 2 * TA * GX[I]; Gay := Gay + 2 * TA * GY[I]; Gaz := Gaz + 2 * TA
     * GZ[I];
   If Not Closed Then Begin
    Rb := Rb + TB * Phi[I];
    Gbx := Gbx + 2 * TB * GX[I]; Gby := Gby + 2 * TB * GY[I]; Gbz := Gbz + 2 *
      TB * GZ[I];
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
  If Closed Then EXdft := EXdft + W * ExPartC(Ra, Saa) Else EXdft := EXdft + W *
    ExPart(Ra, Rb, Saa, Sbb);
  { vecteurs de gradient effectifs }
  QAx := 2 * VSaa * Gax + VSab * Gbx; QAy := 2 * VSaa * Gay + VSab * Gby; QAz :=
    2 * VSaa * Gaz + VSab * Gbz;
  QBx := 2 * VSbb * Gbx + VSab * Gax; QBy := 2 * VSbb * Gby + VSab * Gay;
  QBz := 2 * VSbb * Gbz + VSab * Gaz;
  For I := 1 To NL Do Begin
   AA[I] := W * (0.5 * VRa * Phi[I] + QAx * GX[I] + QAy * GY[I] + QAz * GZ[I]);
   If Not Closed Then AB[I] := W * (0.5 * VRb * Phi[I] + QBx * GX[I] + QBy *
     GY[I] + QBz * GZ[I]);
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
 OrbMullPct : Array[1..MAX_BF, 1..MAX_AT] Of Real;   { contribution des atomes
   aux MOs }

Procedure PopAnalysis;
Var A, B, I, J, K: Integer;
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
   OrbMullPct[K, BFAtom[I]] := OrbMullPct[K, BFAtom[I]] + MCa^[I, K] * MT1^[I,
     K];
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
{ Evaluation de l'energie pour la geometrie courante                     }
{ --------------------------------------------------------------------- }

Function EvalEnergy(Verbose, UseGuess: Boolean; Var E: Real; Var Err: String):
  Boolean;
Begin
 EvalEnergy := False;
 Err := '';
 InitBoysTable;
 BuildPairs;
 BuildOneElectron;
 If Not ComputeERI Then Begin Err :=
   'Memoire insuffisante pour les integrales bielectroniques'; Exit; End;
 If Not BuildOrtho Then Begin Err := 'Dependance lineaire dans la base'; Exit;
   End;
 Enuc := NucRepulsion;
 If Meth = M_DFT Then Begin
  BuildGrid;
  If GridFull Then Begin Err :=
    'Grille DFT trop grande pour cette molecule (utilisez la grille normale)';
    Exit; End;
 End;
 If Not ScfRun(Verbose, UseGuess) Then Begin
  If UserAbort Then Err := 'Calcul interrompu' Else Err :=
    'Le SCF n''a pas converge';
  Exit;
 End;
 ETot := EScf;
 E := ETot;
 EvalEnergy := True;
End;

{ --------------------------------------------------------------------- }
{ Coordonnees internes (contraintes)                                     }
{ --------------------------------------------------------------------- }

{ --------------------------------------------------------------------- }
{ Unites                                                                }
{ --------------------------------------------------------------------- }

Const
 NUNIT = 50;
 UName : Array[1..NUNIT] Of String[14] = (
  'hartree', 'ha', 'au_e', 'ev', 'kcalmol', 'kjmol', 'ry', 'rydberg',
  'k_e', 'j', 'joule', 'wavenumber_e', 'bohr', 'au_l', 'angstrom', 'nm',
  'pm', 'm', 'rad', 'radian', 'deg', 'degree', 'au_t', 'fs', 'ps', 's',
  'ns', 'au_m', 'amu', 'kg', 'g', 'bar', 'atm', 'pa', 'kpa', 'mpa', 'gpa',
  'kbar', 'au_p', 'k', 'e', 'au_q', 'au', 'femtosecond', 'picosecond',
  'second', 'mhartree', 'mev', 'mbar', 'kelvin');
 UFac : Array[1..NUNIT] Of Real = (
  1.0, 1.0, 1.0, 1.0 / HARTREE_EV, 1.0 / HARTREE_KCAL, 1.0 / HARTREE_KJ,
  0.5, 0.5, KB_AU, 2.2937122783963E17, 2.2937122783963E17,
  1.0 / HARTREE_CM, 1.0, 1.0, ANG2BOHR, 10.0 * ANG2BOHR, 0.01 * ANG2BOHR,
  1.0E10 * ANG2BOHR, 1.0, 1.0, 0.017453292519943, 0.017453292519943, 1.0, FS_AU,
  1000.0 * FS_AU, 1.0E15 * FS_AU, 1.0E6 * FS_AU, 1.0, AMU_AU,
  1.0977691057577E30, 1.0977691057577E27, BAR_AU, 1.01325 * BAR_AU,
  1.0 / PRESS_AU_PA, 1000.0 / PRESS_AU_PA, 1.0E6 / PRESS_AU_PA, GPA_AU,
  1000.0 * BAR_AU, 1.0, 1.0, 1.0, 1.0, 1.0, FS_AU, 1000.0 * FS_AU,
  1.0E15 * FS_AU, 0.001, 0.001 / HARTREE_EV, 1.0E-3 * BAR_AU, 1.0);

Function UnitAtom(Const A: String; Var Ok: Boolean): Real;
Var I: Integer;
Begin
 For I := 1 To NUNIT Do
  If UName[I] = A Then Begin UnitAtom := UFac[I]; Exit; End;
 Ok := False;
 UnitAtom := 1;
End;

{ Facteur de conversion d'une unite composee (ex. kcalmol*angstrom^-2)
  vers les unites atomiques }
Function UnitToAU(Const U0: String; Var Ok: Boolean): Real;
Var
 U, Atom: String;
 I, Code: Integer;
 R, V: Real;
 Mul: Boolean;
 Function CloseAtom: Real;
 Var P, E, K: Integer; Nm: String; X: Real;
 Begin
  Nm := Atom; E := 1;
  P := Pos('^', Nm);
  If P > 0 Then Begin
   Val(Copy(Nm, P + 1, 10), E, Code);
   If Code <> 0 Then Begin Ok := False; E := 1; End;
   Nm := Copy(Nm, 1, P - 1);
  End;
  If Nm = '' Then X := 1 Else X := UnitAtom(Nm, Ok);
  If Not Mul Then E := -E;
  { exponent entier }
  V := 1;
  If E >= 0 Then Begin
   For K := 1 To E Do V := V * X;
  End
  Else Begin
   For K := 1 To -E Do V := V / X;
  End;
  CloseAtom := V;
 End;
Begin
 Ok := True;
 U := LoStr(TrimStr(U0));
 If U = '' Then Begin UnitToAU := 1; Exit; End;
 { cas particulier : unites simples avec exposant dans la table }
 For I := 1 To NUNIT Do
  If UName[I] = U Then Begin UnitToAU := UFac[I]; Exit; End;
 R := 1; Atom := ''; Mul := True;
 For I := 1 To Length(U) Do Begin
  If (U[I] = '*') Or (U[I] = '/') Then Begin
   R := R * CloseAtom;
   Atom := '';
   Mul := (U[I] = '*');
  End
  Else Atom := Atom + U[I];
 End;
 R := R * CloseAtom;
 UnitToAU := R;
End;

{ --------------------------------------------------------------------- }
{ Analyse lexicale                                                      }
{ --------------------------------------------------------------------- }

{ Extrait le premier mot de S (les guillemets sont retires) }
Function GetTok(Var S: String): String;
Var I, N: Integer; Q: Char; T: String;
Begin
 N := Length(S);
 I := 1;
 While (I <= N) And (S[I] <= ' ') Do Inc(I);
 If I > N Then Begin GetTok := ''; S := ''; Exit; End;
 T := '';
 If (S[I] = '"') Or (S[I] = '''') Then Begin
  Q := S[I];
  Inc(I);
  While (I <= N) And (S[I] <> Q) Do Begin T := T + S[I]; Inc(I); End;
  Inc(I);
 End
 Else Begin
  While (I <= N) And (S[I] > ' ') Do Begin T := T + S[I]; Inc(I); End;
 End;
 S := Copy(S, I, 255);
 GetTok := T;
End;

{ Retire un commentaire (# ou !) hors des guillemets }
Function StripComment(Const S: String): String;
Var I: Integer; Q: Char;
Begin
 Q := #0;
 For I := 1 To Length(S) Do Begin
  If Q <> #0 Then Begin
   If S[I] = Q Then Q := #0;
  End
  Else If ((S[I] = '"') Or (S[I] = '''')) And ((I = 1) Or (S[I - 1] <= ' '))
    Then Q := S[I]
  Else If (S[I] = '#') Or ((S[I] = '!') And Not ((I < Length(S)) And (S[I + 1] =
    '='))) Then Begin
   StripComment := TrimStr(Copy(S, 1, I - 1));
   Exit;
  End;
 End;
 StripComment := TrimStr(S);
End;

{ Lecture d'un reel : accepte l'exposant Fortran D }
Function ReadReal(Const S: String; Var V: Real): Boolean;
Var T: String; I, Code: Integer; R: Real;
Begin
 T := TrimStr(S);
 For I := 1 To Length(T) Do Begin
  If T[I] = 'D' Then T[I] := 'E';
  If T[I] = 'd' Then T[I] := 'e';
 End;
 ReadReal := False;
 If T = '' Then Exit;
 If (T[1] = '.') Then T := '0' + T;
 If (Length(T) > 1) And (T[1] = '-') And (T[2] = '.') Then T := '-0' + Copy(T,
   2, 255);
 If (Length(T) > 1) And (T[1] = '+') Then T := Copy(T, 2, 255);
 { refuse nan, inf et les valeurs hors limites }
 If Not ((T[1] In ['0'..'9', '.']) Or ((Length(T) > 1) And (T[1] = '-') And
   (T[2] In ['0'..'9', '.']))) Then Exit;
 Val(T, R, Code);
 If (Code = 0) And (R = R) And (Abs(R) < 1E300) Then Begin V := R; ReadReal :=
   True; End;
End;

Function ReadInt(Const S: String; Var V: LongInt): Boolean;
Var T: String; Code: Integer; R: LongInt; X: Real;
Begin
 T := TrimStr(S);
 ReadInt := False;
 If T = '' Then Exit;
 Val(T, R, Code);
 If Code = 0 Then Begin V := R; ReadInt := True; Exit; End;
 If ReadReal(T, X) Then Begin
  If (Abs(X) < 2E9) And (Abs(X - Round(X)) < 1E-9) Then Begin V := Round(X);
    ReadInt := True; End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Preprocesseur                                                         }
{ --------------------------------------------------------------------- }

Const
 MAX_PVAR = 80;

Var
 PVName  : Array[1..MAX_PVAR] Of String[40];
 PVVal   : Array[1..MAX_PVAR] Of String[120];
 PVLock  : Array[1..MAX_PVAR] Of Boolean;
 NPV     : Integer;
 LnP     : Array[1..MAX_LINES] Of PStr;
 LnNo    : Array[1..MAX_LINES] Of Integer;
 NLn     : Integer;
 NInErr  : Integer;
 CheckMode : Boolean;
 RunPhase : Boolean;

Procedure InErr(Line: Integer; Const Msg: String);
Begin
 If RunPhase Then Abort(Msg);
 Inc(NInErr);
 If Line > 0 Then ConMsg('ERREUR (ligne ' + IStr(Line) + ') : ' + Msg)
 Else ConMsg('ERREUR : ' + Msg);
 If OutOpen Then Begin
  If Line > 0 Then Emit(' *** ERROR in line ' + IStr(Line) + ' :: ' + Msg +
    ' ***')
  Else Emit(' *** ERROR :: ' + Msg + ' ***');
 End;
End;

Procedure SetPVar(Const N, V: String);
Var I: Integer; U: String;
Begin
 U := UpStr(N);
 For I := 1 To NPV Do
  If PVName[I] = U Then Begin
   If Not PVLock[I] Then PVVal[I] := V;
   Exit;
  End;
 If NPV < MAX_PVAR Then Begin
  Inc(NPV);
  PVName[NPV] := U;
  PVVal[NPV] := V;
  PVLock[NPV] := False;
 End;
End;

{ Variable fixee par la ligne de commande (-E) : les @SET ne la modifient pas }
Procedure SetPVarLocked(Const N, V: String);
Begin
 SetPVar(N, V);
 PVLock[NPV] := True;
End;

Function GetPVar(Const N: String; Var V: String): Boolean;
Var I: Integer; U: String;
Begin
 U := UpStr(N);
 For I := 1 To NPV Do
  If PVName[I] = U Then Begin V := PVVal[I]; GetPVar := True; Exit; End;
 GetPVar := False;
End;

Function IsIdChar(C: Char): Boolean;
Begin
 IsIdChar := (C In ['A'..'Z', 'a'..'z', '0'..'9', '_']);
End;

{ Substitution des variables du preprocesseur (dollar-accolade et dollar) }
Function Subst(Const S: String; Line: Integer): String;
Var I, J: Integer; R, N, V: String;
Begin
 R := '';
 I := 1;
 While I <= Length(S) Do Begin
  If (S[I] = '$') And (I < Length(S)) Then Begin
   If S[I + 1] = '{' Then Begin
    J := I + 2;
    N := '';
    While (J <= Length(S)) And (S[J] <> '}') Do Begin N := N + S[J]; Inc(J);
      End;
    I := J + 1;
   End
   Else Begin
    J := I + 1;
    N := '';
    While (J <= Length(S)) And IsIdChar(S[J]) Do Begin N := N + S[J]; Inc(J);
      End;
    I := J;
   End;
   If GetPVar(TrimStr(N), V) Then R := R + V
   Else Begin
    InErr(Line, 'Variable du preprocesseur non definie : ' + N);
   End;
  End
  Else Begin
   R := R + S[I];
   Inc(I);
  End;
 End;
 Subst := R;
End;

{ Evalue une condition @IF }
Function EvalCond(Const C0: String): Boolean;
Var C, A, B, Op: String; P: Integer; X, Y: Real; OkA, OkB: Boolean;
 Function Cmp(Const O: String): Boolean;
 Begin
  P := Pos(O, C);
  If P > 0 Then Begin
   A := TrimStr(Copy(C, 1, P - 1));
   B := TrimStr(Copy(C, P + Length(O), 255));
   Op := O;
   Cmp := True;
  End
  Else Cmp := False;
 End;
Begin
 C := TrimStr(C0);
 While (Length(C) > 1) And (C[1] = '(') And (C[Length(C)] = ')') Do
  C := TrimStr(Copy(C, 2, Length(C) - 2));
 If Not (Cmp('==') Or Cmp('/=') Or Cmp('!=') Or Cmp('<=') Or Cmp('>=') Or
         Cmp('<') Or Cmp('>')) Then Begin
  EvalCond := Not ((C = '') Or (C = '0') Or (UpStr(C) = 'FALSE') Or
                   (UpStr(C) = '.FALSE.') Or (UpStr(C) = 'NO'));
  Exit;
 End;
 If (Length(A) > 1) And ((A[1] = '"') Or (A[1] = '''')) Then A := Copy(A, 2,
   Length(A) - 2);
 If (Length(B) > 1) And ((B[1] = '"') Or (B[1] = '''')) Then B := Copy(B, 2,
   Length(B) - 2);
 OkA := ReadReal(A, X);
 OkB := ReadReal(B, Y);
 If OkA And OkB Then Begin
  If Op = '==' Then EvalCond := Abs(X - Y) < 1E-12
  Else If (Op = '/=') Or (Op = '!=') Then EvalCond := Abs(X - Y) >= 1E-12
  Else If Op = '<' Then EvalCond := X < Y
  Else If Op = '>' Then EvalCond := X > Y
  Else If Op = '<=' Then EvalCond := X <= Y + 1E-12
  Else EvalCond := X >= Y - 1E-12;
 End
 Else Begin
  If Op = '==' Then EvalCond := UpStr(A) = UpStr(B)
  Else If (Op = '/=') Or (Op = '!=') Then EvalCond := UpStr(A) <> UpStr(B)
  Else If Op = '<' Then EvalCond := UpStr(A) < UpStr(B)
  Else If Op = '>' Then EvalCond := UpStr(A) > UpStr(B)
  Else If Op = '<=' Then EvalCond := UpStr(A) <= UpStr(B)
  Else EvalCond := UpStr(A) >= UpStr(B);
 End;
End;

Procedure AddLine(Const S: String; N: Integer);
Var P: PStr;
Begin
 If NLn >= MAX_LINES Then Begin
  InErr(N, 'Fichier d''entree trop long');
  Exit;
 End;
 Inc(NLn);
 GetMem(P, Length(S) + 1);
 P^ := S;
 LnP[NLn] := P;
 LnNo[NLn] := N;
End;

Procedure ClearLines;
Var I: Integer;
Begin
 For I := 1 To NLn Do FreeMem(LnP[I], Length(LnP[I]^) + 1);
 NLn := 0;
End;

{ Lit un fichier source en appliquant le preprocesseur }
Procedure ReadSource(Const FN: String; Depth: Integer);
Var
 Fi: Text;
 Raw, Pend, T, Kw, Rest: String;
 LNo, PendNo, StartNo, IORes: Integer;
 Cond: Array[0..30] Of Boolean;     { etat de chaque niveau @IF }
 Taken: Array[0..30] Of Boolean;
 Top: Integer;
 Active: Boolean;
 Function CondActive: Boolean;
 Var K: Integer;
 Begin
  CondActive := True;
  For K := 1 To Top Do If Not Cond[K] Then CondActive := False;
 End;
 Procedure Process(Const L: String; No: Integer);
 Var V, Nm, Val2: String;
 Begin
  T := TrimStr(L);
  If T = '' Then Exit;
  If T[1] = '@' Then Begin
   Rest := Copy(T, 2, 255);
   Kw := UpStr(GetTok(Rest));
   Rest := TrimStr(Rest);
   If Kw = 'IF' Then Begin
    Active := CondActive;
    If Top < 30 Then Inc(Top);
    If Active Then Begin
     Cond[Top] := EvalCond(Subst(Rest, No));
     Taken[Top] := Cond[Top];
    End
    Else Begin
     Cond[Top] := False;
     Taken[Top] := True;
    End;
    Exit;
   End;
   If Kw = 'ELSE' Then Begin
    If Top > 0 Then Begin
     Cond[Top] := Not Taken[Top];
     Taken[Top] := True;
    End
    Else InErr(No, '@ELSE sans @IF');
    Exit;
   End;
   If Kw = 'ENDIF' Then Begin
    If Top > 0 Then Dec(Top) Else InErr(No, '@ENDIF sans @IF');
    Exit;
   End;
   If Not CondActive Then Exit;
   If Kw = 'SET' Then Begin
    Rest := Subst(Rest, No);
    Nm := GetTok(Rest);
    Val2 := TrimStr(Rest);
    If (Length(Val2) > 1) And ((Val2[1] = '"') Or (Val2[1] = '''')) Then Val2 :=
      GetTok(Val2);
    SetPVar(Nm, Val2);
    Exit;
   End;
   If Kw = 'INCLUDE' Then Begin
    Rest := Subst(Rest, No);
    V := GetTok(Rest);
    If Depth < 8 Then ReadSource(V, Depth + 1)
    Else InErr(No, '@INCLUDE trop imbrique');
    Exit;
   End;
   If Kw = 'PRINT' Then Begin
    ConMsg(Subst(Rest, No));
    Exit;
   End;
   InErr(No, 'Directive inconnue du preprocesseur : @' + Kw);
   Exit;
  End;
  If Not CondActive Then Exit;
  AddLine(Subst(T, No), No);
 End;
Begin
 Top := 0;
 Cond[0] := True;
 Taken[0] := True;
 If FN = '' Then Begin
  InErr(0, 'Nom de fichier vide');
  Exit;
 End;
 {$I-}
 If FN = '-' Then Assign(Fi, '') Else Assign(Fi, FN);
 Reset(Fi);
 {$I+}
 IORes := IOResult;
 If IORes <> 0 Then Begin
  InErr(0, 'Impossible d''ouvrir le fichier ' + FN);
  Exit;
 End;
 LNo := 0;
 Pend := '';
 PendNo := 0;
 While Not Eof(Fi) Do Begin
  ReadLn(Fi, Raw);
  Inc(LNo);
  If Length(Raw) >= 255 Then InErr(LNo,
    'Ligne trop longue (255 caracteres au plus)');
  Raw := StripComment(Raw);
  If Pend <> '' Then Begin
   If Length(Pend) + Length(Raw) > 245 Then InErr(PendNo,
     'Ligne logique trop longue (245 caracteres au plus)');
   Raw := Copy(Pend + ' ' + Raw, 1, 250);
   StartNo := PendNo;
  End
  Else StartNo := LNo;
  If (Raw <> '') And (Raw[Length(Raw)] = '\') Then Begin
   Pend := TrimStr(Copy(Raw, 1, Length(Raw) - 1));
   If Pend = '' Then Pend := ' ';
   PendNo := StartNo;
   Continue;
  End;
  Pend := '';
  Process(Raw, StartNo);
 End;
 If Pend <> '' Then Process(Pend, PendNo);
 Close(Fi);
 If Top <> 0 Then InErr(LNo, '@IF sans @ENDIF dans ' + FN);
End;

{ --------------------------------------------------------------------- }
{ Arbre des sections                                                    }
{ --------------------------------------------------------------------- }

Var
 SecName   : Array[0..MAX_SEC] Of String[48];
 SecPar    : Array[0..MAX_SEC] Of String[100];
 SecParent : Array[0..MAX_SEC] Of Integer;
 SecLine   : Array[0..MAX_SEC] Of Integer;
 NSec      : Integer;
 KwTab     : Array[1..MAX_KW] Of PKw;
 NKw       : Integer;

Procedure ClearTree;
Var I: Integer;
Begin
 For I := 1 To NKw Do Dispose(KwTab[I]);
 NKw := 0;
 NSec := 0;
 SecName[0] := 'ROOT';
 SecPar[0] := '';
 SecParent[0] := -1;
 SecLine[0] := 0;
End;

Function AddKw(Sec: Integer; Const O, V: String; Line: Integer): Boolean;
Var P: PKw;
Begin
 AddKw := False;
 If NKw >= MAX_KW Then Begin InErr(Line, 'Trop de mots-cles'); Exit; End;
 New(P);
 P^.O := O;
 P^.U := UpStr(O);
 P^.V := V;
 P^.Sec := Sec;
 P^.Line := Line;
 Inc(NKw);
 KwTab[NKw] := P;
 AddKw := True;
End;

Function AddSec(Parent: Integer; Const Nm, Par: String; Line: Integer): Integer;
Begin
 If NSec >= MAX_SEC Then Begin
  InErr(Line, 'Trop de sections');
  AddSec := Parent;
  Exit;
 End;
 Inc(NSec);
 SecName[NSec] := UpStr(Nm);
 SecPar[NSec] := Par;
 SecParent[NSec] := Parent;
 SecLine[NSec] := Line;
 AddSec := NSec;
End;

{ Analyse les lignes A a B (pretraitees) sous la section Root }
Procedure ParseLines(A, B, Root: Integer);
Var
 I, Cur: Integer;
 S, T, Nm, Rest, T2, Nm2: String;
Begin
 Cur := Root;
 For I := A To B Do Begin
  S := LnP[I]^;
  If S = '' Then Continue;
  If S[1] = '&' Then Begin
   T := Copy(S, 2, 255);
   Nm := GetTok(T);
   Rest := TrimStr(T);
   If UpStr(Nm) = 'END' Then Begin
    If Cur = Root Then InErr(LnNo[I], '&END sans section ouverte')
    Else Begin
     T2 := Rest;
     Nm2 := '';
     If T2 <> '' Then Nm2 := UpStr(GetTok(T2));
     If (Nm2 <> '') And (Nm2 <> SecName[Cur]) Then
      InErr(LnNo[I], '&END ' + Nm2 + ' ne correspond pas a &' + SecName[Cur] +
            ' (ligne ' + IStr(SecLine[Cur]) + ')');
     Cur := SecParent[Cur];
    End;
   End
   Else Cur := AddSec(Cur, Nm, Rest, LnNo[I]);
  End
  Else Begin
   T := S;
   Nm := GetTok(T);
   If (Nm <> '') Then AddKw(Cur, Nm, TrimStr(T), LnNo[I]);
  End;
 End;
 If Cur <> Root Then
  InErr(SecLine[Cur], 'Section &' + SecName[Cur] +
    ' non fermee (&END manquant)');
End;

Procedure ParseTree;
Begin
 ClearTree;
 ParseLines(1, NLn, 0);
End;

{ Section fille Nth (1..) de nom donne ; 0 si absente }
Function FindSec(Parent: Integer; Const Nm: String; Nth: Integer): Integer;
Var I, C: Integer; U: String;
Begin
 U := UpStr(Nm);
 C := 0;
 For I := 1 To NSec Do
  If (SecParent[I] = Parent) And (SecName[I] = U) Then Begin
   Inc(C);
   If C = Nth Then Begin FindSec := I; Exit; End;
  End;
 FindSec := 0;
End;

Function CountSec(Parent: Integer; Const Nm: String): Integer;
Var I, C: Integer; U: String;
Begin
 U := UpStr(Nm);
 C := 0;
 For I := 1 To NSec Do
  If (SecParent[I] = Parent) And (SecName[I] = U) Then Inc(C);
 CountSec := C;
End;

{ Section designee par un chemin NOM/NOM/... depuis la racine }
Function SecPath(Const Path: String): Integer;
Var P, Cur: Integer; Rest, Nm: String;
Begin
 Cur := 0;
 Rest := Path;
 While (Rest <> '') And (Cur >= 0) Do Begin
  P := Pos('/', Rest);
  If P > 0 Then Begin Nm := Copy(Rest, 1, P - 1); Rest := Copy(Rest, P + 1,
    255); End
  Else Begin Nm := Rest; Rest := ''; End;
  Cur := FindSec(Cur, Nm, 1);
  If Cur = 0 Then Begin SecPath := 0; Exit; End;
 End;
 SecPath := Cur;
End;

Function SecPathName(Sec: Integer): String;
Var S: String;
Begin
 S := '';
 While Sec > 0 Do Begin
  If S = '' Then S := SecName[Sec] Else S := SecName[Sec] + '/' + S;
  Sec := SecParent[Sec];
 End;
 SecPathName := S;
End;

{ Index du premier mot-cle de nom donne dans la section ; 0 si absent }
Function KwIdx(Sec: Integer; Const Nm: String): Integer;
Var I: Integer; U: String;
Begin
 U := UpStr(Nm);
 For I := 1 To NKw Do
  If (KwTab[I]^.Sec = Sec) And (KwTab[I]^.U = U) Then Begin KwIdx := I; Exit;
    End;
 KwIdx := 0;
End;

Function KwHas(Sec: Integer; Const Nm: String): Boolean;
Begin
 If Sec <= 0 Then KwHas := False Else KwHas := KwIdx(Sec, Nm) > 0;
End;

{ Valeur brute (sans l'unite) ; Un recoit l'unite entre crochets }
Function KwRaw(Sec: Integer; Const Nm: String; Var Un: String; Var Found:
  Boolean): String;
Var I, P: Integer; V: String;
Begin
 Un := '';
 Found := False;
 KwRaw := '';
 If Sec < 0 Then Exit;
 I := KwIdx(Sec, Nm);
 If I = 0 Then Exit;
 Found := True;
 V := TrimStr(KwTab[I]^.V);
 If (V <> '') And (V[1] = '[') Then Begin
  P := Pos(']', V);
  If P > 0 Then Begin
   Un := Copy(V, 2, P - 2);
   V := TrimStr(Copy(V, P + 1, 255));
  End;
 End;
 KwRaw := V;
End;

Function KwStr(Sec: Integer; Const Nm, Def: String): String;
Var Un, V: String; Fd: Boolean;
Begin
 V := KwRaw(Sec, Nm, Un, Fd);
 If Not Fd Then KwStr := Def
 Else KwStr := GetTok(V);
End;

{ Mot en majuscules }
Function KwWord(Sec: Integer; Const Nm, Def: String): String;
Begin
 KwWord := UpStr(KwStr(Sec, Nm, Def));
End;

{ Le N-ieme mot (1..) de la valeur }
Function KwTok(Sec: Integer; Const Nm: String; N: Integer): String;
Var Un, V, T: String; Fd: Boolean; I: Integer;
Begin
 V := KwRaw(Sec, Nm, Un, Fd);
 T := '';
 For I := 1 To N Do T := GetTok(V);
 KwTok := T;
End;

Function KwBool(Sec: Integer; Const Nm: String; Def: Boolean): Boolean;
Var Un, V, W: String; Fd: Boolean;
Begin
 W := KwRaw(Sec, Nm, Un, Fd);
 V := UpStr(GetTok(W));
 If Not Fd Then KwBool := Def
 Else If V = '' Then KwBool := True
 Else If (V = 'T') Or (V = '.TRUE.') Or (V = 'TRUE') Or (V = 'YES') Or
         (V = 'ON') Or (V = '1') Or (V = '.T.') Then KwBool := True
 Else Begin
  If Not ((V = 'F') Or (V = '.FALSE.') Or (V = 'FALSE') Or (V = 'NO') Or
          (V = 'OFF') Or (V = '0') Or (V = '.F.')) Then
   InErr(KwTab[KwIdx(Sec, Nm)]^.Line, 'Valeur logique invalide pour ' +
     UpStr(Nm) + ' : ' + V);
  KwBool := False;
 End;
End;

Function KwInt(Sec: Integer; Const Nm: String; Def: LongInt): LongInt;
Var Un, V, W: String; Fd: Boolean; R: LongInt;
Begin
 W := KwRaw(Sec, Nm, Un, Fd);
 V := GetTok(W);
 If Fd And ReadInt(V, R) Then KwInt := R
 Else Begin
  If Fd Then InErr(KwTab[KwIdx(Sec, Nm)]^.Line, 'Entier invalide pour ' +
    UpStr(Nm) + ' : ' + V);
  KwInt := Def;
 End;
End;

{ Reel avec unite : Def est exprime dans l'unite DU ; le resultat est en
  unites atomiques (ou brut si DU est vide ou K) }
Function KwReal(Sec: Integer; Const Nm: String; Def: Real; Const DU: String):
  Real;
Var Un, V: String; Fd, Ok: Boolean; R, Fac: Real;
Begin
 V := KwRaw(Sec, Nm, Un, Fd);
 Ok := True;
 If Not Fd Then R := Def
 Else If Not ReadReal(GetTok(V), R) Then Begin
  InErr(KwTab[KwIdx(Sec, Nm)]^.Line, 'Valeur numerique invalide pour ' +
    UpStr(Nm));
  R := Def;
 End
 Else If Un <> '' Then Begin
  Fac := UnitToAU(Un, Ok);
  KwReal := R * Fac;
  Exit;
 End;
 If (DU = '') Or (UpStr(DU) = 'K') Then KwReal := R
 Else Begin
  Fac := UnitToAU(DU, Ok);
  KwReal := R * Fac;
 End;
End;

{ Le N-ieme reel d'un mot-cle (avec unite) ; Ok = vrai si lu }
Function KwRealN(Sec: Integer; Const Nm: String; N: Integer; Const DU: String;
                 Var V: Real): Boolean;
Var Un, S, T: String; Fd, Ok: Boolean; I: Integer; Fac: Real;
Begin
 S := KwRaw(Sec, Nm, Un, Fd);
 KwRealN := False;
 If Not Fd Then Exit;
 T := '';
 For I := 1 To N Do T := GetTok(S);
 If Not ReadReal(T, V) Then Exit;
 Ok := True;
 If Un <> '' Then Fac := UnitToAU(Un, Ok)
 Else If (DU <> '') And (UpStr(DU) <> 'K') Then Fac := UnitToAU(DU, Ok)
 Else Fac := 1;
 V := V * Fac;
 KwRealN := True;
End;

{ --------------------------------------------------------------------- }
{ Ecriture de l'arbre (option --echo et fichiers de redemarrage)         }
{ --------------------------------------------------------------------- }

Var
 RestartMode : Boolean;
 RestartStep : LongInt;
 RestartTime : Real;
 RestartVelOK : Boolean;
 UseResCounters : Boolean;
 ResStep : LongInt;
 ResTime : Real;

Procedure DumpCoordSections(Ind: Integer); Forward;
Function RestartCellKnown: Boolean; Forward;
Procedure DumpCellVectors(Const Pad: String); Forward;

Procedure DumpSec(S, Ind: Integer);
Var I: Integer; Pad: String; Skip: Boolean;
Begin
 Pad := Rep(' ', Ind);
 If S > 0 Then Begin
  If SecPar[S] <> '' Then Emit(Pad + '&' + SecName[S] + ' ' + SecPar[S])
  Else Emit(Pad + '&' + SecName[S]);
 End;
 For I := 1 To NKw Do
  If KwTab[I]^.Sec = S Then Begin
   Skip := False;
   If RestartMode Then Begin
    If (SecName[S] = 'TOPOLOGY') And ((KwTab[I]^.U = 'COORD_FILE_NAME') Or
       (KwTab[I]^.U = 'COORD_FILE_FORMAT')) Then Skip := True;
    If (SecName[S] = 'MD') And ((KwTab[I]^.U = 'STEP_START_VAL') Or
       (KwTab[I]^.U = 'TIME_START_VAL')) Then Skip := True;
    If (SecName[S] = 'GEO_OPT') And (KwTab[I]^.U = 'STEP_START_VAL') Then Skip
      := True;
    If SecName[S] = 'COORD' Then Skip := True;
    If (SecName[S] = 'CELL') And ((KwTab[I]^.U = 'ABC') Or (KwTab[I]^.U =
      'ALPHA_BETA_GAMMA') Or
       (KwTab[I]^.U = 'A') Or (KwTab[I]^.U = 'B') Or (KwTab[I]^.U = 'C')) And
         RestartCellKnown Then Skip := True;
   End;
   If Not Skip Then Begin
    If KwTab[I]^.V <> '' Then Emit(Pad + '   ' + KwTab[I]^.O + ' ' +
      KwTab[I]^.V)
    Else Emit(Pad + '   ' + KwTab[I]^.O);
   End;
  End;
 If RestartMode And RestartCellKnown And (SecName[S] = 'CELL') Then
   DumpCellVectors(Pad + '   ');
 If RestartMode And (SecName[S] = 'MD') Then Begin
  Emit(Pad + '   STEP_START_VAL ' + IStr(RestartStep));
  Emit(Pad + '   TIME_START_VAL ' + RStr(RestartTime / FS_AU, 8));
 End;
 If RestartMode And (SecName[S] = 'GEO_OPT') Then
  Emit(Pad + '   STEP_START_VAL ' + IStr(RestartStep));
 For I := 1 To NSec Do
  If SecParent[I] = S Then Begin
   If RestartMode And ((SecName[I] = 'COORD') Or (SecName[I] = 'VELOCITY')) And
      (SecName[S] = 'SUBSYS') Then Continue;
   If RestartMode And ((SecName[I] = 'EXT_RESTART') Or (SecName[I] =
     'RESTART_DATA')) Then Continue;
   If S = 0 Then DumpSec(I, 0) Else DumpSec(I, Ind + 3);
  End;
 If RestartMode And (SecName[S] = 'SUBSYS') Then DumpCoordSections(Ind + 3);
 If S > 0 Then Emit(Pad + '&END ' + SecName[S]);
End;


{ --------------------------------------------------------------------- }
{ Tableau periodique (H a Rn)                                           }
{ --------------------------------------------------------------------- }

Const
 NPT = 86;
 PtSym : Array[1..NPT] Of String[2] = (
  'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
  'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
  'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
  'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
  'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',
  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn');
 { masses atomiques standard (g/mol) }
 PtMass : Array[1..NPT] Of Real = (
  1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994,
  18.9984032, 20.1797, 22.98976928, 24.305, 26.9815386, 28.0855, 30.973762,
  32.065, 35.453, 39.948, 39.0983, 40.078, 44.955912, 47.867, 50.9415,
  51.9961, 54.938045, 55.845, 58.933195, 58.6934, 63.546, 65.38, 69.723,
  72.64, 74.9216, 78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585, 91.224,
  92.90638, 95.96, 98.0, 101.07, 102.9055, 106.42, 107.8682, 112.411,
  114.818, 118.71, 121.76, 127.6, 126.90447, 131.293, 132.9054519,
  137.327, 138.90547, 140.116, 140.90765, 144.242, 145.0, 150.36, 151.964,
  157.25, 158.92535, 162.5, 164.93032, 167.259, 168.93421, 173.054,
  174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084,
  196.966569, 200.59, 204.3833, 207.2, 208.9804, 209.0, 210.0, 222.0);
 { rayons covalents (angstrom, Cordero) }
 PtRad : Array[1..NPT] Of Real = (
  0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41,
  1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53, 1.39,
  1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16,
  2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
  1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01,
  1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75,
  1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.40,
  1.50, 1.50);

{ Numero atomique d'un symbole exact (insensible a la casse) ; 0 sinon }
Function ElemOfSym(Const S: String): Integer;
Var I: Integer; U: String;
Begin
 U := UpStr(S);
 For I := 1 To NPT Do
  If UpStr(PtSym[I]) = U Then Begin ElemOfSym := I; Exit; End;
 ElemOfSym := 0;
End;

{ Element deduit d'une etiquette (O, O1, Ow, CL, H_1, ...) }
Function ElemFromLabel(Const L: String): Integer;
Var A: String; I, Z1, Z2: Integer;
Begin
 A := '';
 I := 1;
 While (I <= Length(L)) And (UpStr(L[I]) >= 'A') And (UpStr(L[I]) <= 'Z') Do
   Begin
  A := A + L[I];
  Inc(I);
 End;
 ElemFromLabel := 0;
 If A = '' Then Exit;
 Z1 := ElemOfSym(Copy(A, 1, 1));
 Z2 := 0;
 If Length(A) >= 2 Then Z2 := ElemOfSym(Copy(A, 1, 2));
 If Length(A) = 1 Then ElemFromLabel := Z1
 Else If (Length(A) = 2) And (Z2 > 0) Then ElemFromLabel := Z2
 Else If Z1 > 0 Then ElemFromLabel := Z1
 Else ElemFromLabel := Z2;
End;

{ --------------------------------------------------------------------- }
{ Etat du systeme                                                       }
{ --------------------------------------------------------------------- }

Const
 MAX_NB = 8;

Var
 SecGlobal, SecFE, SecSub, SecMotion, SecMM, SecDFT : Integer;
 FEMethod   : Integer;
 RunType    : Integer;
 RunName    : String[24];

 NAtoms     : Integer;
 SPos, SVel, SFrc : TXYZ;
 SLab       : Array[1..MAX_ATOM] Of String[12];
 SZ, SKind, SMol : Array[1..MAX_ATOM] Of Integer;
 SMass, SQ  : Array[1..MAX_ATOM] Of Real;
 NMol       : Integer;

 NKinds     : Integer;
 KName      : Array[1..MAX_KIND] Of String[12];
 KZ, KBas, KCount : Array[1..MAX_KIND] Of Integer;
 KMassA, KQ : Array[1..MAX_KIND] Of Real;
 KBasName   : Array[1..MAX_KIND] Of String[40];
 KPotName   : Array[1..MAX_KIND] Of String[40];

 CellH, CellHI : TMat3;
 CellVol    : Real;
 HasCell, CellOrtho : Boolean;
 Per        : Array[1..3] Of Boolean;
 PerStr     : String[8];

 { connectivite }
 NBond      : Integer;
 BondA, BondB : Array[1..MAX_BONDS] Of Integer;
 AdjN       : Array[1..MAX_ATOM] Of Integer;
 Adj        : Array[1..MAX_ATOM, 1..MAX_NB] Of Integer;
 ExN        : Array[1..MAX_ATOM] Of Integer;
 ExJ        : Array[1..MAX_ATOM, 1..MAX_EXCL] Of Integer;
 ExD        : Array[1..MAX_ATOM, 1..MAX_EXCL] Of Integer;
 ExclMax    : Integer;      { profondeur maximale d'exclusion : 0 a 3 }
 ExclMaxVdw : Integer;
 Scale14Ei, Scale14Vdw : Real;

{ --------------------------------------------------------------------- }
{ Algebre 3x3 et cellule                                                }
{ --------------------------------------------------------------------- }

Function Inv3(Const M: TMat3; Var R: TMat3): Real;
Var D: Real;
Begin
 D := M[1,1] * (M[2,2] * M[3,3] - M[2,3] * M[3,2])
    - M[1,2] * (M[2,1] * M[3,3] - M[2,3] * M[3,1])
    + M[1,3] * (M[2,1] * M[3,2] - M[2,2] * M[3,1]);
 If Abs(D) < 1E-30 Then Begin
  FillChar(R, SizeOf(R), 0);
  Inv3 := 0;
  Exit;
 End;
 R[1,1] := (M[2,2] * M[3,3] - M[2,3] * M[3,2]) / D;
 R[1,2] := (M[1,3] * M[3,2] - M[1,2] * M[3,3]) / D;
 R[1,3] := (M[1,2] * M[2,3] - M[1,3] * M[2,2]) / D;
 R[2,1] := (M[2,3] * M[3,1] - M[2,1] * M[3,3]) / D;
 R[2,2] := (M[1,1] * M[3,3] - M[1,3] * M[3,1]) / D;
 R[2,3] := (M[1,3] * M[2,1] - M[1,1] * M[2,3]) / D;
 R[3,1] := (M[2,1] * M[3,2] - M[2,2] * M[3,1]) / D;
 R[3,2] := (M[1,2] * M[3,1] - M[1,1] * M[3,2]) / D;
 R[3,3] := (M[1,1] * M[2,2] - M[1,2] * M[2,1]) / D;
 Inv3 := D;
End;

{ Met a jour l'inverse, le volume et l'orthogonalite de la cellule }
Procedure UpdateCell;
Var D: Real;
Begin
 D := Inv3(CellH, CellHI);
 CellVol := Abs(D);
 CellOrtho := (Abs(CellH[1,2]) + Abs(CellH[1,3]) + Abs(CellH[2,1]) +
               Abs(CellH[2,3]) + Abs(CellH[3,1]) + Abs(CellH[3,2])) < 1E-10;
End;

{ Image minimale d'un vecteur de deplacement }
Procedure MinImg(Var X, Y, Z: Real);
Var S1, S2, S3: Real;
Begin
 If Not HasCell Then Exit;
 If CellOrtho Then Begin
  If Per[1] Then X := X - CellH[1,1] * Round(X / CellH[1,1]);
  If Per[2] Then Y := Y - CellH[2,2] * Round(Y / CellH[2,2]);
  If Per[3] Then Z := Z - CellH[3,3] * Round(Z / CellH[3,3]);
  Exit;
 End;
 S1 := X * CellHI[1,1] + Y * CellHI[2,1] + Z * CellHI[3,1];
 S2 := X * CellHI[1,2] + Y * CellHI[2,2] + Z * CellHI[3,2];
 S3 := X * CellHI[1,3] + Y * CellHI[2,3] + Z * CellHI[3,3];
 If Per[1] Then S1 := S1 - Round(S1);
 If Per[2] Then S2 := S2 - Round(S2);
 If Per[3] Then S3 := S3 - Round(S3);
 X := S1 * CellH[1,1] + S2 * CellH[2,1] + S3 * CellH[3,1];
 Y := S1 * CellH[1,2] + S2 * CellH[2,2] + S3 * CellH[3,2];
 Z := S1 * CellH[1,3] + S2 * CellH[2,3] + S3 * CellH[3,3];
End;

Function VecLen(Const V: TVec3): Real;
Begin
 VecLen := Sqrt(V[1] * V[1] + V[2] * V[2] + V[3] * V[3]);
End;

{ Largeur de la cellule perpendiculaire a la direction I }
Function CellWidth(I: Integer): Real;
Var J, K: Integer; C: TVec3;
Begin
 J := I Mod 3 + 1;
 K := J Mod 3 + 1;
 C[1] := CellH[J,2] * CellH[K,3] - CellH[J,3] * CellH[K,2];
 C[2] := CellH[J,3] * CellH[K,1] - CellH[J,1] * CellH[K,3];
 C[3] := CellH[J,1] * CellH[K,2] - CellH[J,2] * CellH[K,1];
 If VecLen(C) < 1E-12 Then CellWidth := 0
 Else CellWidth := CellVol / VecLen(C);
End;

Function CellLen(I: Integer): Real;
Begin
 CellLen := Sqrt(Sqr(CellH[I,1]) + Sqr(CellH[I,2]) + Sqr(CellH[I,3]));
End;

Function CellAngle(I, J: Integer): Real;
Var D: Real;
Begin
 D := (CellH[I,1] * CellH[J,1] + CellH[I,2] * CellH[J,2] + CellH[I,3] *
   CellH[J,3]) /
      (CellLen(I) * CellLen(J));
 CellAngle := ArcCos(D) * 180 / Pi;
End;

{ Construit la cellule a partir de la section &CELL }
Procedure BuildCell(Sub: Integer);
Var
 Sc, I, J, NAbg: Integer;
 A, B, C, Al, Be, Ga, V: Real;
 CosA, CosB, CosG, SinG, Cx, Cy, Cz: Real;
 P: String;
 UseVec: Boolean;
Begin
 FillChar(CellH, SizeOf(CellH), 0);
 HasCell := False;
 Per[1] := True; Per[2] := True; Per[3] := True;
 PerStr := 'XYZ';
 Sc := FindSec(Sub, 'CELL', 1);
 If Sc = 0 Then Begin
  Per[1] := False; Per[2] := False; Per[3] := False;
  PerStr := 'NONE';
  Exit;
 End;
 P := KwWord(Sc, 'PERIODIC', 'XYZ');
 If Not ((P = 'XYZ') Or (P = 'XY') Or (P = 'XZ') Or (P = 'YZ') Or (P = 'X') Or
   (P = 'Y') Or
         (P = 'Z') Or (P = 'NONE')) Then Begin
  InErr(SecLine[Sc], 'PERIODIC invalide : ' + P);
  P := 'NONE';
 End;
 Per[1] := Pos('X', P) > 0; Per[2] := Pos('Y', P) > 0; Per[3] := Pos('Z', P) >
   0;
 If P = 'NONE' Then Begin Per[1] := False; Per[2] := False; Per[3] := False;
   End;
 PerStr := P;
 UseVec := KwHas(Sc, 'A') And KwHas(Sc, 'B') And KwHas(Sc, 'C');
 If UseVec Then Begin
  For I := 1 To 3 Do
   For J := 1 To 3 Do Begin
    V := 0;
    If Not KwRealN(Sc, Chr(Ord('A') + I - 1), J, 'angstrom', V) Then
     InErr(SecLine[Sc], 'Vecteur de cellule ' + Chr(Ord('A') + I - 1) +
       ' incomplet');
    CellH[I,J] := V;
   End;
 End
 Else If KwHas(Sc, 'ABC') Then Begin
  A := 0; B := 0; C := 0;
  If Not (KwRealN(Sc, 'ABC', 1, 'angstrom', A) And KwRealN(Sc, 'ABC', 2,
    'angstrom', B) And
          KwRealN(Sc, 'ABC', 3, 'angstrom', C)) Then Begin
   InErr(SecLine[Sc], 'ABC doit contenir trois valeurs');
   Exit;
  End;
  Al := 90; Be := 90; Ga := 90;
  If KwHas(Sc, 'ALPHA_BETA_GAMMA') Then Begin
   NAbg := 0;
   If KwRealN(Sc, 'ALPHA_BETA_GAMMA', 1, 'deg', Al) Then Inc(NAbg);
   If KwRealN(Sc, 'ALPHA_BETA_GAMMA', 2, 'deg', Be) Then Inc(NAbg);
   If KwRealN(Sc, 'ALPHA_BETA_GAMMA', 3, 'deg', Ga) Then Inc(NAbg);
   If NAbg < 3 Then Begin
    InErr(SecLine[Sc], 'ALPHA_BETA_GAMMA doit contenir trois angles');
    Exit;
   End;
   If (Al < 0.01) Or (Al > Pi - 0.01) Or (Be < 0.01) Or (Be > Pi - 0.01) Or (Ga
     < 0.01) Or (Ga > Pi - 0.01) Then Begin
    InErr(SecLine[Sc],
      'Les angles de la cellule doivent etre compris entre 0 et 180 degres');
    Exit;
   End;
  End;
  KwRealN(Sc, 'ALPHA_BETA_GAMMA', 1, 'deg', Al);
  KwRealN(Sc, 'ALPHA_BETA_GAMMA', 2, 'deg', Be);
  KwRealN(Sc, 'ALPHA_BETA_GAMMA', 3, 'deg', Ga);
  If KwHas(Sc, 'ALPHA_BETA_GAMMA') Then Begin
   { les angles lus sont deja en radians (unite deg convertie) }
  End
  Else Begin
   Al := Pi / 2; Be := Pi / 2; Ga := Pi / 2;
  End;
  CosA := Cos(Al); CosB := Cos(Be); CosG := Cos(Ga); SinG := Sin(Ga);
  If SinG < 1E-6 Then Begin InErr(SecLine[Sc], 'Cellule degeneree (gamma)');
    Exit; End;
  CellH[1,1] := A;
  CellH[2,1] := B * CosG;
  CellH[2,2] := B * SinG;
  Cx := C * CosB;
  Cy := C * (CosA - CosB * CosG) / SinG;
  Cz := C * C - Cx * Cx - Cy * Cy;
  If Cz < 1E-12 * C * C Then Begin InErr(SecLine[Sc],
    'Cellule degeneree (angles incompatibles)'); Exit; End;
  CellH[3,1] := Cx;
  CellH[3,2] := Cy;
  CellH[3,3] := Sqrt(Cz);
  If Abs(CellH[2,1]) < 1E-12 * A Then CellH[2,1] := 0;
  If Abs(CellH[3,1]) < 1E-12 * A Then CellH[3,1] := 0;
  If Abs(CellH[3,2]) < 1E-12 * A Then CellH[3,2] := 0;
 End
 Else Begin
  InErr(SecLine[Sc], 'La section &CELL exige ABC ou A, B et C');
  Exit;
 End;
 UpdateCell;
 If CellVol < 1E-10 Then Begin
  InErr(SecLine[Sc], 'Cellule de volume nul');
  Exit;
 End;
 HasCell := True;
End;

{ --------------------------------------------------------------------- }
{ Lecture des coordonnees                                               }
{ --------------------------------------------------------------------- }

Function FindKind(Const Lbl: String): Integer;
Var I: Integer; U: String;
Begin
 U := UpStr(Lbl);
 For I := 1 To NKinds Do
  If UpStr(KName[I]) = U Then Begin FindKind := I; Exit; End;
 FindKind := 0;
End;

Function AddKind(Const Lbl: String; Z: Integer): Integer;
Begin
 If NKinds >= MAX_KIND Then Begin
  InErr(0, 'Trop de sortes d''atomes (maximum ' + IStr(MAX_KIND) + ')');
  AddKind := NKinds;
  Exit;
 End;
 Inc(NKinds);
 KName[NKinds] := Lbl;
 KZ[NKinds] := Z;
 If Z > 0 Then KMassA[NKinds] := PtMass[Z] * AMU_AU Else KMassA[NKinds] :=
   AMU_AU;
 KQ[NKinds] := 0;
 KBas[NKinds] := -1;
 KBasName[NKinds] := '';
 KPotName[NKinds] := '';
 KCount[NKinds] := 0;
 AddKind := NKinds;
End;

Function BasisFromName(Const N0: String): Integer;
Var N: String;
Begin
 N := UpStr(N0);
 BasisFromName := -1;
 If (Pos('STO-3G', N) > 0) Or (Pos('STO3G', N) > 0) Or (Pos('SZV', N) > 0) Then
  BasisFromName := B_STO3G
 Else If (Pos('3-21G', N) > 0) Or (Pos('321G', N) > 0) Then BasisFromName :=
   B_321G
 Else If (Pos('6-31G', N) > 0) Or (Pos('631G', N) > 0) Or (Pos('DZV', N) > 0) Or
         (Pos('TZV', N) > 0) Or (Pos('QZV', N) > 0) Or (Pos('MOLOPT', N) > 0) Or
         (Pos('DZ', N) > 0) Then BasisFromName := B_631G;
End;

Procedure AddAtom(Const Lbl: String; X, Y, Z: Real);
Var K: Integer;
Begin
 If NAtoms >= MAX_ATOM Then Begin
  InErr(0, 'Trop d''atomes (maximum ' + IStr(MAX_ATOM) + ')');
  Exit;
 End;
 Inc(NAtoms);
 SLab[NAtoms] := Lbl;
 SPos[NAtoms,1] := X; SPos[NAtoms,2] := Y; SPos[NAtoms,3] := Z;
 K := FindKind(Lbl);
 If K = 0 Then Begin
  K := AddKind(Lbl, ElemFromLabel(Lbl));
 End;
 SKind[NAtoms] := K;
End;

{ Lecture d'un fichier XYZ (premier modele) ; FacL convertit en bohr }
Function ReadXYZCoords(Const FN: String; FacL: Real): Boolean;
Var
 Fi: Text; Ln, T, Lbl: String; N, I, IORes: Integer; X, Y, Z: Real; NN: LongInt;
Begin
 ReadXYZCoords := False;
 If FN = '' Then Begin InErr(0, 'COORD_FILE_NAME vide'); Exit; End;
 {$I-}
 Assign(Fi, FN);
 Reset(Fi);
 {$I+}
 IORes := IOResult;
 If IORes <> 0 Then Begin
  InErr(0, 'Fichier de coordonnees introuvable : ' + FN);
  Exit;
 End;
 ReadLn(Fi, Ln);
 If Not ReadInt(Ln, NN) Then Begin
  InErr(0, 'Fichier XYZ invalide (nombre d''atomes) : ' + FN);
  Close(Fi);
  Exit;
 End;
 If (NN < 1) Or (NN > MAX_ATOM) Then Begin
  InErr(0, 'Nombre d''atomes invalide dans ' + FN + ' (1 a ' + IStr(MAX_ATOM) +
    ')');
  Close(Fi);
  Exit;
 End;
 N := NN;
 ReadLn(Fi, Ln);
 For I := 1 To N Do Begin
  If Eof(Fi) Then Begin InErr(0, 'Fichier XYZ tronque : ' + FN); Close(Fi);
    Exit; End;
  ReadLn(Fi, Ln);
  T := TrimStr(Ln);
  Lbl := GetTok(T);
  If Not (ReadReal(GetTok(T), X) And ReadReal(GetTok(T), Y) And
    ReadReal(GetTok(T), Z)) Then Begin
   InErr(0, 'Ligne XYZ invalide : ' + Ln);
   Close(Fi);
   Exit;
  End;
  AddAtom(Lbl, X * FacL, Y * FacL, Z * FacL);
 End;
 Close(Fi);
 ReadXYZCoords := True;
End;

{ Lecture d'un fichier PDB (ATOM / HETATM) }
Function ReadPDBCoords(Const FN: String): Boolean;
Var
 Fi: Text; Ln, Lbl, El: String; IORes: Integer; X, Y, Z: Real;
Begin
 ReadPDBCoords := False;
 If FN = '' Then Begin InErr(0, 'COORD_FILE_NAME vide'); Exit; End;
 {$I-}
 Assign(Fi, FN);
 Reset(Fi);
 {$I+}
 IORes := IOResult;
 If IORes <> 0 Then Begin
  InErr(0, 'Fichier de coordonnees introuvable : ' + FN);
  Exit;
 End;
 While Not Eof(Fi) Do Begin
  ReadLn(Fi, Ln);
  If (Copy(Ln, 1, 4) = 'ATOM') Or (Copy(Ln, 1, 6) = 'HETATM') Then Begin
   If Not (ReadReal(Copy(Ln, 31, 8), X) And ReadReal(Copy(Ln, 39, 8), Y) And
           ReadReal(Copy(Ln, 47, 8), Z)) Then Begin
    InErr(0, 'Ligne PDB invalide : ' + Ln);
    Close(Fi);
    Exit;
   End;
   El := TrimStr(Copy(Ln, 77, 2));
   If (El = '') Or (ElemOfSym(El) = 0) Then Begin
    Lbl := TrimStr(Copy(Ln, 13, 4));
    While (Lbl <> '') And (Lbl[1] >= '0') And (Lbl[1] <= '9') Do Delete(Lbl, 1,
      1);
    El := Lbl;
    If ElemFromLabel(El) > 0 Then El := PtSym[ElemFromLabel(El)];
   End;
   AddAtom(El, X * ANG2BOHR, Y * ANG2BOHR, Z * ANG2BOHR);
  End;
 End;
 Close(Fi);
 ReadPDBCoords := True;
End;

{ Charge les atomes depuis &COORD ou depuis un fichier de coordonnees }
Procedure LoadCoords(SrcSec: Integer);
Var
 Sc, St, I: Integer;
 U, FN, Fmt, T, Lbl: String;
 FacL, X, Y, Z, S1, S2, S3: Real;
 Ok, Scaled: Boolean;
Begin
 NAtoms := 0;
 Sc := FindSec(SrcSec, 'COORD', 1);
 St := FindSec(SrcSec, 'TOPOLOGY', 1);
 If Sc > 0 Then Begin
  U := KwStr(Sc, 'UNIT', 'angstrom');
  FacL := UnitToAU(U, Ok);
  If Not Ok Then InErr(SecLine[Sc], 'Unite inconnue : ' + U);
  Scaled := KwBool(Sc, 'SCALED', False);
  If Scaled And Not HasCell Then InErr(SecLine[Sc], 'SCALED exige une cellule');
  For I := 1 To NKw Do
   If KwTab[I]^.Sec = Sc Then Begin
    If (KwTab[I]^.U = 'UNIT') Or (KwTab[I]^.U = 'SCALED') Then Continue;
    T := KwTab[I]^.V;
    Lbl := KwTab[I]^.O;
    If Not (ReadReal(GetTok(T), X) And ReadReal(GetTok(T), Y) And
      ReadReal(GetTok(T), Z)) Then Begin
     InErr(KwTab[I]^.Line, 'Coordonnees invalides pour l''atome ' + Lbl);
     Continue;
    End;
    If Scaled Then Begin
     S1 := X; S2 := Y; S3 := Z;
     X := S1 * CellH[1,1] + S2 * CellH[2,1] + S3 * CellH[3,1];
     Y := S1 * CellH[1,2] + S2 * CellH[2,2] + S3 * CellH[3,2];
     Z := S1 * CellH[1,3] + S2 * CellH[2,3] + S3 * CellH[3,3];
     AddAtom(Lbl, X, Y, Z);
    End
    Else AddAtom(Lbl, X * FacL, Y * FacL, Z * FacL);
   End;
 End
 Else If (St > 0) And KwHas(St, 'COORD_FILE_NAME') Then Begin
  FN := KwStr(St, 'COORD_FILE_NAME', '');
  Fmt := KwWord(St, 'COORD_FILE_FORMAT', 'XYZ');
  If Fmt = 'XYZ' Then ReadXYZCoords(FN, ANG2BOHR)
  Else If Fmt = 'PDB' Then ReadPDBCoords(FN)
  Else InErr(SecLine[St], 'Format de coordonnees non pris en charge : ' + Fmt);
 End
 Else InErr(SecLine[SrcSec],
   'Aucune coordonnee : &COORD ou COORD_FILE_NAME requis');
End;

{ Charge les vitesses depuis &VELOCITY (une ligne par atome) ; renvoie
  vrai si presentes }
Function LoadVelocities(SrcSec: Integer): Boolean;
Var Sv, I, N: Integer; T: String; X, Y, Z: Real;
Begin
 LoadVelocities := False;
 Sv := FindSec(SrcSec, 'VELOCITY', 1);
 If Sv = 0 Then Exit;
 N := 0;
 For I := 1 To NKw Do
  If KwTab[I]^.Sec = Sv Then Begin
   T := KwTab[I]^.O + ' ' + KwTab[I]^.V;
   If Not (ReadReal(GetTok(T), X) And ReadReal(GetTok(T), Y) And
     ReadReal(GetTok(T), Z)) Then Begin
    InErr(KwTab[I]^.Line, 'Vitesse invalide');
    Continue;
   End;
   Inc(N);
   If N <= NAtoms Then Begin
    SVel[N,1] := X; SVel[N,2] := Y; SVel[N,3] := Z;
   End;
  End;
 If N <> NAtoms Then Begin
  InErr(SecLine[Sv], 'Nombre de vitesses (' + IStr(N) +
    ') different du nombre d''atomes (' +
        IStr(NAtoms) + ')');
  Exit;
 End;
 LoadVelocities := True;
End;

{ Sorties : sections &COORD et &VELOCITY des fichiers de redemarrage }
Function RestartCellKnown: Boolean;
Begin
 RestartCellKnown := HasCell;
End;

Procedure DumpCellVectors(Const Pad: String);
Var K: Integer;
Begin
 For K := 1 To 3 Do
  Emit(Pad + Chr(Ord('A') + K - 1) + ' [angstrom] ' + FR(CellH[K,1] * BOHR_ANG,
    18, 10) +
       FR(CellH[K,2] * BOHR_ANG, 18, 10) + FR(CellH[K,3] * BOHR_ANG, 18, 10));
End;

Procedure DumpCoordSections(Ind: Integer);
Var I: Integer; P1, P2: String;
Begin
 P1 := Rep(' ', Ind);
 P2 := Rep(' ', Ind + 3);
 Emit(P1 + '&COORD');
 Emit(P2 + 'UNIT angstrom');
 For I := 1 To NAtoms Do
  Emit(P2 + PadR(SLab[I], 6) + ' ' + FR(SPos[I,1] * BOHR_ANG, 20, 12) + ' ' +
       FR(SPos[I,2] * BOHR_ANG, 20, 12) + ' ' + FR(SPos[I,3] * BOHR_ANG, 20,
         12));
 Emit(P1 + '&END COORD');
 Emit(P1 + '&VELOCITY');
 For I := 1 To NAtoms Do
  Emit(P2 + PadL(EStr(SVel[I,1], 13), 21) + PadL(EStr(SVel[I,2], 13), 21) +
    PadL(EStr(SVel[I,3], 13), 21));
 Emit(P1 + '&END VELOCITY');
End;

{ --------------------------------------------------------------------- }
{ Sortes d'atomes                                                       }
{ --------------------------------------------------------------------- }

Procedure SetupKinds;
Var I, K, Sk, N: Integer; Nm, El, B, T: String; Z: Integer;
Begin
 N := CountSec(SecSub, 'KIND');
 For I := 1 To N Do Begin
  Sk := FindSec(SecSub, 'KIND', I);
  T := SecPar[Sk];
  Nm := GetTok(T);
  If Nm = '' Then Begin InErr(SecLine[Sk], '&KIND exige un nom'); Continue; End;
  K := FindKind(Nm);
  El := KwStr(Sk, 'ELEMENT', '');
  If El <> '' Then Begin
   Z := ElemOfSym(El);
   If Z = 0 Then InErr(SecLine[Sk], 'Element inconnu : ' + El);
  End
  Else Z := ElemFromLabel(Nm);
  If K = 0 Then K := AddKind(Nm, Z)
  Else KZ[K] := Z;
  If (Z > 0) Then KMassA[K] := PtMass[Z] * AMU_AU;
  If KwHas(Sk, 'MASS') Then KMassA[K] := KwReal(Sk, 'MASS', 1.0, 'amu');
  B := KwStr(Sk, 'BASIS_SET', '');
  KBasName[K] := B;
  KPotName[K] := KwStr(Sk, 'POTENTIAL', '');
  If B <> '' Then Begin
   KBas[K] := BasisFromName(B);
   If KBas[K] < 0 Then Begin
    Warn('Base inconnue ' + B + ' pour la sorte ' + Nm + ' : 6-31G utilisee');
    KBas[K] := B_631G;
   End;
  End;
 End;
End;

Procedure ValidateKinds;
Var I: Integer;
Begin
 For I := 1 To NKinds Do
  If KZ[I] = 0 Then
   InErr(0, 'Element indetermine pour la sorte ' + KName[I] +
     ' (utilisez ELEMENT)');
End;

{ Remplit les tableaux par atome a partir des sortes }
Procedure FillAtomData;
Var I, K: Integer;
Begin
 For K := 1 To NKinds Do KCount[K] := 0;
 For I := 1 To NAtoms Do Begin
  K := SKind[I];
  SZ[I] := KZ[K];
  SMass[I] := KMassA[K];
  SQ[I] := KQ[K];
  Inc(KCount[K]);
  SVel[I,1] := 0; SVel[I,2] := 0; SVel[I,3] := 0;
  SFrc[I,1] := 0; SFrc[I,2] := 0; SFrc[I,3] := 0;
 End;
End;

{ --------------------------------------------------------------------- }
{ Atomes fixes (&MOTION / &CONSTRAINT / &FIXED_ATOMS)                   }
{ --------------------------------------------------------------------- }

Var
 FixM      : Array[1..MAX_ATOM, 1..3] Of Boolean;
 NFixComp  : Integer;
 FixActive : Boolean;

{ Lit les listes d'atomes (entiers et intervalles a..b) }
Procedure SetupFixed;
Var
 I, J, K, N, Sc, Sk, P, A: Integer;
 Comp, T, W: String;
 Cx, Cy, Cz: Boolean;
 LA, LB, AA: LongInt;
Begin
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do FixM[I, K] := False;
 NFixComp := 0;
 FixActive := True;
 If SecMotion = 0 Then Exit;
 Sc := FindSec(SecMotion, 'CONSTRAINT', 1);
 If Sc = 0 Then Exit;
 N := CountSec(Sc, 'FIXED_ATOMS');
 For J := 1 To N Do Begin
  Sk := FindSec(Sc, 'FIXED_ATOMS', J);
  Comp := KwWord(Sk, 'COMPONENTS_TO_FIX', 'XYZ');
  Cx := Pos('X', Comp) > 0; Cy := Pos('Y', Comp) > 0; Cz := Pos('Z', Comp) > 0;
  For I := 1 To NKw Do
   If (KwTab[I]^.Sec = Sk) And (KwTab[I]^.U = 'LIST') Then Begin
    T := KwTab[I]^.V;
    While T <> '' Do Begin
     W := GetTok(T);
     If W = '' Then Break;
     P := Pos('..', W);
     If P > 0 Then Begin
      If Not (ReadInt(Copy(W, 1, P - 1), LA) And ReadInt(Copy(W, P + 2, 20),
        LB)) Then Begin
       InErr(KwTab[I]^.Line, 'Intervalle invalide dans LIST : ' + W);
       Continue;
      End;
     End
     Else If Not ReadInt(W, LA) Then Begin
      InErr(KwTab[I]^.Line, 'Numero d''atome invalide dans LIST : ' + W);
      Continue;
     End
     Else LB := LA;
     If (LA < 1) Or (LB > NAtoms) Or (LA > LB) Then Begin
      InErr(KwTab[I]^.Line, 'Intervalle ' + IStr(LA) + '..' + IStr(LB) +
        ' hors de la plage 1..' + IStr(NAtoms));
      Continue;
     End;
     For AA := LA To LB Do Begin
      If Cx Then FixM[AA, 1] := True;
      If Cy Then FixM[AA, 2] := True;
      If Cz Then FixM[AA, 3] := True;
     End;
    End;
   End;
 End;
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do If FixM[I, K] Then Inc(NFixComp);
End;

Procedure ApplyFixed;
Var I, K: Integer;
Begin
 If (NFixComp = 0) Or Not FixActive Then Exit;
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do
   If FixM[I, K] Then Begin
    SFrc[I, K] := 0;
    SVel[I, K] := 0;
   End;
End;

{ --------------------------------------------------------------------- }
{ Connectivite                                                          }
{ --------------------------------------------------------------------- }

Var NbWarned: Boolean;

Procedure AddBond(A, B: Integer);
Begin
 If NBond >= MAX_BONDS Then Begin
  InErr(0, 'Trop de liaisons (maximum ' + IStr(MAX_BONDS) + ')');
  Exit;
 End;
 If (AdjN[A] >= MAX_NB) Or (AdjN[B] >= MAX_NB) Then Begin
  If Not NbWarned Then Warn('Plus de ' + IStr(MAX_NB) +
    ' liaisons pour un atome : liaisons en trop ignorees');
  NbWarned := True;
  Exit;
 End;
 Inc(NBond);
 BondA[NBond] := A; BondB[NBond] := B;
 Inc(AdjN[A]); Adj[A, AdjN[A]] := B;
 Inc(AdjN[B]); Adj[B, AdjN[B]] := A;
End;

{ Genere les liaisons a partir des rayons covalents }
Procedure GenerateBonds(Factor: Real);
Var I, J: Integer; DX, DY, DZ, R, Lim: Real;
Begin
 NBond := 0;
 For I := 1 To NAtoms Do AdjN[I] := 0;
 For I := 1 To NAtoms - 1 Do
  For J := I + 1 To NAtoms Do Begin
   DX := SPos[I,1] - SPos[J,1]; DY := SPos[I,2] - SPos[J,2]; DZ := SPos[I,3] -
     SPos[J,3];
   MinImg(DX, DY, DZ);
   R := Sqrt(DX * DX + DY * DY + DZ * DZ);
   Lim := (PtRad[SZ[I]] + PtRad[SZ[J]]) * ANG2BOHR * Factor;
   If (R < Lim) And (R > 0.1) Then AddBond(I, J);
  End;
End;

{ Molecules : composantes connexes du graphe }
Procedure FindMolecules;
Var I, K, Top, Cur, J: Integer; Stack: Array[1..MAX_ATOM] Of Integer;
Begin
 NMol := 0;
 For I := 1 To NAtoms Do SMol[I] := 0;
 For I := 1 To NAtoms Do
  If SMol[I] = 0 Then Begin
   Inc(NMol);
   Top := 1; Stack[1] := I; SMol[I] := NMol;
   While Top > 0 Do Begin
    Cur := Stack[Top]; Dec(Top);
    For K := 1 To AdjN[Cur] Do Begin
     J := Adj[Cur, K];
     If SMol[J] = 0 Then Begin
      SMol[J] := NMol;
      Inc(Top); Stack[Top] := J;
     End;
    End;
   End;
  End;
End;

{ Liste des paires a exclure ou a ponderer (profondeur 1 a 3) }
Var ExWarned: Boolean;

Procedure BuildExclusions;
Var
 I, J, K, D, H, Cur, NQ: Integer;
 Dist: Array[1..MAX_ATOM] Of Integer;
 Q: Array[1..MAX_ATOM] Of Integer;
 Maxd: Integer;
Begin
 Maxd := ExclMax;
 If ExclMaxVdw > Maxd Then Maxd := ExclMaxVdw;
 For I := 1 To NAtoms Do ExN[I] := 0;
 If Maxd = 0 Then Exit;
 For I := 1 To NAtoms Do Begin
  If AdjN[I] = 0 Then Continue;
  For J := 1 To NAtoms Do Dist[J] := 0;
  Dist[I] := -1;
  NQ := 1; Q[1] := I;
  H := 1;
  While H <= NQ Do Begin
   Cur := Q[H]; Inc(H);
   If Dist[Cur] < 0 Then D := 0 Else D := Dist[Cur];
   If D >= Maxd Then Continue;
   For K := 1 To AdjN[Cur] Do Begin
    J := Adj[Cur, K];
    If Dist[J] = 0 Then Begin
     Dist[J] := D + 1;
     Inc(NQ); Q[NQ] := J;
    End;
   End;
  End;
  For J := I + 1 To NAtoms Do
   If (Dist[J] > 0) And (ExN[I] >= MAX_EXCL) Then Begin
    If Not ExWarned Then
      Warn('Liste d''exclusions trop longue pour un atome : paires en trop ' +
        'ignorees');
    ExWarned := True;
   End
   Else If Dist[J] > 0 Then Begin
    Inc(ExN[I]);
    ExJ[I, ExN[I]] := J;
    ExD[I, ExN[I]] := Dist[J];
   End;
 End;
End;

{ Cree la topologie selon la section &TOPOLOGY }
Procedure SetupTopology;
Var St, Sg: Integer; Fm: String; Fac: Real; I: Integer;
 Function Depth(Const S: String; Def: Integer): Integer;
 Begin
  If S = '1-1' Then Depth := 0
  Else If S = '1-2' Then Depth := 1
  Else If S = '1-3' Then Depth := 2
  Else If S = '1-4' Then Depth := 3
  Else Depth := Def;
 End;
Begin
 NBond := 0;
 For I := 1 To NAtoms Do AdjN[I] := 0;
 St := FindSec(SecSub, 'TOPOLOGY', 1);
 Fm := '';
 If St > 0 Then Fm := KwWord(St, 'CONN_FILE_FORMAT', '');
 Sg := 0;
 If St > 0 Then Sg := FindSec(St, 'GENERATE', 1);
 If (Fm = 'GENERATE') Or (Sg > 0) Then Begin
  Fac := 1.1;
  If Sg > 0 Then Fac := KwReal(Sg, 'BONDPARM_FACTOR', 1.1, '');
  GenerateBonds(Fac);
 End;
 FindMolecules;
 ExclMax := 3;
 ExclMaxVdw := 3;
 Scale14Ei := 0;
 Scale14Vdw := 0;
 If St > 0 Then Begin
  ExclMax := Depth(KwWord(St, 'EXCLUDE_EI', '1-4'), 3);
  ExclMaxVdw := Depth(KwWord(St, 'EXCLUDE_VDW', '1-4'), 3);
 End;
 BuildExclusions;
End;

{ Recentre les coordonnees si &CENTER_COORDINATES est present }
Procedure CenterCoords;
Var St, Sc, I, K: Integer; C: TVec3; Tg: TVec3;
Begin
 St := FindSec(SecSub, 'TOPOLOGY', 1);
 If St = 0 Then Exit;
 Sc := FindSec(St, 'CENTER_COORDINATES', 1);
 If Sc = 0 Then Exit;
 For K := 1 To 3 Do Begin
  C[K] := 0;
  For I := 1 To NAtoms Do C[K] := C[K] + SPos[I,K];
  C[K] := C[K] / NAtoms;
  Tg[K] := 0;
  If HasCell Then Tg[K] := 0.5 * (CellH[1,K] + CellH[2,K] + CellH[3,K]);
 End;
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do SPos[I,K] := SPos[I,K] - C[K] + Tg[K];
End;

{ --------------------------------------------------------------------- }
{ Impression des informations sur le systeme                            }
{ --------------------------------------------------------------------- }

Procedure PrintCell;
Var I: Integer;
Begin
 EmitBlank;
 Info('CELL', 'Volume [angstrom^3]:', FR(CellVol * Cu(BOHR_ANG), 12, 3));
 For I := 1 To 3 Do
  Emit(' CELL| Vector ' + Chr(Ord('a') + I - 1) + ' [angstrom]:' +
       FR(CellH[I,1] * BOHR_ANG, 12, 3) + FR(CellH[I,2] * BOHR_ANG, 10, 3) +
       FR(CellH[I,3] * BOHR_ANG, 10, 3) + '    |' + Chr(Ord('a') + I - 1) +
         '| = ' +
       FR(CellLen(I) * BOHR_ANG, 13, 3));
 Info('CELL', 'Angle (b,c), alpha [degree]:', FR(CellAngle(2, 3), 20, 3));
 Info('CELL', 'Angle (a,c), beta  [degree]:', FR(CellAngle(1, 3), 20, 3));
 Info('CELL', 'Angle (a,b), gamma [degree]:', FR(CellAngle(1, 2), 20, 3));
 Info('CELL', 'Requested initial symmetry:', 'NONE');
 Info('CELL', 'Numerically orthorhombic:', YesNo(CellOrtho));
 Info('CELL', 'Periodicity:', PerStr);
End;

Procedure PrintKinds;
Var K: Integer;
Begin
 EmitBlank;
 EmitBlank;
 Emit(' ATOMIC KIND INFORMATION');
 EmitBlank;
 For K := 1 To NKinds Do Begin
  Emit(' ' + ISt(K, 2) + '. Atomic kind: ' + PadR(KName[K], 20) +
    ' Number of atoms: ' + ISt(KCount[K], 6));
  EmitBlank;
  Emit('     Atomic weight:' + FR(KMassA[K] / AMU_AU, 40, 4));
  If KZ[K] > 0 Then Emit('     Element symbol: ' + Rep(' ', 50) +
    PadL(PtSym[KZ[K]], 3));
  If KBasName[K] <> '' Then Begin
   Emit('     Orbital basis set:' + PadL(KBasName[K], 55));
   If (KBas[K] >= 0) And (FEMethod = FE_QS) And (UpStr(KBasName[K]) <>
     UpStr(BasisNames[KBas[K]])) Then
    Emit('     Engine basis set:' + PadL(BasisNames[KBas[K]], 56));
  End;
  EmitBlank;
 End;
End;

Function ZEff(I: Integer; Nuc: Boolean): Real;
Begin
 If Nuc Then ZEff := SZ[I] Else ZEff := SQ[I];
End;

Procedure PrintCoords(Const Title: String; NuclearZ: Boolean);
Var I: Integer;
Begin
 EmitBlank;
 Emit(' MODULE ' + Title + ':  ATOMIC COORDINATES IN angstrom');
 EmitBlank;
 Emit('  Atom  Kind  Element       X           Y           Z          Z(eff) ' +
   ' ' +
   '     Mass');
 EmitBlank;
 For I := 1 To NAtoms Do
  Emit(ISt(I, 7) + ISt(SKind[I], 6) + '     ' + PadR(PtSym[SZ[I]], 2) +
    ISt(SZ[I], 3) +
       FR(SPos[I,1] * BOHR_ANG, 12, 6) + FR(SPos[I,2] * BOHR_ANG, 12, 6) +
       FR(SPos[I,3] * BOHR_ANG, 12, 6) + FR(ZEff(I, NuclearZ), 12, 4) +
         FR(SMass[I] / AMU_AU, 12, 4));
 EmitBlank;
End;

Procedure PrintTotals;
Begin
 EmitBlank;
 Emit(' TOTAL NUMBERS AND MAXIMUM NUMBERS');
 EmitBlank;
 Emit('  Total number of            - Atomic kinds:' + ISt(NKinds, 30));
 Emit('                             - Atoms:' + ISt(NAtoms, 37));
 Emit('                             - Molecules:' + ISt(NMol, 33));
 Emit('                             - Bonds:' + ISt(NBond, 37));
 EmitBlank;
End;

{ Fichier XYZ : une image }
Procedure WriteXYZFrame(Var Fo: Text; Const Comment: String);
Var I: Integer;
Begin
 WriteLn(Fo, NAtoms:7);
 WriteLn(Fo, Comment);
 For I := 1 To NAtoms Do
  WriteLn(Fo, ' ', PadR(PtSym[SZ[I]], 2), FR(SPos[I,1] * BOHR_ANG, 20, 10),
          FR(SPos[I,2] * BOHR_ANG, 20, 10), FR(SPos[I,3] * BOHR_ANG, 20, 10));
End;


{ --------------------------------------------------------------------- }
{ Moteur FIST : champ de force classique                                }
{ --------------------------------------------------------------------- }

Const
 {$IFDEF FPC}
 ERFC_N   = 6144;       { points de la table erfc (pas 1/1024) }
 ERFC_H   = 1.0 / 1024.0;
 {$ELSE}
 ERFC_N   = 384;
 ERFC_H   = 1.0 / 64.0;
 {$ENDIF}
 SQRT_PI  = 1.7724538509055160;

Var
 { table erfc }
 ErfcF, ErfcD : Array[0..ERFC_N + 1] Of Real;
 ErfcReady    : Boolean;

 { potentiels non lies }
 LJType  : Array[1..MAX_KIND, 1..MAX_KIND] Of Integer;  { 0 aucun, 1 LJ, 2
   Williams }
 LJEps   : Array[1..MAX_KIND, 1..MAX_KIND] Of Real;
 LJSig   : Array[1..MAX_KIND, 1..MAX_KIND] Of Real;
 LJRc    : Array[1..MAX_KIND, 1..MAX_KIND] Of Real;
 LJSh    : Array[1..MAX_KIND, 1..MAX_KIND] Of Real;
 WilA, WilB, WilC : Array[1..MAX_KIND, 1..MAX_KIND] Of Real;
 ShiftCut : Boolean;
 RcutNB   : Real;
 RcMax    : Real;
 MMCharged : Boolean;
 NPairsLJ : Integer;

 { Ewald }
 EwaldOn  : Boolean;
 EwAlpha, EwRc, EwKc : Real;

 { parametres lies }
 NBP      : Integer;
 BPA, BPB, BPKind : Array[1..MAX_PTAB] Of Integer;
 BPK, BPR0 : Array[1..MAX_PTAB] Of Real;
 NAP      : Integer;
 APA, APB, APC, APKind : Array[1..MAX_PTAB] Of Integer;
 APK, APT0 : Array[1..MAX_PTAB] Of Real;
 NTP      : Integer;
 TPA, TPB, TPC, TPD, TPM : Array[1..MAX_PTAB] Of Integer;
 TPK, TPPhi : Array[1..MAX_PTAB] Of Real;

 { instances }
 BondP    : Array[1..MAX_BONDS] Of Integer;
 NBend    : Integer;
 BendA, BendB, BendC, BendP : Array[1..MAX_BENDS] Of Integer;
 NTors    : Integer;
 TorA, TorB, TorC, TorD, TorP : Array[1..MAX_TORS] Of Integer;

 { resultats }
 MMEBond, MMEBend, MMETors, MMEVdw, MMECoul : Real;
 MMEEwReal, MMEEwRecip, MMEEwSelf : Real;
 DEps     : TMat3;        { derivee de l'energie par rapport a la deformation }

Procedure InitErfc;
Var I, N: Integer; X, Terme, Sm, T: Real;
Begin
 For I := 0 To ERFC_N + 1 Do Begin
  X := I * ERFC_H;
  ErfcD[I] := -2.0 / SQRT_PI * Exp(-X * X);
 End;
 { erfc(x) par integration depuis l'infini : utilise le developpement
   en fraction continue pour les grands x puis la serie pour les petits }
 For I := 0 To ERFC_N + 1 Do Begin
  X := I * ERFC_H;
  If X < 2.2 Then Begin
   { erf(x) = 2/sqrt(pi) exp(-x^2) sum 2^n x^(2n+1) / (2n+1)!! }
   Terme := X; Sm := X; N := 0;
   While Terme > 1E-18 * Sm Do Begin
    Inc(N);
    Terme := Terme * 2 * X * X / (2 * N + 1);
    Sm := Sm + Terme;
    If N > 200 Then Break;
   End;
   ErfcF[I] := 1 - 2.0 / SQRT_PI * Exp(-X * X) * Sm;
  End
  Else Begin
   T := X;
   For N := 120 DownTo 1 Do T := X + (N / 2.0) / T;
   ErfcF[I] := Exp(-X * X) / (SQRT_PI * T);
  End;
 End;
 ErfcReady := True;
End;

{ erfc(x) par interpolation d'Hermite cubique }
Function ErfcT(X: Real): Real;
Var U, T, T2, T3, H00, H10, H01, H11: Real; I: Integer;
Begin
 If X <= 0 Then Begin ErfcT := 1; Exit; End;
 U := X / ERFC_H;
 If U >= ERFC_N Then Begin ErfcT := 0; Exit; End;
 I := Trunc(U);
 T := U - I;
 T2 := T * T; T3 := T2 * T;
 H00 := 2 * T3 - 3 * T2 + 1;
 H10 := T3 - 2 * T2 + T;
 H01 := -2 * T3 + 3 * T2;
 H11 := T3 - T2;
 ErfcT := H00 * ErfcF[I] + H10 * ERFC_H * ErfcD[I] + H01 * ErfcF[I + 1] + H11 *
   ERFC_H * ErfcD[I + 1];
End;

{ --------------------------------------------------------------------- }
{ Lecture du champ de force                                             }
{ --------------------------------------------------------------------- }

Procedure MMSetup;
Var
 MmL: LongInt;
 Sff, Sp, Sc, Sn, I, J, N, A, B, C, D, Ka, Kb, Kc, Kd, Cnt, Sw, Se: Integer;
 Nm, Kd2, Wrd: String;
 Eps, Sig, Rc, V: Real;
 Mixed: Boolean;
 Names: Array[1..4] Of String;
 Function AtomsOf(Sc: Integer; Cnt: Integer): Boolean;
 Var Q: Integer; U, Un2: String; Fd2: Boolean;
 Begin
  AtomsOf := False;
  For Q := 1 To 4 Do Names[Q] := '';
  U := KwRaw(Sc, 'ATOMS', Un2, Fd2);
  For Q := 1 To Cnt Do Begin
   Names[Q] := GetTok(U);
   If Names[Q] = '' Then Begin
    InErr(SecLine[Sc], 'ATOMS doit contenir ' + IStr(Cnt) +
      ' sortes d''atomes');
    Exit;
   End;
  End;
  AtomsOf := True;
 End;
 Function KindOrWarn(Const S: String; Sc: Integer): Integer;
 Var Z: Integer;
 Begin
  Z := FindKind(S);
  If Z = 0 Then Warn('Sorte d''atome ' + S +
    ' inconnue dans le champ de force (ligne ' +
                     IStr(SecLine[Sc]) + ')');
  KindOrWarn := Z;
 End;
Begin
 FillChar(LJType, SizeOf(LJType), 0);
 FillChar(LJEps, SizeOf(LJEps), 0);
 FillChar(LJSig, SizeOf(LJSig), 0);
 FillChar(LJRc, SizeOf(LJRc), 0);
 FillChar(LJSh, SizeOf(LJSh), 0);
 FillChar(WilA, SizeOf(WilA), 0);
 FillChar(WilB, SizeOf(WilB), 0);
 FillChar(WilC, SizeOf(WilC), 0);
 NBP := 0; NAP := 0; NTP := 0;
 NBend := 0; NTors := 0;
 ShiftCut := True;
 RcutNB := 10.0 * ANG2BOHR;
 Sff := FindSec(SecMM, 'FORCEFIELD', 1);
 If Sff = 0 Then Begin
  InErr(SecLine[SecMM],
    'La section &FORCEFIELD est obligatoire pour METHOD FIST');
  Exit;
 End;
 ShiftCut := KwBool(Sff, 'SHIFT_CUTOFF', True);
 Scale14Ei := KwReal(Sff, 'EI_SCALE14', 0.0, '');
 Scale14Vdw := KwReal(Sff, 'VDW_SCALE14', 0.0, '');
 Sp := FindSec(Sff, 'SPLINE', 1);
 If Sp > 0 Then RcutNB := KwReal(Sp, 'RCUT_NB', 10.0, 'angstrom');
 If Not KwBool(Sff, 'DO_NONBONDED', True) Then RcutNB := 0;
 { charges }
 For I := 1 To NKinds Do KQ[I] := 0;
 N := CountSec(Sff, 'CHARGE');
 For I := 1 To N Do Begin
  Sc := FindSec(Sff, 'CHARGE', I);
  Nm := KwStr(Sc, 'ATOM', '');
  A := FindKind(Nm);
  If A = 0 Then Begin Warn('Charge pour la sorte inconnue ' + Nm); Continue;
    End;
  KQ[A] := KwReal(Sc, 'CHARGE', 0.0, '');
 End;
 { Lennard-Jones et Williams }
 Sn := FindSec(Sff, 'NONBONDED', 1);
 If Sn > 0 Then Begin
  N := CountSec(Sn, 'LENNARD-JONES');
  For I := 1 To N Do Begin
   Sc := FindSec(Sn, 'LENNARD-JONES', I);
   If Not AtomsOf(Sc, 2) Then Continue;
   Ka := KindOrWarn(Names[1], Sc); Kb := KindOrWarn(Names[2], Sc);
   If (Ka = 0) Or (Kb = 0) Then Continue;
   Eps := KwReal(Sc, 'EPSILON', 0.0, 'hartree');
   Sig := KwReal(Sc, 'SIGMA', 0.0, 'bohr');
   Rc := KwReal(Sc, 'RCUT', 0.0, 'bohr');
   If Rc <= 0 Then Rc := RcutNB;
   LJType[Ka, Kb] := 1; LJType[Kb, Ka] := 1;
   LJEps[Ka, Kb] := Eps; LJEps[Kb, Ka] := Eps;
   LJSig[Ka, Kb] := Sig; LJSig[Kb, Ka] := Sig;
   LJRc[Ka, Kb] := Rc; LJRc[Kb, Ka] := Rc;
  End;
  N := CountSec(Sn, 'WILLIAMS');
  For I := 1 To N Do Begin
   Sc := FindSec(Sn, 'WILLIAMS', I);
   If Not AtomsOf(Sc, 2) Then Continue;
   Ka := KindOrWarn(Names[1], Sc); Kb := KindOrWarn(Names[2], Sc);
   If (Ka = 0) Or (Kb = 0) Then Continue;
   Rc := KwReal(Sc, 'RCUT', 0.0, 'bohr');
   If Rc <= 0 Then Rc := RcutNB;
   LJType[Ka, Kb] := 2; LJType[Kb, Ka] := 2;
   WilA[Ka, Kb] := KwReal(Sc, 'A', 0.0, 'hartree'); WilA[Kb, Ka] := WilA[Ka,
     Kb];
   WilB[Ka, Kb] := KwReal(Sc, 'B', 0.0, 'bohr^-1'); WilB[Kb, Ka] := WilB[Ka,
     Kb];
   WilC[Ka, Kb] := KwReal(Sc, 'C', 0.0, 'hartree*bohr^6'); WilC[Kb, Ka] :=
     WilC[Ka, Kb];
   LJRc[Ka, Kb] := Rc; LJRc[Kb, Ka] := Rc;
  End;
 End;
 { regles de melange de Lorentz-Berthelot pour les paires manquantes }
 Mixed := False;
 For A := 1 To NKinds Do
  For B := A + 1 To NKinds Do
   If (LJType[A, B] = 0) And (LJType[A, A] = 1) And (LJType[B, B] = 1) Then
     Begin
    LJType[A, B] := 1; LJType[B, A] := 1;
    LJEps[A, B] := Sqrt(LJEps[A, A] * LJEps[B, B]); LJEps[B, A] := LJEps[A, B];
    LJSig[A, B] := 0.5 * (LJSig[A, A] + LJSig[B, B]); LJSig[B, A] := LJSig[A,
      B];
    LJRc[A, B] := RMin(LJRc[A, A], LJRc[B, B]); LJRc[B, A] := LJRc[A, B];
    Mixed := True;
   End;
 If Mixed Then
   Warn('Parametres Lennard-Jones croises generes (regles de ' +
     'Lorentz-Berthelot)');
 { decalage des potentiels a la coupure }
 RcMax := 0;
 NPairsLJ := 0;
 For A := 1 To NKinds Do
  For B := 1 To NKinds Do Begin
   If LJType[A, B] = 0 Then Continue;
   If B >= A Then Inc(NPairsLJ);
   Rc := LJRc[A, B];
   If Rc > RcMax Then RcMax := Rc;
   If ShiftCut Then Begin
    If LJType[A, B] = 1 Then Begin
     V := Cu(Sqr(LJSig[A, B] / Rc));
     LJSh[A, B] := 4 * LJEps[A, B] * (V * V - V);
    End
    Else
     LJSh[A, B] := WilA[A, B] * Exp(-WilB[A, B] * Rc) - WilC[A, B] / Cu(Rc *
       Rc);
   End;
  End;
 { liaisons }
 N := CountSec(Sff, 'BOND');
 For I := 1 To N Do Begin
  Sc := FindSec(Sff, 'BOND', I);
  If Not AtomsOf(Sc, 2) Then Continue;
  Ka := KindOrWarn(Names[1], Sc); Kb := KindOrWarn(Names[2], Sc);
  If (Ka = 0) Or (Kb = 0) Then Continue;
  If NBP >= MAX_PTAB Then Begin InErr(SecLine[Sc],
    'Trop de parametres de liaison'); Break; End;
  Inc(NBP);
  BPA[NBP] := Ka; BPB[NBP] := Kb;
  Kd2 := KwWord(Sc, 'KIND', 'HARMONIC');
  If Kd2 = 'HARMONIC' Then BPKind[NBP] := 0
  Else If (Kd2 = 'CHARMM') Or (Kd2 = 'AMBER') Then BPKind[NBP] := 1
  Else Begin InErr(SecLine[Sc], 'KIND de liaison non pris en charge : ' + Kd2);
    BPKind[NBP] := 0; End;
  BPK[NBP] := KwReal(Sc, 'K', 0.0, 'hartree*bohr^-2');
  BPR0[NBP] := KwReal(Sc, 'R0', 0.0, 'bohr');
 End;
 { angles }
 N := CountSec(Sff, 'BEND');
 For I := 1 To N Do Begin
  Sc := FindSec(Sff, 'BEND', I);
  If Not AtomsOf(Sc, 3) Then Continue;
  Ka := KindOrWarn(Names[1], Sc); Kb := KindOrWarn(Names[2], Sc); Kc :=
    KindOrWarn(Names[3], Sc);
  If (Ka = 0) Or (Kb = 0) Or (Kc = 0) Then Continue;
  If NAP >= MAX_PTAB Then Begin InErr(SecLine[Sc],
    'Trop de parametres d''angle'); Break; End;
  Inc(NAP);
  APA[NAP] := Ka; APB[NAP] := Kb; APC[NAP] := Kc;
  Kd2 := KwWord(Sc, 'KIND', 'HARMONIC');
  If Kd2 = 'HARMONIC' Then APKind[NAP] := 0
  Else If (Kd2 = 'CHARMM') Or (Kd2 = 'AMBER') Then APKind[NAP] := 1
  Else Begin InErr(SecLine[Sc], 'KIND d''angle non pris en charge : ' + Kd2);
    APKind[NAP] := 0; End;
  APK[NAP] := KwReal(Sc, 'K', 0.0, 'hartree*rad^-2');
  APT0[NAP] := KwReal(Sc, 'THETA0', 0.0, 'rad');
 End;
 { torsions }
 N := CountSec(Sff, 'TORSION');
 For I := 1 To N Do Begin
  Sc := FindSec(Sff, 'TORSION', I);
  If Not AtomsOf(Sc, 4) Then Continue;
  Ka := KindOrWarn(Names[1], Sc); Kb := KindOrWarn(Names[2], Sc);
  Kc := KindOrWarn(Names[3], Sc); Kd := KindOrWarn(Names[4], Sc);
  If (Ka = 0) Or (Kb = 0) Or (Kc = 0) Or (Kd = 0) Then Continue;
  If NTP >= MAX_PTAB Then Begin InErr(SecLine[Sc],
    'Trop de parametres de torsion'); Break; End;
  Inc(NTP);
  TPA[NTP] := Ka; TPB[NTP] := Kb; TPC[NTP] := Kc; TPD[NTP] := Kd;
  TPK[NTP] := KwReal(Sc, 'K', 0.0, 'hartree');
  TPPhi[NTP] := KwReal(Sc, 'PHI0', 0.0, 'rad');
  MmL := KwInt(Sc, 'M', 1);
  If (MmL < 0) Or (MmL > 12) Then Begin
   InErr(SecLine[Sc], 'M doit etre compris entre 0 et 12');
   MmL := 1;
  End;
  TPM[NTP] := MmL;
  Kd2 := KwWord(Sc, 'KIND', 'CHARMM');
  If (Kd2 <> 'CHARMM') And (Kd2 <> 'AMBER') Then InErr(SecLine[Sc],
    'KIND de torsion non pris en charge : ' + Kd2);
 End;
 { charges et Coulomb }
 MMCharged := False;
 For I := 1 To NKinds Do If Abs(KQ[I]) > 1E-12 Then MMCharged := True;
 For I := 1 To NAtoms Do SQ[I] := KQ[SKind[I]];
 { Ewald }
 EwaldOn := False;
 Sw := FindSec(SecMM, 'POISSON', 1);
 Se := 0;
 If Sw > 0 Then Se := FindSec(Sw, 'EWALD', 1);
 Wrd := 'SPME';
 If Se > 0 Then Wrd := KwWord(Se, 'EWALD_TYPE', 'SPME');
 EwAlpha := 0;
 If Se > 0 Then EwAlpha := KwReal(Se, 'ALPHA', 0.0, 'angstrom^-1');
 If (Per[1] Or Per[2] Or Per[3]) And MMCharged And (Wrd <> 'NONE') Then EwaldOn
   := True;
 If EwaldOn And Not (Per[1] And Per[2] And Per[3]) Then
  InErr(SecLine[SecMM],
    'La sommation d''Ewald exige une periodicite XYZ (ou EWALD_TYPE NONE)');
 { verification de la cellule pour les potentiels non lies }
 If HasCell Then
  For I := 1 To 3 Do
   If Per[I] And (RcMax > 0) And (CellWidth(I) < 2 * RcMax - 1E-9) Then
    InErr(SecLine[SecMM], 'La coupure (' + RStr(RcMax * BOHR_ANG, 3) +
      ' angstrom) depasse la moitie de la largeur de la cellule (' +
          RStr(CellWidth(I) * 0.5 * BOHR_ANG, 3) +
            ' angstrom) selon la direction ' + IStr(I));
 If EwaldOn Then Begin
  If RcutNB <= 0 Then RcutNB := 10.0 * ANG2BOHR;
  EwRc := RcutNB;
  For I := 1 To 3 Do
   If Per[I] And (0.5 * CellWidth(I) - 1E-9 < EwRc) Then EwRc := 0.5 *
     CellWidth(I) - 1E-9;
  If EwAlpha = 0 Then EwAlpha := 4.0 / EwRc;
  
  EwKc := 2 * EwAlpha * Sqrt(-Ln(1E-9));
  If RcMax < EwRc Then RcMax := EwRc;
 End;
 { instances d'angles et de torsions a partir de la connectivite }
 For I := 1 To NBond Do Begin
  BondP[I] := 0;
  A := SKind[BondA[I]]; B := SKind[BondB[I]];
  For J := 1 To NBP Do
   If ((BPA[J] = A) And (BPB[J] = B)) Or ((BPA[J] = B) And (BPB[J] = A)) Then
     Begin
    BondP[I] := J;
    Break;
   End;
  If (BondP[I] = 0) And (NBP > 0) Then
   Warn('Pas de parametre de liaison pour ' + KName[A] + '-' + KName[B]);
 End;
 For C := 1 To NAtoms Do
  For I := 1 To AdjN[C] - 1 Do
   For J := I + 1 To AdjN[C] Do Begin
    A := Adj[C, I]; B := Adj[C, J];
    For D := 1 To NAP Do
     If (APB[D] = SKind[C]) And (((APA[D] = SKind[A]) And (APC[D] = SKind[B]))
       Or
        ((APA[D] = SKind[B]) And (APC[D] = SKind[A]))) Then Begin
      If NBend >= MAX_BENDS Then Begin InErr(0, 'Trop d''angles'); Exit; End;
      Inc(NBend);
      BendA[NBend] := A; BendB[NBend] := C; BendC[NBend] := B; BendP[NBend] :=
        D;
      Break;
     End;
   End;
 For I := 1 To NBond Do Begin
  B := BondA[I]; C := BondB[I];
  For J := 1 To AdjN[B] Do Begin
   A := Adj[B, J];
   If A = C Then Continue;
   For Cnt := 1 To AdjN[C] Do Begin
    D := Adj[C, Cnt];
    If (D = B) Or (D = A) Then Continue;
    For N := 1 To NTP Do
     If ((TPA[N] = SKind[A]) And (TPB[N] = SKind[B]) And (TPC[N] = SKind[C]) And
       (TPD[N] = SKind[D])) Or
        ((TPA[N] = SKind[D]) And (TPB[N] = SKind[C]) And (TPC[N] = SKind[B]) And
          (TPD[N] = SKind[A])) Then Begin
      If NTors >= MAX_TORS Then Begin InErr(0, 'Trop de torsions'); Exit; End;
      Inc(NTors);
      TorA[NTors] := A; TorB[NTors] := B; TorC[NTors] := C; TorD[NTors] := D;
        TorP[NTors] := N;
     End;
   End;
  End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Energie et forces                                                     }
{ --------------------------------------------------------------------- }

Procedure AddF(A: Integer; X, Y, Z: Real);
Begin
 SFrc[A,1] := SFrc[A,1] + X;
 SFrc[A,2] := SFrc[A,2] + Y;
 SFrc[A,3] := SFrc[A,3] + Z;
End;

{ virial : dE/d(eps_ab) -= r_a * F_b }
Procedure AddVir(Const R: TVec3; FX, FY, FZ: Real);
Var K: Integer;
Begin
 For K := 1 To 3 Do Begin
  DEps[K,1] := DEps[K,1] - R[K] * FX;
  DEps[K,2] := DEps[K,2] - R[K] * FY;
  DEps[K,3] := DEps[K,3] - R[K] * FZ;
 End;
End;

Procedure PairTerms(NeedF: Boolean);
Var
 I, J, Ki, Kj, E: Integer;
 DX, DY, DZ, R2, R, Fs, Qq, Sv, Se, Rc, X6, X12, Ex, Ec, A2, Er: Real;
 MkEi, MkVdw : Array[1..MAX_ATOM] Of Real;
 RcLim2: Real;
 Rv: TVec3;
 Fx, Fy, Fz: Real;
 Dp: Integer;
Begin
 For I := 1 To NAtoms Do Begin MkEi[I] := 1; MkVdw[I] := 1; End;
 If MMCharged And Not EwaldOn Then RcLim2 := 1E30
 Else RcLim2 := RcMax * RcMax;
 For I := 1 To NAtoms - 1 Do Begin
  Ki := SKind[I];
  For E := 1 To ExN[I] Do Begin
   J := ExJ[I, E];
   Dp := ExD[I, E];
   If Dp <= ExclMax Then Begin
    If Dp = 3 Then MkEi[J] := Scale14Ei Else MkEi[J] := 0;
   End;
   If Dp <= ExclMaxVdw Then Begin
    If Dp = 3 Then MkVdw[J] := Scale14Vdw Else MkVdw[J] := 0;
   End;
  End;
  For J := I + 1 To NAtoms Do Begin
   DX := SPos[I,1] - SPos[J,1]; DY := SPos[I,2] - SPos[J,2]; DZ := SPos[I,3] -
     SPos[J,3];
   MinImg(DX, DY, DZ);
   R2 := DX * DX + DY * DY + DZ * DZ;
   Sv := MkVdw[J]; Se := MkEi[J];
   If (R2 > RcLim2) And (Se = 1) Then Continue;
   If R2 < 1E-12 Then Continue;
   Kj := SKind[J];
   Fs := 0;
   { van der Waals }
   If (Sv <> 0) And (LJType[Ki, Kj] > 0) Then Begin
    Rc := LJRc[Ki, Kj];
    If R2 < Rc * Rc Then Begin
     R := Sqrt(R2);
     If LJType[Ki, Kj] = 1 Then Begin
      X6 := Cu(LJSig[Ki, Kj] * LJSig[Ki, Kj] / R2);
      X12 := X6 * X6;
      MMEVdw := MMEVdw + Sv * (4 * LJEps[Ki, Kj] * (X12 - X6) - LJSh[Ki, Kj]);
      Fs := Fs + Sv * 4 * LJEps[Ki, Kj] * (12 * X12 - 6 * X6) / R2;
     End
     Else Begin
      Ex := WilA[Ki, Kj] * Exp(-WilB[Ki, Kj] * R);
      X6 := WilC[Ki, Kj] / Cu(R2);
      MMEVdw := MMEVdw + Sv * (Ex - X6 - LJSh[Ki, Kj]);
      Fs := Fs + Sv * (WilB[Ki, Kj] * Ex / R - 6 * X6 / R2);
     End;
    End;
   End;
   { Coulomb }
   Qq := SQ[I] * SQ[J];
   If (Abs(Qq) > 1E-14) And (Se <> 0) Then Begin
    R := Sqrt(R2);
    If EwaldOn Then Begin
     If (R2 < EwRc * EwRc) Or (Se <> 1) Then Begin
      Ec := ErfcT(EwAlpha * R);
      Er := Qq * (Ec - (1 - Se)) / R;
      MMEEwReal := MMEEwReal + Er;
      A2 := EwAlpha * EwAlpha;
      { -dE/dr / r }
      Fs := Fs + Qq * ((Ec - (1 - Se)) / (R2 * R) + 2 * EwAlpha / SQRT_PI *
        Exp(-A2 * R2) / R2);
     End;
    End
    Else Begin
     MMECoul := MMECoul + Se * Qq / R;
     Fs := Fs + Se * Qq / (R2 * R);
    End;
   End
   Else If (Abs(Qq) > 1E-14) And (Se = 0) And EwaldOn Then Begin
    { paire exclue : correction -q q erf(alpha r) / r }
    R := Sqrt(R2);
    Ec := ErfcT(EwAlpha * R);
    Er := Qq * (Ec - 1) / R;
    MMEEwReal := MMEEwReal + Er;
    A2 := EwAlpha * EwAlpha;
    Fs := Fs + Qq * ((Ec - 1) / (R2 * R) + 2 * EwAlpha / SQRT_PI * Exp(-A2 * R2)
      / R2);
   End;
   If Fs <> 0 Then Begin
    Fx := Fs * DX; Fy := Fs * DY; Fz := Fs * DZ;
    AddF(I, Fx, Fy, Fz);
    AddF(J, -Fx, -Fy, -Fz);
    Rv[1] := DX; Rv[2] := DY; Rv[3] := DZ;
    AddVir(Rv, Fx, Fy, Fz);
   End;
  End;
  For E := 1 To ExN[I] Do Begin
   MkEi[ExJ[I, E]] := 1;
   MkVdw[ExJ[I, E]] := 1;
  End;
 End;
End;

{ Sommation d'Ewald : partie reciproque, auto-interaction et fond neutralisant }
Procedure EwaldRecip(NeedF: Boolean);
Var
 N, A, J, N1, N2, N3, M1, M2, M3: Integer;
 Len: Array[1..3] Of Real;
 NMax: Array[1..3] Of Integer;
 Cs, Sn: Array[1..3] Of PRealArr;
 Sz, I: LongInt;
 S: Array[1..3] Of Real;
 X, Y, Z, TwoPi, K2, Af, Sc, Ss, Pref, EK, Gk, C1, S1, C2, S2, C3, S3, Cp, Sp,
   Th, Dd, Alp2: Real;
 PC, PS: PRealArr;
 Sgn2, Sgn3: Real;
 HK: TVec3;
 Qtot, Qsq, Ebg: Real;
Begin
 N := NAtoms;
 TwoPi := 2 * Pi;
 Alp2 := EwAlpha * EwAlpha;
 Qtot := 0; Qsq := 0;
 For I := 1 To N Do Begin Qtot := Qtot + SQ[I]; Qsq := Qsq + SQ[I] * SQ[I]; End;
 MMEEwSelf := -EwAlpha / SQRT_PI * Qsq;
 Ebg := -Pi * Qtot * Qtot / (2 * CellVol * Alp2);
 MMEEwSelf := MMEEwSelf + Ebg;
 { le fond neutralisant depend de 1/V : dE/d(eps) = -Ebg delta }
 DEps[1,1] := DEps[1,1] - Ebg; DEps[2,2] := DEps[2,2] - Ebg; DEps[3,3] :=
   DEps[3,3] - Ebg;
 For A := 1 To 3 Do Begin
  If Per[A] Then Begin
   Len[A] := CellLen(A);
   NMax[A] := Trunc(EwKc * Len[A] / TwoPi) + 1;
  End
  Else NMax[A] := 0;
 End;
 For A := 1 To 3 Do Begin
  Sz := LongInt(NMax[A] + 1) * N * SizeOf(Real);
  GetMem(Cs[A], Sz);
  GetMem(Sn[A], Sz);
 End;
 GetMem(PC, N * SizeOf(Real));
 GetMem(PS, N * SizeOf(Real));
 For J := 1 To N Do Begin
  X := SPos[J,1]; Y := SPos[J,2]; Z := SPos[J,3];
  For A := 1 To 3 Do Begin
   S[A] := X * CellHI[1,A] + Y * CellHI[2,A] + Z * CellHI[3,A];
   Th := TwoPi * S[A];
   Cs[A]^[J - 1] := 1; Sn[A]^[J - 1] := 0;
   If NMax[A] >= 1 Then Begin
    Cs[A]^[LongInt(N) + J - 1] := Cos(Th);
    Sn[A]^[LongInt(N) + J - 1] := Sin(Th);
   End;
   For N1 := 2 To NMax[A] Do Begin
    C1 := Cs[A]^[LongInt(N) + J - 1]; S1 := Sn[A]^[LongInt(N) + J - 1];
    Cp := Cs[A]^[LongInt(N1 - 1) * N + J - 1]; Sp := Sn[A]^[LongInt(N1 - 1) * N
      + J - 1];
    Cs[A]^[LongInt(N1) * N + J - 1] := Cp * C1 - Sp * S1;
    Sn[A]^[LongInt(N1) * N + J - 1] := Sp * C1 + Cp * S1;
   End;
  End;
 End;
 Pref := 4 * Pi / CellVol;
 MMEEwRecip := 0;
 For N1 := 0 To NMax[1] Do
  For N2 := -NMax[2] To NMax[2] Do Begin
   If (N1 = 0) And (N2 < 0) Then Continue;
   For N3 := -NMax[3] To NMax[3] Do Begin
    If (N1 = 0) And (N2 = 0) And (N3 <= 0) Then Continue;
    For A := 1 To 3 Do Begin
     HK[A] := TwoPi * (N1 * CellHI[A,1] + N2 * CellHI[A,2] + N3 * CellHI[A,3]);
    End;
    K2 := HK[1] * HK[1] + HK[2] * HK[2] + HK[3] * HK[3];
    If K2 > EwKc * EwKc Then Continue;
    Sgn2 := 1; If N2 < 0 Then Sgn2 := -1;
    Sgn3 := 1; If N3 < 0 Then Sgn3 := -1;
    M2 := Abs(N2); M3 := Abs(N3);
    Sc := 0; Ss := 0;
    For J := 1 To N Do Begin
     C1 := Cs[1]^[LongInt(N1) * N + J - 1]; S1 := Sn[1]^[LongInt(N1) * N + J -
       1];
     C2 := Cs[2]^[LongInt(M2) * N + J - 1]; S2 := Sgn2 * Sn[2]^[LongInt(M2) * N
       + J - 1];
     C3 := Cs[3]^[LongInt(M3) * N + J - 1]; S3 := Sgn3 * Sn[3]^[LongInt(M3) * N
       + J - 1];
     Cp := C1 * C2 - S1 * S2;
     Sp := S1 * C2 + C1 * S2;
     PC^[J - 1] := Cp * C3 - Sp * S3;
     PS^[J - 1] := Sp * C3 + Cp * S3;
     Sc := Sc + SQ[J] * PC^[J - 1];
     Ss := Ss + SQ[J] * PS^[J - 1];
    End;
    Af := Exp(-K2 / (4 * Alp2)) / K2;
    EK := Pref * Af * (Sc * Sc + Ss * Ss);
    MMEEwRecip := MMEEwRecip + EK;
    { virial reciproque : E_k [ -delta + 2 (1/k2 + 1/(4 a2)) k_a k_b ] }
    Gk := 2 * (1 / K2 + 1 / (4 * Alp2));
    For A := 1 To 3 Do Begin
     DEps[A,A] := DEps[A,A] - EK;
     For M1 := 1 To 3 Do DEps[A,M1] := DEps[A,M1] + EK * Gk * HK[A] * HK[M1];
    End;
    If NeedF Then
     For J := 1 To N Do Begin
      Dd := 2 * Pref * Af * SQ[J] * (Ss * PC^[J - 1] - Sc * PS^[J - 1]);
      AddF(J, -Dd * HK[1], -Dd * HK[2], -Dd * HK[3]);
     End;
   End;
  End;
 For A := 1 To 3 Do Begin
  FreeMem(Cs[A], LongInt(NMax[A] + 1) * N * SizeOf(Real));
  FreeMem(Sn[A], LongInt(NMax[A] + 1) * N * SizeOf(Real));
 End;
 FreeMem(PC, N * SizeOf(Real));
 FreeMem(PS, N * SizeOf(Real));
End;

Procedure Cross3(Const A, B: TVec3; Var C: TVec3);
Begin
 C[1] := A[2] * B[3] - A[3] * B[2];
 C[2] := A[3] * B[1] - A[1] * B[3];
 C[3] := A[1] * B[2] - A[2] * B[1];
End;

Function Dot3(Const A, B: TVec3): Real;
Begin
 Dot3 := A[1] * B[1] + A[2] * B[2] + A[3] * B[3];
End;

Procedure BondedTerms(NeedF: Boolean);
Var
 I, P, Q: Integer;
 U, V, Zero: TVec3;
 R, Dv, Ct, St, Th, Dedth, Lu, Lv: Real;
 Fi, Fk: TVec3;
 B1, B2, B3, Fv, G, H, A, B, Fa, Fb, Fc, Fd: TVec3;
 Ga, A2, B2n, Phi, CosP, SinP, Dedp, FG, HG: Real;
 Rel: TVec3;
Begin
 Zero[1] := 0; Zero[2] := 0; Zero[3] := 0;
 { liaisons }
 For I := 1 To NBond Do Begin
  P := BondP[I];
  If P = 0 Then Continue;
  U[1] := SPos[BondA[I],1] - SPos[BondB[I],1];
  U[2] := SPos[BondA[I],2] - SPos[BondB[I],2];
  U[3] := SPos[BondA[I],3] - SPos[BondB[I],3];
  MinImg(U[1], U[2], U[3]);
  R := Sqrt(Dot3(U, U));
  If R < 1E-10 Then Continue;
  If BPKind[P] = 0 Then Begin
   MMEBond := MMEBond + 0.5 * BPK[P] * Sqr(R - BPR0[P]);
   Dv := BPK[P] * (R - BPR0[P]);
  End
  Else Begin
   MMEBond := MMEBond + BPK[P] * Sqr(R - BPR0[P]);
   Dv := 2 * BPK[P] * (R - BPR0[P]);
  End;
  Fi[1] := -Dv * U[1] / R; Fi[2] := -Dv * U[2] / R; Fi[3] := -Dv * U[3] / R;
  AddF(BondA[I], Fi[1], Fi[2], Fi[3]);
  AddF(BondB[I], -Fi[1], -Fi[2], -Fi[3]);
  AddVir(U, Fi[1], Fi[2], Fi[3]);
 End;
 { angles }
 For I := 1 To NBend Do Begin
  P := BendP[I];
  U[1] := SPos[BendA[I],1] - SPos[BendB[I],1];
  U[2] := SPos[BendA[I],2] - SPos[BendB[I],2];
  U[3] := SPos[BendA[I],3] - SPos[BendB[I],3];
  MinImg(U[1], U[2], U[3]);
  V[1] := SPos[BendC[I],1] - SPos[BendB[I],1];
  V[2] := SPos[BendC[I],2] - SPos[BendB[I],2];
  V[3] := SPos[BendC[I],3] - SPos[BendB[I],3];
  MinImg(V[1], V[2], V[3]);
  Lu := Sqrt(Dot3(U, U)); Lv := Sqrt(Dot3(V, V));
  If (Lu < 1E-10) Or (Lv < 1E-10) Then Continue;
  Ct := Dot3(U, V) / (Lu * Lv);
  If Ct > 1 Then Ct := 1;
  If Ct < -1 Then Ct := -1;
  Th := ArcCos(Ct);
  St := Sqrt(1 - Ct * Ct);
  If St < 1E-8 Then St := 1E-8;
  If APKind[P] = 0 Then Begin
   MMEBend := MMEBend + 0.5 * APK[P] * Sqr(Th - APT0[P]);
   Dedth := APK[P] * (Th - APT0[P]);
  End
  Else Begin
   MMEBend := MMEBend + APK[P] * Sqr(Th - APT0[P]);
   Dedth := 2 * APK[P] * (Th - APT0[P]);
  End;
  { F_i = dE/dth / sin(th) * (v^ - cos u^) / |u| }
  Fi[1] := Dedth / St * (V[1] / Lv - Ct * U[1] / Lu) / Lu;
  Fi[2] := Dedth / St * (V[2] / Lv - Ct * U[2] / Lu) / Lu;
  Fi[3] := Dedth / St * (V[3] / Lv - Ct * U[3] / Lu) / Lu;
  Fk[1] := Dedth / St * (U[1] / Lu - Ct * V[1] / Lv) / Lv;
  Fk[2] := Dedth / St * (U[2] / Lu - Ct * V[2] / Lv) / Lv;
  Fk[3] := Dedth / St * (U[3] / Lu - Ct * V[3] / Lv) / Lv;
  AddF(BendA[I], Fi[1], Fi[2], Fi[3]);
  AddF(BendC[I], Fk[1], Fk[2], Fk[3]);
  AddF(BendB[I], -Fi[1] - Fk[1], -Fi[2] - Fk[2], -Fi[3] - Fk[3]);
  AddVir(U, Fi[1], Fi[2], Fi[3]);
  AddVir(V, Fk[1], Fk[2], Fk[3]);
 End;
 { torsions }
 For I := 1 To NTors Do Begin
  P := TorP[I];
  B1[1] := SPos[TorB[I],1] - SPos[TorA[I],1]; B1[2] := SPos[TorB[I],2] -
    SPos[TorA[I],2];
  B1[3] := SPos[TorB[I],3] - SPos[TorA[I],3];
  MinImg(B1[1], B1[2], B1[3]);
  B2[1] := SPos[TorC[I],1] - SPos[TorB[I],1]; B2[2] := SPos[TorC[I],2] -
    SPos[TorB[I],2];
  B2[3] := SPos[TorC[I],3] - SPos[TorB[I],3];
  MinImg(B2[1], B2[2], B2[3]);
  B3[1] := SPos[TorD[I],1] - SPos[TorC[I],1]; B3[2] := SPos[TorD[I],2] -
    SPos[TorC[I],2];
  B3[3] := SPos[TorD[I],3] - SPos[TorC[I],3];
  MinImg(B3[1], B3[2], B3[3]);
  { F = r_i - r_j, G = r_j - r_k, H = r_l - r_k }
  Fv[1] := -B1[1]; Fv[2] := -B1[2]; Fv[3] := -B1[3];
  G[1] := -B2[1]; G[2] := -B2[2]; G[3] := -B2[3];
  H[1] := B3[1]; H[2] := B3[2]; H[3] := B3[3];
  Cross3(Fv, G, A);
  Cross3(H, G, B);
  A2 := Dot3(A, A); B2n := Dot3(B, B);
  Ga := Sqrt(Dot3(G, G));
  If (A2 < 1E-16) Or (B2n < 1E-16) Or (Ga < 1E-8) Then Continue;
  CosP := Dot3(A, B) / Sqrt(A2 * B2n);
  Cross3(B, A, Fa);
  SinP := Dot3(Fa, G) / (Sqrt(A2 * B2n) * Ga);
  Phi := ATan2(SinP, CosP);
  MMETors := MMETors + TPK[P] * (1 + Cos(TPM[P] * Phi - TPPhi[P]));
  Dedp := -TPK[P] * TPM[P] * Sin(TPM[P] * Phi - TPPhi[P]);
  FG := Dot3(Fv, G); HG := Dot3(H, G);
  { derivees de l'angle dihedre (Blondel et Karplus) }
  For Q := 1 To 3 Do Begin
   Fa[Q] := -Ga / A2 * A[Q];
   Fd[Q] := Ga / B2n * B[Q];
   Fb[Q] := Ga / A2 * A[Q] + FG / (A2 * Ga) * A[Q] - HG / (B2n * Ga) * B[Q];
   Fc[Q] := -Ga / B2n * B[Q] - FG / (A2 * Ga) * A[Q] + HG / (B2n * Ga) * B[Q];
  End;
  AddF(TorA[I], -Dedp * Fa[1], -Dedp * Fa[2], -Dedp * Fa[3]);
  AddF(TorB[I], -Dedp * Fb[1], -Dedp * Fb[2], -Dedp * Fb[3]);
  AddF(TorC[I], -Dedp * Fc[1], -Dedp * Fc[2], -Dedp * Fc[3]);
  AddF(TorD[I], -Dedp * Fd[1], -Dedp * Fd[2], -Dedp * Fd[3]);
  { virial : positions relatives a l'atome A }
  Rel := Zero;
  AddVir(Rel, -Dedp * Fa[1], -Dedp * Fa[2], -Dedp * Fa[3]);
  Rel := B1;
  AddVir(Rel, -Dedp * Fb[1], -Dedp * Fb[2], -Dedp * Fb[3]);
  Rel[1] := B1[1] + B2[1]; Rel[2] := B1[2] + B2[2]; Rel[3] := B1[3] + B2[3];
  AddVir(Rel, -Dedp * Fc[1], -Dedp * Fc[2], -Dedp * Fc[3]);
  Rel[1] := B1[1] + B2[1] + B3[1]; Rel[2] := B1[2] + B2[2] + B3[2];
  Rel[3] := B1[3] + B2[3] + B3[3];
  AddVir(Rel, -Dedp * Fd[1], -Dedp * Fd[2], -Dedp * Fd[3]);
 End;
End;

{ Energie totale FIST ; forces dans SFrc et deformation dans DEps }
Function MMEval(NeedF: Boolean): Real;
Var I, K: Integer;
Begin
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do SFrc[I,K] := 0;
 FillChar(DEps, SizeOf(DEps), 0);
 MMEBond := 0; MMEBend := 0; MMETors := 0; MMEVdw := 0; MMECoul := 0;
 MMEEwReal := 0; MMEEwRecip := 0; MMEEwSelf := 0;
 If Not ErfcReady Then InitErfc;
 PairTerms(NeedF);
 If EwaldOn Then EwaldRecip(NeedF);
 BondedTerms(NeedF);
 MMEval := MMEBond + MMEBend + MMETors + MMEVdw + MMECoul + MMEEwReal +
   MMEEwRecip + MMEEwSelf;
End;


{ --------------------------------------------------------------------- }
{ Interface avec le moteur quantique (Quickstep)                        }
{ --------------------------------------------------------------------- }

Var
 QSPotIgnored : Boolean;
 QSFuncName : String[24];
 QSMethName : String[24];
 QSNumScf   : LongInt;
 QSEps      : Real;
 QSMaxScf   : Integer;

{ Fonctionnelle d'echange-correlation selon la section &XC }
Procedure QSReadXC;
Var
 Sx, Sf, Sh: Integer;
 Short, T: String;
 Frac: Real;
 HasPBE, HasB88, HasLYP, HasPade, HasHF: Boolean;
Begin
 Meth := M_DFT; Func := F_PBE; HFExch := 0.0;
 QSFuncName := 'PBE';
 Sx := FindSec(SecDFT, 'XC', 1);
 If Sx = 0 Then Begin
  Warn('Aucune section &XC : fonctionnelle PBE utilisee');
  Exit;
 End;
 Sf := FindSec(Sx, 'XC_FUNCTIONAL', 1);
 Sh := FindSec(Sx, 'HF', 1);
 Frac := 1.0;
 HasHF := Sh > 0;
 If HasHF Then Frac := KwReal(Sh, 'FRACTION', 1.0, '');
 Short := '';
 If Sf > 0 Then Begin
  T := SecPar[Sf];
  Short := UpStr(GetTok(T));
 End;
 HasPBE := (Sf > 0) And (FindSec(Sf, 'PBE', 1) > 0);
 HasB88 := (Sf > 0) And (FindSec(Sf, 'BECKE88', 1) > 0);
 HasLYP := (Sf > 0) And (FindSec(Sf, 'LYP', 1) > 0);
 HasPade := (Sf > 0) And ((FindSec(Sf, 'PADE', 1) > 0) Or (FindSec(Sf, 'LDA_X',
   1) > 0) Or
                         (FindSec(Sf, 'SLATER', 1) > 0));
 If (Short = 'NONE') Or ((Sf = 0) And HasHF) Then Begin
  If HasHF Then Begin
   Meth := M_HF; Func := F_NONE; HFExch := 1.0; QSFuncName := 'HARTREE-FOCK';
   If Abs(Frac - 1.0) > 1E-9 Then
     Warn('FRACTION de &HF ignoree sans fonctionnelle : Hartree-Fock complet');
  End
  Else Begin
   Warn('Aucune fonctionnelle d''echange-correlation : PBE utilisee');
  End;
  Exit;
 End;
 If (Short = 'PBE') Or ((Short = '') And HasPBE) Or ((Short = 'NO_SHORTCUT') And
   HasPBE) Then Begin
  Func := F_PBE; HFExch := 0.0; QSFuncName := 'PBE';
  If HasHF Then Begin
   If Abs(Frac - 0.25) < 0.01 Then Begin Func := F_PBE0; HFExch := 0.25;
     QSFuncName := 'PBE0'; End
   Else Warn('Fraction d''echange HF ' + RStr(Frac, 3) +
     ' non prise en charge avec PBE : ignoree');
  End;
 End
 Else If Short = 'PBE0' Then Begin Func := F_PBE0; HFExch := 0.25; QSFuncName :=
   'PBE0'; End
 Else If Short = 'BLYP' Then Begin Func := F_BLYP; HFExch := 0.0; QSFuncName :=
   'BLYP'; End
 Else If Short = 'B3LYP' Then Begin Func := F_B3LYP; HFExch := 0.2; QSFuncName
   := 'B3LYP'; End
 Else If (Short = 'LDA') Or (Short = 'PADE') Or (Short = 'SVWN') Or HasPade Then
   Begin
  Func := F_LDA; HFExch := 0.0; QSFuncName := 'LDA (PADE)';
 End
 Else If HasB88 And HasLYP Then Begin
  Func := F_BLYP; HFExch := 0.0; QSFuncName := 'BLYP';
  If HasHF Then Begin
   If Abs(Frac - 0.2) < 0.01 Then Begin Func := F_B3LYP; HFExch := 0.2;
     QSFuncName := 'B3LYP'; End
   Else Warn('Fraction d''echange HF ' + RStr(Frac, 3) +
     ' non prise en charge avec BLYP : ignoree');
  End;
 End
 Else Begin
  Warn('Fonctionnelle ' + Short + ' non prise en charge : PBE utilisee');
  Func := F_PBE; HFExch := 0.0; QSFuncName := 'PBE';
 End;
End;

{ Parametres du calcul quantique selon la section &DFT }
Procedure QSReadParams;
Var Sq, Sc, Sm, I, Nel: Integer;
Begin
 DefaultParams;
 QSReadXC;
 Sq := FindSec(SecDFT, 'QS', 1);
 QSMethName := 'GPW';
 If Sq > 0 Then QSMethName := KwWord(Sq, 'METHOD', 'GPW');
 If (QSMethName <> 'GPW') And (QSMethName <> 'GAPW') And (QSMethName <>
   'GAPW_XC') Then
  InErr(SecLine[Sq], 'METHOD ' + QSMethName +
    ' non pris en charge (seules GPW et GAPW sont emulees)');
 Charge := KwInt(SecDFT, 'CHARGE', 0);
 Nel := 0;
 For I := 1 To NAtoms Do Nel := Nel + SZ[I];
 Nel := Nel - Charge;
 If KwHas(SecDFT, 'MULTIPLICITY') And (KwInt(SecDFT, 'MULTIPLICITY', 0) > 0)
   Then
  Mult := KwInt(SecDFT, 'MULTIPLICITY', 1)
 Else If Odd(Nel) Then Mult := 2 Else Mult := 1;
 UnrestrictedFlag := KwBool(SecDFT, 'UKS', False) Or KwBool(SecDFT, 'LSD',
   False) Or (Mult > 1);
 QSEps := 1E-5;
 QSMaxScf := 50;
 Sc := FindSec(SecDFT, 'SCF', 1);
 If Sc > 0 Then Begin
  QSEps := KwReal(Sc, 'EPS_SCF', 1E-5, '');
  QSMaxScf := KwInt(Sc, 'MAX_SCF', 50);
 End;
 ScfTolErr := QSEps;
 ScfTolE := QSEps * 1E-3;
 ScfTolD := QSEps * 0.1;
 If QSMaxScf < 100 Then ScfMaxIter := 100 Else ScfMaxIter := QSMaxScf;
 GridScale := 1;
 Sm := FindSec(SecDFT, 'MGRID', 1);
 If (Sm > 0) And (KwReal(Sm, 'CUTOFF', 280.0, '') >= 500.0) Then GridScale := 2;
End;

{ Verifie le systeme et prepare les bases ; Faux en cas d'erreur }
Function QSSetupSystem(Var Err: String): Boolean;
Var I: Integer;
Begin
 QSSetupSystem := False;
 Err := '';
 If NAtoms > MAX_AT Then Begin
  Err := 'Le moteur Quickstep est limite a ' + IStr(MAX_AT) + ' atomes (' +
    IStr(NAtoms) + ' fournis)';
  Exit;
 End;
 NAt := NAtoms;
 If Per[1] Or Per[2] Or Per[3] Then
  Warn('Le moteur Quickstep ne gere pas la periodicite : systeme traite ' +
    'comme ' +
    'isole (utilisez PERIODIC NONE)');
 QSPotIgnored := False;
 For I := 1 To NKinds Do If KPotName[I] <> '' Then QSPotIgnored := True;
 For I := 1 To NAtoms Do Begin
  If (SZ[I] < 1) Or (SZ[I] > 18) Then Begin
   Err := 'L''element ' + PtSym[SZ[I]] +
     ' n''est pas pris en charge par le moteur Quickstep (H a Ar)';
   Exit;
  End;
  AtZ[I] := SZ[I];
  If KBas[SKind[I]] >= 0 Then AtBas[I] := KBas[SKind[I]]
  Else Begin
   AtBas[I] := B_STO3G;
   If KBasName[SKind[I]] = '' Then Begin
    Warn('BASIS_SET manquant pour la sorte ' + KName[SKind[I]] +
      ' : STO-3G utilisee');
    KBasName[SKind[I]] := 'STO-3G';
    KBas[SKind[I]] := B_STO3G;
   End;
  End;
  If Not BasisSupports(AtBas[I], AtZ[I]) Then Begin
   Err := 'La base ' + BasisNames[AtBas[I]] +
     ' n''est pas definie pour l''element ' + PtSym[SZ[I]];
   Exit;
  End;
 End;
 If Not BuildBasis(Err) Then Exit;
 If Not SetupElectrons(Err) Then Exit;
 AllocMats;
 ERIThresh := 1E-11;
 OrbsValid := False;
 UserAbort := False;
 QSSetupSystem := True;
End;

Procedure QSTight(On: Boolean);
Begin
 If On Then Begin
  ScfTolE := 1E-10; ScfTolD := 1E-8; ScfTolErr := 1E-8;
 End
 Else Begin
  ScfTolErr := QSEps;
  ScfTolE := QSEps * 1E-3;
  ScfTolD := QSEps * 0.1;
 End;
End;

Var
 QSErr : String;

{ Energie quantique pour la geometrie courante (SPos) }
Function QSCalc(Verbose, UseGuess: Boolean; Var E: Real): Boolean;
Var I, K: Integer;
Begin
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do AtX[I, K] := SPos[I, K];
 QSCalc := EvalEnergy(Verbose, UseGuess, E, QSErr);
 Inc(QSNumScf);
End;

{ Forces par differences finies centrees ; l'energie au point central est
  calculee avec Verbose (tableau SCF) puis les orbitales sont restaurees }
Function QSEnergyForces(Verbose: Boolean; Var E: Real): Boolean;
Var
 A, K: Integer; X0, Ep, Em, H: Real;
Begin
 QSEnergyForces := False;
 QSTight(True);
 If Not QSCalc(Verbose, False, E) Then Begin QSTight(False); Exit; End;
 H := 0.005;
 For A := 1 To NAtoms Do
  For K := 1 To 3 Do Begin
   X0 := SPos[A, K];
   SPos[A, K] := X0 + H;
   If Not QSCalc(False, True, Ep) Then Begin SPos[A, K] := X0; QSTight(False);
     Exit; End;
   SPos[A, K] := X0 - H;
   If Not QSCalc(False, True, Em) Then Begin SPos[A, K] := X0; QSTight(False);
     Exit; End;
   SPos[A, K] := X0;
   SFrc[A, K] := -(Ep - Em) / (2 * H);
  End;
 { retour au point central pour retrouver l'etat du SCF }
 If Not QSCalc(False, True, E) Then Begin QSTight(False); Exit; End;
 QSTight(False);
 FillChar(DEps, SizeOf(DEps), 0);
 QSEnergyForces := True;
End;

{ Resultats du calcul quantique au format Quickstep }
Procedure QSPrintResults;
Var I: Integer; Pop, Gap: Real;
Begin
 PopAnalysis;
 ComputeDipole;
 EmitBlank;
 Emit('  Core Hamiltonian energy:' + FR(EOne, 49, 10));
 Emit('  Hartree energy:' + FR(ECoul, 58, 10));
 If Abs(HFExch) > 1E-12 Then Emit('  Hartree-Fock exchange energy:' + FR(EExch,
   44, 10));
 If Meth = M_DFT Then Emit('  Exchange-correlation energy:' + FR(EXCdft, 45,
   10));
 Emit('  Nuclear repulsion energy:' + FR(Enuc, 48, 10));
 EmitBlank;
 Emit('  Total energy:' + FR(EScf, 60, 10));
 EmitBlank;
 Emit(' Mulliken Population Analysis');
 EmitBlank;
 Emit(' #  Atom  Element  Kind  Atomic population                           ' +
   'Net charge');
 For I := 1 To NAtoms Do Begin
  Pop := SZ[I] - QMull[I];
  Emit(' ' + ISt(I, 1) + PadL(IStr(I), 6 - Length(IStr(I)) + Length(IStr(I))) +
    '  ' +
       PadR(PtSym[SZ[I]], 8) + PadR(IStr(SKind[I]), 4) + FR(Pop, 20, 8) +
         FR(QMull[I], 32, 8));
 End;
 Emit(' # Total charge' + FR(Charge, 61, 8));
 EmitBlank;
 Emit('  Dipole moment [Debye]');
 Emit('     X=' + FR(Dip[1] * BOHR_DEBYE, 15, 8) + ' Y=' + FR(Dip[2] *
   BOHR_DEBYE, 15, 8) +
      ' Z=' + FR(Dip[3] * BOHR_DEBYE, 15, 8) + '     Total=' + FR(DipTot *
        BOHR_DEBYE, 15, 8));
 EmitBlank;
 If PrintLevel >= PL_HIGH Then Begin
  Emit(' MO EIGENVALUES AND MO OCCUPATION NUMBERS');
  EmitBlank;
  If NOccA <> NOccB Then Emit('  Spin alpha');
  Emit('     Index   Eigenvalue [a.u.]   Eigenvalue [eV]   Occupation');
  For I := 1 To NBF Do Begin
   If I <= NOccA Then Pop := 1 Else Pop := 0;
   If NOccA = NOccB Then Pop := 2 * Pop;
   Emit(ISt(I, 10) + FR(EpsA[I], 20, 8) + FR(EpsA[I] * HARTREE_EV, 18, 6) +
     FR(Pop, 14, 6));
  End;
  If NOccA <> NOccB Then Begin
   Emit('  Spin beta');
   For I := 1 To NBF Do Begin
    If I <= NOccB Then Pop := 1 Else Pop := 0;
    Emit(ISt(I, 10) + FR(EpsB[I], 20, 8) + FR(EpsB[I] * HARTREE_EV, 18, 6) +
      FR(Pop, 14, 6));
   End;
  End;
  EmitBlank;
 End;
 If NOccA > 0 Then Begin
  If NOccA < NBF Then Begin
   Gap := (EpsA[NOccA + 1] - EpsA[NOccA]) * HARTREE_EV;
   Emit('  HOMO - LUMO gap [eV]:' + FR(Gap, 52, 6));
  End;
  Emit('  HOMO energy [eV]:' + FR(EpsA[NOccA] * HARTREE_EV, 56, 6));
 End;
End;


{ --------------------------------------------------------------------- }
{ Evaluation de l'energie et des forces                                 }
{ --------------------------------------------------------------------- }

Var
 ECur       : Real;
 NEval      : LongInt;
 TEval      : LongInt;       { temps cumule des evaluations, en ms }
 FEName     : String[8];

Function EvalAll(NeedF, Verb: Boolean): Boolean;
Var T1: LongInt; Ok: Boolean;
Begin
 T1 := NowMs;
 Inc(NEval);
 Ok := True;
 Case FEMethod Of
  FE_FIST: ECur := MMEval(NeedF);
  FE_QS: Begin
   If NeedF Then Ok := QSEnergyForces(Verb, ECur)
   Else Ok := QSCalc(Verb, OrbsValid, ECur);
   If Not Ok Then Abort(QSErr);
  End;
 End;
 If NeedF Then ApplyFixed;
 TEval := TEval + ElapsedMs(T1);
 EvalAll := Ok;
End;

{ Composantes de l'energie FIST }
Procedure PrintFistEnergies;
Begin
 EmitBlank;
 Emit(' FIST ENERGIES [hartree]');
 EmitBlank;
 Info('FIST', 'Bond energy:', FR(MMEBond, 20, 10));
 Info('FIST', 'Bend energy:', FR(MMEBend, 20, 10));
 Info('FIST', 'Torsion energy:', FR(MMETors, 20, 10));
 Info('FIST', 'Van der Waals energy:', FR(MMEVdw, 20, 10));
 If EwaldOn Then Begin
  Info('FIST', 'Ewald real space energy:', FR(MMEEwReal, 20, 10));
  Info('FIST', 'Ewald reciprocal energy:', FR(MMEEwRecip, 20, 10));
  Info('FIST', 'Ewald self and background energy:', FR(MMEEwSelf, 20, 10));
 End
 Else Info('FIST', 'Coulomb energy:', FR(MMECoul, 20, 10));
 EmitBlank;
End;

Procedure PrintEnergyLine;
Begin
 Info('ENERGY', 'Total FORCE_EVAL ( ' + FEName + ' ) energy [a.u.]:', FR(ECur,
   20, 10));
End;

Procedure PrintForces(Const Title: String);
Var I: Integer; SX, SY, SZ2: Real;
Begin
 EmitBlank;
 Emit(' ' + Title);
 EmitBlank;
 Emit(' # Atom   Kind   Element          X              Y              Z');
 SX := 0; SY := 0; SZ2 := 0;
 For I := 1 To NAtoms Do Begin
  Emit(ISt(I, 6) + ISt(SKind[I], 7) + PadL(PtSym[SZ[I]], 8) + FR(SFrc[I,1], 19,
    8) +
       FR(SFrc[I,2], 15, 8) + FR(SFrc[I,3], 15, 8));
  SX := SX + SFrc[I,1]; SY := SY + SFrc[I,2]; SZ2 := SZ2 + SFrc[I,3];
 End;
 Emit(' SUM OF ATOMIC FORCES' + FR(SX, 17, 8) + FR(SY, 15, 8) + FR(SZ2, 15, 8) +
      FR(Sqrt(SX * SX + SY * SY + SZ2 * SZ2), 15, 8));
 EmitBlank;
End;

{ Tenseur des contraintes : dE/d(eps) en [GPa] (P = -dE/dV) }
Procedure PrintStress;
Var A, B: Integer; Pr: TMat3; Tr: Real;
Begin
 For A := 1 To 3 Do
  For B := 1 To 3 Do Pr[A,B] := -0.5 * (DEps[A,B] + DEps[B,A]) / CellVol /
    GPA_AU;
 EmitBlank;
 Emit(' STRESS| Analytical stress tensor [GPa]');
 Emit(' STRESS|                        x                   y                 ' +
   ' ' +
   ' z');
 Tr := 0;
 For A := 1 To 3 Do Begin
  Emit(' STRESS|      ' + Chr(Ord('x') + A - 1) + FR(Pr[A,1], 20, 9) +
    FR(Pr[A,2], 20, 9) + FR(Pr[A,3], 20, 9));
  Tr := Tr + Pr[A,A];
 End;
 Emit(' STRESS| 1/3 Trace' + FR(Tr / 3, 31, 9));
 EmitBlank;
End;

{ --------------------------------------------------------------------- }
{ Utilitaires dynamique                                                 }
{ --------------------------------------------------------------------- }

Var
 MDDof    : Integer;
 GaussSpare : Real;
 GaussHave : Boolean;

Function Gauss: Real;
Var U1, U2, R: Real;
Begin
 If GaussHave Then Begin GaussHave := False; Gauss := GaussSpare; Exit; End;
 Repeat
  U1 := 2 * Random - 1;
  U2 := 2 * Random - 1;
  R := U1 * U1 + U2 * U2;
 Until (R < 1) And (R > 1E-12);
 R := Sqrt(-2 * Ln(R) / R);
 GaussSpare := U2 * R;
 GaussHave := True;
 Gauss := U1 * R;
End;

Function KinEnergy: Real;
Var I: Integer; K: Real;
Begin
 K := 0;
 For I := 1 To NAtoms Do
  K := K + 0.5 * SMass[I] * (Sqr(SVel[I,1]) + Sqr(SVel[I,2]) + Sqr(SVel[I,3]));
 KinEnergy := K;
End;

Function InstTemp: Real;
Begin
 If MDDof > 0 Then InstTemp := 2 * KinEnergy / (MDDof * KB_AU) Else InstTemp :=
   0;
End;

Procedure RemoveComVel;
Var I, K: Integer; P: TVec3; M: Real;
Begin
 M := 0;
 For K := 1 To 3 Do P[K] := 0;
 For I := 1 To NAtoms Do Begin
  M := M + SMass[I];
  For K := 1 To 3 Do P[K] := P[K] + SMass[I] * SVel[I,K];
 End;
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do SVel[I,K] := SVel[I,K] - P[K] / M;
End;

{ Vitesses de Maxwell-Boltzmann a la temperature T (en kelvin) }
Procedure InitVelocities(T: Real);
Var I, K: Integer; Sc: Real;
Begin
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do SVel[I,K] := Gauss * Sqrt(KB_AU * T / SMass[I]);
 ApplyFixed;
 If (NAtoms > 1) And (NFixComp = 0) Then RemoveComVel;
 Sc := InstTemp;
 If Sc > 1E-12 Then Begin
  Sc := Sqrt(T / Sc);
  For I := 1 To NAtoms Do
   For K := 1 To 3 Do SVel[I,K] := SVel[I,K] * Sc;
 End;
End;

{ Pression instantanee [au] ; Pi = (2K_ab + virial) / V }
Function InstPressure: Real;
Var I: Integer; K, W: Real;
Begin
 K := 0;
 For I := 1 To NAtoms Do
  K := K + SMass[I] * (Sqr(SVel[I,1]) + Sqr(SVel[I,2]) + Sqr(SVel[I,3]));
 W := DEps[1,1] + DEps[2,2] + DEps[3,3];
 InstPressure := (K - W) / (3 * CellVol);
End;

{ Frequence d'ecriture : &EACH / cle sous &MOTION/&PRINT/&NOM }
Function EachOf(Const Nm, Key: String; Def: LongInt): LongInt;
Var Sp, Sx, Se: Integer;
Begin
 EachOf := Def;
 If SecMotion = 0 Then Exit;
 Sp := FindSec(SecMotion, 'PRINT', 1);
 If Sp = 0 Then Exit;
 Sx := FindSec(Sp, Nm, 1);
 If Sx = 0 Then Exit;
 Se := FindSec(Sx, 'EACH', 1);
 If Se = 0 Then Exit;
 EachOf := KwInt(Se, Key, Def);
End;

{ Vrai si la duree maximale (WALLTIME) est depassee }
Function WallExceeded: Boolean;
Begin
 WallExceeded := (WallTime > 0) And (RunSec > WallTime);
End;

{ Fichiers de sortie }
Function ProjFile(Const Suffix: String): String;
Begin
 ProjFile := ProjectName + Suffix;
End;

{ Ouvre un fichier de sortie (en ajout lors d'une reprise) ; Fresh vrai si le
  fichier a ete (re)cree }
Function OpenOut(Var F: Text; Const Nm: String; App: Boolean; Var Fresh:
  Boolean): Boolean;
Begin
 Assign(F, Nm);
 Fresh := True;
 {$I-}
 If App Then Begin
  Append(F);
  If IOResult = 0 Then Fresh := False Else Rewrite(F);
 End
 Else Rewrite(F);
 {$I+}
 OpenOut := IOResult = 0;
End;

{ Ecrit le fichier de redemarrage : l'entree avec coordonnees et vitesses }
Procedure WriteRestart(Step: LongInt; Time: Real);
Var Ok: Integer;
Begin
 {$I-}
 Assign(RedirF, ProjFile('-1.restart'));
 Rewrite(RedirF);
 {$I+}
 Ok := IOResult;
 If Ok <> 0 Then Begin Warn('Impossible d''ecrire ' + ProjFile('-1.restart'));
   Exit; End;
 RedirOn := True;
 RestartMode := True;
 RestartStep := Step;
 RestartTime := Time;
 DumpSec(0, 0);
 RestartMode := False;
 RedirOn := False;
 Close(RedirF);
End;


{ --------------------------------------------------------------------- }
{ Vecteurs de travail                                                   }
{ --------------------------------------------------------------------- }

Function NewVec(N: LongInt): PRealArr;
Var P: PRealArr;
Begin
 GetMem(P, N * SizeOf(Real));
 FillChar(P^, N * SizeOf(Real), 0);
 NewVec := P;
End;

Procedure FreeVec(P: PRealArr; N: LongInt);
Begin
 FreeMem(P, N * SizeOf(Real));
End;

Function DotV(N: LongInt; A, B: PRealArr): Real;
Var I: LongInt; S: Real;
Begin
 S := 0;
 For I := 0 To N - 1 Do S := S + A^[I] * B^[I];
 DotV := S;
End;

{ --------------------------------------------------------------------- }
{ Energie simple et energie + forces                                    }
{ --------------------------------------------------------------------- }

Procedure PrintFinalSystem;
Begin
 If FEMethod = FE_FIST Then PrintFistEnergies Else QSPrintResults;
End;

Procedure RunEnergy(WithForces: Boolean);
Begin
 If PrintLevel >= PL_LOW Then Begin
  EvalAll(WithForces, PrintLevel >= PL_MEDIUM);
  If PrintLevel >= PL_MEDIUM Then PrintFinalSystem;
  EmitBlank;
  PrintEnergyLine;
  If WithForces And (PrintLevel >= PL_MEDIUM) Then Begin
   PrintForces('ATOMIC FORCES in [a.u.]');
   If (FEMethod = FE_FIST) And HasCell And (Per[1] Or Per[2] Or Per[3]) Then
     PrintStress;
  End;
 End
 Else EvalAll(WithForces, False);
End;

{ --------------------------------------------------------------------- }
{ Optimisation de geometrie et de cellule (L-BFGS)                       }
{ --------------------------------------------------------------------- }

Const
 LBM = 24;

Var
 HS, HY : Array[1..LBM] Of PRealArr;
 HRho   : Array[1..LBM] Of Real;
 HN     : Integer;

Procedure LbfgsDir(N: LongInt; G, D: PRealArr; H0: Real);
Var
 Q: PRealArr; Al: Array[1..LBM] Of Real; I: Integer; J: LongInt; Be, Gam: Real;
Begin
 Q := NewVec(N);
 For J := 0 To N - 1 Do Q^[J] := G^[J];
 For I := HN DownTo 1 Do Begin
  Al[I] := HRho[I] * DotV(N, HS[I], Q);
  For J := 0 To N - 1 Do Q^[J] := Q^[J] - Al[I] * HY[I]^[J];
 End;
 If HN > 0 Then Gam := DotV(N, HS[HN], HY[HN]) / DotV(N, HY[HN], HY[HN])
 Else Gam := H0;
 For J := 0 To N - 1 Do Q^[J] := Q^[J] * Gam;
 For I := 1 To HN Do Begin
  Be := HRho[I] * DotV(N, HY[I], Q);
  For J := 0 To N - 1 Do Q^[J] := Q^[J] + HS[I]^[J] * (Al[I] - Be);
 End;
 For J := 0 To N - 1 Do D^[J] := -Q^[J];
 FreeVec(Q, N);
End;

Procedure LbfgsPush(N: LongInt; S, Y: PRealArr);
Var SY: Real; I: Integer; J: LongInt; P1, P2: PRealArr;
Begin
 SY := DotV(N, S, Y);
 If SY <= 1E-12 * Sqrt(DotV(N, S, S) * DotV(N, Y, Y)) Then Exit;
 If HN = LBM Then Begin
  P1 := HS[1]; P2 := HY[1];
  For I := 1 To LBM - 1 Do Begin HS[I] := HS[I + 1]; HY[I] := HY[I + 1]; HRho[I]
    := HRho[I + 1]; End;
  HS[LBM] := P1; HY[LBM] := P2;
  HN := LBM - 1;
 End;
 Inc(HN);
 For J := 0 To N - 1 Do Begin HS[HN]^[J] := S^[J]; HY[HN]^[J] := Y^[J]; End;
 HRho[HN] := 1 / SY;
End;

Procedure OptGradient(G: PRealArr; CellToo: Boolean; Ls, Pext: Real);
Var I, K, Bs: Integer;
Begin
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do G^[(I - 1) * 3 + K - 1] := -SFrc[I, K];
 If CellToo Then Begin
  Bs := 3 * NAtoms;
  G^[Bs] := (DEps[1,1] + Pext * CellVol) / Ls;
  G^[Bs + 1] := (DEps[2,2] + Pext * CellVol) / Ls;
  G^[Bs + 2] := (DEps[3,3] + Pext * CellVol) / Ls;
  G^[Bs + 3] := (DEps[1,2] + DEps[2,1]) / Ls;
  G^[Bs + 4] := (DEps[1,3] + DEps[3,1]) / Ls;
  G^[Bs + 5] := (DEps[2,3] + DEps[3,2]) / Ls;
 End;
End;

{ Deplace le systeme de Alpha * D (positions et deformation de la cellule) }
Procedure OptApply(D: PRealArr; Alpha: Real; CellToo: Boolean; Ls: Real);
Var I, K, L, Bs: Integer; Eps: TMat3; R: TVec3; Hn: TMat3;
Begin
 If CellToo Then Begin
  Bs := 3 * NAtoms;
  FillChar(Eps, SizeOf(Eps), 0);
  Eps[1,1] := Alpha * D^[Bs] / Ls;
  Eps[2,2] := Alpha * D^[Bs + 1] / Ls;
  Eps[3,3] := Alpha * D^[Bs + 2] / Ls;
  Eps[1,2] := Alpha * D^[Bs + 3] / Ls; Eps[2,1] := Eps[1,2];
  Eps[1,3] := Alpha * D^[Bs + 4] / Ls; Eps[3,1] := Eps[1,3];
  Eps[2,3] := Alpha * D^[Bs + 5] / Ls; Eps[3,2] := Eps[2,3];
  For I := 1 To NAtoms Do Begin
   For K := 1 To 3 Do Begin
    R[K] := SPos[I, K];
    For L := 1 To 3 Do R[K] := R[K] + SPos[I, L] * Eps[L, K];
   End;
   For K := 1 To 3 Do SPos[I, K] := R[K];
  End;
  For I := 1 To 3 Do
   For K := 1 To 3 Do Begin
    Hn[I, K] := CellH[I, K];
    For L := 1 To 3 Do Hn[I, K] := Hn[I, K] + CellH[I, L] * Eps[L, K];
   End;
  CellH := Hn;
  UpdateCell;
 End;
 For I := 1 To NAtoms Do
  For K := 1 To 3 Do SPos[I, K] := SPos[I, K] + Alpha * D^[(I - 1) * 3 + K - 1];
End;

Procedure PrintOptBlock(Step: Integer; Const MethodName: String; E, DE, Pred,
  Scale, StepSize,
                        Trust: Real; Decr: Boolean; UsedT, MaxStep, RmsStep,
                          MaxG, RmsG,
                        LimDr, LimRms, LimF, LimRmsF: Real; Cs, Cr, Cg, Crg:
                          Boolean);
 Procedure L(Const Lbl: String; V: Real; Dec: Integer);
 Begin
  Emit('  ' + PadR(Lbl, 27) + '= ' + FR(V, 22, Dec));
 End;
 Procedure LB(Const Lbl: String; B: Boolean);
 Begin
  Emit('  ' + PadR(Lbl, 27) + '= ' + PadL(YesNo(B), 22));
 End;
Begin
 EmitBlank;
 Emit(' --------  Informations at step =' + ISt(Step, 6) + ' ------------');
 Emit('  ' + PadR('Optimization Method', 27) + '= ' + PadL(MethodName, 22));
 L('Total Energy', E, 10);
 L('Real energy change', DE, 10);
 L('Predicted change in energy', Pred, 10);
 L('Scaling factor', Scale, 7);
 L('Step size', StepSize, 7);
 L('Trust radius', Trust, 7);
 LB('Decrease in energy', Decr);
 L('Used time', UsedT, 3);
 EmitBlank;
 Emit('  Convergence check :');
 L('Max. step size', MaxStep, 7);
 L('Conv. limit for step size', LimDr, 7);
 LB('Convergence in step size', Cs);
 L('RMS step size', RmsStep, 7);
 L('Conv. limit for RMS step', LimRms, 7);
 LB('Convergence in RMS step', Cr);
 L('Max. gradient', MaxG, 7);
 L('Conv. limit for gradients', LimF, 7);
 LB('Conv. for gradients', Cg);
 L('RMS gradient', RmsG, 7);
 L('Conv. limit for RMS grad.', LimRmsF, 7);
 LB('Conv. in RMS gradients', Crg);
 Emit(' ' + Rep('-', 51));
End;

Procedure RunOpt(CellToo: Boolean);
Var
 Smo, Sgo, Sb: Integer;
 MaxIter, It, Tries, I, K, N3, Nv, Step0, EachGeo: LongInt;
 MaxDr, RmsDr, MaxF, RmsF, TrustR, PressTol, Pext, Ls: Real;
 Opt, Title: String;
 G, Gn, D, Sv, Yv: PRealArr;
 E, En, Alpha, Dmax, StepMax, StepRms, GMax, GRms, Pred, T1, PMaxDev, Rr: Real;
 Cs, Cr, Cg, Crg, Cp, Conv, Decr: Boolean;
 SavePos: TXYZ; SaveCell: TMat3;
 Fx: Text; Ok: Integer; A, B: Integer;
Begin
 Smo := SecMotion;
 If CellToo Then Sgo := FindSec(Smo, 'CELL_OPT', 1) Else Sgo := FindSec(Smo,
   'GEO_OPT', 1);
 If Smo = 0 Then Sgo := 0;
 Opt := 'BFGS';
 MaxIter := 200; MaxDr := 0.003; RmsDr := 0.0015; MaxF := 0.00045; RmsF :=
   0.0003;
 TrustR := 0.3; PressTol := 100.0 * BAR_AU; Pext := 100.0 * BAR_AU;
 Step0 := 0;
 If Sgo > 0 Then Begin
  Opt := KwWord(Sgo, 'OPTIMIZER', 'BFGS');
  MaxIter := KwInt(Sgo, 'MAX_ITER', 200);
  MaxDr := KwReal(Sgo, 'MAX_DR', 0.003, 'bohr');
  RmsDr := KwReal(Sgo, 'RMS_DR', 0.0015, 'bohr');
  MaxF := KwReal(Sgo, 'MAX_FORCE', 0.00045, 'hartree*bohr^-1');
  RmsF := KwReal(Sgo, 'RMS_FORCE', 0.0003, 'hartree*bohr^-1');
  Step0 := KwInt(Sgo, 'STEP_START_VAL', 0);
  If CellToo Then Begin
   PressTol := KwReal(Sgo, 'PRESSURE_TOLERANCE', 100.0, 'bar');
   Pext := KwReal(Sgo, 'EXTERNAL_PRESSURE', 100.0, 'bar');
  End;
  Sb := FindSec(Sgo, 'BFGS', 1);
  If Sb > 0 Then TrustR := KwReal(Sb, 'TRUST_RADIUS', 0.3, 'bohr');
 End;
 If (Opt <> 'BFGS') And (Opt <> 'CG') And (Opt <> 'LBFGS') Then
  Abort('OPTIMIZER ' + Opt + ' non pris en charge (BFGS, CG, LBFGS)');
 If (Opt <> 'BFGS') Then Warn('OPTIMIZER ' + Opt +
   ' remplace par BFGS a memoire limitee');
 If CellToo Then Begin
  If FEMethod <> FE_FIST Then
    Abort('CELL_OPT n''est disponible qu''avec METHOD FIST (tenseur des ' +
      'contraintes)');
  If Not HasCell Then Abort('CELL_OPT exige une section &CELL');
  Title := 'CELL';
 End
 Else Title := 'GEOMETRY';
 N3 := 3 * NAtoms;
 Nv := N3;
 If CellToo Then Nv := N3 + 6;
 Ls := Exp(Ln(CellVol + 1E-30) / 3);
 If Not CellToo Then Ls := 1;
 G := NewVec(Nv); Gn := NewVec(Nv); D := NewVec(Nv); Sv := NewVec(Nv); Yv :=
   NewVec(Nv);
 For I := 1 To LBM Do Begin HS[I] := NewVec(Nv); HY[I] := NewVec(Nv); End;
 HN := 0;
 EmitBlank;
 BoxBanner('STARTING ' + Title + ' OPTIMIZATION', Opt);
 EvalAll(True, PrintLevel >= PL_HIGH);
 E := ECur;
 If CellToo Then E := E + Pext * CellVol;
 OptGradient(G, CellToo, Ls, Pext);
 Assign(Fx, ProjFile('-pos-1.xyz'));
 {$I-}
 Rewrite(Fx);
 {$I+}
 Ok := IOResult;
 If Ok = 0 Then Begin
  WriteXYZFrame(Fx, ' i =' + ISt(Step0, 9) + ', E =' + FR(E, 20, 10));
 End;
 EachGeo := EachOf('TRAJECTORY', 'GEO_OPT', 1);
 If CellToo Then EachGeo := EachOf('TRAJECTORY', 'CELL_OPT', EachGeo);
 If EachGeo < 1 Then EachGeo := 1;
 Conv := False;
 It := 0;
 While (It < MaxIter) And Not Conv Do Begin
  If WallExceeded Then Begin
   Emit(' *** WALLTIME exceeded : optimization stopped ***');
   Warn('WALLTIME depasse : optimisation interrompue');
   Break;
  End;
  Inc(It);
  T1 := RunSec;
  LbfgsDir(Nv, G, D, 2.5);
  Dmax := 0;
  For I := 0 To Nv - 1 Do If Abs(D^[I]) > Dmax Then Dmax := Abs(D^[I]);
  If Dmax < 1E-14 Then Begin Conv := True; Break; End;
  Alpha := 1;
  If Dmax > TrustR Then Alpha := TrustR / Dmax;
  SavePos := SPos; SaveCell := CellH;
  Tries := 0;
  Decr := True;
  Pred := 0.5 * Alpha * DotV(Nv, G, D);
  Repeat
   OptApply(D, Alpha, CellToo, Ls);
   EvalAll(True, False);
   En := ECur;
   If CellToo Then En := En + Pext * CellVol;
   If En <= E + 1E-9 Then Break;
   { energie en hausse : on restaure et on reduit le pas }
   SPos := SavePos; CellH := SaveCell;
   If CellToo Then UpdateCell;
   Alpha := Alpha * 0.5;
   Inc(Tries);
   HN := 0;
   If Tries >= 8 Then Break;
  Until False;
  If Tries >= 8 Then Begin
   { aucun pas ne fait baisser l'energie : on revient au point precedent }
   EvalAll(True, False);
   Emit(' *** WARNING :: la recherche de pas n''a pas reduit l''energie : ' +
     'optimisation arretee ***');
   Warn('Optimisation arretee : aucun pas ne reduit l''energie');
   Break;
  End;
  OptGradient(Gn, CellToo, Ls, Pext);
  For I := 0 To Nv - 1 Do Begin
   Sv^[I] := Alpha * D^[I];
   Yv^[I] := Gn^[I] - G^[I];
  End;
  LbfgsPush(Nv, Sv, Yv);
  StepMax := 0; StepRms := 0;
  For I := 0 To Nv - 1 Do Begin
   If Abs(Sv^[I]) > StepMax Then StepMax := Abs(Sv^[I]);
   StepRms := StepRms + Sv^[I] * Sv^[I];
  End;
  StepRms := Sqrt(StepRms / Nv);
  GMax := 0; GRms := 0;
  For I := 0 To N3 - 1 Do Begin
   If Abs(Gn^[I]) > GMax Then GMax := Abs(Gn^[I]);
   GRms := GRms + Gn^[I] * Gn^[I];
  End;
  GRms := Sqrt(GRms / N3);
  Cs := StepMax <= MaxDr; Cr := StepRms <= RmsDr;
  Cg := GMax <= MaxF; Crg := GRms <= RmsF;
  Conv := Cs And Cr And Cg And Crg;
  If CellToo Then Begin
   PMaxDev := 0;
   For A := 1 To 3 Do
    For B := 1 To 3 Do Begin
     Rr := -0.5 * (DEps[A,B] + DEps[B,A]) / CellVol;
     If A = B Then Rr := Rr - Pext;
     If Abs(Rr) > PMaxDev Then PMaxDev := Abs(Rr);
    End;
   Cp := PMaxDev <= PressTol;
   Conv := Conv And Cp;
  End;
  If PrintLevel >= PL_MEDIUM Then Begin
   PrintOptBlock(It, Opt, ECur, En - E, Pred, 0.0, StepMax, TrustR, Decr,
                 RunSec - T1, StepMax, StepRms, GMax, GRms, MaxDr, RmsDr, MaxF,
                   RmsF,
                 Cs, Cr, Cg, Crg);
   If CellToo Then Begin
    Emit('  ' + PadR('Pressure deviation [bar]', 27) + '= ' + FR(PMaxDev /
      BAR_AU, 22, 7));
    Emit('  ' + PadR('Pressure tolerance [bar]', 27) + '= ' + FR(PressTol /
      BAR_AU, 22, 7));
    Emit('  ' + PadR('Convergence in pressure', 27) + '= ' + PadL(YesNo(Cp),
      22));
   End;
   If PrintLevel >= PL_HIGH Then PrintForces('ATOMIC FORCES in [a.u.]');
  End;
  If (Ok = 0) And (It Mod EachGeo = 0) Then
   WriteXYZFrame(Fx, ' i =' + ISt(Step0 + It, 9) + ', E =' + FR(ECur, 20, 10));
  E := En;
  For I := 0 To Nv - 1 Do G^[I] := Gn^[I];
 End;
 If Ok = 0 Then Close(Fx);
 EmitBlank;
 If Conv Then Begin
  BoxBanner(Title + ' OPTIMIZATION COMPLETED', '');
  EmitBlank;
  Emit('                    Reevaluating energy at the minimum');
 End
 Else Begin
  BoxBanner('MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED',
    'OPTIMIZATION NOT CONVERGED');
  Warn('Optimisation non convergee apres ' + IStr(It) + ' iterations');
 End;
 EvalAll(True, PrintLevel >= PL_MEDIUM);
 If PrintLevel >= PL_MEDIUM Then Begin
  If FEMethod = FE_QS Then PrintCoords('QUICKSTEP', True) Else
    PrintCoords('FIST', False);
  PrintFinalSystem;
  If CellToo Then PrintCell;
 End;
 EmitBlank;
 PrintEnergyLine;
 If PrintLevel >= PL_MEDIUM Then Begin
  PrintForces('ATOMIC FORCES in [a.u.]');
  If (FEMethod = FE_FIST) And HasCell And (Per[1] Or Per[2] Or Per[3]) Then
    PrintStress;
 End;
 For I := 1 To NAtoms Do For K := 1 To 3 Do SVel[I, K] := 0;
 WriteRestart(Step0 + It, 0);
 For I := 1 To LBM Do Begin FreeVec(HS[I], Nv); FreeVec(HY[I], Nv); End;
 FreeVec(G, Nv); FreeVec(Gn, Nv); FreeVec(D, Nv); FreeVec(Sv, Nv); FreeVec(Yv,
   Nv);
End;

{ --------------------------------------------------------------------- }
{ Dynamique moleculaire                                                 }
{ --------------------------------------------------------------------- }

Var
 NhXi, NhV, NhQ : Array[1..10] Of Real;
 NhM            : Integer;
 NhTemp         : Real;
 NhOn           : Boolean;
 CsvrE          : Real;

{ Chaine de Nose-Hoover : demi-pas de duree DtFull / 2 pour la chaine }
Procedure NhcHalf(DtFull: Real);
Var
 W2, W4, W8, K2, Sc, Aa, Kt: Real; J, I: Integer; G: Array[1..10] Of Real;
Begin
 W2 := DtFull / 2; W4 := DtFull / 4; W8 := DtFull / 8;
 Kt := KB_AU * NhTemp;
 K2 := 2 * KinEnergy;
 G[1] := (K2 - MDDof * Kt) / NhQ[1];
 If NhM > 1 Then Begin
  G[NhM] := (NhQ[NhM - 1] * NhV[NhM - 1] * NhV[NhM - 1] - Kt) / NhQ[NhM];
  NhV[NhM] := NhV[NhM] + G[NhM] * W4;
  For J := NhM - 1 DownTo 1 Do Begin
   Aa := Exp(-W8 * NhV[J + 1]);
   NhV[J] := NhV[J] * Aa;
   If J = 1 Then G[J] := (K2 - MDDof * Kt) / NhQ[1]
   Else G[J] := (NhQ[J - 1] * NhV[J - 1] * NhV[J - 1] - Kt) / NhQ[J];
   NhV[J] := NhV[J] + G[J] * W4;
   NhV[J] := NhV[J] * Aa;
  End;
 End
 Else Begin
  G[1] := (K2 - MDDof * Kt) / NhQ[1];
  NhV[1] := NhV[1] + G[1] * W4;
 End;
 Sc := Exp(-W2 * NhV[1]);
 For I := 1 To NAtoms Do Begin
  SVel[I,1] := SVel[I,1] * Sc; SVel[I,2] := SVel[I,2] * Sc; SVel[I,3] :=
    SVel[I,3] * Sc;
 End;
 K2 := K2 * Sc * Sc;
 For J := 1 To NhM Do NhXi[J] := NhXi[J] + NhV[J] * W2;
 G[1] := (K2 - MDDof * Kt) / NhQ[1];
 If NhM > 1 Then Begin
  Aa := Exp(-W8 * NhV[2]);
  NhV[1] := NhV[1] * Aa;
  NhV[1] := NhV[1] + G[1] * W4;
  NhV[1] := NhV[1] * Aa;
  For J := 2 To NhM - 1 Do Begin
   Aa := Exp(-W8 * NhV[J + 1]);
   NhV[J] := NhV[J] * Aa;
   G[J] := (NhQ[J - 1] * NhV[J - 1] * NhV[J - 1] - Kt) / NhQ[J];
   NhV[J] := NhV[J] + G[J] * W4;
   NhV[J] := NhV[J] * Aa;
  End;
  G[NhM] := (NhQ[NhM - 1] * NhV[NhM - 1] * NhV[NhM - 1] - Kt) / NhQ[NhM];
  NhV[NhM] := NhV[NhM] + G[NhM] * W4;
 End
 Else NhV[1] := NhV[1] + G[1] * W4;
End;

Function NhEnergy: Real;
Var J: Integer; E, Kt: Real;
Begin
 Kt := KB_AU * NhTemp;
 E := MDDof * Kt * NhXi[1];
 For J := 1 To NhM Do E := E + 0.5 * NhQ[J] * NhV[J] * NhV[J];
 For J := 2 To NhM Do E := E + Kt * NhXi[J];
 NhEnergy := E;
End;

{ Thermostat a echantillonnage canonique (Bussi) : DtH est la duree }
Procedure CsvrApply(DtH, Tau, T0: Real);
Var
 K, Kt, C, R1, S, Kn, Alpha: Real; J: Integer; I: Integer;
Begin
 K := KinEnergy;
 If K < 1E-30 Then Exit;
 Kt := 0.5 * MDDof * KB_AU * T0;
 C := Exp(-DtH / Tau);
 R1 := Gauss;
 S := 0;
 For J := 2 To MDDof Do S := S + Sqr(Gauss);
 Kn := K + (1 - C) * (Kt / MDDof * (R1 * R1 + S) - K) + 2 * R1 * Sqrt(K * Kt /
   MDDof * (1 - C) * C);
 If Kn < 0 Then Kn := 0;
 Alpha := Sqrt(Kn / K);
 For I := 1 To NAtoms Do Begin
  SVel[I,1] := SVel[I,1] * Alpha; SVel[I,2] := SVel[I,2] * Alpha; SVel[I,3] :=
    SVel[I,3] * Alpha;
 End;
 CsvrE := CsvrE + (K - Kn);
End;

Procedure RunMD;
Var
 Smd, Sth, Sn, Sbar, Sp: Integer;
 Ens, ThType: String;
 NSteps, Step, Step0, EachRes: LongInt;
 Dt, Temp0, TimeNow, Time0, TauT, TauP, PrefP, Cons, Cons0, KinE, PotE, Tmp,
   Pres, Mu, T1: Real;
 SumPot, SumKin, SumTmp, SumPres, SumCpu, DriftAt: Real;
 NVT, NPT, DoCsvr, HaveVel: Boolean;
 I, K, J: Integer;
 Jn, NhL: LongInt;
 FrE, FrP, FrC, FrV, FrF: Boolean;
 Fe, Fp, Fc, Fv, Ff: Text;
 OkE, OkP, OkC, OkV, OkF: Boolean;
 WantVel, WantFrc: Boolean;
 BetaP: Real;
 EachTraj, EachVelo, EachFrc, LastStep: LongInt;
Begin
 Smd := 0;
 If SecMotion > 0 Then Smd := FindSec(SecMotion, 'MD', 1);
 If Smd = 0 Then Abort('RUN_TYPE MD exige une section &MOTION / &MD');
 Ens := KwWord(Smd, 'ENSEMBLE', 'NVE');
 NSteps := KwInt(Smd, 'STEPS', 3);
 Dt := KwReal(Smd, 'TIMESTEP', 0.5, 'fs');
 Temp0 := KwReal(Smd, 'TEMPERATURE', 300.0, 'K');
 Step0 := KwInt(Smd, 'STEP_START_VAL', 0);
 Time0 := KwReal(Smd, 'TIME_START_VAL', 0.0, 'fs');
 If UseResCounters Then Begin Step0 := ResStep; Time0 := ResTime / FS_AU; End;
 NVT := (Ens = 'NVT');
 NPT := (Ens = 'NPT_I') Or (Ens = 'NPT_F');
 If (Ens <> 'NVE') And Not NVT And Not NPT Then
  Abort('ENSEMBLE ' + Ens + ' non pris en charge (NVE, NVT, NPT_I)');
 If NPT And (FEMethod <> FE_FIST) Then
   Abort('L''ensemble NPT exige METHOD FIST');
 If NPT And Not (HasCell And Per[1] And Per[2] And Per[3]) Then
   Abort('L''ensemble NPT exige une cellule periodique (PERIODIC XYZ)');
 If Dt <= 0 Then Abort('TIMESTEP doit etre positif');
 If NSteps < 0 Then Abort('STEPS doit etre positif ou nul');
 If Temp0 < 0 Then Abort('TEMPERATURE doit etre positive');
 If (NVT Or NPT) And (Temp0 <= 0) Then
   Abort('TEMPERATURE doit etre strictement positive pour un thermostat');
 If NFixComp > 0 Then MDDof := 3 * NAtoms - NFixComp
 Else If NAtoms > 1 Then MDDof := 3 * NAtoms - 3
 Else MDDof := 3;
 If MDDof < 1 Then Abort('Tous les degres de liberte sont fixes');
 ThType := 'NONE';
 TauT := 1000.0 * FS_AU;
 NhM := 3;
 If NVT Or NPT Then Begin
  ThType := 'NOSE';
  Sth := FindSec(Smd, 'THERMOSTAT', 1);
  If Sth > 0 Then Begin
   ThType := KwWord(Sth, 'TYPE', 'NOSE');
   Sn := FindSec(Sth, 'NOSE', 1);
   If Sn > 0 Then Begin
    NhL := KwInt(Sn, 'LENGTH', 3);
    If (NhL < 1) Or (NhL > 10) Then
      Abort('LENGTH du thermostat de Nose doit etre compris entre 1 et 10');
    NhM := NhL;
    TauT := KwReal(Sn, 'TIMECON', 1000.0, 'fs');
   End;
   Sn := FindSec(Sth, 'CSVR', 1);
   If (Sn > 0) And (ThType = 'CSVR') Then TauT := KwReal(Sn, 'TIMECON', 1000.0,
     'fs');
  End;
  If (ThType <> 'NOSE') And (ThType <> 'CSVR') And (ThType <> 'NONE') Then
   Abort('Thermostat ' + ThType + ' non pris en charge (NOSE, CSVR)');
 End;
 If TauT <= 0 Then Abort('TIMECON du thermostat doit etre positif');
 DoCsvr := (ThType = 'CSVR');
 NhOn := (ThType = 'NOSE');
 NhTemp := Temp0;
 TauP := 1000.0 * FS_AU;
 PrefP := 1.0 * BAR_AU;
 If NPT Then Begin
  Sbar := FindSec(Smd, 'BAROSTAT', 1);
  If Sbar > 0 Then Begin
   PrefP := KwReal(Sbar, 'PRESSURE', 1.0, 'bar');
   TauP := KwReal(Sbar, 'TIMECON', 1000.0, 'fs');
  End;
  If TauP <= 0 Then Abort('TIMECON du barostat doit etre positif');
 End;
 BetaP := 4.5E-5 / BAR_AU;
 Sp := 0;
 If SecMotion > 0 Then Sp := FindSec(SecMotion, 'PRINT', 1);
 WantVel := (Sp > 0) And (FindSec(Sp, 'VELOCITIES', 1) > 0);
 WantFrc := (Sp > 0) And (FindSec(Sp, 'FORCES', 1) > 0);
 EachRes := EachOf('RESTART', 'MD', 20);
 If EachRes < 1 Then EachRes := 1;
 EachTraj := EachOf('TRAJECTORY', 'MD', 1);
 EachVelo := EachOf('VELOCITIES', 'MD', 1);
 EachFrc := EachOf('FORCES', 'MD', 1);
 If EachTraj < 1 Then EachTraj := 1;
 If EachVelo < 1 Then EachVelo := 1;
 If EachFrc < 1 Then EachFrc := 1;
 { chaine de thermostats }
 For J := 1 To 10 Do Begin NhXi[J] := 0; NhV[J] := 0; End;
 For J := 1 To NhM Do Begin
  NhQ[J] := KB_AU * Temp0 * TauT * TauT;
  If J = 1 Then NhQ[J] := NhQ[J] * MDDof;
 End;
 CsvrE := 0;
 { vitesses initiales }
 HaveVel := False;
 If SecSub > 0 Then HaveVel := LoadVelocities(SecSub);
 If RestartVelOK Then HaveVel := True;
 If Not HaveVel Then InitVelocities(Temp0) Else ApplyFixed;
 { energie initiale }
 EvalAll(True, PrintLevel >= PL_HIGH);
 PotE := ECur;
 KinE := KinEnergy;
 Cons0 := PotE + KinE;
 TimeNow := Time0 * FS_AU;
 If PrintLevel >= PL_LOW Then Begin
  EmitBlank;
  Emit(' ' + Rep('*', 31) + ' GO CP2K GO! ' + Rep('*', 33));
  Emit(' INITIAL POTENTIAL ENERGY[hartree]     =' + FR(PotE, 20, 10));
  Emit(' INITIAL KINETIC ENERGY[hartree]       =' + FR(KinE, 20, 10));
  Emit(' INITIAL TEMPERATURE[K]                =' + FR(InstTemp, 20, 4));
  If NPT Then Emit(' INITIAL PRESSURE[bar]                 =' + FR(InstPressure
    / BAR_AU, 20, 4));
  Emit(' ' + Rep('*', 79));
  EmitBlank;
 End;
 { fichiers }
 OkE := OpenOut(Fe, ProjFile('-1.ener'), Step0 > 0, FrE);
 OkP := OpenOut(Fp, ProjFile('-pos-1.xyz'), Step0 > 0, FrP);
 OkC := False; OkV := False; OkF := False; FrC := True; FrV := True; FrF :=
   True;
 If NPT Then OkC := OpenOut(Fc, ProjFile('-1.cell'), Step0 > 0, FrC);
 If WantVel Then OkV := OpenOut(Fv, ProjFile('-vel-1.xyz'), Step0 > 0, FrV);
 If WantFrc Then OkF := OpenOut(Ff, ProjFile('-frc-1.xyz'), Step0 > 0, FrF);
 If OkE And FrE Then WriteLn(Fe,
   '#     Step Nr.          Time[fs]        Kin.[a.u.]          Temp[K]      ' +
   ' ' +
   '     Pot.[a.u.]        Cons Qty[a.u.]        UsedTime[s]');
 If OkC And FrC Then WriteLn(Fc,
   '#   Step Nr.   Time[fs]           Ax [Angstrom]      Ay [Angstrom]      ' +
   'Az [Angstrom]      Bx [Angstrom]      By [Angstrom]      Bz [Angstrom]   ' +
   ' ' +
   '  Cx [Angstrom]      Cy [Angstrom]      Cz [Angstrom]   ' +
   'Volume [Angstrom^3]');
 { image initiale }
 If OkE And FrE Then WriteLn(Fe, Step0:10, FR(TimeNow / FS_AU, 18, 6), FR(KinE,
   20, 9), FR(InstTemp, 20, 9),
                      FR(PotE, 20, 9), FR(Cons0, 20, 9), FR(0, 20, 9));
 If OkP And FrP Then WriteXYZFrame(Fp, ' i =' + ISt(Step0, 9) + ', time =' +
   FR(TimeNow / FS_AU, 13, 3) + ', E =' + FR(PotE, 20, 10));
 SumPot := 0; SumKin := 0; SumTmp := 0; SumPres := 0; SumCpu := 0;
 LastStep := Step0;
 For Step := Step0 + 1 To Step0 + NSteps Do Begin
  If WallExceeded Then Begin
   Emit(' *** WALLTIME exceeded : MD stopped ***');
   Warn('WALLTIME depasse : dynamique interrompue');
   Break;
  End;
  LastStep := Step;
  T1 := RunSec;
  If NhOn Then NhcHalf(Dt);
  If DoCsvr Then CsvrApply(Dt / 2, TauT, Temp0);
  For I := 1 To NAtoms Do
   For K := 1 To 3 Do Begin
    SVel[I, K] := SVel[I, K] + 0.5 * Dt * SFrc[I, K] / SMass[I];
    SPos[I, K] := SPos[I, K] + Dt * SVel[I, K];
   End;
  EvalAll(True, False);
  PotE := ECur;
  For I := 1 To NAtoms Do
   For K := 1 To 3 Do SVel[I, K] := SVel[I, K] + 0.5 * Dt * SFrc[I, K] /
     SMass[I];
  If NhOn Then NhcHalf(Dt);
  If DoCsvr Then CsvrApply(Dt / 2, TauT, Temp0);
  If NPT Then Begin
   Pres := InstPressure;
   Mu := 1 - BetaP * Dt / TauP * (PrefP - Pres);
   If Mu < 0.9 Then Mu := 0.9;
   If Mu > 1.1 Then Mu := 1.1;
   Mu := Exp(Ln(Mu) / 3);
   For I := 1 To NAtoms Do
    For K := 1 To 3 Do
     If Not (FixActive And FixM[I, K]) Then SPos[I, K] := SPos[I, K] * Mu;
   For I := 1 To 3 Do
    For K := 1 To 3 Do CellH[I, K] := CellH[I, K] * Mu;
   UpdateCell;
   For I := 1 To NAtoms Do
    For K := 1 To 3 Do SVel[I, K] := SVel[I, K] / Mu;
   EvalAll(True, False);
   PotE := ECur;
  End;
  TimeNow := TimeNow + Dt;
  KinE := KinEnergy;
  Tmp := InstTemp;
  Pres := 0;
  If FEMethod = FE_FIST Then Pres := InstPressure;
  Cons := PotE + KinE;
  If NhOn Then Cons := Cons + NhEnergy;
  If DoCsvr Then Cons := Cons + CsvrE;
  SumPot := SumPot + PotE; SumKin := SumKin + KinE; SumTmp := SumTmp + Tmp;
    SumPres := SumPres + Pres;
  SumCpu := SumCpu + (RunSec - T1);
  Jn := Step - Step0;
  DriftAt := (Cons - Cons0) / (NAtoms * KB_AU);
  If PrintLevel >= PL_MEDIUM Then Begin
   Emit(' ' + Rep('*', 79));
   Emit(' ENSEMBLE TYPE                =' + PadL(Ens, 38));
   Emit(' STEP NUMBER                  =' + ISt(Step, 38));
   Emit(' TIME [fs]                    =' + FR(TimeNow / FS_AU, 38, 6));
   Emit(' CONSERVED QUANTITY [hartree] =' + FR(Cons, 38, 9));
   EmitBlank;
   Emit(Rep(' ', 47) + 'INSTANTANEOUS             AVERAGES');
   Emit(' CPU TIME [s]                 =' + FR(RunSec - T1, 22, 2) + FR(SumCpu /
     Jn, 22, 2));
   Emit(' ENERGY DRIFT PER ATOM [K]    =' + PadL(EStr(DriftAt, 9), 22) +
     PadL(EStr(DriftAt, 9), 22));
   Emit(' POTENTIAL ENERGY[hartree]    =' + PadL(EStr(PotE, 9), 22) +
     PadL(EStr(SumPot / Jn, 9), 22));
   Emit(' KINETIC ENERGY [hartree]     =' + PadL(EStr(KinE, 9), 22) +
     PadL(EStr(SumKin / Jn, 9), 22));
   Emit(' TEMPERATURE [K]              =' + FR(Tmp, 22, 2) + FR(SumTmp / Jn, 22,
     2));
   If FEMethod = FE_FIST Then
    Emit(' PRESSURE [bar]               =' + FR(Pres / BAR_AU, 22, 3) +
      FR(SumPres / Jn / BAR_AU, 22, 3));
   If NPT Then Emit(' VOLUME [angstrom^3]          =' + FR(CellVol *
     Cu(BOHR_ANG), 38, 4));
   Emit(' ' + Rep('*', 79));
   EmitBlank;
  End;
  If OkE Then WriteLn(Fe, Step:10, FR(TimeNow / FS_AU, 18, 6), FR(KinE, 20, 9),
    FR(Tmp, 20, 9),
                       FR(PotE, 20, 9), FR(Cons, 20, 9), FR(SumCpu, 20, 9));
  If OkP And (Jn Mod EachTraj = 0) Then
   WriteXYZFrame(Fp, ' i =' + ISt(Step, 9) + ', time =' + FR(TimeNow / FS_AU,
     13, 3) + ', E =' + FR(PotE, 20, 10));
  If OkC Then Begin
   Write(Fc, Step:8, FR(TimeNow / FS_AU, 12, 3));
   For I := 1 To 3 Do For K := 1 To 3 Do Write(Fc, FR(CellH[I, K] * BOHR_ANG,
     19, 9));
   WriteLn(Fc, FR(CellVol * Cu(BOHR_ANG), 20, 9));
  End;
  If OkV And (Jn Mod EachVelo = 0) Then Begin
   WriteLn(Fv, NAtoms:7);
   WriteLn(Fv, ' i =' + ISt(Step, 9) + ', time =' + FR(TimeNow / FS_AU, 13, 3));
   For I := 1 To NAtoms Do
    WriteLn(Fv, ' ', PadR(PtSym[SZ[I]], 2), FR(SVel[I,1], 20, 10), FR(SVel[I,2],
      20, 10), FR(SVel[I,3], 20, 10));
  End;
  If OkF And (Jn Mod EachFrc = 0) Then Begin
   WriteLn(Ff, NAtoms:7);
   WriteLn(Ff, ' i =' + ISt(Step, 9) + ', time =' + FR(TimeNow / FS_AU, 13, 3) +
     ', E =' + FR(PotE, 20, 10));
   For I := 1 To NAtoms Do
    WriteLn(Ff, ' ', PadR(PtSym[SZ[I]], 2), FR(SFrc[I,1], 20, 10), FR(SFrc[I,2],
      20, 10), FR(SFrc[I,3], 20, 10));
  End;
  If (Jn Mod EachRes = 0) Then WriteRestart(Step, TimeNow);
 End;
 WriteRestart(LastStep, TimeNow);
 If OkE Then Close(Fe);
 If OkP Then Close(Fp);
 If OkC Then Close(Fc);
 If OkV Then Close(Fv);
 If OkF Then Close(Ff);
 EmitBlank;
 BoxBanner('MD run completed', '');
 EmitBlank;
 PrintEnergyLine;
End;

{ --------------------------------------------------------------------- }
{ Analyse vibrationnelle                                                }
{ --------------------------------------------------------------------- }

Procedure RunVib;
Var
 Sv, N3, I, J, K, L, NP, NM, Ok: Integer;
 Dx, Mi, Mj, Lam, Nrm, Dd, Zpe, Fmax: Real;
 H, V, W, A, Pj: PRealArr;
 Fp, Fm: PRealArr;
 Vec, Bas: Array[1..6] Of PRealArr;
 Rcm: TVec3; Mt: Real;
 Idx: Array[1..MAX_VIB * 3] Of Integer;
 Freqs: Array[1..MAX_VIB * 3] Of Real;
 X0: Real;
 Fmol: Text;
 Full: Boolean;
 Tmp: Real;
 Cnt: Integer;
Begin
 If NAtoms > MAX_VIB Then Abort('Analyse vibrationnelle limitee a ' +
   IStr(MAX_VIB) + ' atomes');
 Sv := FindSec(0, 'VIBRATIONAL_ANALYSIS', 1);
 Dx := 0.01;
 Full := False;
 If Sv > 0 Then Begin
  Dx := KwReal(Sv, 'DX', 0.01, 'bohr');
  Full := KwBool(Sv, 'FULLY_PERIODIC', False);
 End;
 N3 := 3 * NAtoms;
 If NFixComp > 0 Then
   Warn('Les atomes fixes sont ignores par VIBRATIONAL_ANALYSIS');
 FixActive := False;
 EmitBlank;
 BoxBanner('VIBRATIONAL ANALYSIS', 'Finite difference of the forces');
 EvalAll(True, PrintLevel >= PL_MEDIUM);
 EmitBlank;
 PrintEnergyLine;
 Fmax := 0;
 For I := 1 To NAtoms Do For K := 1 To 3 Do If Abs(SFrc[I, K]) > Fmax Then Fmax
   := Abs(SFrc[I, K]);
 If Fmax > 1E-3 Then
  Warn('Les forces ne sont pas nulles (max ' + EStr(Fmax, 3) +
    ' Ha/bohr) : la geometrie n''est pas un minimum');
 H := NewVec(N3 * N3); V := NewVec(N3 * N3); W := NewVec(N3); A := NewVec(N3 *
   N3); Pj := NewVec(N3 * N3);
 Fp := NewVec(N3); Fm := NewVec(N3);
 For I := 1 To N3 Do Begin
  K := (I - 1) Mod 3 + 1; J := (I - 1) Div 3 + 1;
  X0 := SPos[J, K];
  SPos[J, K] := X0 + Dx;
  EvalAll(True, False);
  For L := 0 To N3 - 1 Do Fp^[L] := SFrc[L Div 3 + 1, L Mod 3 + 1];
  SPos[J, K] := X0 - Dx;
  EvalAll(True, False);
  For L := 0 To N3 - 1 Do Fm^[L] := SFrc[L Div 3 + 1, L Mod 3 + 1];
  SPos[J, K] := X0;
  For L := 0 To N3 - 1 Do H^[(I - 1) * N3 + L] := -(Fp^[L] - Fm^[L]) / (2 * Dx);
  If PrintLevel >= PL_HIGH Then Emit(' VIB| Displacement ' + IStr(I) + '/' +
    IStr(N3) + ' done');
 End;
 EvalAll(True, False);
 { symetrisation et ponderation par les masses }
 For I := 1 To N3 Do
  For J := I To N3 Do Begin
   Tmp := 0.5 * (H^[(I - 1) * N3 + J - 1] + H^[(J - 1) * N3 + I - 1]);
   Mi := SMass[(I - 1) Div 3 + 1]; Mj := SMass[(J - 1) Div 3 + 1];
   Tmp := Tmp / Sqrt(Mi * Mj);
   A^[(I - 1) * N3 + J - 1] := Tmp;
   A^[(J - 1) * N3 + I - 1] := Tmp;
  End;
 { projecteur des translations et des rotations }
 NP := 0;
 Mt := 0;
 For K := 1 To 3 Do Rcm[K] := 0;
 For I := 1 To NAtoms Do Begin
  Mt := Mt + SMass[I];
  For K := 1 To 3 Do Rcm[K] := Rcm[K] + SMass[I] * SPos[I, K];
 End;
 For K := 1 To 3 Do Rcm[K] := Rcm[K] / Mt;
 For K := 1 To 6 Do Vec[K] := NewVec(N3);
 For K := 1 To 3 Do
  For I := 1 To NAtoms Do Vec[K]^[(I - 1) * 3 + K - 1] := Sqrt(SMass[I]);
 For I := 1 To NAtoms Do Begin
  Vec[4]^[(I - 1) * 3 + 1] := -Sqrt(SMass[I]) * (SPos[I, 3] - Rcm[3]);
  Vec[4]^[(I - 1) * 3 + 2] := Sqrt(SMass[I]) * (SPos[I, 2] - Rcm[2]);
  Vec[5]^[(I - 1) * 3] := Sqrt(SMass[I]) * (SPos[I, 3] - Rcm[3]);
  Vec[5]^[(I - 1) * 3 + 2] := -Sqrt(SMass[I]) * (SPos[I, 1] - Rcm[1]);
  Vec[6]^[(I - 1) * 3] := -Sqrt(SMass[I]) * (SPos[I, 2] - Rcm[2]);
  Vec[6]^[(I - 1) * 3 + 1] := Sqrt(SMass[I]) * (SPos[I, 1] - Rcm[1]);
 End;
 { Gram-Schmidt }
 Cnt := 0;
 For K := 1 To 6 Do Begin
  If Full And (K > 3) Then Break;
  For L := 1 To Cnt Do Begin
   Dd := DotV(N3, Vec[K], Bas[L]);
   For I := 0 To N3 - 1 Do Vec[K]^[I] := Vec[K]^[I] - Dd * Bas[L]^[I];
  End;
  Nrm := Sqrt(DotV(N3, Vec[K], Vec[K]));
  If Nrm > 1E-6 Then Begin
   Inc(Cnt);
   Bas[Cnt] := NewVec(N3);
   For I := 0 To N3 - 1 Do Bas[Cnt]^[I] := Vec[K]^[I] / Nrm;
  End;
 End;
 NP := Cnt;
 { P = 1 - sum v v^T ; A := P A P }
 For I := 0 To N3 - 1 Do
  For J := 0 To N3 - 1 Do Begin
   Pj^[I * N3 + J] := 0;
   If I = J Then Pj^[I * N3 + J] := 1;
   For L := 1 To NP Do Pj^[I * N3 + J] := Pj^[I * N3 + J] - Bas[L]^[I] *
     Bas[L]^[J];
  End;
 { H := P A ; A := H P }
 For I := 0 To N3 - 1 Do
  For J := 0 To N3 - 1 Do Begin
   Tmp := 0;
   For L := 0 To N3 - 1 Do Tmp := Tmp + Pj^[I * N3 + L] * A^[L * N3 + J];
   H^[I * N3 + J] := Tmp;
  End;
 For I := 0 To N3 - 1 Do
  For J := 0 To N3 - 1 Do Begin
   Tmp := 0;
   For L := 0 To N3 - 1 Do Tmp := Tmp + H^[I * N3 + L] * Pj^[L * N3 + J];
   A^[I * N3 + J] := Tmp;
  End;
 Jacobi(N3, N3, A, V, W);
 { les NP valeurs propres de plus petit module sont les modes projetes }
 For I := 1 To N3 Do Idx[I] := I - 1;
 For I := 1 To N3 - 1 Do
  For J := I + 1 To N3 Do
   If Abs(W^[Idx[J]]) < Abs(W^[Idx[I]]) Then Begin
    K := Idx[I]; Idx[I] := Idx[J]; Idx[J] := K;
   End;
 { modes restants tries par valeur propre croissante }
 NM := N3 - NP;
 For I := 1 To NM Do Idx[I] := Idx[I + NP];
 For I := 1 To NM - 1 Do
  For J := I + 1 To NM Do
   If W^[Idx[J]] < W^[Idx[I]] Then Begin K := Idx[I]; Idx[I] := Idx[J]; Idx[J]
     := K; End;
 Zpe := 0;
 For I := 1 To NM Do Begin
  Lam := W^[Idx[I]];
  If Lam >= 0 Then Freqs[I] := Sqrt(Lam) * FREQ_FAC Else Freqs[I] := -Sqrt(-Lam)
    * FREQ_FAC;
  If Freqs[I] > 0 Then Zpe := Zpe + 0.5 * Freqs[I];
 End;
 EmitBlank;
 Emit(' VIB| Vibrational Analysis Info');
 Info('VIB', 'Number of atoms:', ISt(NAtoms, 20));
 Info('VIB', 'Degrees of freedom:', ISt(N3, 20));
 Info('VIB', 'Projected translations and rotations:', ISt(NP, 20));
 Info('VIB', 'Finite difference step DX [bohr]:', FR(Dx, 20, 6));
 EmitBlank;
 Emit(' VIB| Frequencies after mass weighting and projection');
 Emit(' VIB| Mode        Frequency [cm^-1]');
 For I := 1 To NM Do
  Emit(' VIB|' + ISt(I, 5) + FR(Freqs[I], 25, 6));
 EmitBlank;
 Info('VIB', 'Zero point energy [cm^-1]:', FR(Zpe, 20, 4));
 Info('VIB', 'Zero point energy [hartree]:', FR(Zpe / FREQ_FAC, 20, 10));
 Info('VIB', 'Zero point energy [kcal/mol]:', FR(Zpe / FREQ_FAC * HARTREE_KCAL,
   20, 4));
 EmitBlank;
 { fichier Molden }
 {$I-}
 Assign(Fmol, ProjFile('-VIBRATIONS-1.mol'));
 Rewrite(Fmol);
 {$I+}
 Ok := IOResult;
 If Ok = 0 Then Begin
  WriteLn(Fmol, ' [Molden Format]');
  WriteLn(Fmol, ' [FREQ]');
  For I := 1 To NM Do WriteLn(Fmol, FR(Freqs[I], 18, 8));
  WriteLn(Fmol, ' [FR-COORD]');
  For I := 1 To NAtoms Do
   WriteLn(Fmol, ' ', PadR(PtSym[SZ[I]], 2), FR(SPos[I,1], 18, 8), FR(SPos[I,2],
     18, 8), FR(SPos[I,3], 18, 8));
  WriteLn(Fmol, ' [FR-NORM-COORD]');
  For I := 1 To NM Do Begin
   WriteLn(Fmol, ' vibration', I:4);
   For J := 1 To NAtoms Do
    WriteLn(Fmol, FR(V^[((J - 1) * 3) * N3 + Idx[I]] / Sqrt(SMass[J]), 18, 8),
            FR(V^[((J - 1) * 3 + 1) * N3 + Idx[I]] / Sqrt(SMass[J]), 18, 8),
            FR(V^[((J - 1) * 3 + 2) * N3 + Idx[I]] / Sqrt(SMass[J]), 18, 8));
  End;
  Close(Fmol);
 End;
 For K := 1 To 6 Do FreeVec(Vec[K], N3);
 For K := 1 To NP Do FreeVec(Bas[K], N3);
 FreeVec(H, N3 * N3); FreeVec(V, N3 * N3); FreeVec(W, N3); FreeVec(A, N3 * N3);
   FreeVec(Pj, N3 * N3);
 FreeVec(Fp, N3); FreeVec(Fm, N3);
End;

{ --------------------------------------------------------------------- }
{ Verification des derivees                                             }
{ --------------------------------------------------------------------- }

Procedure RunDebug;
Var
 Sd, I, K: Integer;
 Dx, E1, E2, Fnum, Fan, Diff, Err, SumD, X0, Delta: Real;
 FA: TXYZ; DA: TMat3; Hb: TMat3; Pos0: TXYZ; SMax: Real;
 DoF, DoS: Boolean; A, B: Integer;
Begin
 If FEMethod <> FE_FIST Then
   Abort('RUN_TYPE DEBUG n''est disponible qu''avec METHOD FIST');
 Sd := FindSec(0, 'DEBUG', 1);
 Dx := 0.001;
 DoF := True; DoS := False;
 If Sd > 0 Then Begin
  Dx := KwReal(Sd, 'DX', 0.001, 'bohr');
  DoF := KwBool(Sd, 'DEBUG_FORCES', True);
  DoS := KwBool(Sd, 'DEBUG_STRESS_TENSOR', False);
 End;
 EvalAll(True, False);
 FA := SFrc; DA := DEps;
 EmitBlank;
 PrintEnergyLine;
 If DoF Then Begin
  EmitBlank;
  Emit(' DEBUG| Atom   Coordinate   f(numerical)   f(analytical)   ' +
    'Difference ' +
    '  Error [%]');
  SumD := 0;
  For I := 1 To NAtoms Do
   For K := 1 To 3 Do Begin
    X0 := SPos[I, K];
    SPos[I, K] := X0 + Dx; E1 := MMEval(False);
    SPos[I, K] := X0 - Dx; E2 := MMEval(False);
    SPos[I, K] := X0;
    Fnum := -(E1 - E2) / (2 * Dx);
    Fan := FA[I, K];
    Diff := Fnum - Fan;
    If Abs(Fan) > 1E-8 Then Err := 100 * Abs(Diff) / Abs(Fan) Else Err := 0;
    SumD := SumD + Abs(Diff);
    Emit(' DEBUG|' + ISt(I, 5) + '      ' + Chr(Ord('x') + K - 1) + FR(Fnum, 17,
      10) + FR(Fan, 17, 10) +
         FR(Diff, 17, 10) + FR(Err, 10, 3));
   End;
  Emit(' DEBUG| Sum of differences:' + FR(SumD, 20, 10));
  EmitBlank;
 End;
 If DoS And HasCell Then Begin
  Pos0 := SPos; Hb := CellH;
  Emit(' DEBUG| Cell strain  d(E)/d(eps) numerical  analytical   Difference');
  SMax := 0;
  For A := 1 To 3 Do
   For B := 1 To 3 Do Begin
    Delta := 1E-4;
    SPos := Pos0; CellH := Hb;
    For I := 1 To NAtoms Do SPos[I, B] := SPos[I, B] + Delta * Pos0[I, A];
    For I := 1 To 3 Do CellH[I, B] := CellH[I, B] + Delta * Hb[I, A];
    UpdateCell;
    E1 := MMEval(False);
    SPos := Pos0; CellH := Hb;
    For I := 1 To NAtoms Do SPos[I, B] := SPos[I, B] - Delta * Pos0[I, A];
    For I := 1 To 3 Do CellH[I, B] := CellH[I, B] - Delta * Hb[I, A];
    UpdateCell;
    E2 := MMEval(False);
    Fnum := (E1 - E2) / (2 * Delta);
    Diff := Fnum - DA[A, B];
    If Abs(Diff) > SMax Then SMax := Abs(Diff);
    Emit(' DEBUG|' + ISt(A, 4) + ISt(B, 4) + FR(Fnum, 20, 10) + FR(DA[A, B], 20,
      10) + FR(Diff, 20, 10));
   End;
  SPos := Pos0; CellH := Hb;
  UpdateCell;
  Emit(' DEBUG| Maximum stress difference:' + FR(SMax, 20, 10));
  EmitBlank;
 End;
End;


Const
 NSCHEMA = 82;
 Schema : Array[1..NSCHEMA] Of String[240] = (
  'GLOBAL:PROJECT PROJECT_NAME RUN_TYPE PRINT_LEVEL SEED WALLTIME ' +
  'PREFERRED_FFT_LIBRARY FFTW_PLAN_TYPE FFT_POOL_SCRATCH_LIMIT ' +
  'EXTENDED_FFT_LENGTHS ECHO_INPUT ECHO_ALL_HOSTS FLUSH_SHOULD_FLUSH SAVE_MEM',
  'GLOBAL:PROGRAM_NAME TRACE TRACE_MASTER TRACE_MAX TRACE_ROUTINES ' +
  'ENABLE_MPI_MAPPING PREFERRED_DIAG_LIBRARY SPHERICAL_HARMONICS CALLGRAPH',
  'GLOBAL/TIMINGS:THRESHOLD SORT_BY_SELF_TIME TIME_MPI REPORT_MAXLOCAL',
  'GLOBAL/DBCSR:*',
  'GLOBAL/FM:*',
  'GLOBAL/FM_DIAG_SETTINGS:*',
  'GLOBAL/REFERENCES:*',
  'FORCE_EVAL:METHOD STRESS_TENSOR NCOPIES',
  'FORCE_EVAL/SUBSYS:*',
  'FORCE_EVAL/SUBSYS/CELL:A B C ABC ALPHA_BETA_GAMMA PERIODIC ' +
  'MULTIPLE_UNIT_CELL SYMMETRY CELL_FILE_NAME CELL_FILE_FORMAT',
  'FORCE_EVAL/SUBSYS/CELL/CELL_REF:A B C ABC ALPHA_BETA_GAMMA PERIODIC ' +
  'MULTIPLE_UNIT_CELL SYMMETRY CELL_FILE_NAME CELL_FILE_FORMAT',
  'FORCE_EVAL/SUBSYS/COORD:*',
  'FORCE_EVAL/SUBSYS/VELOCITY:*',
  'FORCE_EVAL/SUBSYS/KIND:ELEMENT BASIS_SET POTENTIAL MASS GHOST AUTO_BASIS ' +
  'NAME CORE_CORRECTION LEBEDEV_GRID RADIAL_GRID ATOMIC_RADIUS DFT_PLUS_U',
  'FORCE_EVAL/SUBSYS/TOPOLOGY:COORD_FILE_NAME COORD_FILE_FORMAT ' +
  'CONN_FILE_NAME CONN_FILE_FORMAT MOL_CHECK USE_ELEMENT_AS_KIND CHARGE_BETA ' +
  'CHARGE_EXTENDED CHARGE_OCCUP EXCLUDE_VDW EXCLUDE_EI MULTIPLE_UNIT_CELL ' +
  'NUMBER_OF_ATOMS',
  'FORCE_EVAL/SUBSYS/TOPOLOGY:DISABLE_EXCLUSION_LISTS AUTOGEN_EXCLUDE_LISTS',
  'FORCE_EVAL/SUBSYS/TOPOLOGY/GENERATE:BONDPARM BONDPARM_FACTOR ' +
  'CREATE_MOLECULES REORDER_ATOMS BONDLENGTH_MAX BONDLENGTH_MIN',
  'FORCE_EVAL/SUBSYS/TOPOLOGY/CENTER_COORDINATES:CENTER_POINT',
  'FORCE_EVAL/SUBSYS/TOPOLOGY/DUMP_PDB:*',
  'FORCE_EVAL/SUBSYS/TOPOLOGY/DUMP_PSF:*',
  'FORCE_EVAL/SUBSYS/COLVAR:*',
  'FORCE_EVAL/SUBSYS/SHELL_COORD:*',
  'FORCE_EVAL/SUBSYS/CORE_COORD:*',
  'FORCE_EVAL/MM:*',
  'FORCE_EVAL/MM/FORCEFIELD:PARMTYPE PARM_FILE_NAME DO_NONBONDED ' +
  'SHIFT_CUTOFF ' +
  'EI_SCALE14 VDW_SCALE14 IGNORE_MISSING_CRITICAL_PARAMS DO_ELECTROSTATICS',
  'FORCE_EVAL/MM/FORCEFIELD/CHARGE:ATOM CHARGE',
  'FORCE_EVAL/MM/FORCEFIELD/NONBONDED:*',
  'FORCE_EVAL/MM/FORCEFIELD/NONBONDED/LENNARD-JONES:ATOMS EPSILON SIGMA RCUT ' +
  'RMIN RMAX',
  'FORCE_EVAL/MM/FORCEFIELD/NONBONDED/WILLIAMS:ATOMS A B C RCUT RMIN RMAX',
  'FORCE_EVAL/MM/FORCEFIELD/BOND:ATOMS K R0 KIND CS',
  'FORCE_EVAL/MM/FORCEFIELD/BEND:ATOMS K THETA0 KIND CB R012 R032 KBS12 ' +
  'KBS32 ' +
  'KSS THETA0',
  'FORCE_EVAL/MM/FORCEFIELD/TORSION:ATOMS K PHI0 M KIND',
  'FORCE_EVAL/MM/FORCEFIELD/SPLINE:RCUT_NB EMAX_SPLINE EPS_SPLINE ' +
  'EMAX_ACCURACY NPOINTS UNIQUE_SPLINE',
  'FORCE_EVAL/MM/POISSON:POISSON_SOLVER PERIODIC',
  'FORCE_EVAL/MM/POISSON/EWALD:EWALD_TYPE ALPHA GMAX O_SPLINE NS_MAX EPSILON ' +
  'EWALD_ACCURACY RCUT',
  'FORCE_EVAL/MM/NEIGHBOR_LISTS:*',
  'FORCE_EVAL/DFT:BASIS_SET_FILE_NAME POTENTIAL_FILE_NAME CHARGE ' +
  'MULTIPLICITY ' +
  'UKS LSD ROKS WFN_RESTART_FILE_NAME RESTART_WFN PLUS_U_METHOD AUTO_BASIS ' +
  'RELAX_MULTIPLICITY SUBCELLS',
  'FORCE_EVAL/DFT/QS:METHOD EPS_DEFAULT EPS_PGF_ORB EPS_FILTER_MATRIX ' +
  'EXTRAPOLATION EXTRAPOLATION_ORDER MAP_CONSISTENT REF_EMBED_SUBSYS CORE_PPL',
  'FORCE_EVAL/DFT/MGRID:CUTOFF REL_CUTOFF NGRIDS COMMENSURATE ' +
  'PROGRESSION_FACTOR MULTIGRID_SET SKIP_LOAD_BALANCE_DISTRIBUTED',
  'FORCE_EVAL/DFT/XC:DENSITY_CUTOFF GRADIENT_CUTOFF TAU_CUTOFF',
  'FORCE_EVAL/DFT/XC/XC_FUNCTIONAL:*',
  'FORCE_EVAL/DFT/XC/XC_GRID:*',
  'FORCE_EVAL/DFT/XC/HF:FRACTION',
  'FORCE_EVAL/DFT/XC/HF/SCREENING:*',
  'FORCE_EVAL/DFT/XC/HF/MEMORY:*',
  'FORCE_EVAL/DFT/XC/VDW_POTENTIAL:*',
  'FORCE_EVAL/DFT/SCF:EPS_SCF MAX_SCF SCF_GUESS EPS_EIGVAL ADDED_MOS ' +
  'EPS_DIIS ' +
  'MAX_DIIS CHOLESKY IGNORE_CONVERGENCE_FAILURE MAX_ITER_LUMOS EPS_LUMOS',
  'FORCE_EVAL/DFT/SCF/OT:*',
  'FORCE_EVAL/DFT/SCF/DIAGONALIZATION:*',
  'FORCE_EVAL/DFT/SCF/MIXING:*',
  'FORCE_EVAL/DFT/SCF/SMEAR:*',
  'FORCE_EVAL/DFT/SCF/OUTER_SCF:*',
  'FORCE_EVAL/DFT/POISSON:*',
  'FORCE_EVAL/DFT/KPOINTS:*',
  'FORCE_EVAL/DFT/LOCALIZE:*',
  'FORCE_EVAL/DFT/SCCS:*',
  'FORCE_EVAL/PROPERTIES:*',
  'MOTION:*',
  'MOTION/GEO_OPT:TYPE OPTIMIZER MAX_ITER MAX_DR RMS_DR MAX_FORCE RMS_FORCE ' +
  'STEP_START_VAL',
  'MOTION/GEO_OPT/BFGS:TRUST_RADIUS RESTART_HESSIAN USE_MODEL_HESSIAN ' +
  'USE_RAT_FUN_OPT',
  'MOTION/GEO_OPT/CG:MAX_STEEP_STEPS RESTART_LIMIT',
  'MOTION/GEO_OPT/CG/LINE_SEARCH:*',
  'MOTION/GEO_OPT/LBFGS:MAX_H_RANK MAX_F_PER_ITER TRUST_RADIUS',
  'MOTION/CELL_OPT:TYPE OPTIMIZER MAX_ITER MAX_DR RMS_DR MAX_FORCE RMS_FORCE ' +
  'PRESSURE_TOLERANCE EXTERNAL_PRESSURE KEEP_ANGLES KEEP_SYMMETRY CONSTRAINT ' +
  'STEP_START_VAL',
  'MOTION/CELL_OPT/BFGS:TRUST_RADIUS RESTART_HESSIAN USE_MODEL_HESSIAN',
  'MOTION/CELL_OPT/CG:MAX_STEEP_STEPS RESTART_LIMIT',
  'MOTION/CELL_OPT/CG/LINE_SEARCH:*',
  'MOTION/CELL_OPT/LBFGS:MAX_H_RANK MAX_F_PER_ITER TRUST_RADIUS',
  'MOTION/MD:ENSEMBLE STEPS TIMESTEP TEMPERATURE TEMP_TOL TEMP_KIND ' +
  'STEP_START_VAL TIME_START_VAL COMVEL_TOL ANGVEL_TOL ANGVEL_ZERO ' +
  'INITIALIZATION_METHOD SCALE_TEMP_KIND',
  'MOTION/MD/THERMOSTAT:TYPE REGION',
  'MOTION/MD/THERMOSTAT/NOSE:LENGTH YOSHIDA TIMECON MTS',
  'MOTION/MD/THERMOSTAT/CSVR:TIMECON',
  'MOTION/MD/THERMOSTAT/DEFINE_REGION:*',
  'MOTION/MD/BAROSTAT:PRESSURE TIMECON TEMPERATURE TEMP_TOL VIRIAL',
  'MOTION/MD/AVERAGES:*',
  'MOTION/CONSTRAINT:SHAKE_TOLERANCE ROLL_TOLERANCE',
  'MOTION/CONSTRAINT/FIXED_ATOMS:LIST COMPONENTS_TO_FIX MM_SUBSYS QM_SUBSYS ' +
  'MOLNAME',
  'MOTION/PRINT:*',
  'VIBRATIONAL_ANALYSIS:DX FULLY_PERIODIC INTENSITIES NPROC_REP TC_PRESSURE ' +
  'TC_TEMPERATURE THERMOCHEMISTRY ROTATIONAL_SYMMETRY',
  'DEBUG:DX DEBUG_FORCES DEBUG_STRESS_TENSOR DEBUG_DIPOLE ' +
  'DEBUG_POLARIZABILITY CHECK_ATOM_FORCE EPS_NO_ERROR_CHECK ' +
  'MAX_RELATIVE_ERROR STOP_ON_MISMATCH',
  'EXT_RESTART:RESTART_FILE_NAME RESTART_DEFAULT RESTART_POS RESTART_VEL ' +
  'RESTART_RANDOMG RESTART_SCF RESTART_HOMO RESTART_COUNTERS ' +
  'RESTART_THERMOSTAT RESTART_CELL RESTART_BAROSTAT RESTART_CONSTRAINT ' +
  'RESTART_METADYNAMICS',
  'EXT_RESTART:RESTART_QMMM RESTART_BAND RESTART_OPTIMIZE_INPUT_VARIABLES ' +
  'RESTART_AVERAGES RESTART_GEOOPT SCF_GUESS_RESTART_FILE_NAME ' +
  'BINARY_RESTART_FILE_NAME');


{ --------------------------------------------------------------------- }
{ Verification de l'entree (--check)                                    }
{ --------------------------------------------------------------------- }

Const
 PLNames : Array[0..4] Of String[8] = ('SILENT', 'LOW', 'MEDIUM', 'HIGH',
   'DEBUG');

Var
 SecUnknown : Array[0..MAX_SEC] Of Boolean;
 NCheckWarn : Integer;

Function InWordList(Const List, W: String): Boolean;
Begin
 InWordList := Pos(' ' + W + ' ', ' ' + List + ' ') > 0;
End;

{ Vrai si le chemin figure dans le schema }
Function SchemaHasPath(Const Path: String): Boolean;
Var I, P: Integer;
Begin
 SchemaHasPath := True;
 For I := 1 To NSCHEMA Do Begin
  P := Pos(':', Schema[I]);
  If Copy(Schema[I], 1, P - 1) = Path Then Exit;
 End;
 SchemaHasPath := False;
End;

{ Vrai si le mot-cle W (ou *) est admis pour ce chemin }
Function SchemaAccepts(Const Path, W: String): Boolean;
Var I, P: Integer;
Begin
 SchemaAccepts := True;
 For I := 1 To NSCHEMA Do Begin
  P := Pos(':', Schema[I]);
  If Copy(Schema[I], 1, P - 1) = Path Then
   If InWordList(Copy(Schema[I], P + 1, 255), W) Or InWordList(Copy(Schema[I], P
     + 1, 255), '*') Then Exit;
 End;
 SchemaAccepts := False;
End;

{ Chemin du parent }
Function ParentPath(Const Path: String): String;
Var I: Integer;
Begin
 I := Length(Path);
 While (I > 0) And (Path[I] <> '/') Do Dec(I);
 If I = 0 Then ParentPath := '' Else ParentPath := Copy(Path, 1, I - 1);
End;

Function PathUnderPrint(Const Path: String): Boolean;
Begin
 PathUnderPrint := (Path = 'PRINT') Or (Copy(Path, Length(Path) - 5, 6) =
   '/PRINT') Or
                   (Pos('/PRINT/', Path) > 0) Or (Copy(Path, 1, 6) = 'PRINT/');
End;

Procedure CheckWarn(Line: Integer; Const Msg: String);
Begin
 Inc(NCheckWarn);
 If Line > 0 Then ConMsg('AVERTISSEMENT (ligne ' + IStr(Line) + ') : ' + Msg)
 Else ConMsg('AVERTISSEMENT : ' + Msg);
End;

Procedure CheckTree;
Var I, J, P: Integer; Path, Un: String; Ok: Boolean;
Begin
 NCheckWarn := 0;
 For I := 1 To NSec Do SecUnknown[I] := False;
 For I := 1 To NSec Do Begin
  Path := SecPathName(I);
  If PathUnderPrint(Path) Then Continue;
  If (SecParent[I] > 0) And SecUnknown[SecParent[I]] Then Begin
   SecUnknown[I] := True;
   Continue;
  End;
  If Not SchemaHasPath(Path) Then Begin
   SecUnknown[I] := True;
   If Not SchemaAccepts(ParentPath(Path), '*') Then
    CheckWarn(SecLine[I], 'Section &' + SecName[I] + ' (' + Path +
      ') non prise en charge : ignoree');
   Continue;
  End;
  For J := 1 To NKw Do
   If KwTab[J]^.Sec = I Then
    If Not SchemaAccepts(Path, KwTab[J]^.U) Then
     CheckWarn(KwTab[J]^.Line, 'Mot-cle ' + KwTab[J]^.U +
       ' inconnu dans la section &' + SecName[I] + ' : ignore');
 End;
 { unites entre crochets }
 For J := 1 To NKw Do Begin
  Un := TrimStr(KwTab[J]^.V);
  If (Un <> '') And (Un[1] = '[') Then Begin
   P := Pos(']', Un);
   If P = 0 Then InErr(KwTab[J]^.Line, 'Crochet fermant manquant apres ' +
     KwTab[J]^.U)
   Else Begin
    UnitToAU(Copy(Un, 2, P - 2), Ok);
    If Not Ok Then InErr(KwTab[J]^.Line, 'Unite inconnue : ' + Copy(Un, 1, P));
   End;
  End;
 End;
 { mots-cles hors de toute section }
 For J := 1 To NKw Do
  If KwTab[J]^.Sec = 0 Then
   InErr(KwTab[J]^.Line, 'Mot-cle ' + KwTab[J]^.U + ' hors de toute section');
End;

{ Description XML des sections et mots-cles reconnus }
Procedure WriteXml;
Var
 I, P, Top: Integer; E, Path, Kws, W, Last: String;
 Stack: Array[0..15] Of String[120];
 Procedure PopTo(Const Pth: String);
 Begin
  While (Top > 0) And (Stack[Top] <> Pth) And (Copy(Pth, 1, Length(Stack[Top]) +
    1) <> Stack[Top] + '/') Do Begin
   WriteLn(Rep(' ', 2 * Top), '</SECTION>');
   Dec(Top);
  End;
 End;
Begin
 WriteLn('<?xml version="1.0" encoding="UTF-8"?>');
 WriteLn('<CP2K_INPUT>');
 WriteLn('  <CP2K_VERSION>', CP2K_VERSION, '</CP2K_VERSION>');
 Top := 0;
 For I := 1 To NSCHEMA Do Begin
  E := Schema[I];
  P := Pos(':', E);
  Path := Copy(E, 1, P - 1);
  Kws := Copy(E, P + 1, 255);
  PopTo(Path);
  If (Top = 0) Or (Stack[Top] <> Path) Then Begin
   Last := Path;
   While Pos('/', Last) > 0 Do Delete(Last, 1, Pos('/', Last));
   Inc(Top);
   Stack[Top] := Path;
   WriteLn(Rep(' ', 2 * Top), '<SECTION name="', Last, '">');
  End;
  W := Kws;
  While W <> '' Do Begin
   P := Pos(' ', W);
   If P = 0 Then Begin
    If W <> '*' Then WriteLn(Rep(' ', 2 * Top + 2), '<KEYWORD name="', W,
      '"/>');
    W := '';
   End
   Else Begin
    If Copy(W, 1, P - 1) <> '*' Then WriteLn(Rep(' ', 2 * Top + 2),
      '<KEYWORD name="', Copy(W, 1, P - 1), '"/>');
    W := Copy(W, P + 1, 255);
   End;
  End;
 End;
 PopTo('');
 WriteLn('</CP2K_INPUT>');
End;

{ --------------------------------------------------------------------- }
{ En-tete et pied de page                                               }
{ --------------------------------------------------------------------- }

Procedure ArtLine(N: Integer; Const Lbl, Val: String);
Const
 Art : Array[1..5] Of String[24] = (
  ' **** **** ******  **  ', ' ***** ** ***  *** **   ',
    ' **    ****   ******    ',
  ' ***** ** ** ** ** **   ', '  **** **  *******  **  ');
Var S: String; K: Integer;
Begin
 S := Art[N] + Lbl;
 K := 80 - Length(S) - Length(Val);
 If K < 1 Then K := 1;
 Emit(S + Rep(' ', K) + Val);
End;

Function HostName: String;
Var S: String;
Begin
 S := GetEnv('HOSTNAME');
 If S = '' Then S := GetEnv('COMPUTERNAME');
 If S = '' Then S := 'localhost';
 HostName := S;
End;

Function UserName: String;
Var S: String;
Begin
 S := GetEnv('USER');
 If S = '' Then S := GetEnv('USERNAME');
 If S = '' Then S := 'user';
 UserName := S;
End;

Procedure PrintStartBanner;
Begin
 ArtLine(1, 'PROGRAM STARTED AT', StartDate);
 ArtLine(2, 'PROGRAM STARTED ON', HostName);
 ArtLine(3, 'PROGRAM STARTED BY', UserName);
 ArtLine(4, 'PROGRAM PROCESS ID', '0');
 ArtLine(5, 'PROGRAM STARTED IN', CurDir);
 EmitBlank;
 Info('CP2K', 'version string:', 'CP2K version ' + CP2K_VERSION);
 Info('CP2K', 'source code revision number:', 'QChimie Pascal ' + PROG_VERSION);
 Info('CP2K', 'cp2kflags:', 'pascal clone');
 Info('CP2K', 'is freely available from', 'https://www.cp2k.org/');
 Info('CP2K', 'Program compiled on', {$IFDEF
   FPC}'Free Pascal'{$ELSE}'Turbo Pascal'{$ENDIF});
 Info('CP2K', 'Input file name', BaseName(InputName));
 EmitBlank;
End;

Procedure PrintGlobalInfo;
Begin
 Info('GLOBAL', 'Force Environment number', '1');
 Info('GLOBAL', 'Basis set file name', 'BASIS_SET');
 Info('GLOBAL', 'Potential file name', 'POTENTIAL');
 Info('GLOBAL', 'MM Potential file name', 'MM_POTENTIAL');
 Info('GLOBAL', 'Coordinate file name', '__STD_INPUT__');
 Info('GLOBAL', 'Method name', 'CP2K');
 Info('GLOBAL', 'Project name', ProjectName);
 Info('GLOBAL', 'Run type', RunName);
 Info('GLOBAL', 'Global print level', PLNames[PrintLevel]);
 Info('GLOBAL', 'Total number of message passing processes', '1');
 Info('GLOBAL', 'Number of threads for this process', '1');
 Info('GLOBAL', 'This output is from process', '0');
 EmitBlank;
End;

Procedure PrintTiming;
Var Tot, Ev: Real;
 Procedure Row(Const Nm: String; Calls: LongInt; Self, Tt: Real);
 Begin
  Emit(' ' + PadR(Nm, 30) + ISt(Calls, 9) + '  1.0' + FR(Self, 9, 3) + FR(Self,
    9, 3) +
       FR(Tt, 9, 3) + FR(Tt, 9, 3));
 End;
Begin
 Tot := RunSec;
 Ev := TEval / 1000.0;
 EmitBlank;
 Emit(' ' + Rep('-', 79));
 Emit(' -' + Rep(' ', 77) + '-');
 Emit(' -' + Rep(' ', 32) + 'T I M I N G' + Rep(' ', 34) + '-');
 Emit(' -' + Rep(' ', 77) + '-');
 Emit(' ' + Rep('-', 79));
 Emit(' SUBROUTINE                       CALLS  ASD         SELF TIME        ' +
   'TOTAL TIME');
 Emit('                                MAXIMUM       AVERAGE  MAXIMUM  ' +
   'AVERAGE  MAXIMUM');
 Row('CP2K', 1, Tot - Ev, Tot);
 Row('force_env_calc_energy_force', NEval, 0.0, Ev);
 If FEMethod = FE_FIST Then Row('fist_calc_energy_force', NEval, Ev, Ev)
 Else Begin
  Row('qs_energies', NEval, 0.0, Ev);
  Row('scf_env_do_scf', QSNumScf, Ev, Ev);
 End;
 Emit(' ' + Rep('-', 79));
End;

Procedure PrintEndBanner;
Begin
 EmitBlank;
 Emit(' The number of warnings for this run is : ' + IStr(NWarn));
 EmitBlank;
 Emit(' ' + Rep('-', 79));
 ArtLine(1, 'PROGRAM ENDED AT', DateTimeStr);
 ArtLine(2, 'PROGRAM RAN ON', HostName);
 ArtLine(3, 'PROGRAM RAN BY', UserName);
 ArtLine(4, 'PROGRAM PROCESS ID', '0');
 ArtLine(5, 'PROGRAM STOPPED IN', CurDir);
End;

{ --------------------------------------------------------------------- }
{ Mise en place du systeme                                              }
{ --------------------------------------------------------------------- }

Function RunTypeOf(Const S: String): Integer;
Begin
 RunTypeOf := -1;
 If (S = 'NONE') Then RunTypeOf := RT_NONE
 Else If (S = 'ENERGY') Or (S = 'WFN_OPT') Or (S = 'WAVEFUNCTION_OPTIMIZATION')
   Then RunTypeOf := RT_ENERGY
 Else If (S = 'ENERGY_FORCE') Then RunTypeOf := RT_FORCE
 Else If (S = 'GEO_OPT') Or (S = 'GEOMETRY_OPTIMIZATION') Then RunTypeOf :=
   RT_GEOOPT
 Else If (S = 'CELL_OPT') Then RunTypeOf := RT_CELLOPT
 Else If (S = 'MD') Or (S = 'MOLECULAR_DYNAMICS') Then RunTypeOf := RT_MD
 Else If (S = 'VIBRATIONAL_ANALYSIS') Or (S = 'NORMAL_MODES') Then RunTypeOf :=
   RT_VIB
 Else If (S = 'DEBUG') Then RunTypeOf := RT_DEBUG
 Else If (S = 'MC') Or (S = 'BAND') Or (S = 'PINT') Or (S = 'BSSE') Or (S =
   'LR') Or (S = 'TAMC') Or
         (S = 'SPECTRA') Or (S = 'ELECTRONIC_SPECTRA') Or (S = 'ROT') Or (S =
           'RT_PROPAGATION') Or
         (S = 'EHRENFEST_DYN') Or (S = 'DRIVER') Or (S = 'LINEAR_RESPONSE') Or
           (S = 'TMC') Then RunTypeOf := -2;
End;

Function PrintLevelOf(Const S: String): Integer;
Begin
 PrintLevelOf := -1;
 If S = 'SILENT' Then PrintLevelOf := PL_SILENT
 Else If S = 'LOW' Then PrintLevelOf := PL_LOW
 Else If S = 'MEDIUM' Then PrintLevelOf := PL_MEDIUM
 Else If S = 'HIGH' Then PrintLevelOf := PL_HIGH
 Else If S = 'DEBUG' Then PrintLevelOf := PL_DEBUG;
End;

{ Lit &GLOBAL }
Procedure SetupGlobal;
Var S: String; N: Integer; Seed: LongInt; V: Real;
Begin
 WallTime := 0;
 SecGlobal := FindSec(0, 'GLOBAL', 1);
 ProjectName := 'PROJECT';
 RunName := 'ENERGY';
 PrintLevel := PL_MEDIUM;
 If SecGlobal > 0 Then Begin
  ProjectName := KwStr(SecGlobal, 'PROJECT_NAME', KwStr(SecGlobal, 'PROJECT',
    'PROJECT'));
  RunName := KwWord(SecGlobal, 'RUN_TYPE', 'ENERGY');
  S := KwWord(SecGlobal, 'PRINT_LEVEL', 'MEDIUM');
  N := PrintLevelOf(S);
  If N < 0 Then InErr(SecLine[SecGlobal], 'PRINT_LEVEL invalide : ' + S) Else
    PrintLevel := N;
  Seed := KwInt(SecGlobal, 'SEED', 2000);
  S := KwStr(SecGlobal, 'WALLTIME', '0');
  WallTime := 0;
  While S <> '' Do Begin
   N := Pos(':', S);
   If N = 0 Then Begin
    If ReadReal(S, V) Then WallTime := WallTime * 60 + V;
    S := '';
   End
   Else Begin
    If ReadReal(Copy(S, 1, N - 1), V) Then WallTime := WallTime * 60 + V;
    S := Copy(S, N + 1, 255);
   End;
  End;
 End
 Else Seed := 2000;
 RandSeed := Seed;
 RunType := RunTypeOf(RunName);
 If RunType = -1 Then InErr(SecLine[SecGlobal], 'RUN_TYPE invalide : ' +
   RunName)
 Else If RunType = -2 Then
  InErr(SecLine[SecGlobal], 'RUN_TYPE ' + RunName +
    ' n''est pas pris en charge par ce clone');
End;

{ Avertit des mots-cles acceptes par CP2K mais non pris en charge ici }
Procedure WarnUnused(Sec: Integer; Const Names: String);
Var L, W: String; P: Integer;
Begin
 If Sec <= 0 Then Exit;
 L := Names;
 While L <> '' Do Begin
  P := Pos(' ', L);
  If P = 0 Then Begin W := L; L := ''; End
  Else Begin W := Copy(L, 1, P - 1); L := Copy(L, P + 1, 255); End;
  If KwHas(Sec, W) Then Warn('Mot-cle ' + W + ' de &' + SecName[Sec] +
    ' non pris en charge : ignore');
 End;
End;

{ Construit tout le systeme a partir de l'arbre d'entree }
Procedure SetupSystem;
Var
 Mth: String;
 Sr, Rs, Rsub, Rmd, Src: Integer;
 RFile: String;
 OldN: Integer;
 UseRPos, UseRVel, UseRCnt: Boolean;
Begin
 SecFE := FindSec(0, 'FORCE_EVAL', 1);
 If SecFE = 0 Then Begin
  InErr(0, 'La section &FORCE_EVAL est obligatoire');
  Exit;
 End;
 Mth := KwWord(SecFE, 'METHOD', 'QS');
 If (Mth = 'QS') Or (Mth = 'QUICKSTEP') Then Begin FEMethod := FE_QS; FEName :=
   'QS'; End
 Else If Mth = 'FIST' Then Begin FEMethod := FE_FIST; FEName := 'FIST'; End
 Else Begin
  InErr(SecLine[SecFE], 'METHOD ' + Mth +
    ' n''est pas pris en charge (QS ou FIST)');
  Exit;
 End;
 SecSub := FindSec(SecFE, 'SUBSYS', 1);
 If SecSub = 0 Then Begin
  InErr(SecLine[SecFE], 'La section &SUBSYS est obligatoire');
  Exit;
 End;
 SecMM := FindSec(SecFE, 'MM', 1);
 SecDFT := FindSec(SecFE, 'DFT', 1);
 SecMotion := FindSec(0, 'MOTION', 1);
 If (FEMethod = FE_FIST) And (SecMM = 0) Then InErr(SecLine[SecFE],
   'METHOD FIST exige une section &MM');
 If (FEMethod = FE_QS) And (SecDFT = 0) Then InErr(SecLine[SecFE],
   'METHOD QS exige une section &DFT');
 If NInErr > 0 Then Exit;
 WarnUnused(FindSec(SecSub, 'CELL', 1),
   'MULTIPLE_UNIT_CELL SYMMETRY CELL_FILE_NAME CELL_FILE_FORMAT');
 WarnUnused(FindSec(SecSub, 'TOPOLOGY', 1),
   'MULTIPLE_UNIT_CELL CONN_FILE_NAME');
 If SecMM > 0 Then
  WarnUnused(FindSec(SecMM, 'FORCEFIELD', 1),
    'PARMTYPE PARM_FILE_NAME DO_ELECTROSTATICS');
 { redemarrage externe }
 UseRPos := False; UseRVel := False; UseRCnt := False;
 Src := SecSub;
 Sr := FindSec(0, 'EXT_RESTART', 1);
 RestartVelOK := False;
 UseResCounters := False;
 If (Sr > 0) And KwHas(Sr, 'RESTART_FILE_NAME') Then Begin
  RFile := KwStr(Sr, 'RESTART_FILE_NAME', '');
  OldN := NLn;
  ReadSource(RFile, 1);
  If NInErr > 0 Then Exit;
  Rs := AddSec(0, 'RESTART_DATA', '', 0);
  ParseLines(OldN + 1, NLn, Rs);
  Rsub := FindSec(FindSec(Rs, 'FORCE_EVAL', 1), 'SUBSYS', 1);
  If Rsub = 0 Then InErr(0,
    'Le fichier de redemarrage ne contient pas de &SUBSYS');
  UseRPos := KwBool(Sr, 'RESTART_POS', KwBool(Sr, 'RESTART_DEFAULT', True));
  UseRVel := KwBool(Sr, 'RESTART_VEL', KwBool(Sr, 'RESTART_DEFAULT', True));
  UseRCnt := KwBool(Sr, 'RESTART_COUNTERS', KwBool(Sr, 'RESTART_DEFAULT',
    True));
  Rmd := FindSec(FindSec(Rs, 'MOTION', 1), 'MD', 1);
  If UseRCnt And (Rmd > 0) Then Begin
   UseResCounters := True;
   ResStep := KwInt(Rmd, 'STEP_START_VAL', 0);
   ResTime := KwReal(Rmd, 'TIME_START_VAL', 0.0, 'fs');
  End;
  If UseRPos And (Rsub > 0) Then Src := Rsub;
  If NInErr > 0 Then Exit;
 End;
 BuildCell(Src);
 NKinds := 0;
 SetupKinds;
 LoadCoords(Src);
 If NInErr > 0 Then Exit;
 If NAtoms = 0 Then Begin InErr(SecLine[SecSub], 'Aucun atome dans le systeme');
   Exit; End;
 ValidateKinds;
 If NInErr > 0 Then Exit;
 FillAtomData;
 SetupFixed;
 If UseRVel And (Rsub > 0) Then RestartVelOK := LoadVelocities(Rsub);
 If Src = SecSub Then CenterCoords;
 If (Per[1] Or Per[2] Or Per[3]) And Not HasCell Then
  InErr(SecLine[SecSub], 'Systeme periodique sans cellule');
 SetupTopology;
 If NInErr > 0 Then Exit;
 If FEMethod = FE_FIST Then MMSetup
 Else Begin
  QSReadParams;
  If NInErr = 0 Then
   If Not QSSetupSystem(QSErr) Then InErr(SecLine[SecDFT], QSErr);
 End;
End;

{ Informations sur le systeme (section de sortie) }
Procedure PrintSystemInfo;
Begin
 If PrintLevel < PL_MEDIUM Then Exit;
 If HasCell Then PrintCell
 Else Info('CELL', 'No cell defined (isolated system)', '');
 EmitBlank;
 Emit(' TOPOLOGY INFO - GENERATE');
 Info('TOPOLOGY', 'Number of atoms:', ISt(NAtoms, 15));
 Info('TOPOLOGY', 'Number of bonds:', ISt(NBond, 15));
 Info('TOPOLOGY', 'Number of molecules:', ISt(NMol, 15));
 PrintKinds;
 If FEMethod = FE_FIST Then PrintCoords('FIST', False) Else
   PrintCoords('QUICKSTEP', True);
 PrintTotals;
 If FEMethod = FE_FIST Then Begin
  Info('FIST', 'Cutoff radius [angstrom]:', FR(RcMax * BOHR_ANG, 15, 4));
  Info('FIST', 'Number of Lennard-Jones/Williams kind pairs:', ISt(NPairsLJ,
    6));
  If EwaldOn Then Begin
   Info('FIST', 'Ewald summation, alpha [1/angstrom]:', FR(EwAlpha / BOHR_ANG,
     15, 6));
   Info('FIST', 'Ewald real space cutoff [angstrom]:', FR(EwRc * BOHR_ANG, 15,
     4));
   Info('FIST', 'Ewald reciprocal cutoff [1/angstrom]:', FR(EwKc / BOHR_ANG, 15,
     4));
  End
  Else Info('FIST', 'Coulomb:', PadL(BoolStr(MMCharged), 15));
  Info('FIST', 'Bonds / bends / torsions with parameters:',
       IStr(NBond) + ' / ' + IStr(NBend) + ' / ' + IStr(NTors));
 End
 Else Begin
  EmitBlank;
  If UnrestrictedFlag Then Info('DFT',
    'Spin unrestricted Kohn-Sham (UKS) calculation', 'UKS')
  Else Info('DFT', 'Spin restricted Kohn-Sham (RKS) calculation', 'RKS');
  Info('DFT', 'Multiplicity', ISt(Mult, 15));
  Info('DFT', 'Number of spin states', ISt(1 + Ord(NOccA <> NOccB), 15));
  Info('DFT', 'Charge', ISt(Charge, 15));
  Info('DFT', 'Exchange-correlation functional', PadL(QSFuncName, 20));
  Info('QS', 'Method:', PadL('GPW (all-electron Gaussian)', 40));
  Info('QS', 'Number of basis functions:', ISt(NBF, 15));
  Info('QS', 'Number of electrons:', ISt(NElec, 15));
  Info('QS', 'EPS_SCF:', EStr(QSEps, 2));
  If QSPotIgnored Then Info('QS', 'Pseudopotentials (POTENTIAL):',
    PadL('ignored (all-electron)', 30));
 End;
 EmitBlank;
End;

{ Fonctions et procedures d'aide en ligne de commande }
Procedure ShowHelp;
Begin
 WriteLn('CP2K : Clone de la ligne de commande de CP2K (chimie quantique et ' +
   'physique');
 WriteLn('       de l''etat solide, simulations atomistiques).');
 WriteLn;
 WriteLn('Syntaxe : CP2K [options] [-i] fichier.inp [-o sortie.out]');
 WriteLn;
 WriteLn(' -i fichier          Fichier d''entree (ou premier argument sans ' +
   'option ;');
 WriteLn('                     - lit l''entree standard)');
 WriteLn(' -o fichier          Fichier de sortie (par defaut : sortie ' +
   'standard)');
 WriteLn(' -c, --check         Verifie seulement la syntaxe du fichier ' +
   'd''entree');
 WriteLn(' -e, --echo          Affiche l''entree apres pretraitement et ' +
   'analyse');
 WriteLn(' -E, --set nom=val   Definit une variable du preprocesseur (@SET)');
 WriteLn(' -r, --run           Lance le calcul (comportement par defaut)');
 WriteLn(' -v, --version       Affiche la version');
 WriteLn('     --license       Affiche la licence');
 WriteLn('     --xml           Decrit les sections et mots-cles reconnus ' +
   '(XML)');
 WriteLn(' -h, --help          Affiche cette aide');
 WriteLn;
 WriteLn('RUN_TYPE : ENERGY, ENERGY_FORCE, GEO_OPT, CELL_OPT, MD,');
 WriteLn('           VIBRATIONAL_ANALYSIS, DEBUG');
 WriteLn('METHOD   : FIST (champ de force), QS (Quickstep : HF et DFT, bases ' +
   'STO-3G');
 WriteLn('           pour H a Ne, 3-21G et 6-31G pour H, He et C a F)');
End;

Procedure ShowVersion;
Begin
 WriteLn(' CP2K version ', CP2K_VERSION, ' (clone Pascal QChimie ',
   PROG_VERSION, ')');
 WriteLn(' Source code revision number: QChimie');
 WriteLn(' cp2kflags: pascal fist qs');
End;

Procedure ShowLicense;
Begin
 WriteLn(' CP2K est un logiciel libre distribue sous licence GNU GPL v2 ou ' +
   'ulterieure.');
 WriteLn(' Ce programme est un clone ecrit en Pascal : il ne contient aucun ' +
   'code de CP2K.');
 WriteLn(' Voir https://www.cp2k.org/ pour le logiciel original.');
End;

{ --------------------------------------------------------------------- }
{ Programme principal                                                   }
{ --------------------------------------------------------------------- }

Var
 ArgI : Integer;
 Arg, UArg, Vl : String;
 EchoMode, DoXml : Boolean;
 P: Integer;

Begin
 NWarn := 0; NInErr := 0; NLn := 0; NSec := 0; NKw := 0; NPV := 0;
 CheckMode := False; EchoMode := False; DoXml := False; RunPhase := False;
 OutOpen := False; RedirOn := False; RestartMode := False;
 InputName := ''; OutputName := '';
 NEval := 0; TEval := 0; QSNumScf := 0;
 NAtoms := 0; NKinds := 0; NBond := 0;
 ErfcReady := False;
 T0Run := NowMs;
 StartDate := DateTimeStr;
 GetDir(0, CurDir);
 ClearTree;
 ArgI := 1;
 While ArgI <= ParamCount Do Begin
  Arg := ParamStr(ArgI);
  UArg := UpStr(Arg);
  If (Arg = '-h') Or (UArg = '--HELP') Or (Arg = '-?') Or (Arg = '/?') Or (UArg
    = '/H') Then Begin
   ShowHelp; Halt(0);
  End
  Else If (Arg = '-v') Or (UArg = '--VERSION') Then Begin ShowVersion; Halt(0);
    End
  Else If UArg = '--LICENSE' Then Begin ShowLicense; Halt(0); End
  Else If UArg = '--XML' Then DoXml := True
  Else If (Arg = '-c') Or (UArg = '--CHECK') Then CheckMode := True
  Else If (Arg = '-e') Or (UArg = '--ECHO') Then EchoMode := True
  Else If (Arg = '-r') Or (UArg = '--RUN') Then Begin End
  Else If Arg = '-i' Then Begin
   Inc(ArgI);
   If ArgI > ParamCount Then Begin
     ConMsg('ERREUR : -i exige un nom de fichier'); Halt(1); End;
   InputName := ParamStr(ArgI);
  End
  Else If Arg = '-o' Then Begin
   Inc(ArgI);
   If ArgI > ParamCount Then Begin
     ConMsg('ERREUR : -o exige un nom de fichier'); Halt(1); End;
   OutputName := ParamStr(ArgI);
  End
  Else If (Arg = '-E') Or (UArg = '--SET') Then Begin
   Inc(ArgI);
   If ArgI > ParamCount Then Begin ConMsg('ERREUR : -E exige nom=valeur');
     Halt(1); End;
   Vl := ParamStr(ArgI);
   P := Pos('=', Vl);
   If P = 0 Then Begin ConMsg('ERREUR : -E exige nom=valeur'); Halt(1); End;
   SetPVarLocked(Copy(Vl, 1, P - 1), Copy(Vl, P + 1, 255));
  End
  Else If (Arg <> '') And (Arg <> '-') And (Arg[1] = '-') Then Begin
   ConMsg('ERREUR : option inconnue : ' + Arg);
   ConMsg('Utilisez CP2K --help pour la liste des options.');
   Halt(1);
  End
  Else Begin
   If InputName = '' Then InputName := Arg
   Else If OutputName = '' Then OutputName := Arg
   Else Begin ConMsg('ERREUR : argument inattendu : ' + Arg); Halt(1); End;
  End;
  Inc(ArgI);
 End;
 If DoXml Then Begin WriteXml; Halt(0); End;
 If InputName = '' Then Begin
  ConMsg('ERREUR : aucun fichier d''entree specifie');
  ShowHelp;
  Halt(1);
 End;
 If (OutputName <> '') And (FExpand(OutputName) = FExpand(InputName)) And
   (InputName <> '-') Then Begin
  ConMsg('ERREUR : le fichier de sortie ne doit pas etre le fichier d''entree');
  Halt(1);
 End;
 { sortie }
 If (OutputName <> '') And Not CheckMode And Not EchoMode Then Begin
  Assign(OutF, OutputName);
  {$I-}
  Rewrite(OutF);
  {$I+}
  If IOResult <> 0 Then Begin ConMsg('ERREUR : impossible de creer ' +
    OutputName); Halt(1); End;
 End
 Else Begin
  Assign(OutF, '');
  Rewrite(OutF);
 End;
 OutOpen := True;
 If Not CheckMode And Not EchoMode Then PrintStartBanner;
 { lecture de l'entree }
 ReadSource(InputName, 1);
 If NInErr = 0 Then ParseTree;
 If NInErr = 0 Then CheckTree;
 If CheckMode Then Begin
  If NInErr = 0 Then Begin
   WriteLn(OutF, ' Input file ', InputName, ' : syntax OK (', NSec,
     ' sections, ', NKw, ' keywords, ',
           NCheckWarn, ' warnings)');
   Close(OutF);
   Halt(0);
  End;
  Close(OutF);
  Halt(1);
 End;
 If EchoMode Then Begin
  If NInErr = 0 Then DumpSec(0, 0);
  Close(OutF);
  If NInErr > 0 Then Halt(1);
  Halt(0);
 End;
 If NInErr = 0 Then SetupGlobal;
 If NInErr = 0 Then SetupSystem;
 If NInErr > 0 Then Begin
  Abort('Le fichier d''entree contient ' + IStr(NInErr) +
    ' erreur(s) (voir ci-dessus)');
 End;
 PrintGlobalInfo;
 If (SecGlobal > 0) And KwBool(SecGlobal, 'ECHO_INPUT', False) Then Begin
  Emit(' ******************************* ECHO OF INPUT ' +
    '*********************************');
  DumpSec(0, 0);
  Emit(' ' + Rep('*', 79));
  EmitBlank;
 End;
 PrintSystemInfo;
 RunPhase := True;
 Case RunType Of
  RT_NONE: ;
  RT_ENERGY: RunEnergy(False);
  RT_FORCE: RunEnergy(True);
  RT_GEOOPT: RunOpt(False);
  RT_CELLOPT: RunOpt(True);
  RT_MD: RunMD;
  RT_VIB: RunVib;
  RT_DEBUG: RunDebug;
 End;
 PrintTiming;
 PrintEndBanner;
 Close(OutF);
 OutOpen := False;
End.
