{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2026
  @version: 1.00
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal 7, Free Pascal)
  @description: Clone de l'application VCCTL (Virtual Cement and Concrete
                Testing Laboratory - Laboratoire virtuel d'essais du ciment
                et du beton) : construction d'une microstructure 3D de pate
                de ciment (ciment Portland, cendres volantes, laitier,
                calcaire, fumee de silice), simulation de l'hydratation par
                automate cellulaire puis estimation des proprietes.

  Ligne de commande :
    VCCTL [projet.VCT]
    VCCTL /B [projet.VCT] [sortie.CSV]   (mode lot : calcul sans interface)

  Fonctionnalites :
    - Bibliotheque de materiaux : ciments (types I, II, III, IV, V, blanc),
      cendres volantes, laitier, calcaire, fumee de silice, granulats
    - Granulometrie de type Rosin-Rammler, surface specifique, calculateur
      de Bogue (analyse d'oxydes -> phases du clinker)
    - Formulation : pate, mortier, beton (rapport E/L, proportions, volumes)
    - Generation d'une microstructure 3D (empilement aleatoire de spheres
      numerisees en voxels, mur de granulat pour la zone de transition)
    - Empilement 3D de granulats (courbe de Fuller) et zone de transition
    - Hydratation par automate cellulaire (dissolution, diffusion,
      precipitation de C-S-H, portlandite, ettringite, monosulfate,
      hydrogrenat, reactions pouzzolaniques), loi de temps parabolique,
      temperature (Arrhenius), cure saturee, scellee ou adiabatique
    - Calibrage de la cinetique sur des mesures experimentales
    - Proprietes : degre d'hydratation, chaleur, retrait chimique, eau liee,
      porosite, percolation des solides et des pores, distribution des
      tailles de pores, modules elastiques (milieu effectif auto-coherent),
      resistance a la compression (Powers), conductivite / diffusivite
      (differences finies), profils de la zone de transition
    - Visualisation : coupes XY / XZ / YZ, vue 3D, phases masquables
    - Import / export : projet VCT, microstructure binaire MIC et texte,
      resultats CSV, coupe BMP
    - Menus deroulants adaptes de CLONEDRAW\PAINT.PAS

  Aucune dependance externe : tout le code est dans ce fichier.
}

Program VCCTL;

{$IFDEF FPC}
 {$MODE TP}
 {$Q-}
 {$R-}
{$ELSE}
 {$M 40000,0,655360}
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
 PROG_TITLE      = 'VCCTL Clone 1.00 - QChimie';
 PROG_VERSION    = '1.00';

 MENUBAR_H       = 20;
 TOOLBAR_H       = 24;
 STATUS_H        = 18;

 NUM_MENUS       = 8;
 MAX_MENU_ITEMS  = 22;
 MENU_FILE       = 0;
 MENU_MAT        = 1;
 MENU_MIX        = 2;
 MENU_HYD        = 3;
 MENU_PROP       = 4;
 MENU_VIEW       = 5;
 MENU_TOOLS      = 6;
 MENU_HELP       = 7;

 {$IFDEF FPC}
 MAXN            = 100;          { taille maximale du systeme (voxels) }
 MAXLIST         = 400000;
 {$ELSE}
 MAXN            = 32;
 MAXLIST         = 16000;
 {$ENDIF}
 MAXD            = 50;           { diametre maximal d'une particule (voxels) }

 { Phases }
 P_POR  = 0;   P_C3S  = 1;   P_C2S  = 2;   P_C3A  = 3;   P_C4AF = 4;
 P_GYP  = 5;   P_HEM  = 6;   P_ANH  = 7;   P_LIME = 8;   P_ALK  = 9;
 P_FILL = 10;  P_SF   = 11;  P_FA   = 12;  P_SLAG = 13;  P_CSH  = 14;
 P_CH   = 15;  P_AFT  = 16;  P_AFM  = 17;  P_HG   = 18;  P_PCSH = 19;
 P_VOID = 20;  P_AGG  = 21;
 NPH    = 22;

 { Reactions }
 RX_C3S = 0; RX_C2S = 1; RX_AFT = 2; RX_AFM = 3; RX_HG = 4;
 RX_C4AF = 5; RX_LIME = 6; RX_SF = 7; RX_FA = 8; RX_SLAG = 9;
 RX_CONV = 10;
 NRX    = 11;

 { Series temporelles }
 MAXPT   = 320;
 NSER    = 10;
 SER_T   = 0;  { temps (h) }
 SER_CYC = 1;  { cycles }
 SER_ALF = 2;  { degre d'hydratation du ciment }
 SER_Q   = 3;  { chaleur (J/g de liant) }
 SER_TMP = 4;  { temperature (C) }
 SER_POR = 5;  { porosite capillaire (%) }
 SER_WN  = 6;  { eau liee (g/g de liant) }
 SER_SHR = 7;  { retrait chimique (mL/g de liant) }
 SER_FC  = 8;  { resistance (MPa) }
 SER_ALT = 9;  { degre d'hydratation du liant }

 MAX_ROWS = 24;

{ --------------------------------------------------------------------- }
{ Types                                                                 }
{ --------------------------------------------------------------------- }

Type
 TSlice  = Array[0..MAXN*MAXN-1] Of Byte;
 PSlice  = ^TSlice;
 TGrid   = Array[0..MAXN-1] Of PSlice;
 TRSlice = Array[0..MAXN*MAXN-1] Of Real;
 PRSlice = ^TRSlice;
 TRGrid  = Array[0..MAXN-1] Of PRSlice;

 TLongArr = Array[0..MAXLIST] Of LongInt;
 PLongArr = ^TLongArr;

 TPSD = Record
  D63, Nr, DMin, DMax : Real;       { Rosin-Rammler : d' (um), n, bornes (um) }
 End;

 TCement = Record
  Name : String[30];
  W    : Array[1..9] Of Real;       { % massiques : C3S C2S C3A C4AF Gyp Hem Anh Chaux AlcSO4 }
  P    : TPSD;
 End;

 TSCM = Record
  Name  : String[30];
  Rho   : Real;
  P     : TPSD;
  React : Real;                     { reactivite relative }
  Extra : Real;                     { fraction vitreuse (%) / nucleation }
 End;

 TAggr = Record
  Rho, DMin, DMax, Q : Real;        { densite, mm, mm, exposant de Fuller }
 End;

 TReact = Record
  NP   : Integer;
  Cnt  : Array[0..1] Of Real;       { voxels de produit par voxel de reactif }
  Pr   : Array[0..1] Of Byte;
  Wat, CHc, AFc, Sul, SulV, E, Heat : Real;
 End;

 TSim = Record
  Cyc      : LongInt;
  T, Te, Temp : Real;
  Heat, Wn, Shrk, Pool, SulfRes : Real;
  Lost     : LongInt;
  NS       : Integer;
  NextSamp : LongInt;
  TSet     : Real;
  M0Cl, M0Bind, M0Cem, MSol0, MW0 : Real;
  CurN     : Integer;
  FAInert  : LongInt;
  Dummy    : Array[0..7] Of Real;
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

 TR6 = Array[0..5] Of Real;
 TR10 = Array[0..9] Of Real;
 TPhR = Array[0..21] Of Real;

{ --------------------------------------------------------------------- }
{ Variables globales                                                    }
{ --------------------------------------------------------------------- }

Var
 { Materiaux }
 Cem   : TCement;
 FAm, SlgM, FilM, SFM : TSCM;
 Sand, Coarse : TAggr;

 { Formulation }
 MixKind   : Integer;               { 0 pate, 1 mortier, 2 beton }
 WB        : Real;                  { rapport eau / liant }
 Frac      : Array[0..4] Of Real;   { liant : ciment, cendre, laitier, calcaire, fumee (%) }
 NVox      : Integer;
 ResUm     : Real;                  { micrometres par voxel }
 Seed      : LongInt;
 WallOn    : Boolean;
 ConcBind, ConcSand, ConcCoarse, ConcAir : Real;   { kg/m3, kg/m3, kg/m3, % }

 { Hydratation }
 CureTemp  : Real;                  { C }
 CureMode  : Integer;               { 0 satur., 1 scelle, 2 adiabatique }
 Beta      : Real;                  { h par cycle^2 }
 EaCem     : Real;                  { kJ/mol }
 TRefC     : Real;
 EndDays   : Real;
 MaxCycles : LongInt;
 PNucCH    : Real;
 MaxAge    : Integer;
 DormCyc   : Integer;               { cycles de la periode dormante }
 AccNear, AccFar : Real;
 Rate      : Array[0..NPH-1] Of Real;   { probabilite de dissolution par cycle }
 Fc0, FcExp : Real;                 { resistance de Powers }
 SigmaPore : Real;                  { S/m }
 SigTab    : Array[0..NPH-1] Of Real; { conductivite relative des phases }
 ModE, ModNu : Array[0..NPH-1] Of Real;  { modules (GPa), coefficient de Poisson }

 { Phases }
 PhRho     : Array[0..NPH-1] Of Real;
 PhColor   : Array[0..NPH-1] Of Byte;
 PhShow    : Array[0..NPH-1] Of Boolean;

 { Microstructure }
 N, NN     : Integer;
 NTot      : LongInt;
 Grid      : TGrid;
 GridN     : Integer;
 HaveGrid  : Boolean;
 Wr        : Array[-24..MAXN+24] Of Integer;
 PhCnt     : Array[0..NPH-1] Of LongInt;
 SphVol    : Array[1..MAXD] Of LongInt;
 ZLo, ZHi  : Integer;               { plans extremes du paste (mur de granulat) }
 GenInfo   : String;

 { Simulation }
 Sim       : TSim;
 Ser       : Array[0..NSER-1, 0..MAXPT-1] Of Real;
 PhSer     : Array[0..NPH-1, 0..MAXPT-1] Of LongInt;
 Rx        : Array[0..NRX-1] Of TReact;
 Thr       : Array[0..NPH-1] Of LongInt;
 SulfList  : PLongArr;
 SulfN     : LongInt;
 MolVox    : Array[0..NPH-1] Of Real;
 Running   : Boolean;
 SimActive : Boolean;               { une simulation est en cours }
 SimGoalT  : Real;                  { temps cible (h) }
 CcapJ     : Real;

 { Divers }
 RngS      : LongInt;
 BatchMode : Boolean;
 UserAbort : Boolean;
 StatusMsg : String;
 NeedRedraw : Boolean;
 ProjName  : String;
 Modified  : Boolean;

{ Avancement (defini avec l'interface) }
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

{ Horloge en millisecondes depuis minuit }
Function NowMs: LongInt;
Var H, M, S, S100: Word;
Begin
 GetTime(H, M, S, S100);
 NowMs := LongInt(H) * 3600000 + LongInt(M) * 60000 +
          LongInt(S) * 1000 + LongInt(S100) * 10;
End;

{ --------------------------------------------------------------------- }
{ Generateur pseudo-aleatoire (Park-Miller, identique en TP et en FPC)  }
{ --------------------------------------------------------------------- }

Procedure SetSeed(S: LongInt);
Begin
 If S < 0 Then S := -S;
 S := S Mod 2147483646;
 If S = 0 Then S := 12345;
 RngS := S;
End;

Function RndRaw: LongInt;
Var Hi, Lo, T: LongInt;
Begin
 Hi := RngS Div 127773;
 Lo := RngS Mod 127773;
 T := 16807 * Lo - 2836 * Hi;
 If T <= 0 Then T := T + 2147483647;
 RngS := T;
 RndRaw := T;
End;

{ entier dans 0..M-1 }
Function RndN(M: LongInt): LongInt;
Begin
 If M <= 1 Then RndN := 0 Else RndN := RndRaw Mod M;
End;

Function Rnd15: Integer;
Begin
 Rnd15 := RndRaw Mod 32768;
End;

{ reel dans [0,1) }
Function RndF: Real;
Begin
 RndF := RndRaw / 2147483647.0;
End;

{ arrondi stochastique }
Function StochRound(V: Real): Integer;
Var I: Integer;
Begin
 If V <= 0 Then Begin StochRound := 0; Exit; End;
 I := Trunc(V);
 If RndF < V - I Then Inc(I);
 StochRound := I;
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
{ Tables des phases                                                     }
{ --------------------------------------------------------------------- }

Const
 PhName : Array[0..NPH-1] Of String[24] = (
  'Porosite (eau)', 'C3S (alite)', 'C2S (belite)', 'C3A (aluminate)',
  'C4AF (ferrite)', 'Gypse', 'Hemihydrate', 'Anhydrite', 'Chaux libre',
  'Sulfates alcalins', 'Calcaire / inerte', 'Fumee de silice',
  'Cendre volante', 'Laitier', 'C-S-H', 'Portlandite (CH)',
  'Ettringite (AFt)', 'Monosulfate (AFm)', 'Hydrogrenat', 'C-S-H pouzzol.',
  'Pores vides', 'Granulat');

 PhShort : Array[0..NPH-1] Of String[6] = (
  'POR', 'C3S', 'C2S', 'C3A', 'C4AF', 'GYP', 'HEM', 'ANH', 'CaO',
  'K2SO4', 'FILL', 'SF', 'FA', 'SLAG', 'CSH', 'CH', 'AFt', 'AFm', 'HG',
  'PCSH', 'VOID', 'AGG');

Procedure InitPhases;
Var I: Integer;
Begin
 PhRho[P_POR] := 1.00; PhRho[P_C3S] := 3.21; PhRho[P_C2S] := 3.28;
 PhRho[P_C3A] := 3.03; PhRho[P_C4AF] := 3.73; PhRho[P_GYP] := 2.32;
 PhRho[P_HEM] := 2.74; PhRho[P_ANH] := 2.61; PhRho[P_LIME] := 3.32;
 PhRho[P_ALK] := 2.66; PhRho[P_FILL] := 2.71; PhRho[P_SF] := 2.20;
 PhRho[P_FA] := 2.40; PhRho[P_SLAG] := 2.90; PhRho[P_CSH] := 2.00;
 PhRho[P_CH] := 2.24; PhRho[P_AFT] := 1.75; PhRho[P_AFM] := 2.02;
 PhRho[P_HG] := 2.52; PhRho[P_PCSH] := 2.00; PhRho[P_VOID] := 0.0;
 PhRho[P_AGG] := 2.65;
 PhColor[P_POR] := 0;  PhColor[P_C3S] := 7;  PhColor[P_C2S] := 6;
 PhColor[P_C3A] := 4;  PhColor[P_C4AF] := 14; PhColor[P_GYP] := 11;
 PhColor[P_HEM] := 11; PhColor[P_ANH] := 11; PhColor[P_LIME] := 12;
 PhColor[P_ALK] := 5;  PhColor[P_FILL] := 8; PhColor[P_SF] := 5;
 PhColor[P_FA] := 2;   PhColor[P_SLAG] := 10; PhColor[P_CSH] := 9;
 PhColor[P_CH] := 15;  PhColor[P_AFT] := 13; PhColor[P_AFM] := 3;
 PhColor[P_HG] := 5;   PhColor[P_PCSH] := 9; PhColor[P_VOID] := 1;
 PhColor[P_AGG] := 8;
 For I := 0 To NPH - 1 Do PhShow[I] := True;
 { modules elastiques (GPa) et coefficients de Poisson }
 For I := 0 To NPH - 1 Do Begin ModE[I] := 0; ModNu[I] := 0.25; End;
 ModE[P_C3S] := 135; ModNu[P_C3S] := 0.30;
 ModE[P_C2S] := 130; ModNu[P_C2S] := 0.30;
 ModE[P_C3A] := 145; ModNu[P_C3A] := 0.30;
 ModE[P_C4AF] := 125; ModNu[P_C4AF] := 0.30;
 ModE[P_GYP] := 45; ModNu[P_GYP] := 0.33;
 ModE[P_HEM] := 50; ModNu[P_HEM] := 0.33;
 ModE[P_ANH] := 75; ModNu[P_ANH] := 0.33;
 ModE[P_LIME] := 100; ModNu[P_LIME] := 0.30;
 ModE[P_ALK] := 30; ModNu[P_ALK] := 0.30;
 ModE[P_FILL] := 75; ModNu[P_FILL] := 0.30;
 ModE[P_SF] := 72; ModNu[P_SF] := 0.17;
 ModE[P_FA] := 90; ModNu[P_FA] := 0.22;
 ModE[P_SLAG] := 105; ModNu[P_SLAG] := 0.25;
 ModE[P_CSH] := 22; ModNu[P_CSH] := 0.25;
 ModE[P_CH] := 38; ModNu[P_CH] := 0.305;
 ModE[P_AFT] := 24; ModNu[P_AFT] := 0.25;
 ModE[P_AFM] := 42; ModNu[P_AFM] := 0.25;
 ModE[P_HG] := 80; ModNu[P_HG] := 0.30;
 ModE[P_PCSH] := 22; ModNu[P_PCSH] := 0.25;
 ModE[P_AGG] := 70; ModNu[P_AGG] := 0.20;
 { conductivite relative a la solution interstitielle }
 For I := 0 To NPH - 1 Do SigTab[I] := 0;
 SigTab[P_POR] := 1.0;
 SigTab[P_CSH] := 0.004;
 SigTab[P_PCSH] := 0.004;
 SigTab[P_AFT] := 0.0005;
 SigTab[P_AFM] := 0.0005;
 SigTab[P_HG] := 0.0002;
 SigTab[P_CH] := 0.0002;
 { masse molaire / volume molaire des sulfates (mol par voxel) }
 For I := 0 To NPH - 1 Do MolVox[I] := 0;
End;

{ --------------------------------------------------------------------- }
{ Bibliotheque de materiaux                                             }
{ --------------------------------------------------------------------- }

Const
 NUM_CEM_PRESETS = 6;
 CemPresetName : Array[1..NUM_CEM_PRESETS] Of String[30] = (
  'Type I (CEM I 42.5)', 'Type II (moderee)', 'Type III (CEM I 52.5R)',
  'Type IV (faible chaleur)', 'Type V (resist. sulfates)', 'Ciment blanc');

Procedure SetPSD(Var P: TPSD; D63, Nr, DMin, DMax: Real);
Begin
 P.D63 := D63; P.Nr := Nr; P.DMin := DMin; P.DMax := DMax;
End;

Procedure SetCementW(Var C: TCement; A, B, D, E, F, G, H, I, J: Real);
Begin
 C.W[1] := A; C.W[2] := B; C.W[3] := D; C.W[4] := E; C.W[5] := F;
 C.W[6] := G; C.W[7] := H; C.W[8] := I; C.W[9] := J;
End;

Procedure LoadCementPreset(K: Integer);
Var I: Integer; S: Real;
Begin
 Cem.Name := CemPresetName[K];
 Case K Of
  1: Begin SetCementW(Cem, 57.0, 17.0, 9.5, 8.0, 3.2, 0.5, 0.5, 0.8, 0.6); SetPSD(Cem.P, 17, 1.05, 1, 60); End;
  2: Begin SetCementW(Cem, 52.0, 24.0, 6.0, 11.0, 3.2, 0.5, 0.4, 0.6, 0.5); SetPSD(Cem.P, 18, 1.00, 1, 60); End;
  3: Begin SetCementW(Cem, 63.0, 12.0, 10.0, 7.0, 3.6, 0.7, 0.5, 0.8, 0.6); SetPSD(Cem.P, 11, 1.00, 1, 45); End;
  4: Begin SetCementW(Cem, 28.0, 50.0, 4.0, 12.0, 2.8, 0.4, 0.3, 0.4, 0.4); SetPSD(Cem.P, 19, 1.00, 1, 60); End;
  5: Begin SetCementW(Cem, 48.0, 32.0, 3.5, 13.0, 2.8, 0.4, 0.4, 0.4, 0.4); SetPSD(Cem.P, 17, 1.05, 1, 60); End;
  6: Begin SetCementW(Cem, 67.0, 20.0, 9.0, 1.0, 3.2, 0.5, 0.5, 0.3, 0.3); SetPSD(Cem.P, 15, 1.00, 1, 50); End;
 End;
 S := 0;
 For I := 1 To 9 Do S := S + Cem.W[I];
 For I := 1 To 9 Do Cem.W[I] := Cem.W[I] * 100 / S;
End;

Procedure LoadFAPreset(K: Integer);
Begin
 If K = 1 Then Begin
  FAm.Name := 'Cendre volante classe F'; FAm.Rho := 2.35;
  SetPSD(FAm.P, 14, 1.00, 1, 60); FAm.React := 1.0; FAm.Extra := 70;
 End
 Else Begin
  FAm.Name := 'Cendre volante classe C'; FAm.Rho := 2.65;
  SetPSD(FAm.P, 12, 1.00, 1, 60); FAm.React := 2.0; FAm.Extra := 75;
 End;
End;

Procedure LoadSlagPreset(K: Integer);
Begin
 SlgM.Name := 'Laitier granule'; SlgM.Rho := 2.90;
 SetPSD(SlgM.P, 11, 1.10, 1, 50); SlgM.React := 1.0; SlgM.Extra := 0;
 If K = 2 Then Begin SlgM.Name := 'Laitier broye fin'; SetPSD(SlgM.P, 7, 1.10, 1, 40); SlgM.React := 1.4; End;
End;

Procedure LoadFillerPreset(K: Integer);
Begin
 FilM.Name := 'Calcaire broye'; FilM.Rho := 2.71;
 SetPSD(FilM.P, 10, 1.00, 1, 50); FilM.React := 0; FilM.Extra := 3;
 If K = 2 Then Begin FilM.Name := 'Quartz broye (inerte)'; FilM.Rho := 2.65; SetPSD(FilM.P, 8, 1.00, 1, 40); FilM.Extra := 1; End;
End;

Procedure LoadSFPreset;
Begin
 SFM.Name := 'Fumee de silice'; SFM.Rho := 2.20;
 SetPSD(SFM.P, 1.5, 1.00, 1, 6); SFM.React := 1.0; SFM.Extra := 0;
End;

Procedure LoadAggDefaults;
Begin
 Sand.Rho := 2.65; Sand.DMin := 0.15; Sand.DMax := 4.75; Sand.Q := 0.5;
 Coarse.Rho := 2.70; Coarse.DMin := 4.75; Coarse.DMax := 19.0; Coarse.Q := 0.5;
End;

{ --------------------------------------------------------------------- }
{ Granulometrie Rosin-Rammler                                           }
{ --------------------------------------------------------------------- }

Function RRG(Const P: TPSD; D: Real): Real;
Begin
 If D <= 0 Then RRG := 0
 Else RRG := 1 - Exp(-Pw(D / P.D63, P.Nr));
End;

{ fraction massique cumulee tronquee a [DMin, DMaxEff] }
Function PSDCum(Const P: TPSD; D, DMaxEff: Real): Real;
Var G0, G1: Real;
Begin
 G0 := RRG(P, P.DMin); G1 := RRG(P, DMaxEff);
 If D <= P.DMin Then Begin PSDCum := 0; Exit; End;
 If D >= DMaxEff Then Begin PSDCum := 1; Exit; End;
 If G1 - G0 < 1E-9 Then Begin PSDCum := 1; Exit; End;
 PSDCum := (RRG(P, D) - G0) / (G1 - G0);
End;

{ masse volumique du ciment (g/cm3) }
Function CemRho: Real;
Var I: Integer; S, T: Real;
Begin
 S := 0; T := 0;
 For I := 1 To 9 Do Begin
  S := S + Cem.W[I];
  T := T + Cem.W[I] / PhRho[I];
 End;
 If (S <= 0) Or (T <= 0) Then CemRho := 3.15 Else CemRho := S / T;
End;

{ surface specifique geometrique (m2/kg) d'une granulometrie }
Function SpecSurface(Const P: TPSD; Rho: Real): Real;
Var I: Integer; D0, D1, DM, F, S: Real;
Begin
 S := 0;
 For I := 0 To 79 Do Begin
  D0 := P.DMin * Exp(Ln(P.DMax / P.DMin) * I / 80);
  D1 := P.DMin * Exp(Ln(P.DMax / P.DMin) * (I + 1) / 80);
  DM := Sqrt(D0 * D1);
  F := PSDCum(P, D1, P.DMax) - PSDCum(P, D0, P.DMax);
  S := S + F * 6000.0 / (Rho * DM);
 End;
 SpecSurface := S;
End;

Procedure GetPSD(M: Integer; Var P: TPSD);
Begin
 Case M Of
  0: P := Cem.P;
  1: P := FAm.P;
  2: P := SlgM.P;
  3: P := FilM.P;
 Else P := SFM.P;
 End;
End;

Function MatRho(M: Integer): Real;
Begin
 Case M Of
  0: MatRho := CemRho;
  1: MatRho := FAm.Rho;
  2: MatRho := SlgM.Rho;
  3: MatRho := FilM.Rho;
 Else MatRho := SFM.Rho;
 End;
End;

Function MatName(M: Integer): String;
Begin
 Case M Of
  0: MatName := 'Ciment';
  1: MatName := 'Cendre volante';
  2: MatName := 'Laitier';
  3: MatName := 'Calcaire';
 Else MatName := 'Fumee de silice';
 End;
End;

{ --------------------------------------------------------------------- }
{ Formulation                                                           }
{ --------------------------------------------------------------------- }

Procedure SyncMatDensities;
Begin
 PhRho[P_FA] := FAm.Rho;
 PhRho[P_SLAG] := SlgM.Rho;
 PhRho[P_FILL] := FilM.Rho;
 PhRho[P_SF] := SFM.Rho;
End;

{ volumes (L) par kg de liant : ciment, cendre, laitier, calcaire, fumee, eau }
Procedure MixVolumes(Var V: TR6; Var Mas: TR6);
Var I: Integer; S: Real;
Begin
 SyncMatDensities;
 S := 0;
 For I := 0 To 4 Do S := S + Frac[I];
 If S <= 0 Then Begin Frac[0] := 100; S := 100; End;
 For I := 0 To 4 Do Begin
  Mas[I] := Frac[I] / S;
  V[I] := Mas[I] / MatRho(I);
 End;
 Mas[5] := WB; V[5] := WB;
End;

{ fractions volumiques de la pate (ciment, cendre, laitier, calcaire, fumee, eau) }
Procedure PasteFractions(Var F: TR6);
Var V, M: TR6; I: Integer; T: Real;
Begin
 MixVolumes(V, M);
 T := 0;
 For I := 0 To 5 Do T := T + V[I];
 For I := 0 To 5 Do F[I] := V[I] / T;
End;

{ fractions volumiques de chaque phase du ciment (somme = 1) }
Procedure CemVolFracs(Var F: TR10);
Var I: Integer; T: Real;
Begin
 T := 0;
 For I := 1 To 9 Do T := T + Cem.W[I] / PhRho[I];
 If T <= 0 Then Begin F[1] := 1; Exit; End;
 For I := 1 To 9 Do F[I] := Cem.W[I] / PhRho[I] / T;
End;

Procedure DefaultAll;
Var I: Integer;
Begin
 InitPhases;
 LoadCementPreset(1); LoadFAPreset(1); LoadSlagPreset(1);
 LoadFillerPreset(1); LoadSFPreset; LoadAggDefaults;
 MixKind := 0; WB := 0.40;
 Frac[0] := 100; Frac[1] := 0; Frac[2] := 0; Frac[3] := 0; Frac[4] := 0;
 NVox := 50; ResUm := 1.0; Seed := 12345; WallOn := False;
 ConcBind := 350; ConcSand := 750; ConcCoarse := 1050; ConcAir := 2;
 CureTemp := 25; CureMode := 0; Beta := 0.0007; EaCem := 40; TRefC := 25;
 EndDays := 28; MaxCycles := 6000; PNucCH := 0.02; MaxAge := 50; DormCyc := 90;
 AccNear := 0.80; AccFar := 0.80;
 For I := 0 To NPH - 1 Do Rate[I] := 0;
 Rate[P_C3S] := 0.0072; Rate[P_C2S] := 0.0012; Rate[P_C3A] := 0.0120;
 Rate[P_C4AF] := 0.0040; Rate[P_LIME] := 0.0320; Rate[P_SF] := 0.0060;
 Rate[P_FA] := 0.0012; Rate[P_SLAG] := 0.0016;
 Fc0 := 130; FcExp := 3.0; SigmaPore := 5.0;
 ProjName := ''; Modified := False;
End;

{ --------------------------------------------------------------------- }
{ Reseau 3D de voxels                                                   }
{ --------------------------------------------------------------------- }

Procedure FreeGridG(Var G: TGrid; Var GN: Integer);
Var Z: Integer;
Begin
 If GN > 0 Then
  For Z := 0 To GN - 1 Do
   If G[Z] <> Nil Then Begin
    FreeMem(G[Z], GN * GN);
    G[Z] := Nil;
   End;
 GN := 0;
End;

Function AllocGridG(Var G: TGrid; Var GN: Integer; NewN: Integer): Boolean;
Var Z: Integer;
Begin
 AllocGridG := True;
 If GN = NewN Then Exit;
 FreeGridG(G, GN);
 For Z := 0 To MAXN - 1 Do G[Z] := Nil;
 For Z := 0 To NewN - 1 Do Begin
  {$IFNDEF FPC}
  If MaxAvail < LongInt(NewN) * NewN + 1024 Then Begin
   GN := Z;
   FreeGridG(G, GN);
   AllocGridG := False;
   Exit;
  End;
  {$ENDIF}
  GetMem(G[Z], NewN * NewN);
 End;
 GN := NewN;
End;

Procedure ClearGridG(Var G: TGrid; GN: Integer; V: Byte);
Var Z: Integer;
Begin
 For Z := 0 To GN - 1 Do FillChar(G[Z]^, GN * GN, V);
End;

Procedure InitWrap;
Var I: Integer;
Begin
 For I := -24 To N + 24 Do Wr[I] := ((I Mod N) + N) Mod N;
End;

Function GV(X, Y, Z: Integer): Byte;
Begin
 GV := Grid[Wr[Z]]^[Wr[Y] * N + Wr[X]];
End;

Procedure SV(X, Y, Z: Integer; P: Byte);
Var Q: PSlice; Idx: Integer; Old: Byte;
Begin
 Q := Grid[Wr[Z]];
 Idx := Wr[Y] * N + Wr[X];
 Old := Q^[Idx];
 Dec(PhCnt[Old]); Inc(PhCnt[P]);
 Q^[Idx] := P;
End;

Procedure CountPhases;
Var Z, I: Integer; Q: PSlice;
Begin
 FillChar(PhCnt, SizeOf(PhCnt), 0);
 For Z := 0 To N - 1 Do Begin
  Q := Grid[Z];
  For I := 0 To NN - 1 Do Inc(PhCnt[Q^[I]]);
 End;
End;

{ nombre de voxels de la pate (hors granulat) }
Function PasteVox: LongInt;
Begin
 PasteVox := NTot - PhCnt[P_AGG];
End;

Procedure InitSphVol;
Var D, I, J, K, OI, OJ, OK: Integer; C: LongInt;
Begin
 For D := 1 To MAXD Do Begin
  C := 0;
  For I := 0 To D - 1 Do Begin
   OI := 2 * I - (D - 1);
   For J := 0 To D - 1 Do Begin
    OJ := 2 * J - (D - 1);
    For K := 0 To D - 1 Do Begin
     OK := 2 * K - (D - 1);
     If LongInt(OI) * OI + LongInt(OJ) * OJ + LongInt(OK) * OK <= LongInt(D) * D Then Inc(C);
    End;
   End;
  End;
  SphVol[D] := C;
 End;
End;

{ --------------------------------------------------------------------- }
{ Generation de la microstructure initiale                              }
{ --------------------------------------------------------------------- }

Var
 CumCem  : Array[1..9] Of Real;
 SeedX, SeedY, SeedZ : Array[0..63] Of Integer;
 SeedPh  : Array[0..63] Of Byte;

Function DrawCemPhase: Byte;
Var R: Real; I: Integer;
Begin
 R := RndF;
 For I := 1 To 9 Do
  If R < CumCem[I] Then Begin DrawCemPhase := I; Exit; End;
 DrawCemPhase := P_C3S;
End;

Function SphereFree(X0, Y0, Z0, D: Integer): Boolean;
Var I, J, K, OI, OJ, OK, XX, YY, ZZ, D2: Integer; S: LongInt;
Begin
 SphereFree := False;
 D2 := D;
 For I := 0 To D - 1 Do Begin
  OI := 2 * I - (D - 1);
  XX := (X0 + I) Mod N;
  For J := 0 To D - 1 Do Begin
   OJ := 2 * J - (D - 1);
   YY := (Y0 + J) Mod N;
   For K := 0 To D - 1 Do Begin
    OK := 2 * K - (D - 1);
    S := LongInt(OI) * OI + LongInt(OJ) * OJ + LongInt(OK) * OK;
    If S <= LongInt(D2) * D2 Then Begin
     If WallOn Then ZZ := Z0 + K Else ZZ := (Z0 + K) Mod N;
     If Grid[ZZ]^[YY * N + XX] <> P_POR Then Exit;
    End;
   End;
  End;
 End;
 SphereFree := True;
End;

Procedure SpherePlace(X0, Y0, Z0, D, Mat: Integer);
Var
 I, J, K, OI, OJ, OK, XX, YY, ZZ, NS, S, Best : Integer;
 DX, DY, DZ, DD, BD : LongInt;
 Ph, Single : Byte;
 Q : PSlice;
 Idx : Integer;
Begin
 NS := 1;
 Single := P_C3S;
 Case Mat Of
  0: Begin
   If D >= 6 Then NS := Round(Pw(D / 6.0, 3));
   If NS < 1 Then NS := 1;
   If NS > 60 Then NS := 60;
   For S := 0 To NS - 1 Do Begin
    SeedX[S] := RndN(D); SeedY[S] := RndN(D); SeedZ[S] := RndN(D);
    SeedPh[S] := DrawCemPhase;
   End;
   Single := SeedPh[0];
  End;
  1: Begin
   If RndF * 100 < FAm.Extra Then Single := P_FA Else Begin Single := P_FILL; Inc(Sim.FAInert, SphVol[D]); End;
  End;
  2: Single := P_SLAG;
  3: Single := P_FILL;
 Else Single := P_SF;
 End;
 For I := 0 To D - 1 Do Begin
  OI := 2 * I - (D - 1);
  XX := (X0 + I) Mod N;
  For J := 0 To D - 1 Do Begin
   OJ := 2 * J - (D - 1);
   YY := (Y0 + J) Mod N;
   For K := 0 To D - 1 Do Begin
    OK := 2 * K - (D - 1);
    If LongInt(OI) * OI + LongInt(OJ) * OJ + LongInt(OK) * OK <= LongInt(D) * D Then Begin
     If WallOn Then ZZ := Z0 + K Else ZZ := (Z0 + K) Mod N;
     Ph := Single;
     If NS > 1 Then Begin
      BD := MaxLongInt; Best := 0;
      For S := 0 To NS - 1 Do Begin
       DX := I - SeedX[S]; DY := J - SeedY[S]; DZ := K - SeedZ[S];
       DD := DX * DX + DY * DY + DZ * DZ;
       If DD < BD Then Begin BD := DD; Best := S; End;
      End;
      Ph := SeedPh[Best];
     End;
     Q := Grid[ZZ];
     Idx := YY * N + XX;
     Dec(PhCnt[Q^[Idx]]);
     Inc(PhCnt[Ph]);
     Q^[Idx] := Ph;
    End;
   End;
  End;
 End;
End;

Procedure BuildSulfList;
Var Z, I: Integer; Q: PSlice; Ph: Byte;
Begin
 If SulfList = Nil Then GetMem(SulfList, SizeOf(TLongArr));
 SulfN := 0; Sim.SulfRes := 0;
 For Z := 0 To N - 1 Do Begin
  Q := Grid[Z];
  For I := 0 To NN - 1 Do Begin
   Ph := Q^[I];
   If (Ph >= P_GYP) And (Ph <= P_ANH) Or (Ph = P_ALK) Then
    If SulfN < MAXLIST Then Begin
     SulfList^[SulfN] := LongInt(Z) * NN + I;
     Inc(SulfN);
     Sim.SulfRes := Sim.SulfRes + MolVox[Ph];
    End;
  End;
 End;
End;

Procedure InitReactions; Forward;

{ Reinitialise l'etat de la simulation a partir de la microstructure }
Procedure ResetSimState;
Var I, K: Integer; Cl, Bd, Cm, Sol: Real; FAI: LongInt;
Begin
 FAI := Sim.FAInert;
 FillChar(Sim, SizeOf(Sim), 0);
 Sim.FAInert := FAI;
 Sim.CurN := N;
 Sim.Temp := CureTemp;
 Sim.TSet := -1;
 InitReactions;
 CountPhases;
 Cl := 0; Cm := 0; Bd := 0; Sol := 0;
 For I := 1 To 4 Do Cl := Cl + PhRho[I] * PhCnt[I];
 Cl := Cl + PhRho[P_LIME] * PhCnt[P_LIME];
 For I := 1 To 9 Do Cm := Cm + PhRho[I] * PhCnt[I];
 Bd := Cm + PhRho[P_FA] * PhCnt[P_FA] + FAm.Rho * Sim.FAInert +
       PhRho[P_SLAG] * PhCnt[P_SLAG] + PhRho[P_SF] * PhCnt[P_SF];
 For I := 1 To NPH - 1 Do
  If (I <> P_VOID) Then Sol := Sol + PhRho[I] * PhCnt[I];
 Sim.M0Cl := Cl; Sim.M0Cem := Cm; Sim.M0Bind := Bd; Sim.MSol0 := Sol;
 Sim.Dummy[0] := Cl + PhRho[P_FA] * PhCnt[P_FA] + PhRho[P_SF] * PhCnt[P_SF] +
                PhRho[P_SLAG] * PhCnt[P_SLAG];
 If Sim.Dummy[0] <= 0 Then Sim.Dummy[0] := 1;
 Sim.MW0 := PhCnt[P_POR];
 CcapJ := 0.75 * Sol + 4.18 * Sim.MW0;
 If Sim.M0Bind <= 0 Then Sim.M0Bind := 1;
 If Sim.M0Cl <= 0 Then Sim.M0Cl := 1;
 BuildSulfList;
 For K := 0 To NSER - 1 Do
  For I := 0 To MAXPT - 1 Do Ser[K, I] := 0;
End;

Function GenerateMicro: Boolean;
Var
 Fr : TR6;
 Tgt, Def : Array[0..4] Of LongInt;
 NCnt : Array[0..4, 1..MAXD] Of LongInt;
 CemF : TR10;
 M, D, K, Tries, MaxTry, X0, Y0, Z0, DCap, I, Placed : Integer;
 NPaste, TotalWork, Done : LongInt;
 P : TPSD;
 DMaxEff, Lo, Up, F, V, Carry, S : Real;
Begin
 GenerateMicro := False;
 SyncMatDensities;
 InitSphVol;
 If NVox < 16 Then NVox := 16;
 If NVox > MAXN Then NVox := MAXN;
 N := NVox; NN := N * N; NTot := LongInt(N) * NN;
 If Not AllocGridG(Grid, GridN, N) Then Begin
  StatusMsg := 'Memoire insuffisante pour la microstructure';
  Exit;
 End;
 HaveGrid := False;
 ClearGridG(Grid, N, P_POR);
 InitWrap;
 SetSeed(Seed);
 FillChar(PhCnt, SizeOf(PhCnt), 0);
 PhCnt[P_POR] := NTot;
 FillChar(Sim, SizeOf(Sim), 0);
 ZLo := 0; ZHi := N - 1;
 If WallOn Then Begin
  FillChar(Grid[0]^, NN, P_AGG);
  FillChar(Grid[N - 1]^, NN, P_AGG);
  Dec(PhCnt[P_POR], 2 * NN); Inc(PhCnt[P_AGG], 2 * NN);
  ZLo := 1; ZHi := N - 2;
 End;
 NPaste := PasteVox;
 PasteFractions(Fr);
 For M := 0 To 4 Do Tgt[M] := Round(Fr[M] * NPaste);
 CemVolFracs(CemF);
 S := 0;
 For I := 1 To 9 Do Begin S := S + CemF[I]; CumCem[I] := S; End;
 DCap := N Div 2;
 If DCap > MAXD Then DCap := MAXD;
 If WallOn And (DCap > N - 3) Then DCap := N - 3;
 { histogrammes de tailles }
 For M := 0 To 4 Do Begin
  Def[M] := 0;
  For D := 1 To MAXD Do NCnt[M, D] := 0;
  If Tgt[M] > 0 Then Begin
   GetPSD(M, P);
   DMaxEff := RMin(P.DMax, (DCap + 0.5) * ResUm);
   If DMaxEff < P.DMin + 0.5 Then DMaxEff := P.DMin + 0.5;
   Carry := 0;
   For D := DCap Downto 1 Do Begin
    Lo := (D - 0.5) * ResUm; Up := (D + 0.5) * ResUm;
    If D = 1 Then Lo := 0;
    If D = DCap Then Up := DMaxEff + 1;
    F := PSDCum(P, Up, DMaxEff) - PSDCum(P, Lo, DMaxEff);
    V := F * Tgt[M] + Carry;
    K := Trunc(V / SphVol[D]);
    NCnt[M, D] := K;
    Carry := V - K * SphVol[D];
   End;
   Def[M] := Round(Carry);
  End;
 End;
 { placement du plus gros au plus petit }
 TotalWork := 0;
 For M := 0 To 4 Do For D := 1 To DCap Do TotalWork := TotalWork + NCnt[M, D];
 If TotalWork < 1 Then TotalWork := 1;
 Done := 0;
 For D := DCap Downto 1 Do Begin
  For M := 0 To 4 Do Begin
   For K := 1 To NCnt[M, D] Do Begin
    If D <= 1 Then MaxTry := 400 Else If D = 2 Then MaxTry := 250 Else MaxTry := 120;
    Placed := 0;
    For Tries := 1 To MaxTry Do Begin
     X0 := RndN(N); Y0 := RndN(N);
     If WallOn Then Begin
      If N - 1 - D < 1 Then Break;
      Z0 := 1 + RndN(N - 1 - D);
     End
     Else Z0 := RndN(N);
     If SphereFree(X0, Y0, Z0, D) Then Begin
      SpherePlace(X0, Y0, Z0, D, M);
      Placed := 1;
      Break;
     End;
    End;
    If Placed = 0 Then Def[M] := Def[M] + SphVol[D];
    Inc(Done);
    If (Done And 255) = 0 Then
     ShowProgress('Generation de la microstructure...', Round(80.0 * Done / TotalWork));
   End;
  End;
 End;
 { complement : voxels isoles }
 For M := 0 To 4 Do Begin
  K := 0;
  While (Def[M] > 0) And (K < 30000) Do Begin
   X0 := RndN(N); Y0 := RndN(N);
   If WallOn Then Z0 := 1 + RndN(N - 2) Else Z0 := RndN(N);
   Inc(K);
   If Grid[Z0]^[Y0 * N + X0] = P_POR Then Begin
    SpherePlace(X0, Y0, Z0, 1, M);
    Dec(Def[M]);
   End;
  End;
 End;
 CountPhases;
 ResetSimState;
 HaveGrid := True;
 GenInfo := 'N=' + IStr(N) + '  graine=' + IStr(Seed);
 ShowProgress('Microstructure generee', 100);
 GenerateMicro := True;
End;

{ --------------------------------------------------------------------- }
{ Percolation (algorithme du feu)                                       }
{ --------------------------------------------------------------------- }

Const
 QBS = 8192;

Type
 TQBlk = Array[0..QBS-1] Of LongInt;
 PQBlk = ^TQBlk;

Var
 QB      : Array[0..127] Of PQBlk;
 QHead, QTail : LongInt;
 VisG    : TGrid;
 VisN    : Integer;
 PMask   : Array[0..NPH-1] Of Byte;   { 0 non, 1 oui, 2 oui si lie a un hydrate }

Procedure QPush(V: LongInt);
Var B: Integer;
Begin
 B := QTail Div QBS;
 If QB[B] = Nil Then GetMem(QB[B], SizeOf(TQBlk));
 QB[B]^[QTail Mod QBS] := V;
 Inc(QTail);
End;

Function QPop: LongInt;
Begin
 QPop := QB[QHead Div QBS]^[QHead Mod QBS];
 Inc(QHead);
End;

Procedure QFree;
Var I: Integer;
Begin
 For I := 0 To 127 Do
  If QB[I] <> Nil Then Begin FreeMem(QB[I], SizeOf(TQBlk)); QB[I] := Nil; End;
 QHead := 0; QTail := 0;
End;

{ solides lies par les hydrates (grains non hydrates comptes s'ils touchent un produit) }
Procedure MaskSolids;
Var I: Integer;
Begin
 For I := 0 To NPH - 1 Do Begin
  If (I = P_POR) Or (I = P_VOID) Or (I = P_AGG) Then PMask[I] := 0
  Else If (I >= P_CSH) And (I <= P_PCSH) Then PMask[I] := 1
  Else PMask[I] := 2;
 End;
End;

Procedure MaskPores;
Var I: Integer;
Begin
 For I := 0 To NPH - 1 Do
  If (I = P_POR) Or (I = P_VOID) Then PMask[I] := 1 Else PMask[I] := 0;
End;

Const
 DX6 : Array[0..5] Of Integer = (1, -1, 0, 0, 0, 0);
 DY6 : Array[0..5] Of Integer = (0, 0, 1, -1, 0, 0);
 DZ6 : Array[0..5] Of Integer = (0, 0, 0, 0, 1, -1);

Procedure Burn(ZF: Integer; Bit: Byte);
Var Idx, R: LongInt; X, Y, Z : Integer;

 Function Bonded(NX, NY, NZ: Integer): Boolean;
 Var I: Integer; Ph: Byte;
 Begin
  Bonded := True;
  For I := 0 To 5 Do Begin
   Ph := Grid[Wr[NZ + DZ6[I]]]^[Wr[NY + DY6[I]] * N + Wr[NX + DX6[I]]];
   If (Ph >= P_CSH) And (Ph <= P_PCSH) Then Exit;
  End;
  Bonded := False;
 End;

 Procedure Visit(NX, NY, NZ: Integer);
 Var K: Integer; M: Byte;
 Begin
  K := NY * N + NX;
  M := PMask[Grid[NZ]^[K]];
  If (M <> 0) And ((VisG[NZ]^[K] And Bit) = 0) Then Begin
   If (M = 2) And Not Bonded(NX, NY, NZ) Then Exit;
   VisG[NZ]^[K] := VisG[NZ]^[K] Or Bit;
   QPush(LongInt(NZ) * NN + K);
  End;
 End;

Begin
 QHead := 0; QTail := 0;
 For Y := 0 To N - 1 Do
  For X := 0 To N - 1 Do Visit(X, Y, ZF);
 While QHead < QTail Do Begin
  Idx := QPop;
  Z := Idx Div NN; R := Idx - LongInt(Z) * NN; Y := R Div N; X := R - Y * N;
  If X < N - 1 Then Visit(X + 1, Y, Z) Else Visit(0, Y, Z);
  If X > 0 Then Visit(X - 1, Y, Z) Else Visit(N - 1, Y, Z);
  If Y < N - 1 Then Visit(X, Y + 1, Z) Else Visit(X, 0, Z);
  If Y > 0 Then Visit(X, Y - 1, Z) Else Visit(X, N - 1, Z);
  If Z < ZHi Then Visit(X, Y, Z + 1);
  If Z > ZLo Then Visit(X, Y, Z - 1);
 End;
End;

{ Percolation selon z du masque PMask. Full : compte les voxels de l'amas
  traversant ; Vis conserve les marques (bits 1 et 2). }
Function Percolate(Full: Boolean; Var Cnt: LongInt): Boolean;
Var Z, I: Integer; Res: Boolean;
Begin
 Cnt := 0;
 If Not AllocGridG(VisG, VisN, N) Then Begin Percolate := False; Exit; End;
 ClearGridG(VisG, N, 0);
 Burn(ZLo, 1);
 Res := False;
 For I := 0 To NN - 1 Do
  If (VisG[ZHi]^[I] And 1) <> 0 Then Begin Res := True; Break; End;
 If Res And Full Then Begin
  Burn(ZHi, 2);
  For Z := ZLo To ZHi Do
   For I := 0 To NN - 1 Do
    If (VisG[Z]^[I] And 3) = 3 Then Inc(Cnt);
 End;
 QFree;
 Percolate := Res;
End;

{ --------------------------------------------------------------------- }
{ Reactions chimiques                                                   }
{ --------------------------------------------------------------------- }

Procedure InitReactions;
Const VW = 18.015;
Var VCH, VCSH, VCSH11, VAFT, VAFM, VHG, VHGF, VGYP : Real;

 Procedure SetRx(I: Integer; V0, HSpec, RhoR: Real; P1: Byte; Mol1, VP1: Real;
                 P2: Byte; Mol2, VP2: Real; WatMol, CHMol, SulMol: Real);
 Var Sc, Tot: Real;
 Begin
  Sc := 1 / V0;
  Rx[I].NP := 1;
  Rx[I].Pr[0] := P1; Rx[I].Cnt[0] := Mol1 * VP1 * Sc;
  Rx[I].Pr[1] := P2; Rx[I].Cnt[1] := 0;
  Tot := Rx[I].Cnt[0];
  If Mol2 > 0 Then Begin
   Rx[I].NP := 2; Rx[I].Cnt[1] := Mol2 * VP2 * Sc; Tot := Tot + Rx[I].Cnt[1];
  End;
  Rx[I].Wat := WatMol * VW * Sc;
  Rx[I].CHc := CHMol * VCH * Sc;
  Rx[I].AFc := 0;
  Rx[I].Sul := SulMol * Sc;
  Rx[I].SulV := SulMol * VGYP * Sc;
  Rx[I].E := 1 + Rx[I].CHc + Rx[I].SulV + Rx[I].Wat - Tot;
  If Rx[I].E < 0 Then Rx[I].E := 0;
  Rx[I].Heat := RhoR * HSpec;
 End;

Begin
 SyncMatDensities;
 VCH := 74.09 / PhRho[P_CH];
 VCSH := 227.48 / PhRho[P_CSH];
 VCSH11 := 193.8 / PhRho[P_PCSH];
 VAFT := 1255.1 / PhRho[P_AFT];
 VAFM := 622.5 / PhRho[P_AFM];
 VHG := 378.3 / PhRho[P_HG];
 VHGF := 407.15 / PhRho[P_HG];
 VGYP := 172.17 / PhRho[P_GYP];
 { C3S + 5.3 H -> C-S-H + 1.3 CH }
 SetRx(RX_C3S, 228.32 / PhRho[P_C3S], 517, PhRho[P_C3S], P_CSH, 1, VCSH, P_CH, 1.3, VCH, 5.3, 0, 0);
 { C2S + 4.3 H -> C-S-H + 0.3 CH }
 SetRx(RX_C2S, 172.24 / PhRho[P_C2S], 262, PhRho[P_C2S], P_CSH, 1, VCSH, P_CH, 0.3, VCH, 4.3, 0, 0);
 { C3A + 3 gypse + 26 H -> ettringite }
 SetRx(RX_AFT, 270.19 / PhRho[P_C3A], 1672, PhRho[P_C3A], P_AFT, 1, VAFT, 0, 0, 0, 26, 0, 3);
 { C3A + gypse + 10 H -> monosulfate }
 SetRx(RX_AFM, 270.19 / PhRho[P_C3A], 1144, PhRho[P_C3A], P_AFM, 1, VAFM, 0, 0, 0, 10, 0, 1);
 { C3A + 6 H -> hydrogrenat }
 SetRx(RX_HG, 270.19 / PhRho[P_C3A], 908, PhRho[P_C3A], P_HG, 1, VHG, 0, 0, 0, 6, 0, 0);
 { C4AF + 2 CH + 10 H -> 2 hydrogrenat (Fe) }
 SetRx(RX_C4AF, 485.96 / PhRho[P_C4AF], 418, PhRho[P_C4AF], P_HG, 2, VHGF, 0, 0, 0, 10, 2, 0);
 { 2 C3A + AFt + 4 H -> 3 monosulfate (conversion de l'ettringite) }
 SetRx(RX_CONV, 270.19 / PhRho[P_C3A], 1144, PhRho[P_C3A], P_AFM, 1.5, VAFM, 0, 0, 0, 2, 0, 0);
 Rx[RX_CONV].AFc := 0.5 * VAFT * (1 / (270.19 / PhRho[P_C3A]));
 Rx[RX_CONV].E := 1 + Rx[RX_CONV].AFc + Rx[RX_CONV].Wat - Rx[RX_CONV].Cnt[0];
 If Rx[RX_CONV].E < 0 Then Rx[RX_CONV].E := 0;
 { CaO + H -> CH }
 SetRx(RX_LIME, 56.08 / PhRho[P_LIME], 1166, PhRho[P_LIME], P_CH, 1, VCH, 0, 0, 0, 1, 0, 0);
 { SiO2 + 1.1 CH + 2.9 H -> C1.1SH4 }
 SetRx(RX_SF, 60.08 / PhRho[P_SF], 780, PhRho[P_SF], P_PCSH, 1, VCSH11, 0, 0, 0, 2.9, 1.1, 0);
 { verre alumino-silicate (M=100) + CH + 2 H -> gel pouzzolanique }
 SetRx(RX_FA, 100.0 / PhRho[P_FA], 570, PhRho[P_FA], P_PCSH, 1, 210.12 / PhRho[P_PCSH], 0, 0, 0, 2, 1, 0);
 { laitier (M=100) + 1.3 H -> gel }
 SetRx(RX_SLAG, 100.0 / PhRho[P_SLAG], 461, PhRho[P_SLAG], P_PCSH, 1, 123.4 / PhRho[P_PCSH], 0, 0, 0, 1.3, 0, 0);
 MolVox[P_GYP] := PhRho[P_GYP] / 172.17;
 MolVox[P_HEM] := PhRho[P_HEM] / 145.15;
 MolVox[P_ANH] := PhRho[P_ANH] / 136.14;
 MolVox[P_ALK] := PhRho[P_ALK] / 174.26;
End;

{ --------------------------------------------------------------------- }
{ Moteur d'hydratation                                                  }
{ --------------------------------------------------------------------- }

Const
 DXs : Array[0..5] Of Integer = (1, -1, 0, 0, 0, 0);
 DYs : Array[0..5] Of Integer = (0, 0, 1, -1, 0, 0);
 DZs : Array[0..5] Of Integer = (0, 0, 0, 0, 1, -1);

Function IsSol(Ph: Byte): Boolean;
Begin
 IsSol := (Ph <> P_POR) And (Ph <> P_VOID);
End;

{ pore (eau) au hasard dans le voisinage }
Function FindPore(X, Y, Z: Integer; Var PX, PY, PZ: Integer): Boolean;
Const Rad : Array[0..4] Of Integer = (1, 2, 3, 5, 8);
Var R, T, DX, DY, DZ: Integer;
Begin
 FindPore := True;
 For R := 0 To 4 Do
  For T := 1 To 6 + 3 * R Do Begin
   DX := RndN(2 * Rad[R] + 1) - Rad[R];
   DY := RndN(2 * Rad[R] + 1) - Rad[R];
   DZ := RndN(2 * Rad[R] + 1) - Rad[R];
   If GV(X + DX, Y + DY, Z + DZ) = P_POR Then Begin
    PX := X + DX; PY := Y + DY; PZ := Z + DZ;
    Exit;
   End;
  End;
 FindPore := False;
End;

Function Precip(X, Y, Z: Integer; Pr: Byte): Boolean;
Var I, Same, Sol, Fil: Integer; Ph: Byte; P: Real;
Begin
 Same := 0; Sol := 0; Fil := 0;
 For I := 0 To 5 Do Begin
  Ph := GV(X + DXs[I], Y + DYs[I], Z + DZs[I]);
  If Ph = Pr Then Inc(Same);
  If IsSol(Ph) Then Inc(Sol);
  If Ph = P_FILL Then Inc(Fil);
 End;
 Case Pr Of
  P_CSH, P_PCSH: Precip := Sol > 0;
  P_CH: Begin
   If Same > 0 Then Precip := True
   Else If Sol > 0 Then Begin
    P := PNucCH * (1 + Fil * (FilM.Extra - 1));
    Precip := RndF < P;
   End
   Else Precip := False;
  End;
  P_AFT: Precip := (Same > 0) Or ((Sol > 0) And (RndF < 0.08));
  P_AFM: Precip := (Same > 0) Or ((Sol > 0) And (RndF < 0.15));
 Else Precip := (Same > 0) Or ((Sol > 0) And (RndF < 0.3));
 End;
End;

{ pore (eau) n'importe ou dans le reseau : dernier recours quand le voisinage est plein }
Function FindPoreAny(Var PX, PY, PZ: Integer): Boolean;
Var T: Integer;
Begin
 FindPoreAny := True;
 If PhCnt[P_POR] > 0 Then
  For T := 1 To 400 Do Begin
   PX := RndN(N); PY := RndN(N); PZ := RndN(N);
   If Grid[PZ]^[PY * N + PX] = P_POR Then Exit;
  End;
 FindPoreAny := False;
End;

Procedure PlaceProduct(X, Y, Z: Integer; Pr: Byte);
Var SX, SY, SZ, Age, D: Integer;
Begin
 If GV(X, Y, Z) = P_POR Then Begin SX := X; SY := Y; SZ := Z; End
 Else If Not FindPore(X, Y, Z, SX, SY, SZ) Then Begin
  If Not FindPoreAny(SX, SY, SZ) Then Begin
   Inc(Sim.Lost);
   Exit;
  End;
 End;
 Age := 0;
 While Age < MaxAge Do Begin
  If Precip(SX, SY, SZ, Pr) Then Break;
  D := RndN(6);
  If GV(SX + DXs[D], SY + DYs[D], SZ + DZs[D]) = P_POR Then Begin
   SX := SX + DXs[D]; SY := SY + DYs[D]; SZ := SZ + DZs[D];
   SX := Wr[SX]; SY := Wr[SY]; SZ := Wr[SZ];
  End;
  Inc(Age);
 End;
 SV(SX, SY, SZ, Pr);
End;

Procedure PlaceVoid(X, Y, Z: Integer);
Var PX, PY, PZ: Integer;
Begin
 If FindPore(X, Y, Z, PX, PY, PZ) Then SV(PX, PY, PZ, P_VOID);
End;

{ sulfate dissous a la demande }
Procedure SulfTake(Need: Real);
Var K, Idx, R: LongInt; X, Y, Z: Integer; Ph: Byte;
Begin
 While (Sim.Pool < Need) And (SulfN > 0) Do Begin
  K := RndN(SulfN);
  Idx := SulfList^[K];
  SulfList^[K] := SulfList^[SulfN - 1];
  Dec(SulfN);
  Z := Idx Div NN; R := Idx - LongInt(Z) * NN; Y := R Div N; X := R - Y * N;
  Ph := Grid[Z]^[R];
  If (Ph >= P_GYP) And (Ph <= P_ANH) Or (Ph = P_ALK) Then Begin
   SV(X, Y, Z, P_POR);
   Sim.Pool := Sim.Pool + MolVox[Ph];
   Sim.SulfRes := Sim.SulfRes - MolVox[Ph];
  End;
 End;
 Sim.Pool := Sim.Pool - Need;
 If Sim.Pool < 0 Then Sim.Pool := 0;
End;

Function DoReact(X, Y, Z, R: Integer): Boolean;
Var
 NC, NA, K, J, NPj, NE, Tries, PX, PY, PZ : Integer;
 CX, CY, CZ : Array[0..3] Of Integer;
 AX, AY, AZ : Array[0..11] Of Integer;
 Found : Boolean;
 Tot : Real;
Begin
 DoReact := False;
 { assez d'eau pour loger les produits ? }
 Tot := Rx[R].Cnt[0] + Rx[R].Cnt[1] - 1 - Rx[R].CHc - Rx[R].AFc;
 If Tot + 3 > PhCnt[P_POR] Then Exit;
 NA := StochRound(Rx[R].AFc);
 If NA > 12 Then NA := 12;
 If NA > 0 Then Begin
  If PhCnt[P_AFT] < NA + 4 Then Exit;
  For K := 0 To NA - 1 Do Begin
   Found := False;
   For Tries := 1 To 800 Do Begin
    PX := RndN(N); PY := RndN(N); PZ := RndN(N);
    If Grid[PZ]^[PY * N + PX] = P_AFT Then Begin
     Found := True;
     For J := 0 To K - 1 Do
      If (AX[J] = PX) And (AY[J] = PY) And (AZ[J] = PZ) Then Found := False;
     If Found Then Break;
    End;
   End;
   If Not Found Then Exit;
   AX[K] := PX; AY[K] := PY; AZ[K] := PZ;
  End;
 End;
 NC := StochRound(Rx[R].CHc);
 If NC > 4 Then NC := 4;
 If NC > 0 Then Begin
  If PhCnt[P_CH] < NC Then Exit;
  For K := 0 To NC - 1 Do Begin
   Found := False;
   For Tries := 1 To 600 Do Begin
    PX := RndN(N); PY := RndN(N); PZ := RndN(N);
    If Grid[PZ]^[PY * N + PX] = P_CH Then Begin
     Found := True;
     For J := 0 To K - 1 Do
      If (CX[J] = PX) And (CY[J] = PY) And (CZ[J] = PZ) Then Found := False;
     If Found Then Break;
    End;
   End;
   If Not Found Then Exit;
   CX[K] := PX; CY[K] := PY; CZ[K] := PZ;
  End;
  For K := 0 To NC - 1 Do SV(CX[K], CY[K], CZ[K], P_POR);
 End;
 For K := 0 To NA - 1 Do SV(AX[K], AY[K], AZ[K], P_POR);
 SV(X, Y, Z, P_POR);
 For J := 0 To Rx[R].NP - 1 Do Begin
  NPj := StochRound(Rx[R].Cnt[J]);
  For K := 1 To NPj Do PlaceProduct(X, Y, Z, Rx[R].Pr[J]);
 End;
 Sim.Shrk := Sim.Shrk + Rx[R].E;
 If CureMode <> 0 Then Begin
  NE := StochRound(Rx[R].E);
  For K := 1 To NE Do PlaceVoid(X, Y, Z);
 End;
 Sim.Heat := Sim.Heat + Rx[R].Heat;
 Sim.Wn := Sim.Wn + Rx[R].Wat;
 DoReact := True;
End;

{ le voxel est-il atteint par l'eau (diffusion a travers les produits) ? }
Function Accessible(X, Y, Z: Integer): Boolean;
Var DX, DY, DZ: Integer;
Begin
 Accessible := True;
 If GV(X + 1, Y, Z) = P_POR Then Exit;
 If GV(X - 1, Y, Z) = P_POR Then Exit;
 If GV(X, Y + 1, Z) = P_POR Then Exit;
 If GV(X, Y - 1, Z) = P_POR Then Exit;
 If GV(X, Y, Z + 1) = P_POR Then Exit;
 If GV(X, Y, Z - 1) = P_POR Then Exit;
 For DX := -1 To 1 Do
  For DY := -1 To 1 Do
   For DZ := -1 To 1 Do
    If GV(X + DX, Y + DY, Z + DZ) = P_POR Then Begin
     Accessible := RndF < AccNear;
     Exit;
    End;
 For DX := -2 To 2 Do
  For DY := -2 To 2 Do
   For DZ := -2 To 2 Do
    If GV(X + DX, Y + DY, Z + DZ) = P_POR Then Begin
     Accessible := RndF < AccFar;
     Exit;
    End;
 For DX := -3 To 3 Do
  For DY := -3 To 3 Do
   For DZ := -3 To 3 Do
    If GV(X + DX, Y + DY, Z + DZ) = P_POR Then Begin
     Accessible := RndF < AccFar * 0.25;
     Exit;
    End;
 Accessible := False;
End;

Procedure TryDissolve(X, Y, Z: Integer; Ph: Byte);
Var Av: Real; Ok: Boolean;
Begin
 If Not Accessible(X, Y, Z) Then Exit;
 Case Ph Of
  P_C3S: Ok := DoReact(X, Y, Z, RX_C3S);
  P_C2S: Ok := DoReact(X, Y, Z, RX_C2S);
  P_C3A: Begin
   Av := Sim.Pool + Sim.SulfRes;
   If Av >= Rx[RX_AFT].Sul Then Begin
    SulfTake(Rx[RX_AFT].Sul); Ok := DoReact(X, Y, Z, RX_AFT);
   End
   Else If PhCnt[P_AFT] >= 12 Then Ok := DoReact(X, Y, Z, RX_CONV)
   Else If Av >= Rx[RX_AFM].Sul Then Begin
    SulfTake(Rx[RX_AFM].Sul); Ok := DoReact(X, Y, Z, RX_AFM);
   End
   Else Ok := DoReact(X, Y, Z, RX_HG);
  End;
  P_C4AF: Ok := DoReact(X, Y, Z, RX_C4AF);
  P_LIME: Ok := DoReact(X, Y, Z, RX_LIME);
  P_SF: Ok := DoReact(X, Y, Z, RX_SF);
  P_FA: Ok := DoReact(X, Y, Z, RX_FA);
  P_SLAG: Ok := DoReact(X, Y, Z, RX_SLAG);
 End;
End;

Procedure ComputeThr;
Var Ph: Integer; P: Real;
Begin
 For Ph := 0 To NPH - 1 Do Thr[Ph] := 0;
 For Ph := 1 To NPH - 1 Do Begin
  P := Rate[Ph];
  If P > 0 Then Begin
   Case Ph Of
    P_C3A: If Sim.Pool + Sim.SulfRes > 0 Then P := P * 0.3;
    P_FA: P := P * FAm.React;
    P_SF: P := P * SFM.React;
    P_SLAG: Begin
     P := P * SlgM.React;
     If PhCnt[P_CH] <= 0 Then P := P * 0.2;
    End;
   End;
   If Sim.Cyc < DormCyc Then P := P * Sqr((Sim.Cyc + 1.0) / DormCyc);
   If P > 0.95 Then P := 0.95;
   Thr[Ph] := Round(P * 32768);
  End;
 End;
End;

Procedure RefillVoids;
Var Z, I: Integer; Q: PSlice;
Begin
 For Z := 0 To N - 1 Do Begin
  Q := Grid[Z];
  For I := 0 To NN - 1 Do
   If Q^[I] = P_VOID Then Begin
    Q^[I] := P_POR; Dec(PhCnt[P_VOID]); Inc(PhCnt[P_POR]);
   End;
 End;
End;

{ --------------------------------------------------------------------- }
{ Grandeurs calculees a partir de la microstructure                     }
{ --------------------------------------------------------------------- }

Function AlphaCement: Real;
Var I: Integer; Cl: Real;
Begin
 Cl := PhRho[P_LIME] * PhCnt[P_LIME];
 For I := 1 To 4 Do Cl := Cl + PhRho[I] * PhCnt[I];
 AlphaCement := 1 - Cl / Sim.M0Cl;
End;

Function AlphaBinder: Real;
Var Cl: Real; I: Integer;
Begin
 Cl := PhRho[P_LIME] * PhCnt[P_LIME] + PhRho[P_FA] * PhCnt[P_FA] +
       PhRho[P_SF] * PhCnt[P_SF] + PhRho[P_SLAG] * PhCnt[P_SLAG];
 For I := 1 To 4 Do Cl := Cl + PhRho[I] * PhCnt[I];
 AlphaBinder := 1 - Cl / Sim.Dummy[0];
End;

Function CapPorosity: Real;
Begin
 CapPorosity := (PhCnt[P_POR] + PhCnt[P_VOID]) / PasteVox;
End;

{ rapport gel / espace de Powers }
Function GelSpace: Real;
Var H, C: Real;
Begin
 H := PhCnt[P_CSH] + PhCnt[P_CH] + PhCnt[P_AFT] + PhCnt[P_AFM] + PhCnt[P_HG] + PhCnt[P_PCSH];
 C := PhCnt[P_POR] + PhCnt[P_VOID];
 If H + C <= 0 Then GelSpace := 0 Else GelSpace := H / (H + C);
End;

Function PowersStrength: Real;
Begin
 PowersStrength := Fc0 * Pw(GelSpace, FcExp);
End;

Procedure Sample;
Var I, K, J: Integer; Cnt: LongInt;
Begin
 I := Sim.NS;
 If I >= MAXPT Then Begin
  For K := 0 To NSER - 1 Do
   For J := 0 To MAXPT Div 2 - 1 Do Ser[K, J] := Ser[K, 2 * J];
  For K := 0 To NPH - 1 Do
   For J := 0 To MAXPT Div 2 - 1 Do PhSer[K, J] := PhSer[K, 2 * J];
  I := MAXPT Div 2;
 End;
 For K := 0 To NPH - 1 Do PhSer[K, I] := PhCnt[K];
 Ser[SER_T, I] := Sim.T;
 Ser[SER_CYC, I] := Sim.Cyc;
 Ser[SER_ALF, I] := AlphaCement;
 Ser[SER_ALT, I] := AlphaBinder;
 Ser[SER_Q, I] := Sim.Heat / Sim.M0Bind;
 Ser[SER_TMP, I] := Sim.Temp;
 Ser[SER_POR, I] := 100 * CapPorosity;
 Ser[SER_WN, I] := Sim.Wn / Sim.M0Bind;
 Ser[SER_SHR, I] := Sim.Shrk / Sim.M0Bind;
 Ser[SER_FC, I] := PowersStrength;
 Sim.NS := I + 1;
 Sim.NextSamp := Sim.Cyc + 1 + Sim.Cyc Div 40;
 If Sim.TSet < 0 Then Begin
  MaskSolids;
  If Percolate(False, Cnt) Then Sim.TSet := Sim.T;
 End;
End;

Procedure RunCycle;
Var
 Z, ZZ, Z0, I, X, Y : Integer;
 Q : PSlice;
 Ph : Byte;
 DTe, KT, DT, HeatBefore : Real;
Begin
 If (CureMode = 0) And (PhCnt[P_VOID] > 0) Then RefillVoids;
 ComputeThr;
 HeatBefore := Sim.Heat;
 Z0 := RndN(N);
 For ZZ := 0 To N - 1 Do Begin
  Z := (Z0 + ZZ) Mod N;
  Q := Grid[Z];
  For I := 0 To NN - 1 Do Begin
   Ph := Q^[I];
   If Thr[Ph] > 0 Then
    If Rnd15 < Thr[Ph] Then Begin
     Y := I Div N; X := I - Y * N;
     TryDissolve(X, Y, Z, Ph);
    End;
  End;
 End;
 DTe := Beta * (2.0 * Sim.Cyc + 1);
 KT := Exp(EaCem * 1000 / 8.314 * (1 / (TRefC + 273.15) - 1 / (Sim.Temp + 273.15)));
 If KT < 0.05 Then KT := 0.05;
 DT := DTe / KT;
 Sim.Te := Sim.Te + DTe;
 Sim.T := Sim.T + DT;
 If CureMode = 2 Then Begin
  Sim.Temp := Sim.Temp + (Sim.Heat - HeatBefore) / CcapJ;
  If Sim.Temp > 100 Then Sim.Temp := 100;
 End
 Else Sim.Temp := CureTemp;
 Inc(Sim.Cyc);
 If Sim.Cyc >= Sim.NextSamp Then Sample;
End;

{ etat initial de la serie (t = 0) }
Procedure SampleInitial;
Begin
 Sim.NS := 0;
 Sim.NextSamp := 0;
 Sample;
 Sim.NextSamp := 1;
End;

{ --------------------------------------------------------------------- }
{ Distribution des tailles de pores                                     }
{ --------------------------------------------------------------------- }

Const
 MAXR = 40;

Var
 AuxG     : TGrid;
 AuxN     : Integer;
 PoreHist : Array[1..MAXR] Of LongInt;
 PoreTot  : LongInt;

Procedure ChamferPass(Forward: Boolean);
Var Z, Y, X, DZ, DY, DX, W, NV, ZA, ZB, ZS: Integer; Cur: Byte; Pos: Boolean;
Begin
 If Forward Then Begin ZA := 0; ZB := N - 1; ZS := 1; End
 Else Begin ZA := N - 1; ZB := 0; ZS := -1; End;
 Z := ZA;
 While True Do Begin
  For Y := 0 To N - 1 Do
   For X := 0 To N - 1 Do Begin
    Cur := VisG[Z]^[Y * N + X];
    If Cur > 0 Then Begin
     For DZ := -1 To 1 Do
      For DY := -1 To 1 Do
       For DX := -1 To 1 Do Begin
        If (DX <> 0) Or (DY <> 0) Or (DZ <> 0) Then Begin
         If Forward Then Pos := (DZ < 0) Or ((DZ = 0) And (DY < 0)) Or ((DZ = 0) And (DY = 0) And (DX < 0))
         Else Pos := (DZ > 0) Or ((DZ = 0) And (DY > 0)) Or ((DZ = 0) And (DY = 0) And (DX > 0));
         If Pos Then Begin
          W := 3 + Ord(DX <> 0) + Ord(DY <> 0) + Ord(DZ <> 0) - 1;
          NV := VisG[Wr[Z + DZ]]^[Wr[Y + DY] * N + Wr[X + DX]] + W;
          If NV < Cur Then Cur := NV;
         End;
        End;
       End;
     VisG[Z]^[Y * N + X] := Cur;
    End;
   End;
  If Z = ZB Then Break;
  Z := Z + ZS;
 End;
End;

{ Remplit PoreHist (rayon en voxels) pour les pores capillaires }
Function PoreSizeDist: Boolean;
Var
 Z, Y, X, I, R, DZ, DY, DX, Pass, ZZ, YY, XX : Integer;
 D, M : Byte;
 IsMax : Boolean;
Begin
 PoreSizeDist := False;
 FillChar(PoreHist, SizeOf(PoreHist), 0);
 PoreTot := 0;
 If Not AllocGridG(VisG, VisN, N) Then Exit;
 If Not AllocGridG(AuxG, AuxN, N) Then Exit;
 ClearGridG(AuxG, N, 0);
 For Z := 0 To N - 1 Do
  For I := 0 To NN - 1 Do Begin
   If (Grid[Z]^[I] = P_POR) Or (Grid[Z]^[I] = P_VOID) Then Begin
    VisG[Z]^[I] := 250; Inc(PoreTot);
   End
   Else VisG[Z]^[I] := 0;
  End;
 If PoreTot = 0 Then Exit;
 For Pass := 1 To 8 Do Begin
  ChamferPass(True);
  ChamferPass(False);
  ShowProgress('Distribution des pores : transformee de distance...', Pass * 6);
 End;
 { balayage des boules maximales }
 For Z := 0 To N - 1 Do Begin
  For Y := 0 To N - 1 Do
   For X := 0 To N - 1 Do Begin
    D := VisG[Z]^[Y * N + X];
    If D >= 3 Then Begin
     IsMax := True;
     For DZ := -1 To 1 Do
      For DY := -1 To 1 Do
       For DX := -1 To 1 Do
        If VisG[Wr[Z + DZ]]^[Wr[Y + DY] * N + Wr[X + DX]] > D Then IsMax := False;
     If IsMax Then Begin
      R := D Div 3;
      If R > MAXR Then R := MAXR;
      If R > 12 Then R := 12 + (R - 12);
      For DZ := -R To R Do
       For DY := -R To R Do
        For DX := -R To R Do
         If DX * DX + DY * DY + DZ * DZ <= R * R + 1 Then Begin
          ZZ := Z + DZ; YY := Y + DY; XX := X + DX;
          If ZZ < 0 Then ZZ := ZZ + N; If ZZ >= N Then ZZ := ZZ - N;
          If YY < 0 Then YY := YY + N; If YY >= N Then YY := YY - N;
          If XX < 0 Then XX := XX + N; If XX >= N Then XX := XX - N;
          If VisG[ZZ]^[YY * N + XX] > 0 Then
           If AuxG[ZZ]^[YY * N + XX] < R Then AuxG[ZZ]^[YY * N + XX] := R;
         End;
     End;
    End;
   End;
  If (Z And 3) = 0 Then ShowProgress('Distribution des pores : boules inscrites...', 50 + (Z * 50) Div N);
  If UserAbort Then Exit;
 End;

 For Z := 0 To N - 1 Do
  For I := 0 To NN - 1 Do
   If VisG[Z]^[I] > 0 Then Begin
    M := AuxG[Z]^[I];
    If M < 1 Then M := 1;
    If M > MAXR Then M := MAXR;
    Inc(PoreHist[M]);
   End;
 PoreSizeDist := True;
End;

{ --------------------------------------------------------------------- }
{ Modules elastiques : milieu effectif auto-coherent                    }
{ --------------------------------------------------------------------- }

Var
 ElK, ElG, ElE, ElNu, ElVoigt, ElPerc : Real;

Procedure BulkShear(Ph: Integer; Var K, G: Real);
Var E, Nu: Real;
Begin
 E := ModE[Ph]; Nu := ModNu[Ph];
 If E <= 0 Then Begin K := 0; G := 0; Exit; End;
 K := E / (3 * (1 - 2 * Nu));
 G := E / (2 * (1 + Nu));
End;

{ Fractions F[ph] (fraction du volume total) ; resout K, G }
Procedure SelfConsistent(Var F: TPhR);
Var
 Ph, It : Integer;
 K, G, KP, GP, A, GS, S1, S2, Fp, KN, GN, Tot : Real;
Begin
 Tot := 0; K := 0; G := 0;
 For Ph := 1 To NPH - 1 Do Begin
  BulkShear(Ph, KP, GP);
  K := K + F[Ph] * KP; G := G + F[Ph] * GP;
  Tot := Tot + F[Ph];
 End;
 ElVoigt := 0;
 For Ph := 1 To NPH - 1 Do ElVoigt := ElVoigt + F[Ph] * ModE[Ph];
 Fp := 1 - Tot;
 If Fp < 0 Then Fp := 0;
 For It := 1 To 400 Do Begin
  If G < 0.0001 Then Begin K := 0; G := 0; Break; End;
  A := 4.0 / 3.0 * G;
  S1 := Fp / A;
  For Ph := 1 To NPH - 1 Do
   If F[Ph] > 0 Then Begin
    BulkShear(Ph, KP, GP);
    S1 := S1 + F[Ph] / (KP + A);
   End;
  KN := 1 / S1 - A;
  GS := G * (9 * K + 8 * G) / (6 * (K + 2 * G));
  S2 := Fp / GS;
  For Ph := 1 To NPH - 1 Do
   If F[Ph] > 0 Then Begin
    BulkShear(Ph, KP, GP);
    S2 := S2 + F[Ph] / (GP + GS);
   End;
  GN := 1 / S2 - GS;
  If KN < 0 Then KN := 0;
  If GN < 0 Then GN := 0;
  If (Abs(KN - K) < 0.00001) And (Abs(GN - G) < 0.00001) Then Begin K := KN; G := GN; Break; End;
  K := 0.5 * K + 0.5 * KN;
  G := 0.5 * G + 0.5 * GN;
 End;
 ElK := K; ElG := G;
 If K + G > 0 Then Begin
  ElE := 9 * K * G / (3 * K + G);
  ElNu := (3 * K - 2 * G) / (2 * (3 * K + G));
 End
 Else Begin ElE := 0; ElNu := 0; End;
End;

{ Modules du reseau de solides lies traversant le systeme }
Function ElasticModuli: Boolean;
Var
 F : TPhR;
 Z, I : Integer;
 Cnt : LongInt;
 Cl : Array[0..NPH-1] Of LongInt;
 Ph : Byte;
Begin
 ElasticModuli := False;
 ElE := 0; ElK := 0; ElG := 0; ElNu := 0; ElVoigt := 0; ElPerc := 0;
 MaskSolids;
 If Not Percolate(True, Cnt) Then Begin
  QFree;
  Exit;
 End;
 FillChar(Cl, SizeOf(Cl), 0);
 For Z := ZLo To ZHi Do
  For I := 0 To NN - 1 Do
   If (VisG[Z]^[I] And 3) = 3 Then Inc(Cl[Grid[Z]^[I]]);
 For Ph := 0 To NPH - 1 Do F[Ph] := Cl[Ph] / PasteVox;
 ElPerc := Cnt / PasteVox;
 SelfConsistent(F);
 ElasticModuli := True;
End;

{ --------------------------------------------------------------------- }
{ Conductivite / diffusivite (differences finies)                       }
{ --------------------------------------------------------------------- }

Var
 CondBeta   : Real;
 CondIters  : Integer;
 CondAxis   : Integer;
 CondConv   : Boolean;

Function Conductivity: Boolean;
Var
 Phi : TRGrid;
 GT  : Array[0..NPH-1, 0..NPH-1] Of Real;
 A, B, Z, Y, X, I, It, XP, XM, YP, YM, Axis : Integer;
 Sg  : ^TGrid;
 Sum, WSum, G, V, NewV, Dm, Om, Flux : Real;
 Q, QL, QU : PSlice;
 P, PL, PU : PRSlice;
Begin
 Conductivity := False;
 CondBeta := 0; CondIters := 0; CondConv := False;
 Axis := 2;
 If WallOn Then Axis := 0;
 CondAxis := Axis;
 For A := 0 To NPH - 1 Do
  For B := 0 To NPH - 1 Do
   If SigTab[A] + SigTab[B] > 0 Then GT[A, B] := 2 * SigTab[A] * SigTab[B] / (SigTab[A] + SigTab[B])
   Else GT[A, B] := 0;

 For Z := 0 To MAXN - 1 Do Phi[Z] := Nil;
 If Axis = 2 Then Sg := @Grid
 Else Begin
  If Not AllocGridG(AuxG, AuxN, N) Then Exit;
  For Z := 0 To N - 1 Do
   For Y := 0 To N - 1 Do
    For X := 0 To N - 1 Do AuxG[X]^[Y * N + Z] := Grid[Z]^[Y * N + X];
  Sg := @AuxG;
 End;
 For Z := 0 To N - 1 Do Begin
  {$IFNDEF FPC}
  If MaxAvail < NN * SizeOf(Real) + 2048 Then Begin
   For I := 0 To Z - 1 Do FreeMem(Phi[I], NN * SizeOf(Real));
   Exit;
  End;
  {$ENDIF}
  GetMem(Phi[Z], NN * SizeOf(Real));
  For I := 0 To NN - 1 Do Phi[Z]^[I] := 1 - (Z + 1.0) / (N + 1.0);
 End;
 Om := 1.8;
 If N <= 32 Then Om := 1.7;
 It := 0;
 Repeat
  Inc(It);
  Dm := 0;
  For Z := 0 To N - 1 Do Begin
   Q := Sg^[Z]; P := Phi[Z];
   If Z > 0 Then Begin QL := Sg^[Z - 1]; PL := Phi[Z - 1]; End;
   If Z < N - 1 Then Begin QU := Sg^[Z + 1]; PU := Phi[Z + 1]; End;
   For Y := 0 To N - 1 Do Begin
    If Y > 0 Then YM := Y - 1 Else YM := N - 1;
    If Y < N - 1 Then YP := Y + 1 Else YP := 0;
    For X := 0 To N - 1 Do Begin
     I := Y * N + X;
     A := Q^[I];
     If SigTab[A] > 0 Then Begin
      If X > 0 Then XM := X - 1 Else XM := N - 1;
      If X < N - 1 Then XP := X + 1 Else XP := 0;
      Sum := 0; WSum := 0;
      B := Q^[Y * N + XP]; G := GT[A, B]; Sum := Sum + G * P^[Y * N + XP]; WSum := WSum + G;
      B := Q^[Y * N + XM]; G := GT[A, B]; Sum := Sum + G * P^[Y * N + XM]; WSum := WSum + G;
      B := Q^[YP * N + X]; G := GT[A, B]; Sum := Sum + G * P^[YP * N + X]; WSum := WSum + G;
      B := Q^[YM * N + X]; G := GT[A, B]; Sum := Sum + G * P^[YM * N + X]; WSum := WSum + G;
      If Z > 0 Then Begin
       B := QL^[I]; G := GT[A, B]; Sum := Sum + G * PL^[I]; WSum := WSum + G;
      End
      Else Begin G := SigTab[A]; Sum := Sum + G; WSum := WSum + G; End;
      If Z < N - 1 Then Begin
       B := QU^[I]; G := GT[A, B]; Sum := Sum + G * PU^[I]; WSum := WSum + G;
      End
      Else Begin G := SigTab[A]; WSum := WSum + G; End;
      If WSum > 0 Then Begin
       V := P^[I];
       NewV := V + Om * (Sum / WSum - V);
       If Abs(NewV - V) > Dm Then Dm := Abs(NewV - V);
       P^[I] := NewV;
      End;
     End;
    End;
   End;
  End;
  If (It Mod 10) = 0 Then Begin
   ShowProgress('Conductivite : iteration ' + IStr(It) + ' (Echap = arret)', IMin(99, It Div 8));
  End;
 Until (Dm < 0.000005) Or (It >= 800) Or UserAbort;
 CondConv := Dm < 0.000005;
 CondIters := It;
 Flux := 0;
 Q := Sg^[0]; P := Phi[0];
 For I := 0 To NN - 1 Do Flux := Flux + SigTab[Q^[I]] * (1 - P^[I]);
 CondBeta := Flux * (N + 1.0) / NN;
 For Z := 0 To N - 1 Do
  If Phi[Z] <> Nil Then FreeMem(Phi[Z], NN * SizeOf(Real));
 Conductivity := True;
End;

{ --------------------------------------------------------------------- }
{ Profil de la zone de transition (mur de granulat)                     }
{ --------------------------------------------------------------------- }

Const
 ITZ_MAX = 50;

Var
 ItzN  : Integer;
 ItzD  : Array[0..ITZ_MAX-1] Of Real;      { distance au granulat (um) }
 ItzP  : Array[0..ITZ_MAX-1] Of Real;      { porosite (%) }
 ItzU  : Array[0..ITZ_MAX-1] Of Real;      { anhydre (%) }
 ItzH  : Array[0..ITZ_MAX-1] Of Real;      { portlandite (%) }
 ItzG  : Array[0..ITZ_MAX-1] Of Real;      { C-S-H et autres hydrates (%) }

Function ItzProfile: Boolean;
Var K, Half : Integer; Cnt : Array[0..NPH-1] Of LongInt; Ph: Integer;
 Tot, Un, Hy: LongInt;
 Procedure AddSlab(Z: Integer);
 Var II: Integer;
 Begin
  For II := 0 To NN - 1 Do Begin Inc(Cnt[Grid[Z]^[II]]); Inc(Tot); End;
 End;
Begin
 ItzProfile := False;
 If Not (HaveGrid And WallOn) Then Exit;
 Half := (N - 2) Div 2;
 If Half > ITZ_MAX Then Half := ITZ_MAX;
 ItzN := Half;
 For K := 0 To Half - 1 Do Begin
  FillChar(Cnt, SizeOf(Cnt), 0);
  Tot := 0;
  AddSlab(1 + K); AddSlab(N - 2 - K);
  Un := 0; Hy := 0;
  For Ph := 1 To 13 Do If (Ph <> P_FILL) Then Un := Un + Cnt[Ph];
  Un := Un + Cnt[P_FILL];
  For Ph := P_CSH To P_PCSH Do Hy := Hy + Cnt[Ph];
  ItzD[K] := (K + 0.5) * ResUm;
  ItzP[K] := 100.0 * (Cnt[P_POR] + Cnt[P_VOID]) / Tot;
  ItzU[K] := 100.0 * Un / Tot;
  ItzH[K] := 100.0 * Cnt[P_CH] / Tot;
  ItzG[K] := 100.0 * (Hy - Cnt[P_CH]) / Tot;
 End;
 ItzProfile := True;
End;

{ --------------------------------------------------------------------- }
{ Ecran, menus et composants d'interface                                }
{ --------------------------------------------------------------------- }

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
 MAX_TL = 300;

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
 FD_MAX_FILES = 256;

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
 MAXSER = 5;
 SerColors : Array[0..MAXSER-1] Of Byte = (1, 4, 2, 5, 6);

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
{ Fenetre principale : vue, panneau, barre d'outils, barre d'etat      }
{ --------------------------------------------------------------------- }

Const
 VX0 = 6;
 VY0 = 48;
 VSZ = 400;
 PNX = 412;
 PNY = 48;

 TB_GEN = 1;  TB_RUN = 2;  TB_PAUSE = 3; TB_STOP = 4;
 TB_XY = 10;  TB_XZ = 11;  TB_YZ = 12;   TB_3D = 13;
 TB_PREV = 14; TB_NEXT = 15;
 TB_CURV = 20; TB_PROP = 21; TB_FORM = 22; TB_HELP = 23;
 MAX_TB = 24;

Var
 ViewPlane : Integer;               { 0 XY, 1 XZ, 2 YZ, 3 vue 3D }
 SliceIdx  : Integer;
 SimPaused : Boolean;
 NTB       : Integer;
 TBId      : Array[0..MAX_TB-1] Of Integer;
 TBX1, TBX2 : Array[0..MAX_TB-1] Of Integer;
 TBCap     : Array[0..MAX_TB-1] Of String[5];
 TBHint    : Array[0..MAX_TB-1] Of String[50];
 LegPh     : Array[0..NPH-1] Of Integer;
 LegN      : Integer;
 LegY0     : Integer;
 LastShowT : LongInt;

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
 AddTB(TB_GEN, 'Gener', 'Generer la microstructure (Melange)', X);
 AddTB(TB_RUN, 'Hydr.', 'Demarrer / reprendre l''hydratation', X);
 AddTB(TB_PAUSE, 'Pause', 'Suspendre l''hydratation', X);
 AddTB(TB_STOP, 'Arret', 'Arreter l''hydratation', X);
 X := X + 6;
 AddTB(TB_XY, 'XY', 'Coupe dans le plan XY', X);
 AddTB(TB_XZ, 'XZ', 'Coupe dans le plan XZ', X);
 AddTB(TB_YZ, 'YZ', 'Coupe dans le plan YZ', X);
 AddTB(TB_3D, '3D', 'Vue 3D du cube', X);
 X := X + 6;
 AddTB(TB_PREV, '<', 'Coupe precedente (PgUp)', X);
 AddTB(TB_NEXT, '>', 'Coupe suivante (PgDn)', X);
 X := X + 6;
 AddTB(TB_CURV, 'Courb', 'Courbes de l''hydratation', X);
 AddTB(TB_PROP, 'Prop.', 'Resume des proprietes', X);
 AddTB(TB_FORM, 'Form.', 'Formulation du melange', X);
 AddTB(TB_HELP, '?', 'Aide', X);
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
  TB_RUN: TBEnabled := HaveGrid And Not SimActive;
  TB_PAUSE, TB_STOP: TBEnabled := SimActive;
  TB_XY..TB_NEXT, TB_CURV, TB_PROP: TBEnabled := HaveGrid;
 Else TBEnabled := True;
 End;
End;

Function TBPressed(Id: Integer): Boolean;
Begin
 TBPressed := (Id >= TB_XY) And (Id <= TB_3D) And (ViewPlane = Id - TB_XY);
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
  If TBPressed(TBId[I]) Then SetFillStyle(SolidFill, LightCyan);
  Bar(TBX1[I], Y1 + 2, TBX2[I], Y2 - 2);
  Bevel(TBX1[I], Y1 + 2, TBX2[I], Y2 - 2, TBPressed(TBId[I]));
  If TBEnabled(TBId[I]) Then SetColor(Black) Else SetColor(DarkGray);
  OutTextXY(TBX1[I] + 4, Y1 + 8, TBCap[I]);
 End;
End;

{ --------------------------------------------------------------------- }
{ Vue des coupes                                                        }
{ --------------------------------------------------------------------- }

Function PhDispColor(Ph: Byte): Byte;
Begin
 If PhShow[Ph] Then PhDispColor := PhColor[Ph] Else PhDispColor := Black;
End;

Function VoxAt(Plane, I, J: Integer): Byte;
Begin
 Case Plane Of
  0: VoxAt := Grid[SliceIdx]^[J * N + I];
  1: VoxAt := Grid[N - 1 - J]^[SliceIdx * N + I];
 Else VoxAt := Grid[N - 1 - J]^[I * N + SliceIdx];
 End;
End;

Procedure DrawSlice;
Var S, OX, OY, I, J, I0: Integer; C, C0: Byte;
Begin
 S := VSZ Div N;
 If S < 1 Then S := 1;
 OX := VX0 + (VSZ - S * N) Div 2; OY := VY0 + (VSZ - S * N) Div 2;
 For J := 0 To N - 1 Do Begin
  I0 := 0; C0 := PhDispColor(VoxAt(ViewPlane, 0, J));
  For I := 1 To N Do Begin
   If I < N Then C := PhDispColor(VoxAt(ViewPlane, I, J)) Else C := 255;
   If C <> C0 Then Begin
    SetFillStyle(SolidFill, C0);
    Bar(OX + I0 * S, OY + J * S, OX + I * S - 1, OY + J * S + S - 1);
    I0 := I; C0 := C;
   End;
  End;
 End;
End;

Procedure Draw3D;
Var S, H, X, Y, Z, FX, FY, TotW, TotH, X1, Y1: Integer; C: Byte;
Begin
 S := 8;
 While (S > 1) And (N * S + N * (S Div 2) > VSZ) Do Dec(S);
 H := S Div 2;
 If H < 1 Then H := 1;
 TotW := N * S + N * H; TotH := N * S + N * H;
 FX := VX0 + (VSZ - TotW) Div 2;
 FY := VY0 + (VSZ - TotH) Div 2 + N * H + N * S - N * S;
 FY := VY0 + (VSZ - TotH) Div 2 + N * H;
 { face du dessus : plan XZ en y = 0 }
 For Z := N - 1 Downto 0 Do
  For X := 0 To N - 1 Do Begin
   C := PhDispColor(Grid[Z]^[X]);
   X1 := FX + X * S + Z * H; Y1 := FY - (Z + 1) * H;
   SetFillStyle(SolidFill, C);
   Bar(X1, Y1, X1 + S - 1, Y1 + H - 1);
  End;
 { face de droite : plan YZ en x = N-1 }
 For Z := N - 1 Downto 0 Do
  For Y := 0 To N - 1 Do Begin
   C := PhDispColor(Grid[Z]^[Y * N + N - 1]);
   X1 := FX + N * S + Z * H; Y1 := FY + Y * S - (Z + 1) * H;
   SetFillStyle(SolidFill, C);
   Bar(X1, Y1, X1 + H - 1, Y1 + S - 1);
  End;
 { face avant : plan XY en z = 0 }
 For Y := 0 To N - 1 Do
  For X := 0 To N - 1 Do Begin
   C := PhDispColor(Grid[0]^[Y * N + X]);
   SetFillStyle(SolidFill, C);
   Bar(FX + X * S, FY + Y * S, FX + X * S + S - 1, FY + Y * S + S - 1);
  End;
 SetColor(White);
 Rectangle(FX, FY, FX + N * S, FY + N * S);
 Line(FX, FY, FX + N * H, FY - N * H);
 Line(FX + N * S, FY, FX + N * S + N * H, FY - N * H);
 Line(FX + N * H, FY - N * H, FX + N * S + N * H, FY - N * H);
 Line(FX + N * S, FY + N * S, FX + N * S + N * H, FY + N * S - N * H);
 Line(FX + N * S + N * H, FY - N * H, FX + N * S + N * H, FY + N * S - N * H);
End;

Procedure DrawWelcome;
Var Y: Integer;

 Procedure L(Const S: String; Ok: Boolean; Col: Byte);
 Begin
  SetColor(Col);
  If Ok Then OutTextXY(VX0 + 14, Y, '[x] ' + S) Else OutTextXY(VX0 + 14, Y, '[ ] ' + S);
  Y := Y + 14;
 End;

Begin
 SetFillStyle(SolidFill, White);
 Bar(VX0, VY0, VX0 + VSZ - 1, VY0 + VSZ - 1);
 SetColor(Blue);
 SetTextStyle(DefaultFont, HorizDir, 1);
 OutTextXY(VX0 + 14, VY0 + 14, 'VCCTL - Laboratoire virtuel d''essais du');
 OutTextXY(VX0 + 14, VY0 + 26, 'ciment et du beton (clone)');
 SetColor(DarkGray);
 Line(VX0 + 14, VY0 + 40, VX0 + VSZ - 14, VY0 + 40);
 Y := VY0 + 56;
 SetColor(Black);
 OutTextXY(VX0 + 14, Y, 'Marche a suivre :'); Y := Y + 20;
 L('1. Materiaux : ciment, ajouts, granulats', True, Black);
 L('2. Melange : E/L, proportions, taille', True, Black);
 L('3. Melange > Generer la microstructure', False, Black);
 L('4. Hydratation > Demarrer', False, Black);
 L('5. Proprietes : porosite, modules,', False, Black);
 Y := Y - 0;
 OutTextXY(VX0 + 14 + 32, Y, 'resistance, conductivite...'); Y := Y + 20;
 SetColor(DarkGray);
 OutTextXY(VX0 + 14, Y, 'Raccourcis :'); Y := Y + 14;
 OutTextXY(VX0 + 14, Y, ' Ctrl+G   generer la microstructure'); Y := Y + 12;
 OutTextXY(VX0 + 14, Y, ' F9       demarrer l''hydratation'); Y := Y + 12;
 OutTextXY(VX0 + 14, Y, ' Espace   pause / reprise'); Y := Y + 12;
 OutTextXY(VX0 + 14, Y, ' PgUp/Dn  coupe precedente / suivante'); Y := Y + 12;
 OutTextXY(VX0 + 14, Y, ' F10      menus     F1  aide'); Y := Y + 12;
 SetColor(Blue);
 OutTextXY(VX0 + 14, VY0 + VSZ - 22, 'Ciment : ' + Copy(Cem.Name, 1, 34));
 OutTextXY(VX0 + 14, VY0 + VSZ - 12, 'E/L = ' + RTrim(WB, 3) + '   N = ' + IStr(NVox) + '   graine = ' + IStr(Seed));
End;

Procedure DrawViewer;
Var S: String;
Begin
 SetFillStyle(SolidFill, LightGray);
 Bar(VX0 - 3, VY0 - 3, VX0 + VSZ + 2, VY0 + VSZ + 12);
 If Not HaveGrid Then DrawWelcome
 Else Begin
  SetFillStyle(SolidFill, Black);
  Bar(VX0, VY0, VX0 + VSZ - 1, VY0 + VSZ - 1);
  If SliceIdx >= N Then SliceIdx := N - 1;
  If ViewPlane = 3 Then Draw3D Else DrawSlice;
 End;
 Bevel(VX0 - 2, VY0 - 2, VX0 + VSZ + 1, VY0 + VSZ + 1, True);
 SetColor(Black);
 SetTextStyle(DefaultFont, HorizDir, 1);
 If HaveGrid Then Begin
  Case ViewPlane Of
   0: S := 'Coupe XY  z = ' + IStr(SliceIdx + 1) + ' / ' + IStr(N);
   1: S := 'Coupe XZ  y = ' + IStr(SliceIdx + 1) + ' / ' + IStr(N);
   2: S := 'Coupe YZ  x = ' + IStr(SliceIdx + 1) + ' / ' + IStr(N);
  Else S := 'Vue 3D : faces z=1, y=1, x=' + IStr(N);
  End;
  OutTextXY(VX0, VY0 + VSZ + 4, S);
 End;
End;

{ --------------------------------------------------------------------- }
{ Panneau d'information                                                 }
{ --------------------------------------------------------------------- }

Function StateText: String;
Begin
 If Not HaveGrid Then StateText := 'Aucune'
 Else If SimActive Then StateText := 'En cours'
 Else If Sim.Cyc = 0 Then StateText := 'Pret (t = 0)'
 Else If SimPaused Then StateText := 'En pause'
 Else StateText := 'Termine';
End;

Function TimeText(H: Real): String;
Begin
 If H < 48 Then TimeText := RStr(H, 1) + ' h'
 Else TimeText := RStr(H / 24, 1) + ' j';
End;

Procedure DrawPanel;
Var Y, I: Integer; Pv: Real; NP: LongInt;

 Procedure PLn(Const A, B: String);
 Begin
  SetColor(DarkGray); OutTextXY(PNX + 6, Y, A);
  SetColor(Black); OutTextXY(PNX + 6 + 8 * 13, Y, Copy(B, 1, 12));
  Y := Y + 11;
 End;

Begin
 SetFillStyle(SolidFill, LightGray);
 Bar(PNX - 2, PNY - 2, ScreenWidth - 1, ScreenHeight - STATUS_H - 1);
 SetFillStyle(SolidFill, White);
 Bar(PNX, PNY, ScreenWidth - 4, ScreenHeight - STATUS_H - 8);
 Bevel(PNX - 1, PNY - 1, ScreenWidth - 3, ScreenHeight - STATUS_H - 7, True);
 SetFillStyle(SolidFill, Blue);
 Bar(PNX, PNY, ScreenWidth - 4, PNY + 11);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(White);
 OutTextXY(PNX + 4, PNY + 2, 'Informations');
 Y := PNY + 16;
 PLn('Etat', StateText);
 If HaveGrid Then Begin
  PLn('Systeme', IStr(N) + 'x' + IStr(N) + 'x' + IStr(N));
  PLn('Taille', RTrim(N * ResUm, 1) + ' um');
  PLn('Temps', TimeText(Sim.T));
  PLn('Cycles', IStr(Sim.Cyc));
  PLn('Temperature', RStr(Sim.Temp, 1) + ' C');
  If Sim.NS > 0 Then Begin
   PLn('alpha ciment', RStr(AlphaCement, 3));
   PLn('Chaleur', RStr(Sim.Heat / Sim.M0Bind, 1) + ' J/g');
   PLn('Porosite cap.', RStr(100 * CapPorosity, 1) + ' %');
   If Sim.TSet >= 0 Then PLn('Prise', TimeText(Sim.TSet)) Else PLn('Prise', 'non atteinte');
  End;
  Y := Y + 4;
  SetColor(Blue);
  OutTextXY(PNX + 6, Y, 'Phases (% volume)'); Y := Y + 12;
  LegY0 := Y; LegN := 0;
  NP := PasteVox;
  For I := 0 To NPH - 1 Do
   If (PhCnt[I] > 0) And (I <> P_AGG) Then Begin
    LegPh[LegN] := I;
    Inc(LegN);
    SetFillStyle(SolidFill, PhColor[I]);
    Bar(PNX + 6, Y, PNX + 16, Y + 8);
    SetColor(Black);
    Rectangle(PNX + 6, Y, PNX + 16, Y + 8);
    If Not PhShow[I] Then Begin SetColor(Red); Line(PNX + 6, Y, PNX + 16, Y + 8); Line(PNX + 16, Y, PNX + 6, Y + 8); End;
    SetColor(Black);
    OutTextXY(PNX + 22, Y + 1, Copy(PhName[I], 1, 18));
    Pv := 100.0 * PhCnt[I] / NP;
    OutTextXY(PNX + 22 + 8 * 18, Y + 1, PadL(RStr(Pv, 1), 5));
    Y := Y + 10;
   End;
 End;
End;

Function LegendHit(MX, MY: Integer): Integer;
Var K: Integer;
Begin
 LegendHit := -1;
 If (MX < PNX) Or (MX > PNX + 130) Or (MY < LegY0) Then Exit;
 K := (MY - LegY0) Div 10;
 If (K >= 0) And (K < LegN) Then LegendHit := LegPh[K];
End;

{ --------------------------------------------------------------------- }
{ Barre d'etat                                                          }
{ --------------------------------------------------------------------- }

Var ProgPct : Integer;

Function BaseName(Const FN: String): String;
Var I: Integer;
Begin
 I := Length(FN);
 While (I > 0) And (FN[I] <> '\') And (FN[I] <> '/') Do Dec(I);
 BaseName := Copy(FN, I + 1, 255);
End;

Procedure DrawStatusBar;
Var Y, W, MaxC: Integer; S, PN: String;
Begin
 Y := ScreenHeight - STATUS_H;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, Y, ScreenWidth - 1, ScreenHeight - 1);
 SetColor(DarkGray);
 Line(0, Y, ScreenWidth - 1, Y);
 SetTextStyle(DefaultFont, HorizDir, 1);
 SetColor(Black);
 PN := BaseName(ProjName);
 If PN = '' Then PN := '(sans nom)';
 If Modified Then PN := PN + ' *';
 S := StatusMsg;
 If ProgPct >= 0 Then MaxC := (ScreenWidth - 170) Div 8
 Else MaxC := (ScreenWidth - 8 * Length(PN) - 24) Div 8;
 OutTextXY(6, Y + 5, Copy(S, 1, MaxC));
 If ProgPct >= 0 Then Begin
  W := 150;
  SetFillStyle(SolidFill, White);
  Bar(ScreenWidth - W - 8, Y + 4, ScreenWidth - 8, Y + STATUS_H - 4);
  SetFillStyle(SolidFill, Blue);
  If ProgPct > 100 Then ProgPct := 100;
  Bar(ScreenWidth - W - 8, Y + 4, ScreenWidth - W - 8 + (W * ProgPct) Div 100, Y + STATUS_H - 4);
  SetColor(Black);
  Rectangle(ScreenWidth - W - 8, Y + 4, ScreenWidth - 8, Y + STATUS_H - 4);
 End
 Else Begin
  SetColor(Blue);
  OutTextXY(ScreenWidth - 8 * Length(PN) - 8, Y + 5, PN);
 End;
End;

Procedure ShowProgress(Const Msg: String; Pct: Integer);
Var K: Char; T: LongInt;
Begin
 If BatchMode Then Exit;
 If KeyPressed Then Begin
  K := ReadKey;
  If K = #27 Then UserAbort := True
  Else If K = #0 Then If KeyPressed Then K := ReadKey;
 End;
 T := NowMs;
 If (Pct < 100) And (Abs(T - LastShowT) < 60) Then Exit;
 LastShowT := T;
 StatusMsg := Msg;
 ProgPct := Pct;
 HideMouseCursor;
 DrawStatusBar;
 ShowMouseCursor;
End;



Procedure RedrawAll;
Begin
 HideMouseCursor;
 SetFillStyle(SolidFill, LightGray);
 Bar(0, MENUBAR_H + TOOLBAR_H, ScreenWidth - 1, ScreenHeight - 1);
 DrawViewer;
 DrawPanel;
 DrawToolbar;
 DrawStatusBar;
 DrawMenuBar;
 If MenuActive Then DrawMenuPopup(ActiveMenu);
 ShowMouseCursor;
 NeedRedraw := False;
End;

{ remise a jour rapide pendant la simulation }
Procedure RedrawSim;
Begin
 HideMouseCursor;
 DrawViewer;
 DrawPanel;
 DrawToolbar;
 DrawStatusBar;
 ShowMouseCursor;
End;

{ --------------------------------------------------------------------- }
{ Dialogues : materiaux                                                 }
{ --------------------------------------------------------------------- }

Procedure CementDialog;
Var I: Integer; Sel: Integer; Sum: Real; Done: Boolean;
Begin
 Done := False;
 While Not Done Do Begin
  Sum := 0;
  For I := 1 To 9 Do Sum := Sum + Cem.W[I];
  PReset;
  PAddE('Modele de la bibliotheque', 0, '(valeurs actuelles)|' + CemPresetName[1] + '|' + CemPresetName[2] +
        '|' + CemPresetName[3] + '|' + CemPresetName[4] + '|' + CemPresetName[5] + '|' + CemPresetName[6]);
  PAddH('Composition (% massique)');
  PAddR('C3S (alite)', Cem.W[1], 0, 100, 2);
  PAddR('C2S (belite)', Cem.W[2], 0, 100, 2);
  PAddR('C3A (aluminate)', Cem.W[3], 0, 100, 2);
  PAddR('C4AF (ferrite)', Cem.W[4], 0, 100, 2);
  PAddR('Gypse', Cem.W[5], 0, 100, 2);
  PAddR('Hemihydrate', Cem.W[6], 0, 100, 2);
  PAddR('Anhydrite', Cem.W[7], 0, 100, 2);
  PAddR('Chaux libre', Cem.W[8], 0, 100, 2);
  PAddR('Sulfates alcalins (K2SO4)', Cem.W[9], 0, 100, 2);
  PAddT('Somme (normalisee a 100 %)', RStr(Sum, 2) + ' %');
  PAddH('Granulometrie de Rosin-Rammler');
  PAddR('d'' - diametre a 63,2 % (um)', Cem.P.D63, 0.5, 200, 2);
  PAddR('n - pente', Cem.P.Nr, 0.3, 4, 2);
  PAddR('Diametre minimal (um)', Cem.P.DMin, 0.5, 20, 1);
  PAddR('Diametre maximal (um)', Cem.P.DMax, 5, 200, 1);
  PAddH('Valeurs deduites');
  PAddT('Masse volumique (g/cm3)', RStr(CemRho, 3));
  PAddT('Surface specifique geometrique (m2/kg)', RStr(SpecSurface(Cem.P, CemRho), 0));
  If Not ParamDialog('Materiaux : ciment - ' + Cem.Name, 'Choisir un modele puis OK pour le charger') Then Exit;
  Sel := PInt(0);
  If Sel > 0 Then Begin
   LoadCementPreset(Sel);
   StatusMsg := 'Ciment charge : ' + Cem.Name;
  End
  Else Begin
   For I := 1 To 9 Do Cem.W[I] := PVal[1 + I];
   Cem.P.D63 := PVal[13]; Cem.P.Nr := PVal[14]; Cem.P.DMin := PVal[15]; Cem.P.DMax := PVal[16];
   If Cem.P.DMax <= Cem.P.DMin Then Cem.P.DMax := Cem.P.DMin + 1;
   Sum := 0;
   For I := 1 To 9 Do Sum := Sum + Cem.W[I];
   If Sum <= 0 Then LoadCementPreset(1)
   Else For I := 1 To 9 Do Cem.W[I] := Cem.W[I] * 100 / Sum;
   Cem.Name := 'Ciment personnalise';
   Modified := True;
   Done := True;
  End;
 End;
End;

Procedure SCMDialog(Kind: Integer; Var Mt: TSCM);
Var Sel, NP: Integer; Done: Boolean; T: String; XL: String;
Begin
 Done := False;
 Case Kind Of
  1: Begin T := 'Cendre volante'; NP := 2; XL := 'Fraction vitreuse reactive (%)'; End;
  2: Begin T := 'Laitier'; NP := 2; XL := 'Parametre libre (non utilise)'; End;
  3: Begin T := 'Calcaire / filler inerte'; NP := 2; XL := 'Facteur de nucleation de la portlandite'; End;
 Else Begin T := 'Fumee de silice'; NP := 1; XL := 'Parametre libre (non utilise)'; End;
 End;
 While Not Done Do Begin
  PReset;
  If NP = 2 Then PAddE('Modele de la bibliotheque', 0, '(valeurs actuelles)|Modele 1|Modele 2')
  Else PAddE('Modele de la bibliotheque', 0, '(valeurs actuelles)|Modele par defaut');
  PAddR('Masse volumique (g/cm3)', Mt.Rho, 1.0, 5.0, 3);
  PAddR('d'' - diametre a 63,2 % (um)', Mt.P.D63, 0.2, 200, 2);
  PAddR('n - pente', Mt.P.Nr, 0.3, 4, 2);
  PAddR('Diametre minimal (um)', Mt.P.DMin, 0.2, 20, 1);
  PAddR('Diametre maximal (um)', Mt.P.DMax, 1, 200, 1);
  PAddR('Reactivite relative', Mt.React, 0, 10, 2);
  PAddR(XL, Mt.Extra, 0, 100, 1);
  PAddT('Surface specifique (m2/kg)', RStr(SpecSurface(Mt.P, Mt.Rho), 0));
  If Not ParamDialog('Materiaux : ' + T + ' - ' + Mt.Name, 'Entree : modifier une valeur') Then Exit;
  Sel := PInt(0);
  If Sel > 0 Then Begin
   Case Kind Of
    1: LoadFAPreset(Sel);
    2: LoadSlagPreset(Sel);
    3: LoadFillerPreset(Sel);
   Else LoadSFPreset;
   End;
  End
  Else Begin
   Mt.Rho := PVal[1]; Mt.P.D63 := PVal[2]; Mt.P.Nr := PVal[3];
   Mt.P.DMin := PVal[4]; Mt.P.DMax := PVal[5]; Mt.React := PVal[6]; Mt.Extra := PVal[7];
   If Mt.P.DMax <= Mt.P.DMin Then Mt.P.DMax := Mt.P.DMin + 1;
   Modified := True;
   Done := True;
  End;
 End;
End;

Procedure AggDialog(Const T: String; Var Ag: TAggr);
Begin
 PReset;
 PAddR('Masse volumique (g/cm3)', Ag.Rho, 1.5, 4.0, 3);
 PAddR('Diametre minimal (mm)', Ag.DMin, 0.01, 50, 3);
 PAddR('Diametre maximal (mm)', Ag.DMax, 0.1, 100, 2);
 PAddR('Exposant de Fuller q', Ag.Q, 0.2, 1.0, 2);
 If ParamDialog('Materiaux : ' + T, 'Courbe de Fuller : passant = (d / Dmax) ^ q') Then Begin
  Ag.Rho := PVal[0]; Ag.DMin := PVal[1]; Ag.DMax := PVal[2]; Ag.Q := PVal[3];
  If Ag.DMax <= Ag.DMin Then Ag.DMax := Ag.DMin * 2;
  Modified := True;
 End;
End;

Procedure ShowPSDTable;
Var M: Integer; P: TPSD; Rho, S: Real;
Begin
 TClear;
 TAdd('Materiau            Masse vol.  d''(um)    n    Surface (m2/kg)');
 TAdd('-------------------------------------------------------------');
 For M := 0 To 4 Do Begin
  GetPSD(M, P); Rho := MatRho(M);
  S := SpecSurface(P, Rho);
  TAdd(PadR(MatName(M), 18) + PadL(RStr(Rho, 3), 8) + PadL(RStr(P.D63, 1), 10) +
       PadL(RStr(P.Nr, 2), 7) + PadL(RStr(S, 0), 12));
 End;
 TAdd('');
 TAdd('Surface specifique geometrique des particules spheriques (valeurs');
 TAdd('indicatives : la surface Blaine est en general un peu plus elevee).');
 TView('Caracteristiques des poudres');
End;

Procedure ShowPSDChart;
Var M, I, S: Integer; P: TPSD; D: Real;
Begin
 ChartReset('Granulometrie des materiaux (Rosin-Rammler)', 'Diametre (um)', 'Passant cumule (% masse)', True);
 For M := 0 To 4 Do Begin
  If (Frac[M] > 0) Or (M = 0) Then Begin
   GetPSD(M, P);
   S := ChartAddSeries(MatName(M), False);
   For I := 0 To 59 Do Begin
    D := P.DMin * Exp(Ln(P.DMax / P.DMin) * I / 59);
    ChartPoint(S, D, 100 * PSDCum(P, D, P.DMax));
   End;
  End;
 End;
 ChFixY0 := True;
 ChartDialog;
End;

{ --------------------------------------------------------------------- }
{ Dialogues : formulation                                               }
{ --------------------------------------------------------------------- }

Procedure MixDialog;
Var I: Integer; S: Real;
Begin
 PReset;
 PAddE('Type de melange', MixKind, 'Pate de ciment|Mortier|Beton');
 PAddR('Rapport eau / liant (E/L)', WB, 0.15, 1.2, 3);
 PAddH('Composition du liant (% massique)');
 PAddR('Ciment', Frac[0], 0, 100, 1);
 PAddR('Cendre volante', Frac[1], 0, 100, 1);
 PAddR('Laitier', Frac[2], 0, 100, 1);
 PAddR('Calcaire / filler', Frac[3], 0, 100, 1);
 PAddR('Fumee de silice', Frac[4], 0, 100, 1);
 PAddH('Microstructure 3D');
 PAddI('Cote du systeme (voxels)', NVox, 16, MAXN);
 PAddR('Resolution (um par voxel)', ResUm, 0.25, 10, 2);
 PAddI('Graine aleatoire', Seed, 1, 2000000000);
 PAddB('Mur de granulat (zone de transition)', WallOn);
 If ParamDialog('Melange : formulation', 'Les fractions du liant sont normalisees a 100 %') Then Begin
  MixKind := PInt(0); WB := PVal[1];
  For I := 0 To 4 Do Frac[I] := PVal[3 + I];
  S := 0;
  For I := 0 To 4 Do S := S + Frac[I];
  If S <= 0 Then Begin Frac[0] := 100; S := 100; End;
  For I := 0 To 4 Do Frac[I] := Frac[I] * 100 / S;
  NVox := PInt(9); ResUm := PVal[10]; Seed := PInt(11); WallOn := PBool(12);
  Modified := True;
  StatusMsg := 'Formulation mise a jour : generer la microstructure';
 End;
End;

Procedure ConcreteDialog;
Begin
 PReset;
 PAddR('Dosage en liant (kg/m3 de beton)', ConcBind, 100, 700, 0);
 PAddR('Sable (kg/m3)', ConcSand, 0, 1500, 0);
 PAddR('Gravier (kg/m3)', ConcCoarse, 0, 1500, 0);
 PAddR('Air occlus (%)', ConcAir, 0, 10, 1);
 If ParamDialog('Melange : proportions du beton', 'Le rapport E/L provient de la formulation') Then Begin
  ConcBind := PVal[0]; ConcSand := PVal[1]; ConcCoarse := PVal[2]; ConcAir := PVal[3];
  Modified := True;
 End;
End;

Procedure ShowProportions;
Var V, M: TR6; I: Integer; T, Wkg, VPaste, VS, VG, VAir, Tot: Real;
Begin
 MixVolumes(V, M);
 T := 0;
 For I := 0 To 5 Do T := T + V[I];
 TClear;
 TAdd('Pate : pour 100 g de liant');
 TAdd('Constituant          Masse (g)   Volume (cm3)   % volume');
 TAdd('---------------------------------------------------------');
 For I := 0 To 4 Do
  If Frac[I] > 0 Then
   TAdd(PadR(MatName(I), 20) + PadL(RStr(100 * M[I], 2), 8) + PadL(RStr(100 * V[I], 2), 14) + PadL(RStr(100 * V[I] / T, 1), 11));
 TAdd(PadR('Eau', 20) + PadL(RStr(100 * M[5], 2), 8) + PadL(RStr(100 * V[5], 2), 14) + PadL(RStr(100 * V[5] / T, 1), 11));
 TAdd('---------------------------------------------------------');
 TAdd('Masse volumique de la pate : ' + RStr((100 + 100 * WB) / (100 * T), 3) + ' g/cm3');
 TAdd('Fraction volumique de solides : ' + RStr(100 * (T - V[5]) / T, 1) + ' %');
 TAdd('Densite du ciment : ' + RStr(CemRho, 3) + ' g/cm3');
 If MixKind > 0 Then Begin
  TAdd('');
  TAdd('Beton / mortier : pour 1 m3');
  Wkg := ConcBind * WB;
  VPaste := ConcBind * T;
  VS := ConcSand / Sand.Rho;
  If MixKind = 2 Then VG := ConcCoarse / Coarse.Rho Else VG := 0;
  VAir := 10 * ConcAir;
  Tot := VPaste + VS + VG + VAir;
  TAdd('Liant        : ' + RStr(ConcBind, 0) + ' kg');
  TAdd('Eau          : ' + RStr(Wkg, 0) + ' kg');
  TAdd('Sable        : ' + RStr(ConcSand, 0) + ' kg  (' + RStr(VS, 0) + ' L)');
  If MixKind = 2 Then TAdd('Gravier      : ' + RStr(ConcCoarse, 0) + ' kg  (' + RStr(VG, 0) + ' L)');
  TAdd('Air occlus   : ' + RStr(VAir, 0) + ' L');
  TAdd('Volume de pate : ' + RStr(VPaste, 0) + ' L');
  TAdd('Volume total   : ' + RStr(Tot, 0) + ' L');
  If Abs(Tot - 1000) > 15 Then TAdd('Attention : le volume total devrait etre voisin de 1000 L.');
  TAdd('Masse volumique du beton : ' + RStr((ConcBind + Wkg + ConcSand + ConcCoarse) / 1000, 3) + ' kg/L');
 End;
 TAdd('');
 TAdd('Microstructure : N = ' + IStr(NVox) + ', ' + RStr(ResUm, 2) + ' um/voxel, ' +
      IStr(LongInt(NVox) * NVox * NVox) + ' voxels');
 TView('Proportions du melange');
End;

{ --------------------------------------------------------------------- }
{ Dialogues : hydratation                                               }
{ --------------------------------------------------------------------- }

Procedure HydParamsDialog;
Begin
 PReset;
 PAddH('Conditions de cure');
 PAddR('Temperature (C)', CureTemp, 0, 90, 1);
 PAddE('Mode de cure', CureMode, 'Isotherme sature (eau externe)|Isotherme scelle|Adiabatique');
 PAddR('Duree de la simulation (jours)', EndDays, 0.05, 3650, 2);
 PAddI('Nombre maximal de cycles', MaxCycles, 100, 200000);
 PAddH('Loi de temps  t = beta x cycle^2');
 PAddR('beta (h par cycle^2)', Beta, 0.000001, 1, 6);
 PAddR('Energie d''activation (kJ/mol)', EaCem, 0, 80, 1);
 PAddR('Temperature de reference (C)', TRefC, 5, 60, 1);
 PAddH('Automate cellulaire');
 PAddI('Cycles de periode dormante', DormCyc, 0, 500);
 PAddR('Nucleation de la portlandite', PNucCH, 0, 1, 3);
 PAddI('Pas de diffusion maximal', MaxAge, 5, 300);
 PAddR('Acces de l''eau : couche proche', AccNear, 0, 1, 2);
 PAddR('Acces de l''eau : couche eloignee', AccFar, 0, 1, 2);
 If ParamDialog('Hydratation : parametres', 'beta calibre la correspondance cycles / temps (voir Calibrer)') Then Begin
  CureTemp := PVal[1]; CureMode := PInt(2); EndDays := PVal[3]; MaxCycles := PInt(4);
  Beta := PVal[6]; EaCem := PVal[7]; TRefC := PVal[8];
  DormCyc := PInt(10); PNucCH := PVal[11]; MaxAge := PInt(12);
  AccNear := PVal[13]; AccFar := PVal[14];
  If HaveGrid Then Begin
   If CureMode = 0 Then RefillVoids;
   If (Sim.Cyc = 0) And (CureMode <> 2) Then Sim.Temp := CureTemp;
  End;
  Modified := True;
 End;
End;

Procedure KineticsDialog;
Begin
 PReset;
 PAddH('Probabilite de dissolution par cycle (surface accessible)');
 PAddR('C3S (alite)', Rate[P_C3S], 0, 0.9, 5);
 PAddR('C2S (belite)', Rate[P_C2S], 0, 0.9, 5);
 PAddR('C3A (aluminate)', Rate[P_C3A], 0, 0.9, 5);
 PAddR('C4AF (ferrite)', Rate[P_C4AF], 0, 0.9, 5);
 PAddR('Chaux libre', Rate[P_LIME], 0, 0.9, 5);
 PAddR('Fumee de silice', Rate[P_SF], 0, 0.9, 5);
 PAddR('Cendre volante (verre)', Rate[P_FA], 0, 0.9, 5);
 PAddR('Laitier', Rate[P_SLAG], 0, 0.9, 5);
 PAddH('Les sulfates sont dissous a la demande (reaction du C3A)');
 If ParamDialog('Hydratation : cinetique des phases', 'Multipliee par la reactivite relative des ajouts') Then Begin
  Rate[P_C3S] := PVal[1]; Rate[P_C2S] := PVal[2]; Rate[P_C3A] := PVal[3];
  Rate[P_C4AF] := PVal[4]; Rate[P_LIME] := PVal[5]; Rate[P_SF] := PVal[6];
  Rate[P_FA] := PVal[7]; Rate[P_SLAG] := PVal[8];
  Modified := True;
 End;
End;

{ --------------------------------------------------------------------- }
{ Generation et simulation                                              }
{ --------------------------------------------------------------------- }

Procedure DoGenerate;
Var Ok: Boolean;
Begin
 If SimActive Then Begin
  ShowMessage('Simulation en cours', 'Arretez d''abord l''hydratation.');
  Exit;
 End;
 If HaveGrid And (Sim.Cyc > 0) Then
  If Not ConfirmBox('Generer', 'La microstructure actuelle (hydratee) sera|remplacee. Continuer ?') Then Exit;
 UserAbort := False;
 ProgPct := 0;
 Ok := GenerateMicro;
 ProgPct := -1;
 If Not Ok Then Begin
  ShowMessage('Generation impossible', StatusMsg);
  Exit;
 End;
 SampleInitial;
 SliceIdx := N Div 2;
 SimPaused := False;
 StatusMsg := 'Microstructure generee : ' + IStr(N) + '^3 voxels';
 Modified := True;
 NeedRedraw := True;
End;

Procedure StartSim(GoalHours: Real);
Begin
 If Not HaveGrid Then Begin
  ShowMessage('Hydratation', 'Generez d''abord une microstructure|(menu Melange).');
  Exit;
 End;
 If Sim.T >= GoalHours Then Begin
  ShowMessage('Hydratation', 'La microstructure a deja atteint ' + RStr(GoalHours / 24, 1) + ' j.|Augmentez la duree dans Parametres.');
  Exit;
 End;
 SimGoalT := GoalHours;
 If Sim.Cyc = 0 Then Begin
  Sim.Temp := CureTemp;
  If Sim.NS = 0 Then SampleInitial;
 End;
 SimActive := True; SimPaused := False; UserAbort := False;
 StatusMsg := 'Hydratation en cours...';
 NeedRedraw := True;
End;

Procedure StopSim(Pause: Boolean);
Begin
 If SimActive Then Begin
  SimActive := False; SimPaused := Pause;
  If Sim.NS > 0 Then Sample;
  If Pause Then StatusMsg := 'Hydratation suspendue' Else StatusMsg := 'Hydratation arretee';
  NeedRedraw := True;
 End;
End;

{ avance la simulation pendant MaxMs millisecondes }
Procedure StepSim(MaxMs: LongInt);
Var T0: LongInt;
Begin
 T0 := NowMs;
 Repeat
  RunCycle;
  If (Sim.T >= SimGoalT) Or (Sim.Cyc >= MaxCycles) Then Begin
   SimActive := False; SimPaused := False;
   If Sim.NS > 0 Then Sample;
   StatusMsg := 'Hydratation terminee : t = ' + TimeText(Sim.T) + ', alpha = ' + RStr(AlphaCement, 3);
   ProgPct := -1;
   Modified := True;
   Exit;
  End;
 Until (NowMs - T0 >= MaxMs) Or (NowMs < T0);
 If SimGoalT > 0 Then ProgPct := Round(100 * Sim.T / SimGoalT);
 StatusMsg := 'Hydratation : t = ' + TimeText(Sim.T) + '  alpha = ' + RStr(AlphaCement, 3) +
              '  (Espace = pause)';
End;

Procedure AdvanceDialog;
Var Idx: Integer; Hrs: Real;
Begin
 TClear;
 TAdd('1 heure'); TAdd('3 heures'); TAdd('6 heures'); TAdd('12 heures');
 TAdd('1 jour'); TAdd('3 jours'); TAdd('7 jours'); TAdd('14 jours');
 TAdd('28 jours'); TAdd('56 jours'); TAdd('90 jours'); TAdd('180 jours'); TAdd('365 jours');
 Idx := ListDialog('Hydratation : avancer jusqu''a', True);
 If Idx < 0 Then Exit;
 Case Idx Of
  0: Hrs := 1; 1: Hrs := 3; 2: Hrs := 6; 3: Hrs := 12; 4: Hrs := 24;
  5: Hrs := 72; 6: Hrs := 168; 7: Hrs := 336; 8: Hrs := 672; 9: Hrs := 1344;
  10: Hrs := 2160; 11: Hrs := 4320;
 Else Hrs := 8760;
 End;
 EndDays := Hrs / 24;
 If MaxCycles < 100 Then MaxCycles := 100;
 StartSim(Hrs);
End;

{ --------------------------------------------------------------------- }
{ Proprietes : rapports et graphiques                                   }
{ --------------------------------------------------------------------- }

Function NeedGrid: Boolean;
Begin
 NeedGrid := HaveGrid;
 If Not HaveGrid Then ShowMessage('Proprietes', 'Aucune microstructure : utilisez Melange > Generer.');
End;

Procedure ShowPhaseTable;
Var I: Integer; NP: LongInt; MTot, Mi: Real;
Begin
 If Not NeedGrid Then Exit;
 NP := PasteVox;
 MTot := 0;
 For I := 0 To NPH - 1 Do If (I <> P_VOID) And (I <> P_AGG) And (I <> P_POR) Then MTot := MTot + PhRho[I] * PhCnt[I];
 TClear;
 TAdd('Phase                      Voxels   % volume  % masse solide');
 TAdd('------------------------------------------------------------');
 For I := 0 To NPH - 1 Do
  If (PhCnt[I] > 0) And (I <> P_AGG) Then Begin
   If (I = P_POR) Or (I = P_VOID) Then Mi := 0 Else Mi := 100 * PhRho[I] * PhCnt[I] / MTot;
   TAdd(PadR(PhName[I], 24) + PadL(IStr(PhCnt[I]), 9) + PadL(RStr(100.0 * PhCnt[I] / NP, 2), 10) +
        PadL(RStr(Mi, 2), 12));
  End;
 TAdd('------------------------------------------------------------');
 TAdd(PadR('Total (pate)', 24) + PadL(IStr(NP), 9) + PadL('100.00', 10));
 TAdd('');
 TAdd('t = ' + TimeText(Sim.T) + ' ; cycles = ' + IStr(Sim.Cyc));
 TView('Fractions volumiques des phases');
End;

Procedure ShowSummary;
Var Cnt, CntP: LongInt; Pc, PPc: Boolean; X: Real;
Begin
 If Not NeedGrid Then Exit;
 MaskSolids;
 Pc := Percolate(True, Cnt);
 QFree;
 MaskPores;
 PPc := Percolate(True, CntP);
 QFree;
 TClear;
 TAdd('Etat : ' + StateText);
 TAdd('Temps : ' + RStr(Sim.T, 2) + ' h (' + RStr(Sim.T / 24, 2) + ' j)   cycles : ' + IStr(Sim.Cyc));
 TAdd('Temperature : ' + RStr(Sim.Temp, 1) + ' C   mode : ' + WordLine('satur|scelle|adiabatique', CureMode + 1));
 TAdd('');
 TAdd('Degre d''hydratation du ciment (clinker) : ' + RStr(AlphaCement, 4));
 TAdd('Degre de reaction du liant                : ' + RStr(AlphaBinder, 4));
 TAdd('Chaleur degagee : ' + RStr(Sim.Heat / Sim.M0Bind, 1) + ' J/g de liant');
 TAdd('Eau liee (non evaporable) : ' + RStr(Sim.Wn / Sim.M0Bind, 4) + ' g/g de liant');
 TAdd('Retrait chimique : ' + RStr(Sim.Shrk / Sim.M0Bind, 4) + ' mL/g de liant');
 TAdd('');
 TAdd('Porosite capillaire (eau) : ' + RStr(100.0 * PhCnt[P_POR] / PasteVox, 2) + ' %');
 TAdd('Pores vides (auto-dessiccation) : ' + RStr(100.0 * PhCnt[P_VOID] / PasteVox, 2) + ' %');
 If PPc Then TAdd('Pores capillaires : reseau connecte (' + RStr(100.0 * CntP / IMax(1, PhCnt[P_POR] + PhCnt[P_VOID]), 1) + ' % du volume poreux)')
 Else TAdd('Pores capillaires : reseau deconnecte');
 If Sim.TSet >= 0 Then TAdd('Prise (percolation des solides) : ' + RStr(Sim.TSet, 2) + ' h')
 Else TAdd('Prise (percolation des solides) : non atteinte');
 If Pc Then TAdd('Solides lies : percolation (' + RStr(100.0 * Cnt / PasteVox, 1) + ' % du volume)')
 Else TAdd('Solides lies : pas de percolation');
 X := GelSpace;
 TAdd('Rapport gel / espace de Powers : ' + RStr(X, 3));
 TAdd('Resistance a la compression (Powers) : ' + RStr(PowersStrength, 1) + ' MPa');
 TAdd('');
 TAdd('Phases perdues (volume non place) : ' + IStr(Sim.Lost) + ' voxels');
 TAdd('Modules, pores et conductivite : menu Proprietes.');
 TView('Resume des proprietes');
End;

Procedure ShowPercolation;
Var Cnt: LongInt; Pc: Boolean;
Begin
 If Not NeedGrid Then Exit;
 TClear;
 MaskSolids;
 Pc := Percolate(True, Cnt);
 QFree;
 TAdd('Solides lies par les hydrates (direction z) :');
 If Pc Then TAdd('   percolation : OUI  (' + RStr(100.0 * Cnt / PasteVox, 2) + ' % du volume de pate)')
 Else TAdd('   percolation : NON');
 If Sim.TSet >= 0 Then TAdd('   temps de prise : ' + RStr(Sim.TSet, 2) + ' h')
 Else TAdd('   temps de prise : non atteint');
 TAdd('');
 MaskPores;
 Pc := Percolate(True, Cnt);
 QFree;
 TAdd('Pores capillaires (eau + vides) :');
 If Pc Then TAdd('   percolation : OUI  (' + RStr(100.0 * Cnt / IMax(1, PhCnt[P_POR] + PhCnt[P_VOID]), 2) + ' % du volume poreux)')
 Else TAdd('   percolation : NON (reseau capillaire deconnecte)');
 TAdd('   porosite capillaire : ' + RStr(100.0 * (PhCnt[P_POR] + PhCnt[P_VOID]) / PasteVox, 2) + ' %');
 TView('Percolation');
End;

Procedure ShowPoreSizes;
Var I, S: Integer; Cum: LongInt; Tot: LongInt;
Begin
 If Not NeedGrid Then Exit;
 UserAbort := False; ProgPct := 0;
 If Not PoreSizeDist Then Begin
  ProgPct := -1;
  ShowMessage('Pores', 'Calcul impossible ou interrompu.');
  Exit;
 End;
 ProgPct := -1;
 Tot := 0;
 For I := 1 To MAXR Do Tot := Tot + PoreHist[I];
 If Tot = 0 Then Tot := 1;
 ChartReset('Distribution des tailles de pores capillaires', 'Diametre des pores (um)', 'Volume poreux cumule (%)', False);
 S := ChartAddSeries('Cumul croissant (boules inscrites)', False);
 Cum := 0;
 ChartPoint(S, 0, 0);
 For I := 1 To MAXR Do Begin
  Cum := Cum + PoreHist[I];
  ChartPoint(S, (2 * I - 1) * ResUm, 100.0 * Cum / Tot);
  If Cum >= Tot Then Break;
 End;
 ChFixY0 := True;
 ChartDialog;
End;

Procedure ShowElastic;
Begin
 If Not NeedGrid Then Exit;
 TClear;
 If ElasticModuli Then Begin
  TAdd('Reseau de solides lies traversant : ' + RStr(100 * ElPerc, 1) + ' % du volume');
  TAdd('Modele : milieu effectif auto-coherent (isotrope) applique aux');
  TAdd('phases percolantes ; pores et solides isoles : rigidite nulle.');
  TAdd('');
  TAdd('Module d''Young E       : ' + RStr(ElE, 2) + ' GPa');
  TAdd('Coefficient de Poisson : ' + RStr(ElNu, 3));
  TAdd('Module de compression K : ' + RStr(ElK, 2) + ' GPa');
  TAdd('Module de cisaillement G : ' + RStr(ElG, 2) + ' GPa');
  TAdd('Borne de Voigt (superieure) : ' + RStr(ElVoigt, 2) + ' GPa');
 End
 Else Begin
  TAdd('Les solides ne forment pas encore un reseau traversant :');
  TAdd('la rigidite de la pate est nulle (avant la prise).');
 End;
 TAdd('');
 TAdd('Modules des phases (GPa) : C3S 135, CH 38, C-S-H 22, AFt 24 ...');
 TView('Modules elastiques');
End;

Procedure StrengthDialog;
Begin
 If Not NeedGrid Then Exit;
 PReset;
 PAddH('Modele de Powers : fc = fc0 x X^n');
 PAddR('Resistance intrinseque fc0 (MPa)', Fc0, 10, 500, 1);
 PAddR('Exposant n', FcExp, 1, 6, 2);
 PAddH('Resultats');
 PAddT('Rapport gel / espace X', RStr(GelSpace, 4));
 PAddT('Resistance a la compression (MPa)', RStr(Fc0 * Pw(GelSpace, FcExp), 1));
 PAddT('Porosite capillaire (%)', RStr(100 * CapPorosity, 2));
 PAddT('Degre d''hydratation du ciment', RStr(AlphaCement, 3));
 If ParamDialog('Resistance a la compression', 'Modifier fc0 et n puis OK pour recalculer (Powers, ajustable)') Then Begin
  Fc0 := PVal[1]; FcExp := PVal[2];
  Modified := True;
  ShowMessage('Resistance a la compression', 'fc = ' + RStr(PowersStrength, 1) + ' MPa|X = ' + RStr(GelSpace, 3));
 End;
End;

Procedure ConductDialog;
Var D0: Real; Err: Boolean;
Begin
 If Not NeedGrid Then Exit;
 PReset;
 PAddR('Conductivite de la solution (S/m)', SigmaPore, 0.1, 30, 2);
 PAddR('Conductivite relative du C-S-H', SigTab[P_CSH], 0, 0.5, 4);
 PAddR('Conductivite relative des autres hydrates', SigTab[P_AFT], 0, 0.5, 5);
 PAddR('D0 de l''ion dans l''eau (1e-9 m2/s)', 2.03, 0.1, 10, 2);
 PAddH('Reseau de conducteurs : un potentiel est applique');
 If Not ParamDialog('Conductivite et diffusivite', 'Differences finies (SOR) ; Echap interrompt le calcul') Then Exit;
 SigmaPore := PVal[0];
 SigTab[P_CSH] := PVal[1]; SigTab[P_PCSH] := PVal[1];
 SigTab[P_AFT] := PVal[2]; SigTab[P_AFM] := PVal[2];
 D0 := PVal[3];
 UserAbort := False; ProgPct := 0;
 Err := Not Conductivity;
 ProgPct := -1;
 If Err Then Begin ShowMessage('Conductivite', 'Memoire insuffisante.'); Exit; End;
 TClear;
 TAdd('Direction du champ : ' + WordLine('x|y|z', CondAxis + 1) + '   iterations : ' + IStr(CondIters));
 If CondConv Then TAdd('Convergence : atteinte') Else TAdd('Convergence : partielle (arret anticipe)');
 TAdd('');
 TAdd('Conductivite relative (sigma / sigma_solution) : ' + RStr(CondBeta, 5));
 TAdd('Facteur de formation F = 1 / beta : ' + RStr(1 / RMax(CondBeta, 0.0000001), 1));
 TAdd('Conductivite electrique : ' + RStr(1000 * SigmaPore * CondBeta, 3) + ' mS/m');
 TAdd('Resistivite electrique : ' + RStr(1 / RMax(SigmaPore * CondBeta, 0.0000001), 1) + ' ohm.m');
 TAdd('Diffusivite relative D / D0 : ' + RStr(CondBeta, 5));
 TAdd('Coefficient de diffusion effectif : ' + RStr(D0 * CondBeta, 4) + ' e-9 m2/s');
 TAdd('');
 TAdd('Porosite capillaire : ' + RStr(100 * CapPorosity, 2) + ' %');
 TView('Conductivite / diffusivite');
End;

Procedure ShowITZ;
Var S1, S2, S3, S4, I: Integer;
Begin
 If Not NeedGrid Then Exit;
 If Not ItzProfile Then Begin
  ShowMessage('Zone de transition', 'Activez le mur de granulat dans Melange > Formulation|puis regenerez la microstructure.');
  Exit;
 End;
 ChartReset('Profil de la zone de transition interfaciale', 'Distance au granulat (um)', 'Fraction volumique (%)', False);
 S1 := ChartAddSeries('Porosite', False);
 S2 := ChartAddSeries('Particules non hydratees', False);
 S3 := ChartAddSeries('Portlandite', False);
 S4 := ChartAddSeries('C-S-H et autres hydrates', False);
 For I := 0 To ItzN - 1 Do Begin
  ChartPoint(S1, ItzD[I], ItzP[I]); ChartPoint(S2, ItzD[I], ItzU[I]);
  ChartPoint(S3, ItzD[I], ItzH[I]); ChartPoint(S4, ItzD[I], ItzG[I]);
 End;
 ChartDialog;
End;

{ courbes en fonction du temps }
Procedure TimeChart(Kind: Integer);
Var I, J1, J2, S, S2, S3, S4, S5: Integer; T, Rt: Real;
Begin
 If Not NeedGrid Then Exit;
 If Sim.NS < 2 Then Begin ShowMessage('Courbes', 'Lancez l''hydratation pour obtenir des courbes.'); Exit; End;
 Case Kind Of
  0: Begin
   ChartReset('Degre d''hydratation', 'Temps (h)', 'Degre d''hydratation', True);
   S := ChartAddSeries('Ciment (clinker)', False);
   S2 := ChartAddSeries('Liant (ciment + ajouts)', False);
   For I := 1 To Sim.NS - 1 Do Begin
    ChartPoint(S, Ser[SER_T, I], Ser[SER_ALF, I]);
    ChartPoint(S2, Ser[SER_T, I], Ser[SER_ALT, I]);
   End;
  End;
  1: Begin
   ChartReset('Chaleur d''hydratation cumulee', 'Temps (h)', 'Chaleur (J/g de liant)', True);
   S := ChartAddSeries('Chaleur cumulee', False);
   For I := 1 To Sim.NS - 1 Do ChartPoint(S, Ser[SER_T, I], Ser[SER_Q, I]);
  End;
  2: Begin
   ChartReset('Flux de chaleur (calorimetrie isotherme)', 'Temps (h)', 'Flux (mW/g de liant)', True);
   S := ChartAddSeries('Flux', False);
   For I := 2 To Sim.NS - 1 Do Begin
    J1 := I; J2 := I;
    While (J1 > 0) And (Ser[SER_T, J1] > Ser[SER_T, I] / 1.3) Do Dec(J1);
    While (J2 < Sim.NS - 1) And (Ser[SER_T, J2] < Ser[SER_T, I] * 1.3) Do Inc(J2);
    If (J2 > J1) And (Ser[SER_T, J2] > Ser[SER_T, J1]) And (Ser[SER_T, I] >= 0.1) Then Begin
     Rt := (Ser[SER_Q, J2] - Ser[SER_Q, J1]) / (Ser[SER_T, J2] - Ser[SER_T, J1]) / 3.6;
     ChartPoint(S, Ser[SER_T, I], Rt);
    End;
   End;
  End;
  3: Begin
   ChartReset('Porosite capillaire', 'Temps (h)', 'Porosite (%)', True);
   ChFixY0 := True;
   S := ChartAddSeries('Porosite capillaire', False);
   For I := 1 To Sim.NS - 1 Do ChartPoint(S, Ser[SER_T, I], Ser[SER_POR, I]);
  End;
  4: Begin
   ChartReset('Eau liee et retrait chimique', 'Temps (h)', 'g/g de liant  -  mL/g de liant', True);
   S := ChartAddSeries('Eau liee (g/g)', False);
   S2 := ChartAddSeries('Retrait chimique (mL/g)', False);
   For I := 1 To Sim.NS - 1 Do Begin
    ChartPoint(S, Ser[SER_T, I], Ser[SER_WN, I]);
    ChartPoint(S2, Ser[SER_T, I], Ser[SER_SHR, I]);
   End;
  End;
  5: Begin
   ChartReset('Resistance a la compression (Powers)', 'Temps (h)', 'Resistance (MPa)', True);
   S := ChartAddSeries('fc', False);
   For I := 1 To Sim.NS - 1 Do ChartPoint(S, Ser[SER_T, I], Ser[SER_FC, I]);
  End;
  6: Begin
   ChartReset('Temperature', 'Temps (h)', 'Temperature (C)', True);
   ChFixY0 := False;
   S := ChartAddSeries('Temperature', False);
   For I := 1 To Sim.NS - 1 Do ChartPoint(S, Ser[SER_T, I], Ser[SER_TMP, I]);
  End;
 Else Begin
  ChartReset('Evolution des phases', 'Temps (h)', 'Fraction volumique de pate (%)', True);
  S := ChartAddSeries('Ciment non hydrate', False);
  S2 := ChartAddSeries('C-S-H', False);
  S3 := ChartAddSeries('Portlandite', False);
  S4 := ChartAddSeries('Ettringite + monosulfate', False);
  S5 := ChartAddSeries('Porosite capillaire', False);
  For I := 1 To Sim.NS - 1 Do Begin
   T := Ser[SER_T, I];
   ChartPoint(S, T, 100.0 * (PhSer[P_C3S, I] + PhSer[P_C2S, I] + PhSer[P_C3A, I] + PhSer[P_C4AF, I]) / PasteVox);
   ChartPoint(S2, T, 100.0 * (PhSer[P_CSH, I] + PhSer[P_PCSH, I]) / PasteVox);
   ChartPoint(S3, T, 100.0 * PhSer[P_CH, I] / PasteVox);
   ChartPoint(S4, T, 100.0 * (PhSer[P_AFT, I] + PhSer[P_AFM, I]) / PasteVox);
   ChartPoint(S5, T, Ser[SER_POR, I]);
  End;
 End;
 End;
 ChartDialog;
End;

Procedure CurvesMenu;
Var Sel: Integer;
Begin
 If Not NeedGrid Then Exit;
 TClear;
 TAdd('Degre d''hydratation');
 TAdd('Chaleur d''hydratation cumulee');
 TAdd('Flux de chaleur');
 TAdd('Porosite capillaire');
 TAdd('Eau liee et retrait chimique');
 TAdd('Resistance a la compression');
 TAdd('Temperature');
 TAdd('Evolution des phases');
 Sel := ListDialog('Courbes en fonction du temps', True);
 If Sel >= 0 Then TimeChart(Sel);
End;

{ --------------------------------------------------------------------- }
{ Calibrage de la cinetique                                             }
{ --------------------------------------------------------------------- }

Const
 MAXEXP = 24;

Var
 ExpT, ExpY : Array[0..MAXEXP-1] Of Real;
 ExpN       : Integer;

Function ModelAt(Kind: Integer; Cyc: Real): Real;
Var I: Integer; F, Y0, Y1: Real; Col: Integer;
Begin
 If Kind = 0 Then Col := SER_ALF Else Col := SER_Q;
 If Sim.NS < 2 Then Begin ModelAt := 0; Exit; End;
 If Cyc >= Ser[SER_CYC, Sim.NS - 1] Then Begin ModelAt := Ser[Col, Sim.NS - 1]; Exit; End;
 If Cyc <= Ser[SER_CYC, 0] Then Begin ModelAt := Ser[Col, 0]; Exit; End;
 I := 1;
 While (I < Sim.NS - 1) And (Ser[SER_CYC, I] < Cyc) Do Inc(I);
 Y0 := Ser[Col, I - 1]; Y1 := Ser[Col, I];
 If Ser[SER_CYC, I] - Ser[SER_CYC, I - 1] < 0.5 Then ModelAt := Y1
 Else Begin
  F := (Cyc - Ser[SER_CYC, I - 1]) / (Ser[SER_CYC, I] - Ser[SER_CYC, I - 1]);
  ModelAt := Y0 + F * (Y1 - Y0);
 End;
End;

Function ParseExpData(Const S: String): Boolean;
Var Rest, Tok: String; P, Q: Integer; T, Y: Real;
Begin
 ParseExpData := False;
 ExpN := 0;
 Rest := S;
 While (Length(Rest) > 0) And (ExpN < MAXEXP) Do Begin
  P := Pos(';', Rest);
  If P > 0 Then Begin Tok := Copy(Rest, 1, P - 1); Delete(Rest, 1, P); End
  Else Begin Tok := Rest; Rest := ''; End;
  Tok := TrimStr(Tok);
  If Tok <> '' Then Begin
   Q := Pos(':', Tok);
   If Q = 0 Then Exit;
   If Not ParseReal(Copy(Tok, 1, Q - 1), T) Then Exit;
   If Not ParseReal(Copy(Tok, Q + 1, 255), Y) Then Exit;
   If T <= 0 Then Exit;
   ExpT[ExpN] := T; ExpY[ExpN] := Y;
   Inc(ExpN);
  End;
 End;
 ParseExpData := ExpN >= 2;
End;

Function CalibError(Kind: Integer; B: Real): Real;
Var I: Integer; E, S: Real;
Begin
 S := 0;
 For I := 0 To ExpN - 1 Do Begin
  E := ModelAt(Kind, Sqrt(ExpT[I] / B)) - ExpY[I];
  If Kind = 1 Then E := E / 100;
  S := S + E * E;
 End;
 CalibError := Sqrt(S / ExpN);
End;

Procedure CalibrateDialog;
Var
 Kind, I, S, S2 : Integer;
 Data : String;
 B, Best, BestB, E, Ratio : Real;
Begin
 If Not NeedGrid Then Exit;
 If Sim.NS < 10 Then Begin
  ShowMessage('Calibrage', 'Simulez d''abord l''hydratation (au moins quelques|jours) : le calibrage utilise la courbe simulee.');
  Exit;
 End;
 PReset;
 PAddE('Type de mesures', 0, 'Degre d''hydratation|Chaleur cumulee (J/g de liant)');
 PAddT('Format des donnees', 'temps(h):valeur ; temps(h):valeur ...');
 If Not ParamDialog('Calibrage de la cinetique', 'Ajuste beta pour que la courbe simulee rejoigne les mesures') Then Exit;
 Kind := PInt(0);
 If Kind = 0 Then Data := '6:0.12;24:0.45;72:0.62;168:0.72'
 Else Data := '6:60;24:240;72:320;168:360';
 If Not InputString('Calibrage : mesures', 'temps:valeur ; temps:valeur ... (temps en heures)', Data) Then Exit;
 If Not ParseExpData(Data) Then Begin
  ShowMessage('Calibrage', 'Donnees invalides : saisir au moins deux couples|temps:valeur separes par des points-virgules.');
  Exit;
 End;
 Best := 1E30; BestB := Beta;
 For I := 0 To 400 Do Begin
  B := 0.00001 * Exp(Ln(10000.0) * I / 400);
  E := CalibError(Kind, B);
  If E < Best Then Begin Best := E; BestB := B; End;
 End;
 ChartReset('Calibrage : simulation (beta ajuste) et mesures', 'Temps (h)', 'Valeur', True);
 S := ChartAddSeries('Simulation avec beta = ' + RStr(BestB, 6), False);
 S2 := ChartAddSeries('Mesures', True);

 For I := 1 To Sim.NS - 1 Do Begin
  If Kind = 0 Then ChartPoint(S, BestB * Sqr(Ser[SER_CYC, I]), Ser[SER_ALF, I])
  Else ChartPoint(S, BestB * Sqr(Ser[SER_CYC, I]), Ser[SER_Q, I]);
 End;
 For I := 0 To ExpN - 1 Do ChartPoint(S2, ExpT[I], ExpY[I]);
 ChartDialog;
 If ConfirmBox('Calibrage', 'beta ajuste = ' + RStr(BestB, 6) + ' (ecart ' + RStr(Best, 4) + ')|Appliquer cette valeur ?') Then Begin
  Ratio := BestB / Beta;
  Beta := BestB;
  For I := 0 To Sim.NS - 1 Do Ser[SER_T, I] := Ser[SER_T, I] * Ratio;
  Sim.T := Sim.T * Ratio; Sim.Te := Sim.Te * Ratio;
  If Sim.TSet >= 0 Then Sim.TSet := Sim.TSet * Ratio;
  Modified := True;
  StatusMsg := 'beta = ' + RStr(Beta, 6) + ' h/cycle^2';
 End;
End;

{ --------------------------------------------------------------------- }
{ Outils : calculateur de Bogue                                         }
{ --------------------------------------------------------------------- }

Procedure BogueDialog;
Var
 C, S, A, F, SO3, FL, C3S, C2S, C3A, C4AF, Gyp : Real;
Begin
 PReset;
 PAddH('Analyse chimique du ciment (% massique)');
 PAddR('CaO total', 64.5, 0, 100, 2);
 PAddR('SiO2', 21.0, 0, 100, 2);
 PAddR('Al2O3', 5.2, 0, 100, 2);
 PAddR('Fe2O3', 3.0, 0, 100, 2);
 PAddR('SO3', 2.8, 0, 20, 2);
 PAddR('CaO libre', 1.0, 0, 20, 2);
 If Not ParamDialog('Calculateur de Bogue', 'Calcul des phases potentielles du clinker (ASTM C150)') Then Exit;
 C := PVal[1] - PVal[6]; S := PVal[2]; A := PVal[3]; F := PVal[4]; SO3 := PVal[5]; FL := PVal[6];
 If (F > 0) And (A / F < 0.64) Then Begin
  C3A := 0;
  C4AF := 2.1 * A + 1.702 * F;
  C3S := 4.071 * C - 7.6 * S - 4.479 * A - 2.859 * F - 2.852 * SO3;
 End
 Else Begin
  C3A := 2.65 * A - 1.692 * F;
  C4AF := 3.043 * F;
  C3S := 4.071 * C - 7.6 * S - 6.718 * A - 1.43 * F - 2.852 * SO3;
 End;
 C2S := 2.867 * S - 0.7544 * C3S;
 Gyp := 2.151 * SO3;
 TClear;
 TAdd('C3S  (alite)     : ' + RStr(C3S, 2) + ' %');
 TAdd('C2S  (belite)    : ' + RStr(C2S, 2) + ' %');
 TAdd('C3A  (aluminate) : ' + RStr(C3A, 2) + ' %');
 TAdd('C4AF (ferrite)   : ' + RStr(C4AF, 2) + ' %');
 TAdd('Sulfate de calcium (gypse equivalent) : ' + RStr(Gyp, 2) + ' %');
 TAdd('Chaux libre : ' + RStr(FL, 2) + ' %');
 If (C3S < 0) Or (C2S < 0) Or (C3A < 0) Then
  TAdd('Attention : valeur negative, analyse incoherente avec Bogue.');
 TView('Resultat de Bogue');
 If (C3S >= 0) And (C2S >= 0) And (C3A >= 0) And (C4AF >= 0) Then
  If ConfirmBox('Bogue', 'Appliquer ces phases au ciment actuel ?') Then Begin
   Cem.W[1] := C3S; Cem.W[2] := C2S; Cem.W[3] := C3A; Cem.W[4] := C4AF;
   Cem.W[5] := 0.70 * Gyp; Cem.W[6] := 0.15 * Gyp; Cem.W[7] := 0.15 * Gyp;
   Cem.W[8] := FL; Cem.W[9] := 0;
   Cem.Name := 'Ciment (Bogue)';
   Modified := True;
   StatusMsg := 'Composition du ciment mise a jour (Bogue)';
  End;
End;

{ --------------------------------------------------------------------- }
{ Outils : empilement 3D de granulats                                   }
{ --------------------------------------------------------------------- }

Const
 MAXSPH  = 6000;
 MAXCELL = 8000;

Var
 SpX, SpY, SpZ, SpR : Array[0..MAXSPH-1] Of Real;
 SpNext   : Array[0..MAXSPH-1] Of Integer;
 CellHead : Array[0..MAXCELL-1] Of Integer;
 NSph     : Integer;
 AgM      : Integer;
 AgCs     : Real;
 AgL, AgDMax, AgDMin, AgQ, AgVa, AgItz : Real;
 AgSeed   : LongInt;
 AgPlaced, AgWanted : Integer;
 AgVolSph, AgVolItz : Real;

Function SphereOverlap(X, Y, Z, R: Real): Boolean;
Var CX, CY, CZ, DX, DY, DZ, I, MX, MY, MZ : Integer; DDX, DDY, DDZ, RS : Real;
Begin
 SphereOverlap := True;
 CX := Trunc(X / AgCs); CY := Trunc(Y / AgCs); CZ := Trunc(Z / AgCs);
 For DX := -1 To 1 Do
  For DY := -1 To 1 Do
   For DZ := -1 To 1 Do Begin
    MX := ((CX + DX) Mod AgM + AgM) Mod AgM;
    MY := ((CY + DY) Mod AgM + AgM) Mod AgM;
    MZ := ((CZ + DZ) Mod AgM + AgM) Mod AgM;
    I := CellHead[(MZ * AgM + MY) * AgM + MX];
    While I >= 0 Do Begin
     DDX := Abs(X - SpX[I]); If DDX > AgL / 2 Then DDX := AgL - DDX;
     DDY := Abs(Y - SpY[I]); If DDY > AgL / 2 Then DDY := AgL - DDY;
     DDZ := Abs(Z - SpZ[I]); If DDZ > AgL / 2 Then DDZ := AgL - DDZ;
     RS := R + SpR[I];
     If DDX * DDX + DDY * DDY + DDZ * DDZ < RS * RS Then Exit;
     I := SpNext[I];
    End;
   End;
 SphereOverlap := False;
End;

Procedure AggPack;
Var
 B : Array[0..40] Of Real;
 NCl, K, I, J, Tries, Cell, N0 : Integer;
 Sum, W, V, DMid, VS, X, Y, Z, T : Real;
 Ok : Boolean;
Begin
 SetSeed(AgSeed);
 NSph := 0; AgPlaced := 0; AgWanted := 0;
 AgM := Trunc(AgL / AgDMax);
 If AgM < 1 Then AgM := 1;
 If AgM > 20 Then AgM := 20;
 AgCs := AgL / AgM;
 For I := 0 To MAXCELL - 1 Do CellHead[I] := -1;
 B[0] := AgDMax; NCl := 0;
 While (B[NCl] > AgDMin * 1.0001) And (NCl < 38) Do Begin
  B[NCl + 1] := B[NCl] / Sqrt(2.0);
  If B[NCl + 1] < AgDMin Then B[NCl + 1] := AgDMin;
  Inc(NCl);
 End;
 Sum := 1 - Pw(AgDMin / AgDMax, AgQ);
 If Sum < 0.0001 Then Sum := 0.0001;
 For K := 0 To NCl - 1 Do Begin
  W := (Pw(B[K] / AgDMax, AgQ) - Pw(B[K + 1] / AgDMax, AgQ)) / Sum;
  V := W * AgVa * AgL * AgL * AgL;
  DMid := Sqrt(B[K] * B[K + 1]);
  VS := Pi / 6 * DMid * DMid * DMid;
  N0 := Trunc(V / VS + RndF);
  For J := 1 To N0 Do Begin
   Inc(AgWanted);
   If NSph >= MAXSPH Then Continue;
   Ok := False;
   For Tries := 1 To 250 Do Begin
    X := RndF * AgL; Y := RndF * AgL; Z := RndF * AgL;
    If Not SphereOverlap(X, Y, Z, DMid / 2) Then Begin Ok := True; Break; End;
   End;
   If Ok Then Begin
    SpX[NSph] := X; SpY[NSph] := Y; SpZ[NSph] := Z; SpR[NSph] := DMid / 2;
    Cell := (Trunc(Z / AgCs) Mod AgM * AgM + Trunc(Y / AgCs) Mod AgM) * AgM + Trunc(X / AgCs) Mod AgM;
    SpNext[NSph] := CellHead[Cell];
    CellHead[Cell] := NSph;
    Inc(NSph); Inc(AgPlaced);
   End;
  End;
  ShowProgress('Empilement des granulats...', (100 * (K + 1)) Div NCl);
 End;
 AgVolSph := 0; AgVolItz := 0;
 T := AgItz / 1000.0;
 For I := 0 To NSph - 1 Do Begin
  AgVolSph := AgVolSph + 4.0 / 3.0 * Pi * Sqr(SpR[I]) * SpR[I];
  AgVolItz := AgVolItz + 4.0 / 3.0 * Pi * (Sqr(SpR[I] + T) * (SpR[I] + T) - Sqr(SpR[I]) * SpR[I]);
 End;
End;

Procedure AggDraw(ZIdx: Integer);
Var CX0, CY0, CS, I, XP, YP, RP : Integer; Z, DZ, R2, Sc : Real;
Begin
 CX0 := 28; CY0 := 40; CS := 400;
 Sc := CS / AgL;
 Z := (ZIdx + 0.5) / 100 * AgL;
 SetFillStyle(SolidFill, White);
 Bar(CX0, CY0, CX0 + CS - 1, CY0 + CS - 1);
 SetViewPort(CX0, CY0, CX0 + CS - 1, CY0 + CS - 1, True);
 For I := 0 To NSph - 1 Do Begin
  DZ := Abs(Z - SpZ[I]); If DZ > AgL / 2 Then DZ := AgL - DZ;
  If DZ < SpR[I] Then Begin
   R2 := Sqr(SpR[I]) - Sqr(DZ);
   RP := Round(Sqrt(R2) * Sc);
   XP := Round(SpX[I] * Sc); YP := Round(SpY[I] * Sc);
   If RP < 1 Then RP := 1;
   SetFillStyle(SolidFill, DarkGray);
   SetColor(DarkGray);
   FillEllipse(XP, YP, RP, RP);
   { copies periodiques }
   If XP - RP < 0 Then FillEllipse(XP + CS, YP, RP, RP);
   If XP + RP >= CS Then FillEllipse(XP - CS, YP, RP, RP);
   If YP - RP < 0 Then FillEllipse(XP, YP + CS, RP, RP);
   If YP + RP >= CS Then FillEllipse(XP, YP - CS, RP, RP);
  End;
 End;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
 SetColor(Black);
 Rectangle(CX0 - 1, CY0 - 1, CX0 + CS, CY0 + CS);
End;

Procedure AggDialog3D;
Const BW = 100; BH = 22;
Var
 X1, Y1, X2, Y2, BY1, BY2, B1, B2, B3, ZIdx, Key, MX, MY : Integer;
 Click, Done, Redraw, Regen : Boolean;
 VT, Paste, Itz : Real;
 Procedure DrawAll;
 Var TY: Integer;
 Begin
  HideMouseCursor;
  DrawDialogBox(X1, Y1, X2, Y2, 'Empilement 3D de granulats (spheres, conditions periodiques)');
  AggDraw(ZIdx);
  SetTextStyle(DefaultFont, HorizDir, 1);
  SetColor(Black);
  TY := Y1 + 40;
  OutTextXY(444, TY, 'Cube : ' + RTrim(AgL, 1) + ' mm'); TY := TY + 14;
  OutTextXY(444, TY, 'd max : ' + RTrim(AgDMax, 2) + ' mm'); TY := TY + 12;
  OutTextXY(444, TY, 'd min : ' + RTrim(AgDMin, 2) + ' mm'); TY := TY + 12;
  OutTextXY(444, TY, 'Fuller q : ' + RTrim(AgQ, 2)); TY := TY + 18;
  OutTextXY(444, TY, 'Spheres : ' + IStr(AgPlaced) + ' / ' + IStr(AgWanted)); TY := TY + 14;
  OutTextXY(444, TY, 'Fraction visee : ' + RStr(100 * AgVa, 1) + ' %'); TY := TY + 12;
  VT := AgVolSph / (AgL * AgL * AgL);
  OutTextXY(444, TY, 'Obtenue : ' + RStr(100 * VT, 1) + ' %'); TY := TY + 18;
  Paste := AgL * AgL * AgL - AgVolSph;
  If Paste > 0 Then Itz := AgVolItz / Paste Else Itz := 0;
  If Itz > 1 Then Itz := 1;
  OutTextXY(444, TY, 'ZTI : ' + RTrim(AgItz, 1) + ' um'); TY := TY + 12;
  OutTextXY(444, TY, 'Volume de ZTI dans'); TY := TY + 12;
  OutTextXY(444, TY, 'la pate : ' + RStr(100 * Itz, 1) + ' %'); TY := TY + 18;
  SetColor(DarkGray);
  OutTextXY(444, TY, 'Coupe ' + IStr(ZIdx + 1) + ' / 100'); TY := TY + 12;
  OutTextXY(444, TY, 'Haut/Bas : +-1'); TY := TY + 12;
  OutTextXY(444, TY, 'PgUp/PgDn : +-10');
  DrawButton(B1, BY1, B1 + BW, BY2, 'Parametres', False);
  DrawButton(B2, BY1, B2 + BW, BY2, 'Regenerer', False);
  DrawButton(B3, BY1, B3 + BW, BY2, 'Fermer', True);
  ShowMouseCursor;
 End;
Begin
 X1 := 12; Y1 := 12; X2 := ScreenWidth - 12; Y2 := ScreenHeight - 12;
 BY2 := Y2 - 8; BY1 := BY2 - BH;
 B3 := X2 - BW - 10; B2 := B3 - BW - 8; B1 := B2 - BW - 8;
 ZIdx := 50; Done := False; Redraw := True;
 UserAbort := False;
 AggPack;
 ProgPct := -1; StatusMsg := '';
 DlgBegin;
 While Not Done Do Begin
  If Redraw Then Begin DrawAll; Redraw := False; End;
  If DlgPoll(Key, MX, MY, Click) Then Begin
   Regen := False;
   Case Key Of
    256 + 72: Begin If ZIdx > 0 Then Dec(ZIdx); Redraw := True; End;
    256 + 80: Begin If ZIdx < 99 Then Inc(ZIdx); Redraw := True; End;
    256 + 73: Begin ZIdx := IMax(0, ZIdx - 10); Redraw := True; End;
    256 + 81: Begin ZIdx := IMin(99, ZIdx + 10); Redraw := True; End;
    27, 13: Done := True;
    Ord('r'), Ord('R'): Regen := True;
   End;
   If Click Then Begin
    If PtInRect(MX, MY, B3, BY1, B3 + BW, BY2) Then Done := True
    Else If PtInRect(MX, MY, B2, BY1, B2 + BW, BY2) Then Regen := True
    Else If PtInRect(MX, MY, B1, BY1, B1 + BW, BY2) Then Begin
     PReset;
     PAddR('Cote du cube (mm)', AgL, 10, 200, 1);
     PAddR('Diametre maximal (mm)', AgDMax, 2, 40, 2);
     PAddR('Diametre minimal (mm)', AgDMin, 0.5, 20, 2);
     PAddR('Exposant de Fuller q', AgQ, 0.2, 1, 2);
     PAddR('Fraction volumique de granulats (%)', AgVa * 100, 5, 70, 1);
     PAddR('Epaisseur de la ZTI (um)', AgItz, 5, 100, 1);
     PAddI('Graine aleatoire', AgSeed, 1, 2000000000);
     If ParamDialog('Empilement de granulats : parametres', 'Les petites classes sont limitees a 6000 spheres') Then Begin
      AgL := PVal[0]; AgDMax := PVal[1]; AgDMin := PVal[2]; AgQ := PVal[3];
      AgVa := PVal[4] / 100; AgItz := PVal[5]; AgSeed := PInt(6);
      If AgDMin >= AgDMax Then AgDMin := AgDMax / 4;
      Regen := True;
     End;
     Redraw := True;
    End;
   End;
   If Regen Then Begin
    HideMouseCursor;
    SetFillStyle(SolidFill, LightGray);
    Bar(X1 + 4, Y1 + 20, X2 - 4, Y2 - 40);
    ShowMouseCursor;
    AggPack;
    ProgPct := -1; StatusMsg := '';
    Redraw := True;
   End;
  End;
  Delay(10);
 End;
 DlgEnd;
End;

Procedure InitAggParams;
Begin
 AgL := 60; AgDMax := Coarse.DMax; AgDMin := Coarse.DMin; AgQ := Coarse.Q;
 AgVa := 0.40; AgItz := 20; AgSeed := 4321;
End;

{ --------------------------------------------------------------------- }
{ Fichiers : projet, microstructure, exports                            }
{ --------------------------------------------------------------------- }

Var IOErr: Boolean;

Procedure IOR(Var F: Text; Wrt: Boolean; Const Key: String; Var V: Real);
Var L: String; P: Integer; R: Real;
Begin
 If Wrt Then WriteLn(F, Key, '=', V:0:8)
 Else Begin
  If Eof(F) Then Begin IOErr := True; Exit; End;
  ReadLn(F, L);
  P := Pos('=', L);
  If (P = 0) Or Not ParseReal(Copy(L, P + 1, 255), R) Then IOErr := True Else V := R;
 End;
End;

Procedure IOL(Var F: Text; Wrt: Boolean; Const Key: String; Var V: LongInt);
Var R: Real;
Begin
 R := V;
 IOR(F, Wrt, Key, R);
 If Not Wrt Then V := Round(R);
End;

Procedure IOI(Var F: Text; Wrt: Boolean; Const Key: String; Var V: Integer);
Var R: Real;
Begin
 R := V;
 IOR(F, Wrt, Key, R);
 If Not Wrt Then V := Round(R);
End;

Procedure IOS(Var F: Text; Wrt: Boolean; Const Key: String; Var V: String);
Var L: String; P: Integer;
Begin
 If Wrt Then WriteLn(F, Key, '=', V)
 Else Begin
  If Eof(F) Then Begin IOErr := True; Exit; End;
  ReadLn(F, L);
  P := Pos('=', L);
  If P = 0 Then IOErr := True Else V := Copy(L, P + 1, 255);
 End;
End;

Procedure IOB(Var F: Text; Wrt: Boolean; Const Key: String; Var V: Boolean);
Var R: Real;
Begin
 If V Then R := 1 Else R := 0;
 IOR(F, Wrt, Key, R);
 If Not Wrt Then V := R > 0.5;
End;

Procedure IOPsd(Var F: Text; Wrt: Boolean; Const Pre: String; Var P: TPSD);
Begin
 IOR(F, Wrt, Pre + '.d63', P.D63); IOR(F, Wrt, Pre + '.n', P.Nr);
 IOR(F, Wrt, Pre + '.dmin', P.DMin); IOR(F, Wrt, Pre + '.dmax', P.DMax);
End;

Procedure IOScm(Var F: Text; Wrt: Boolean; Const Pre: String; Var M: TSCM);
Var NM: String;
Begin
 NM := M.Name;
 IOS(F, Wrt, Pre + '.nom', NM); M.Name := Copy(NM, 1, 30);
 IOR(F, Wrt, Pre + '.rho', M.Rho);
 IOPsd(F, Wrt, Pre, M.P);
 IOR(F, Wrt, Pre + '.reactivite', M.React);
 IOR(F, Wrt, Pre + '.extra', M.Extra);
End;

Procedure IOAgg(Var F: Text; Wrt: Boolean; Const Pre: String; Var A: TAggr);
Begin
 IOR(F, Wrt, Pre + '.rho', A.Rho); IOR(F, Wrt, Pre + '.dmin', A.DMin);
 IOR(F, Wrt, Pre + '.dmax', A.DMax); IOR(F, Wrt, Pre + '.q', A.Q);
End;

Procedure ProjIO(Var F: Text; Wrt: Boolean);
Var I: Integer; NM: String;
Begin
 NM := Cem.Name;
 IOS(F, Wrt, 'ciment.nom', NM); Cem.Name := Copy(NM, 1, 30);
 For I := 1 To 9 Do IOR(F, Wrt, 'ciment.phase' + IStr(I), Cem.W[I]);
 IOPsd(F, Wrt, 'ciment', Cem.P);
 IOScm(F, Wrt, 'cendre', FAm);
 IOScm(F, Wrt, 'laitier', SlgM);
 IOScm(F, Wrt, 'calcaire', FilM);
 IOScm(F, Wrt, 'fumee', SFM);
 IOAgg(F, Wrt, 'sable', Sand);
 IOAgg(F, Wrt, 'gravier', Coarse);
 IOI(F, Wrt, 'melange.type', MixKind);
 IOR(F, Wrt, 'melange.el', WB);
 For I := 0 To 4 Do IOR(F, Wrt, 'melange.liant' + IStr(I), Frac[I]);
 IOI(F, Wrt, 'melange.voxels', NVox);
 IOR(F, Wrt, 'melange.resolution', ResUm);
 IOL(F, Wrt, 'melange.graine', Seed);
 IOB(F, Wrt, 'melange.mur', WallOn);
 IOR(F, Wrt, 'beton.liant', ConcBind); IOR(F, Wrt, 'beton.sable', ConcSand);
 IOR(F, Wrt, 'beton.gravier', ConcCoarse); IOR(F, Wrt, 'beton.air', ConcAir);
 IOR(F, Wrt, 'hydr.temperature', CureTemp);
 IOI(F, Wrt, 'hydr.mode', CureMode);
 IOR(F, Wrt, 'hydr.beta', Beta);
 IOR(F, Wrt, 'hydr.ea', EaCem);
 IOR(F, Wrt, 'hydr.tref', TRefC);
 IOR(F, Wrt, 'hydr.jours', EndDays);
 IOL(F, Wrt, 'hydr.cycles', MaxCycles);
 IOR(F, Wrt, 'hydr.nucch', PNucCH);
 IOI(F, Wrt, 'hydr.age', MaxAge);
 IOI(F, Wrt, 'hydr.dormant', DormCyc);
 IOR(F, Wrt, 'hydr.proche', AccNear);
 IOR(F, Wrt, 'hydr.eloigne', AccFar);
 For I := 1 To 13 Do IOR(F, Wrt, 'vitesse.' + PhShort[I], Rate[I]);
 IOR(F, Wrt, 'resistance.fc0', Fc0);
 IOR(F, Wrt, 'resistance.n', FcExp);
 IOR(F, Wrt, 'conduc.sigma', SigmaPore);
 IOR(F, Wrt, 'conduc.csh', SigTab[P_CSH]);
 IOR(F, Wrt, 'conduc.autres', SigTab[P_AFT]);
End;

Procedure FreeMicro;
Begin
 HaveGrid := False;
 FillChar(Sim, SizeOf(Sim), 0);
 Sim.TSet := -1;
 FillChar(PhCnt, SizeOf(PhCnt), 0);
 SimActive := False; SimPaused := False;
End;

Function SaveProject(Const FN: String): Boolean;
Var F: Text;
Begin
 SaveProject := False;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire ' + FN); Exit; End;
 WriteLn(F, 'VCCTL-PROJET ', PROG_VERSION);
 IOErr := False;
 ProjIO(F, True);
 Close(F);
 ProjName := FN; Modified := False;
 StatusMsg := 'Projet enregistre : ' + FN;
 SaveProject := True;
End;

Function LoadProject(Const FN: String): Boolean;
Var F: Text; L: String;
Begin
 LoadProject := False;
 {$I-}
 Assign(F, FN); Reset(F);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible de lire ' + FN); Exit; End;
 ReadLn(F, L);
 If Copy(L, 1, 12) <> 'VCCTL-PROJET' Then Begin
  Close(F);
  ShowMessage('Erreur', 'Ce fichier n''est pas un projet VCCTL.');
  Exit;
 End;
 IOErr := False;
 ProjIO(F, False);
 Close(F);
 If IOErr Then ShowMessage('Projet', 'Fichier incomplet : certaines valeurs sont inchangees.');
 If NVox < 16 Then NVox := 16;
 If NVox > MAXN Then NVox := MAXN;
 SyncMatDensities;
 FreeMicro;
 ProjName := FN; Modified := False;
 StatusMsg := 'Projet charge : ' + FN;
 LoadProject := True;
End;

Type
 THdr = Record
  Magic : Array[0..3] Of Char;
  Ver, NN0, SzSim, SzSer, SzPh : Integer;
  Wall : Boolean;
  Res : Real;
  Mode : Integer;
 End;

Function SaveMicro(Const FN: String): Boolean;
Var F: File; H: THdr; Z: Integer;
Begin
 SaveMicro := False;
 If Not HaveGrid Then Exit;
 {$I-}
 Assign(F, FN); Rewrite(F, 1);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire ' + FN); Exit; End;
 H.Magic[0] := 'V'; H.Magic[1] := 'M'; H.Magic[2] := 'I'; H.Magic[3] := 'C';
 H.Ver := 1; H.NN0 := N; H.SzSim := SizeOf(Sim); H.SzSer := SizeOf(Ser); H.SzPh := SizeOf(PhSer);
 H.Wall := WallOn; H.Res := ResUm; H.Mode := CureMode;
 BlockWrite(F, H, SizeOf(H));
 BlockWrite(F, Sim, SizeOf(Sim));
 BlockWrite(F, Ser, SizeOf(Ser));
 BlockWrite(F, PhSer, SizeOf(PhSer));
 For Z := 0 To N - 1 Do BlockWrite(F, Grid[Z]^, NN);
 Close(F);
 StatusMsg := 'Microstructure enregistree : ' + FN;
 SaveMicro := True;
End;

Function LoadMicro(Const FN: String): Boolean;
Var F: File; H: THdr; Z: Integer; Got: Word;
Begin
 LoadMicro := False;
 {$I-}
 Assign(F, FN); Reset(F, 1);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible de lire ' + FN); Exit; End;
 BlockRead(F, H, SizeOf(H), Got);
 If (Got <> SizeOf(H)) Or (H.Magic[0] <> 'V') Or (H.Magic[1] <> 'M') Or (H.Magic[2] <> 'I') Or
    (H.Ver <> 1) Or (H.NN0 < 16) Or (H.NN0 > MAXN) Or (H.SzSim <> SizeOf(Sim)) Or
    (H.SzSer <> SizeOf(Ser)) Or (H.SzPh <> SizeOf(PhSer)) Then Begin
  Close(F);
  ShowMessage('Erreur', 'Fichier de microstructure invalide ou incompatible.');
  Exit;
 End;
 N := H.NN0; NN := N * N; NTot := LongInt(N) * NN; NVox := N;
 WallOn := H.Wall; ResUm := H.Res; CureMode := H.Mode;
 If Not AllocGridG(Grid, GridN, N) Then Begin Close(F); ShowMessage('Erreur', 'Memoire insuffisante.'); Exit; End;
 BlockRead(F, Sim, SizeOf(Sim));
 BlockRead(F, Ser, SizeOf(Ser));
 BlockRead(F, PhSer, SizeOf(PhSer));
 For Z := 0 To N - 1 Do BlockRead(F, Grid[Z]^, NN);
 Close(F);
 InitWrap;
 ZLo := 0; ZHi := N - 1;
 If WallOn Then Begin ZLo := 1; ZHi := N - 2; End;
 SyncMatDensities;
 InitReactions;
 CountPhases;
 BuildSulfList;
 CcapJ := 0.75 * Sim.MSol0 + 4.18 * Sim.MW0;
 HaveGrid := True;
 SimActive := False; SimPaused := Sim.Cyc > 0;
 SliceIdx := N Div 2;
 StatusMsg := 'Microstructure chargee : ' + FN;
 LoadMicro := True;
End;

Procedure ExportMicroText(Const FN: String);
Var F: Text; Z, I: Integer;
Begin
 If Not HaveGrid Then Exit;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire ' + FN); Exit; End;
 WriteLn(F, 'VCCTL-CLONE-IMG ', PROG_VERSION);
 WriteLn(F, 'X_Size: ', N);
 WriteLn(F, 'Y_Size: ', N);
 WriteLn(F, 'Z_Size: ', N);
 WriteLn(F, 'Image_Resolution: ', ResUm:0:3);
 WriteLn(F, 'Time_h: ', Sim.T:0:3);
 WriteLn(F, 'Phases: ', NPH);
 For I := 0 To NPH - 1 Do WriteLn(F, '# ', I, ' = ', PhName[I]);
 For Z := 0 To N - 1 Do
  For I := 0 To NN - 1 Do WriteLn(F, Grid[Z]^[I]);
 Close(F);
 StatusMsg := 'Microstructure exportee : ' + FN;
End;

Procedure ExportResultsCSV(Const FN: String);
Var F: Text; I: Integer;
Begin
 If Sim.NS < 1 Then Begin ShowMessage('Resultats', 'Aucun resultat a exporter.'); Exit; End;
 {$I-}
 Assign(F, FN); Rewrite(F);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire ' + FN); Exit; End;
 WriteLn(F, 'temps_h;cycles;alpha_ciment;alpha_liant;chaleur_J_g;temperature_C;porosite_cap_pct;eau_liee_g_g;retrait_mL_g;resistance_MPa;C-S-H_pct;CH_pct;AFt_pct;AFm_pct');
 For I := 0 To Sim.NS - 1 Do
  WriteLn(F, Ser[SER_T, I]:0:4, ';', Ser[SER_CYC, I]:0:0, ';', Ser[SER_ALF, I]:0:5, ';',
          Ser[SER_ALT, I]:0:5, ';', Ser[SER_Q, I]:0:3, ';', Ser[SER_TMP, I]:0:2, ';',
          Ser[SER_POR, I]:0:3, ';', Ser[SER_WN, I]:0:5, ';', Ser[SER_SHR, I]:0:5, ';',
          Ser[SER_FC, I]:0:3, ';',
          100.0 * (PhSer[P_CSH, I] + PhSer[P_PCSH, I]) / PasteVox:0:3, ';',
          100.0 * PhSer[P_CH, I] / PasteVox:0:3, ';',
          100.0 * PhSer[P_AFT, I] / PasteVox:0:3, ';',
          100.0 * PhSer[P_AFM, I] / PasteVox:0:3);
 Close(F);
 StatusMsg := 'Resultats exportes : ' + FN;
End;

Const
 VgaRGB : Array[0..15, 0..2] Of Byte = (
  (0, 0, 0), (0, 0, 170), (0, 170, 0), (0, 170, 170), (170, 0, 0), (170, 0, 170),
  (170, 85, 0), (170, 170, 170), (85, 85, 85), (85, 85, 255), (85, 255, 85),
  (85, 255, 255), (255, 85, 85), (255, 85, 255), (255, 255, 85), (255, 255, 255));

Procedure ExportSliceBMP(Const FN: String);
Var
 F: File; S, W, RowBytes, X, Y, Plane: Integer; Row: Array[0..1599] Of Byte;
 Hd: Array[0..53] Of Byte; Sz: LongInt; C: Byte;
 Procedure PutL(Off: Integer; V: LongInt);
 Begin Hd[Off] := V And 255; Hd[Off + 1] := (V Shr 8) And 255; Hd[Off + 2] := (V Shr 16) And 255; Hd[Off + 3] := (V Shr 24) And 255; End;
Begin
 If Not HaveGrid Then Exit;
 S := 400 Div N;
 If S < 1 Then S := 1;
 If S > 4 Then S := 4;
 W := N * S;
 RowBytes := ((W * 3 + 3) Div 4) * 4;
 Plane := ViewPlane;
 If Plane > 2 Then Plane := 0;
 {$I-}
 Assign(F, FN); Rewrite(F, 1);
 {$I+}
 If IOResult <> 0 Then Begin ShowMessage('Erreur', 'Impossible d''ecrire ' + FN); Exit; End;
 FillChar(Hd, SizeOf(Hd), 0);
 Hd[0] := Ord('B'); Hd[1] := Ord('M');
 Sz := 54 + LongInt(RowBytes) * W;
 PutL(2, Sz); PutL(10, 54); PutL(14, 40); PutL(18, W); PutL(22, W);
 Hd[26] := 1; Hd[28] := 24; PutL(34, LongInt(RowBytes) * W);
 BlockWrite(F, Hd, 54);
 For Y := W - 1 Downto 0 Do Begin
  FillChar(Row, SizeOf(Row), 0);
  For X := 0 To W - 1 Do Begin
   C := PhDispColor(VoxAt(Plane, X Div S, Y Div S));
   Row[X * 3] := VgaRGB[C, 2]; Row[X * 3 + 1] := VgaRGB[C, 1]; Row[X * 3 + 2] := VgaRGB[C, 0];
  End;
  BlockWrite(F, Row, RowBytes);
 End;
 Close(F);
 StatusMsg := 'Coupe exportee : ' + FN;
End;

{ --------------------------------------------------------------------- }
{ Aide                                                                  }
{ --------------------------------------------------------------------- }

Procedure ShowHelp;
Begin
 TClear;
 TAdd('VCCTL Clone - Laboratoire virtuel d''essais du ciment et du beton');
 TAdd('');
 TAdd('Marche a suivre');
 TAdd(' 1. Materiaux : choisir le ciment (bibliotheque ou composition libre),');
 TAdd('    les ajouts (cendre volante, laitier, calcaire, fumee de silice)');
 TAdd('    et les granulats. Granulometrie de type Rosin-Rammler.');
 TAdd(' 2. Melange > Formulation : rapport E/L, proportions du liant, taille');
 TAdd('    du systeme (voxels), resolution, graine, mur de granulat.');
 TAdd(' 3. Melange > Generer la microstructure (Ctrl+G) : empilement aleatoire');
 TAdd('    de spheres numerisees en voxels (periodique en x, y, z).');
 TAdd(' 4. Hydratation > Demarrer (F9) : automate cellulaire. A chaque cycle,');
 TAdd('    les voxels de surface se dissolvent avec une probabilite propre a');
 TAdd('    chaque phase ; les produits (C-S-H, portlandite, ettringite,');
 TAdd('    monosulfate, hydrogrenat, gel pouzzolanique) se deposent par marche');
 TAdd('    aleatoire. Temps reel = beta x cycle^2 (loi parabolique).');
 TAdd(' 5. Proprietes : degre d''hydratation, chaleur, retrait chimique, eau');
 TAdd('    liee, porosite, percolation, pores, modules, resistance,');
 TAdd('    conductivite et diffusivite, zone de transition.');
 TAdd('');
 TAdd('Modele chimique (volumes par voxel de reactif)');
 TAdd(' C3S + 5,3 H -> C-S-H + 1,3 CH        C2S + 4,3 H -> C-S-H + 0,3 CH');
 TAdd(' C3A + 3 gypse + 26 H -> ettringite   C3A + gypse + 10 H -> monosulfate');
 TAdd(' C3A + 6 H -> hydrogrenat             C4AF + 2 CH + 10 H -> hydrogrenat');
 TAdd(' CaO + H -> CH                        SiO2 + 1,1 CH + 2,9 H -> C-S-H');
 TAdd(' Les sulfates sont consommes a la demande ; le retrait chimique cree des');
 TAdd(' pores vides (cure scellee) ou est comble par l''eau externe (cure saturee).');
 TAdd('');
 TAdd('Raccourcis');
 TAdd(' F1 aide  F2 formulation  F3 parametres d''hydratation  F4 resume');
 TAdd(' F5 courbes  F9 demarrer  Espace pause/reprise  Echap pause');
 TAdd(' Ctrl+G generer  Ctrl+N nouveau  Ctrl+O ouvrir  Ctrl+S enregistrer');
 TAdd(' PgUp / PgDn  coupe precedente / suivante  F10 menus  Alt+X quitter');
 TAdd('');
 TAdd('Mode lot : VCCTL /B [projet.VCT] [sortie.CSV]');
 TAdd('');
 TAdd('Limites : modele simplifie et pedagogique ; les valeurs numeriques sont');
 TAdd('indicatives et doivent etre calibrees (Hydratation > Calibrer).');
 TView('Aide de VCCTL');
End;

Procedure ShowPhaseLegend;
Var I: Integer;
Begin
 TClear;
 TAdd('Phase                       Couleur            Densite');
 TAdd('-----------------------------------------------------');
 For I := 0 To NPH - 1 Do
  TAdd(PadR(PhName[I], 26) + PadR('couleur VGA ' + IStr(PhColor[I]), 18) + PadL(RStr(PhRho[I], 2), 6));
 TView('Phases et couleurs');
End;

Procedure ShowAbout;
Begin
 ShowMessage('A propos', PROG_TITLE + '|Clone pedagogique du logiciel VCCTL (NIST)|Microstructure 3D, hydratation par automate cellulaire|Auteur : Sylvain Maltais - www.gladir.com/qchimie');
End;

Procedure PhaseVisDialog;
Var I: Integer;
Begin
 PReset;
 For I := 0 To NPH - 1 Do PAddB(PhName[I], PhShow[I]);
 If ParamDialog('Phases visibles', 'Espace / Entree : afficher ou masquer la phase') Then
  For I := 0 To NPH - 1 Do PhShow[I] := PBool(I);
End;

{ --------------------------------------------------------------------- }
{ Definition des menus                                                  }
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

Procedure MS(Var Mn: TMenu; Idx: Integer);
Begin
 Mn.Items[Idx].Separator := True;
End;

Procedure InitMenus;
Begin
 { Fichier }
 ZeroMenu(Menus[MENU_FILE]);
 Menus[MENU_FILE].Title := 'Fichier';
 Menus[MENU_FILE].ItemCount := 12;
 MI(Menus[MENU_FILE], 0, 'Nouveau projet', 'Ctrl+N');
 MI(Menus[MENU_FILE], 1, 'Ouvrir un projet...', 'Ctrl+O');
 MI(Menus[MENU_FILE], 2, 'Enregistrer le projet', 'Ctrl+S');
 MI(Menus[MENU_FILE], 3, 'Enregistrer sous...', '');
 MS(Menus[MENU_FILE], 4);
 MI(Menus[MENU_FILE], 5, 'Ouvrir une microstructure...', '');
 MI(Menus[MENU_FILE], 6, 'Enreg. la microstructure...', '');
 MI(Menus[MENU_FILE], 7, 'Exporter la microstr. (texte)', '');
 MI(Menus[MENU_FILE], 8, 'Exporter la coupe (BMP)...', '');
 MI(Menus[MENU_FILE], 9, 'Exporter les resultats (CSV)', '');
 MS(Menus[MENU_FILE], 10);
 MI(Menus[MENU_FILE], 11, 'Quitter', 'Alt+X');
 { Materiaux }
 ZeroMenu(Menus[MENU_MAT]);
 Menus[MENU_MAT].Title := 'Materiaux';
 Menus[MENU_MAT].ItemCount := 11;
 MI(Menus[MENU_MAT], 0, 'Ciment...', '');
 MI(Menus[MENU_MAT], 1, 'Cendres volantes...', '');
 MI(Menus[MENU_MAT], 2, 'Laitier...', '');
 MI(Menus[MENU_MAT], 3, 'Calcaire / filler inerte...', '');
 MI(Menus[MENU_MAT], 4, 'Fumee de silice...', '');
 MS(Menus[MENU_MAT], 5);
 MI(Menus[MENU_MAT], 6, 'Granulats fins (sable)...', '');
 MI(Menus[MENU_MAT], 7, 'Granulats grossiers...', '');
 MS(Menus[MENU_MAT], 8);
 MI(Menus[MENU_MAT], 9, 'Caracteristiques des poudres', '');
 MI(Menus[MENU_MAT], 10, 'Courbes granulometriques', '');
 { Melange }
 ZeroMenu(Menus[MENU_MIX]);
 Menus[MENU_MIX].Title := 'Melange';
 Menus[MENU_MIX].ItemCount := 8;
 MI(Menus[MENU_MIX], 0, 'Formulation...', 'F2');
 MI(Menus[MENU_MIX], 1, 'Proportions du beton...', '');
 MI(Menus[MENU_MIX], 2, 'Calcul des proportions', '');
 MS(Menus[MENU_MIX], 3);
 MI(Menus[MENU_MIX], 4, 'Generer la microstructure', 'Ctrl+G');
 MI(Menus[MENU_MIX], 5, 'Empilement de granulats...', '');
 MS(Menus[MENU_MIX], 6);
 MI(Menus[MENU_MIX], 7, 'Fractions volumiques...', '');
 { Hydratation }
 ZeroMenu(Menus[MENU_HYD]);
 Menus[MENU_HYD].Title := 'Hydratation';
 Menus[MENU_HYD].ItemCount := 9;
 MI(Menus[MENU_HYD], 0, 'Parametres...', 'F3');
 MI(Menus[MENU_HYD], 1, 'Cinetique des phases...', '');
 MS(Menus[MENU_HYD], 2);
 MI(Menus[MENU_HYD], 3, 'Demarrer / reprendre', 'F9');
 MI(Menus[MENU_HYD], 4, 'Pause', 'Espace');
 MI(Menus[MENU_HYD], 5, 'Arreter', '');
 MI(Menus[MENU_HYD], 6, 'Avancer jusqu''a...', '');
 MS(Menus[MENU_HYD], 7);
 MI(Menus[MENU_HYD], 8, 'Calibrer la cinetique...', '');
 { Proprietes }
 ZeroMenu(Menus[MENU_PROP]);
 Menus[MENU_PROP].Title := 'Proprietes';
 Menus[MENU_PROP].ItemCount := 8;
 MI(Menus[MENU_PROP], 0, 'Resume des proprietes', 'F4');
 MI(Menus[MENU_PROP], 1, 'Courbes en fonction du temps', 'F5');
 MI(Menus[MENU_PROP], 2, 'Percolation', '');
 MI(Menus[MENU_PROP], 3, 'Distribution des pores', '');
 MI(Menus[MENU_PROP], 4, 'Modules elastiques', '');
 MI(Menus[MENU_PROP], 5, 'Resistance a la compression', '');
 MI(Menus[MENU_PROP], 6, 'Conductivite / diffusivite', '');
 MI(Menus[MENU_PROP], 7, 'Zone de transition (ZTI)', '');
 { Vue }
 ZeroMenu(Menus[MENU_VIEW]);
 Menus[MENU_VIEW].Title := 'Vue';
 Menus[MENU_VIEW].ItemCount := 11;
 MI(Menus[MENU_VIEW], 0, 'Coupe XY', '');
 MI(Menus[MENU_VIEW], 1, 'Coupe XZ', '');
 MI(Menus[MENU_VIEW], 2, 'Coupe YZ', '');
 MI(Menus[MENU_VIEW], 3, 'Vue 3D', '');
 MS(Menus[MENU_VIEW], 4);
 MI(Menus[MENU_VIEW], 5, 'Coupe suivante', 'PgDn');
 MI(Menus[MENU_VIEW], 6, 'Coupe precedente', 'PgUp');
 MS(Menus[MENU_VIEW], 7);
 MI(Menus[MENU_VIEW], 8, 'Phases visibles...', '');
 MI(Menus[MENU_VIEW], 9, 'Afficher toutes les phases', '');
 MI(Menus[MENU_VIEW], 10, 'Masquer l''eau (solides)', '');
 { Outils }
 ZeroMenu(Menus[MENU_TOOLS]);
 Menus[MENU_TOOLS].Title := 'Outils';
 Menus[MENU_TOOLS].ItemCount := 3;
 MI(Menus[MENU_TOOLS], 0, 'Calculateur de Bogue...', '');
 MI(Menus[MENU_TOOLS], 1, 'Valeurs par defaut', '');
 MI(Menus[MENU_TOOLS], 2, 'Phases et couleurs', '');
 { Aide }
 ZeroMenu(Menus[MENU_HELP]);
 Menus[MENU_HELP].Title := 'Aide';
 Menus[MENU_HELP].ItemCount := 3;
 MI(Menus[MENU_HELP], 0, 'Aide', 'F1');
 MI(Menus[MENU_HELP], 1, 'Mode lot (ligne de commande)', '');
 MI(Menus[MENU_HELP], 2, 'A propos de VCCTL', '');
End;

Procedure UpdateMenuChecks;
Var I: Integer;
Begin
 For I := 0 To 3 Do Menus[MENU_VIEW].Items[I].Checked := (ViewPlane = I);
 Menus[MENU_HYD].Items[3].Enabled := HaveGrid And Not SimActive;
 Menus[MENU_HYD].Items[4].Enabled := SimActive;
 Menus[MENU_HYD].Items[5].Enabled := SimActive;
 For I := 0 To 7 Do Menus[MENU_PROP].Items[I].Enabled := HaveGrid;
 Menus[MENU_PROP].Items[7].Enabled := HaveGrid And WallOn;
 Menus[MENU_FILE].Items[6].Enabled := HaveGrid;
 Menus[MENU_FILE].Items[7].Enabled := HaveGrid;
 Menus[MENU_FILE].Items[8].Enabled := HaveGrid;
 Menus[MENU_FILE].Items[9].Enabled := Sim.NS > 0;
End;

{ --------------------------------------------------------------------- }
{ Commandes                                                             }
{ --------------------------------------------------------------------- }

Procedure InitAllDefaults;
Begin
 DefaultAll;
 InitAggParams;
 FreeMicro;
 ViewPlane := 0; SliceIdx := 25;
 StatusMsg := '';
End;

Function ConfirmDiscard: Boolean;
Begin
 ConfirmDiscard := True;
 If Modified Then
  ConfirmDiscard := ConfirmBox('Projet modifie', 'Les modifications non enregistrees seront|perdues. Continuer ?');
End;

Procedure DoNewProject;
Begin
 If SimActive Then Begin ShowMessage('Simulation en cours', 'Arretez d''abord l''hydratation.'); Exit; End;
 If Not ConfirmDiscard Then Exit;
 InitAllDefaults;
 NeedRedraw := True;
End;

Procedure DoOpenProject;
Var FN: String;
Begin
 If SimActive Then Begin ShowMessage('Simulation en cours', 'Arretez d''abord l''hydratation.'); Exit; End;
 If Not ConfirmDiscard Then Exit;
 FN := '';
 If FileOpenDialog('Ouvrir un projet', FN, '*.VCT') Then
  If LoadProject(FN) Then ;
 NeedRedraw := True;
End;

Procedure DoSaveAs;
Var FN: String;
Begin
 FN := ProjName;
 If FN = '' Then FN := 'PROJET.VCT';
 If FileSaveDialog('Enregistrer le projet', FN) Then Begin
  If ExtractExt(FN) = '' Then FN := FN + '.VCT';
  If SaveProject(FN) Then ;
 End;
End;

Procedure DoSave;
Begin
 If ProjName = '' Then DoSaveAs
 Else If SaveProject(ProjName) Then ;
End;

Procedure DoOpenMicro;
Var FN: String;
Begin
 If SimActive Then Begin ShowMessage('Simulation en cours', 'Arretez d''abord l''hydratation.'); Exit; End;
 FN := '';
 If FileOpenDialog('Ouvrir une microstructure', FN, '*.MIC') Then
  If LoadMicro(FN) Then ;
 NeedRedraw := True;
End;

Procedure DoSaveMicro;
Var FN: String;
Begin
 If Not NeedGrid Then Exit;
 FN := 'MICRO.MIC';
 If FileSaveDialog('Enregistrer la microstructure', FN) Then Begin
  If ExtractExt(FN) = '' Then FN := FN + '.MIC';
  If SaveMicro(FN) Then ;
 End;
End;

Procedure DoExportText;
Var FN: String;
Begin
 If Not NeedGrid Then Exit;
 FN := 'MICRO.IMG';
 If FileSaveDialog('Exporter la microstructure (texte)', FN) Then ExportMicroText(FN);
End;

Procedure DoExportBMP;
Var FN: String;
Begin
 If Not NeedGrid Then Exit;
 FN := 'COUPE.BMP';
 If FileSaveDialog('Exporter la coupe (BMP)', FN) Then Begin
  If ExtractExt(FN) = '' Then FN := FN + '.BMP';
  ExportSliceBMP(FN);
 End;
End;

Procedure DoExportCSV;
Var FN: String;
Begin
 FN := 'RESULTAT.CSV';
 If FileSaveDialog('Exporter les resultats (CSV)', FN) Then Begin
  If ExtractExt(FN) = '' Then FN := FN + '.CSV';
  ExportResultsCSV(FN);
 End;
End;

Procedure StepSlice(D: Integer);
Begin
 If Not HaveGrid Then Exit;
 SliceIdx := SliceIdx + D;
 If SliceIdx < 0 Then SliceIdx := 0;
 If SliceIdx > N - 1 Then SliceIdx := N - 1;
 NeedRedraw := True;
End;

Procedure SetPlane(P: Integer);
Begin
 If (P = 3) And Not HaveGrid Then Exit;
 ViewPlane := P;
 NeedRedraw := True;
End;

Procedure DoRun;
Begin
 StartSim(EndDays * 24);
End;

Procedure DoPauseResume;
Begin
 If SimActive Then StopSim(True)
 Else If HaveGrid And SimPaused Then StartSim(EndDays * 24);
End;

Procedure DispatchMenu(M, I: Integer);
Var K: Integer;
Begin
 UpdateMenuChecks;
 RedrawAll;
 NeedRedraw := True;
 Case M Of
  MENU_FILE: Case I Of
   0: DoNewProject;
   1: DoOpenProject;
   2: DoSave;
   3: DoSaveAs;
   5: DoOpenMicro;
   6: DoSaveMicro;
   7: DoExportText;
   8: DoExportBMP;
   9: DoExportCSV;
   11: Begin
    If SimActive Then StopSim(False);
    If ConfirmDiscard Then Running := False;
   End;
  End;
  MENU_MAT: Case I Of
   0: CementDialog;
   1: SCMDialog(1, FAm);
   2: SCMDialog(2, SlgM);
   3: SCMDialog(3, FilM);
   4: SCMDialog(4, SFM);
   6: AggDialog('granulats fins (sable)', Sand);
   7: AggDialog('granulats grossiers', Coarse);
   9: ShowPSDTable;
   10: ShowPSDChart;
  End;
  MENU_MIX: Case I Of
   0: MixDialog;
   1: ConcreteDialog;
   2: ShowProportions;
   4: DoGenerate;
   5: AggDialog3D;
   7: ShowPhaseTable;
  End;
  MENU_HYD: Case I Of
   0: HydParamsDialog;
   1: KineticsDialog;
   3: DoRun;
   4: StopSim(True);
   5: StopSim(False);
   6: AdvanceDialog;
   8: CalibrateDialog;
  End;
  MENU_PROP: Case I Of
   0: ShowSummary;
   1: CurvesMenu;
   2: ShowPercolation;
   3: ShowPoreSizes;
   4: ShowElastic;
   5: StrengthDialog;
   6: ConductDialog;
   7: ShowITZ;
  End;
  MENU_VIEW: Case I Of
   0, 1, 2, 3: SetPlane(I);
   5: StepSlice(1);
   6: StepSlice(-1);
   8: PhaseVisDialog;
   9: For K := 0 To NPH - 1 Do PhShow[K] := True;
   10: Begin PhShow[P_POR] := False; PhShow[P_VOID] := False; End;
  End;
  MENU_TOOLS: Case I Of
   0: BogueDialog;
   1: If ConfirmBox('Valeurs par defaut', 'Restaurer tous les parametres par defaut ?|(la microstructure est conservee)') Then Begin
     DefaultAll; InitAggParams;
    End;
   2: ShowPhaseLegend;
  End;
  MENU_HELP: Case I Of
   0: ShowHelp;
   1: Begin
    TClear;
    TAdd('Mode lot (sans interface) :');
    TAdd('');
    TAdd('   VCCTL /B [projet.VCT] [sortie.CSV]');
    TAdd('');
    TAdd('Charge le projet (ou les valeurs par defaut), genere la');
    TAdd('microstructure, simule l''hydratation jusqu''a la duree du projet');
    TAdd('et ecrit les resultats au format CSV.');
    TView('Mode lot');
   End;
   2: ShowAbout;
  End;
 End;
End;

Procedure DoToolbar(Id: Integer);
Begin
 Case Id Of
  TB_GEN: DoGenerate;
  TB_RUN: DoRun;
  TB_PAUSE: StopSim(True);
  TB_STOP: StopSim(False);
  TB_XY, TB_XZ, TB_YZ, TB_3D: SetPlane(Id - TB_XY);
  TB_PREV: StepSlice(-1);
  TB_NEXT: StepSlice(1);
  TB_CURV: CurvesMenu;
  TB_PROP: ShowSummary;
  TB_FORM: MixDialog;
  TB_HELP: ShowHelp;
 End;
End;

{ --------------------------------------------------------------------- }
{ Clavier et souris                                                     }
{ --------------------------------------------------------------------- }

Procedure HandleKey;
Var K: Char; MZ, MN: Integer;
Begin
 K := ReadKey;
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
   #60: MixDialog;
   #61: HydParamsDialog;
   #62: ShowSummary;
   #63: CurvesMenu;
   #67: DoRun;
   #68: Begin MenuActive := True; ActiveMenu := 0; ActiveMenuItem := 0; MenuNextItem(0); End;
   #45: Begin If SimActive Then StopSim(False); If ConfirmDiscard Then Running := False; End;
   #73: StepSlice(-1);
   #81: StepSlice(1);
   #72: StepSlice(1);
   #80: StepSlice(-1);
  End;
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
  #14: DoNewProject;
  #15: DoOpenProject;
  #19: DoSave;
  #7: DoGenerate;
  #17: Begin If SimActive Then StopSim(False); If ConfirmDiscard Then Running := False; End;
  ' ': DoPauseResume;
  #27: If SimActive Then StopSim(True);
  '1': SetPlane(0);
  '2': SetPlane(1);
  '3': SetPlane(2);
  '4': SetPlane(3);
 End;
 NeedRedraw := True;
End;

Procedure HoverInfo;
Var S, OX, I, J: Integer; Ph: Byte;
Begin
 If Not HaveGrid Or (ViewPlane = 3) Or MenuActive Then Exit;
 S := VSZ Div N;
 If S < 1 Then S := 1;
 OX := VX0 + (VSZ - S * N) Div 2;
 I := (MouseX - OX) Div S; J := (MouseY - (VY0 + (VSZ - S * N) Div 2)) Div S;
 If (MouseX >= OX) And (MouseY >= VY0) And (I >= 0) And (I < N) And (J >= 0) And (J < N) Then Begin
  Ph := VoxAt(ViewPlane, I, J);
  StatusMsg := 'Voxel (' + IStr(I + 1) + ', ' + IStr(J + 1) + ') : ' + PhName[Ph];
  HideMouseCursor;
  DrawStatusBar;
  ShowMouseCursor;
 End;
End;

Procedure HandleMouse;
Var
 Hit, ItemHit, T, MZ, LH : Integer;
 LeftDown, OldLeft, Moved : Boolean;
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
 If Moved And Not MenuActive And Not SimActive Then Begin
  T := ToolbarHit(MouseX, MouseY);
  If T >= 0 Then Begin
   If StatusMsg <> TBHint[T] Then Begin StatusMsg := TBHint[T]; HideMouseCursor; DrawStatusBar; ShowMouseCursor; End;
  End
  Else HoverInfo;
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
  Else If HaveGrid And (MouseX >= PNX) Then Begin
   LH := LegendHit(MouseX, MouseY);
   If LH >= 0 Then Begin PhShow[LH] := Not PhShow[LH]; NeedRedraw := True; End;
  End;
 End;
 If Not DlgHappened Then OldMouseButtons := MouseButtons;
 OldMouseX := MouseX; OldMouseY := MouseY;
End;

{ --------------------------------------------------------------------- }
{ Mode lot                                                              }
{ --------------------------------------------------------------------- }

Procedure BatchRun;
Var Prj, Outf, Hdr: String; Fh: Text;
Begin
 BatchMode := True;
 Prj := ''; Outf := 'VCCTL.CSV';
 If ParamCount >= 2 Then Prj := ParamStr(2);
 If ParamCount >= 3 Then Outf := ParamStr(3);
 If (Prj <> '') And (UpStr(ExtractExt(Prj)) = '.CSV') Then Begin Outf := Prj; Prj := ''; End;
 If Prj <> '' Then Begin
  If Not FileOK(Prj) Then Begin WriteLn('Projet introuvable : ', Prj); Halt(1); End;
  Assign(Fh, Prj); Reset(Fh);
  ReadLn(Fh, Hdr);
  IOErr := False;
  ProjIO(Fh, False);
  Close(Fh);
  If ParamCount >= 3 Then Outf := ParamStr(3) Else Outf := ChangeExt(Prj, '.CSV');
  If NVox > MAXN Then NVox := MAXN;
 End;
 WriteLn('VCCTL - mode lot');
 WriteLn('Ciment : ', Cem.Name, '  E/L = ', WB:0:3, '  N = ', NVox, '  graine = ', Seed);
 If Not GenerateMicro Then Begin WriteLn('Echec de la generation : ', StatusMsg); Halt(1); End;
 WriteLn('Microstructure : ', N, '^3 voxels, ', PhCnt[P_POR], ' voxels d''eau');
 SampleInitial;
 SimGoalT := EndDays * 24;
 If CureMode = 2 Then Sim.Temp := CureTemp;
 While (Sim.T < SimGoalT) And (Sim.Cyc < MaxCycles) Do RunCycle;
 Sample;
 WriteLn('Termine : ', Sim.Cyc, ' cycles, t = ', Sim.T:0:1, ' h, alpha = ', AlphaCement:0:3,
         ', Q = ', Sim.Heat / Sim.M0Bind:0:1, ' J/g, porosite = ', 100 * CapPorosity:0:1, ' %');
 ExportResultsCSV(Outf);
 WriteLn('Resultats : ', Outf);
End;

{ --------------------------------------------------------------------- }
{ Initialisation et boucle principale                                   }
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
 InitAllDefaults;
 InitMenus;
 InitToolbar;
 LayoutMenuBar;
 MenuActive := False; ActiveMenu := -1; ActiveMenuItem := -1;
 OldMouseButtons := 0; OldMouseX := -1; OldMouseY := -1;
 ProgPct := -1;
 InitMouseDriver;
 Running := True;
 NeedRedraw := True;
End;

Procedure MainLoop;
Var LastDraw: LongInt;
Begin
 RedrawAll;
 LastDraw := NowMs;
 While Running Do Begin
  If KeyPressed Then HandleKey;
  HandleMouse;
  If SimActive Then Begin
   StepSim(50);
   If Not SimActive Then NeedRedraw := True
   Else If Not MenuActive And (Abs(NowMs - LastDraw) > 250) Then Begin
    RedrawSim;
    LastDraw := NowMs;
   End;
  End;
  If MenuActive Then UpdateMenuChecks;
  If NeedRedraw Then Begin RedrawAll; LastDraw := NowMs; End;
  If Not SimActive Then Delay(8);
 End;
End;

{ ===================================================================== }
{ Programme principal                                                   }
{ ===================================================================== }

BEGIN
 If (ParamStr(1) = '/?') Or (ParamStr(1) = '--help') Or (ParamStr(1) = '-h') Then Begin
  WriteLn('VCCTL - Cette commande permet de simuler l''hydratation du ciment (laboratoire virtuel).');
  WriteLn;
  WriteLn('Syntaxe : VCCTL [projet.VCT]');
  WriteLn('          VCCTL /B [projet.VCT] [sortie.CSV]');
  WriteLn;
  WriteLn(' projet.VCT   Projet a ouvrir');
  WriteLn(' /B           Mode lot : genere, hydrate et ecrit les resultats sans interface');
  Halt;
 End;
 If (UpStr(ParamStr(1)) = '/B') Or (UpStr(ParamStr(1)) = '-B') Then Begin
  InitAllDefaults;
  BatchRun;
  Halt;
 End;
 InitializeGraphics;
 InitializeAll;
 SetViewPort(0, 0, ScreenWidth - 1, ScreenHeight - 1, True);
 If ParamCount >= 1 Then Begin
  If FileOK(ParamStr(1)) Then Begin
   If LoadProject(ParamStr(1)) Then ;
  End
  Else StatusMsg := 'Fichier introuvable : ' + ParamStr(1);
 End;
 Modified := False;
 MainLoop;
 CloseGraph;
END.
