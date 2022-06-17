{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2022
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal, Free Pascal)
}

Program Formula;

Uses Strings;

Type
 ComposeRec=Record
  Nom,Formula:PChar;
 End;

Const
 Composes:Array[0..145]of ComposeRec=(
  (Nom:'Ac‚tyle';Formula:'CH3CO'),
  (Nom:'Ac‚tylac‚tone';Formula:'C5H8O2'),
  (Nom:'Acide ascorbique';Formula:'C6H8O6'),
  (Nom:'Acide folique';Formula:'C19H19N7O6'),
  (Nom:'Acide pantoth‚nique';Formula:'C9H17NO5'),
  (Nom:'Actinium';Formula:'Ac'),
  (Nom:'Ad‚nine';Formula:'C5H5N5'),
  (Nom:'Aluminium';Formula:'Al'),
  (Nom:'Am‚ricium';Formula:'Am'),
  (Nom:'Antimoine';Formula:'Sb'),
  (Nom:'Argent';Formula:'Ag'),
  (Nom:'Argon';Formula:'Ar'),
  (Nom:'Arsenic';Formula:'As'),
  (Nom:'Astate';Formula:'At'),
  (Nom:'Azote';Formula:'N'),
  (Nom:'B1';Formula:'C12H17N4OS'),
  (Nom:'B2';Formula:'C17H20N4O6'),
  (Nom:'B3';Formula:'C6H5NO2'),
  (Nom:'B4';Formula:'C5H5N5'),
  (Nom:'B5';Formula:'C9H17NO5'),
  (Nom:'B8';Formula:'C10H16N2O3S'),
  (Nom:'B9';Formula:'C19H19N7O6'),
  (Nom:'B12';Formula:'C72H100CoN18O17P'),
  (Nom:'Baryum';Formula:'Ba'),
  (Nom:'Berk‚lium';Formula:'Bk'),
  (Nom:'B‚ryllium';Formula:'Be'),
  (Nom:'Biotine';Formula:'C10H16N2O3S'),
  (Nom:'Bismuth';Formula:'Bi'),
  (Nom:'Bohrium';Formula:'Bh'),
  (Nom:'Bore';Formula:'B'),
  (Nom:'Brome';Formula:'Br'),
  (Nom:'Bromure d''actinium';Formula:'AcBr3'),
  (Nom:'C';Formula:'C6H8O6'),
  (Nom:'Cadmium';Formula:'Cd'),
  (Nom:'Calcium';Formula:'Ca'),
  (Nom:'Californium';Formula:'Cf'),
  (Nom:'Carbone';Formula:'C'),
  (Nom:'Carium';Formula:'Cm'),
  (Nom:'C‚rium';Formula:'Ce'),
  (Nom:'C‚sium';Formula:'Cs'),
  (Nom:'Chlore';Formula:'Cl'),
  (Nom:'Chlorure d''actinium';Formula:'AcCl3'),
  (Nom:'Chrome';Formula:'Cr'),
  (Nom:'Cobalamine';Formula:'C72H100CoN18O17P'),
  (Nom:'Cobalt';Formula:'Co'),
  (Nom:'Copernicium';Formula:'Cn'),
  (Nom:'Cuivre';Formula:'Cu'),
  (Nom:'Darmstaditium';Formula:'Ds'),
  (Nom:'Dihydrure d''actinium';Formula:'AcH2'),
  (Nom:'Dubnium';Formula:'Db'),
  (Nom:'Dysprosium';Formula:'Dy'),
  (Nom:'Eau';Formula:'H2O'),
  (Nom:'Einsteinium';Formula:'Es'),
  (Nom:'Erbium';Formula:'Er'),
  (Nom:'Etain';Formula:'Sn'),
  (Nom:'Europium';Formula:'Eu'),
  (Nom:'Fer';Formula:'Fe'),
  (Nom:'Fermium';Formula:'Fm'),
  (Nom:'Fl‚rovium';Formula:'Fl'),
  (Nom:'Fluor';Formula:'F'),
  (Nom:'Fluorure d''actinium';Formula:'AcF3'),
  (Nom:'Francium';Formula:'Fr'),
  (Nom:'Gadolinium';Formula:'Gd'),
  (Nom:'Gallium';Formula:'Ga'),
  (Nom:'Germanium';Formula:'Ge'),
  (Nom:'Hafnium';Formula:'Hf'),
  (Nom:'Hassium';Formula:'Hs'),
  (Nom:'H‚lium';Formula:'He'),
  (Nom:'Holmium';Formula:'Ho'),
  (Nom:'HydrogŠne';Formula:'H'),
  (Nom:'Indium';Formula:'In'),
  (Nom:'Iode';Formula:'I'),
  (Nom:'Iodure d''actinium';Formula:'AcI3'),
  (Nom:'Iridium';Formula:'Ir'),
  (Nom:'Krypton';Formula:'Kr'),
  (Nom:'Lanthane';Formula:'La'),
  (Nom:'Lawrencium';Formula:'Lr'),
  (Nom:'Lithium';Formula:'Li'),
  (Nom:'Livermorium';Formula:'Lv'),
  (Nom:'Lut‚cium';Formula:'Lu'),
  (Nom:'Magn‚sium';Formula:'Mg'),
  (Nom:'ManganŠse';Formula:'Mn'),
  (Nom:'Meitn‚rium';Formula:'Mt'),
  (Nom:'Mend‚l‚vium';Formula:'Md'),
  (Nom:'Mercure';Formula:'Hg'),
  (Nom:'MolybdŠne';Formula:'Mo'),
  (Nom:'Moscovium';Formula:'Mc'),
  (Nom:'N‚odyme';Formula:'Nd'),
  (Nom:'N‚on';Formula:'Ne'),
  (Nom:'Neptunium';Formula:'Np'),
  (Nom:'Niacine';Formula:'C6H5NO2'),
  (Nom:'Nickel';Formula:'Ni'),
  (Nom:'Nihonium';Formula:'Nh'),
  (Nom:'Niobium';Formula:'Nb'),
  (Nom:'Nob‚lium';Formula:'No'),
  (Nom:'Oganesson';Formula:'Og'),
  (Nom:'Or';Formula:'Au'),
  (Nom:'Osmium';Formula:'Os'),
  (Nom:'OxygŠne';Formula:'O'),
  (Nom:'Palladium';Formula:'Pd'),
  (Nom:'Phosphore';Formula:'P'),
  (Nom:'Platine';Formula:'Pt'),
  (Nom:'Plomb';Formula:'Pb'),
  (Nom:'Plutonium';Formula:'Pu'),
  (Nom:'Polonium';Formula:'Po'),
  (Nom:'Potassium';Formula:'K'),
  (Nom:'Pras‚odyme';Formula:'Pr'),
  (Nom:'Prom‚thium';Formula:'Pm'),
  (Nom:'Protactinium';Formula:'Pa'),
  (Nom:'Radon';Formula:'Rn'),
  (Nom:'Radium';Formula:'Ra'),
  (Nom:'Rh‚nium';Formula:'Re'),
  (Nom:'Rhodium';Formula:'Rh'),
  (Nom:'Riboflavine';Formula:'C17H20N4O6'),
  (Nom:'Rubidium';Formula:'Rb'),
  (Nom:'Roentgenium';Formula:'Rg'),
  (Nom:'Ruth‚nium';Formula:'Ru'),
  (Nom:'Rutherfordium';Formula:'Rf'),
  (Nom:'Samarium';Formula:'Sm'),
  (Nom:'Scandium';Formula:'Sc'),
  (Nom:'Seaborgium';Formula:'Sg'),
  (Nom:'S‚l‚mium';Formula:'Se'),
  (Nom:'Sesquioxyde d''actinium';Formula:'Ac2O3'),
  (Nom:'Silicium';Formula:'Si'),
  (Nom:'Sodium';Formula:'Na'),
  (Nom:'Soufre';Formula:'S'),
  (Nom:'Strontium';Formula:'Sr'),
  (Nom:'Tantale';Formula:'Ta'),
  (Nom:'Techn‚tium';Formula:'Tc'),
  (Nom:'Tellure';Formula:'Te'),
  (Nom:'Tennesse';Formula:'Ts'),
  (Nom:'Terbium';Formula:'Tb'),
  (Nom:'Thallium';Formula:'Tl'),
  (Nom:'Thiamine';Formula:'C12H17N4OS'),
  (Nom:'Thorium';Formula:'Th'),
  (Nom:'Thulium';Formula:'Tm'),
  (Nom:'Titane';Formula:'Ti'),
  (Nom:'Trisulfure de diactinium';Formula:'Ac2S3'),
  (Nom:'TungstŠne';Formula:'W'),
  (Nom:'Uranium';Formula:'U'),
  (Nom:'Vanadium';Formula:'V'),
  (Nom:'X‚non';Formula:'Xe'),
  (Nom:'Ytterbium';Formula:'Yb'),
  (Nom:'Yttrium';Formula:'Y'),
  (Nom:'Zinc';Formula:'Zn'),
  (Nom:'Zirconium';Formula:'Zr')


 );

Function StrToUpper(S:String):String;
Var
 I:Byte;
Begin
 For I:=1 to Length(S)do Begin
  If S[I]in['‚','Š','ˆ',#232,#233,#234]Then S[I]:='E'Else
  If S[I]in['Œ',#238]Then S[I]:='I'Else
  If S[I] in['a'..'z']Then S[I]:=Chr(Ord(S[I])-32);
 End;
 StrToUpper:=S;
End;

Var
 I,J:Integer;

BEGIN
 If(ParamStr(1)='/?')or(ParamStr(1)='--help')or(ParamStr(1)='-h')or
   (ParamStr(1)='/h')or(ParamStr(1)='/H')Then Begin
  WriteLn('FORMULA: Cette commande permet d''afficher ',
                   'la formule d''un compos‚ chimique.');
  WriteLn;
  WriteLn('Syntaxe : FORMULA compose');
  WriteLn;
  WriteLn(' compose  Ce parametre permet d''indiquer la compos‚ a analyser');
 End
  Else
 If ParamCount>0Then Begin
  For I:=1 to ParamCount do Begin
   For J:=Low(Composes) to High(Composes)do Begin
    If StrToUpper(ParamStr(I))=StrToUpper(Composes[J].Nom)Then Begin
     WriteLn('Nom :     ',Composes[J].Nom);
     WriteLn('Formule : ',Composes[J].Formula);
    End;
   End;
  End;
 End;
END.
