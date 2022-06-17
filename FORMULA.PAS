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
 Composes:Array[0..208]of ComposeRec=(
  (Nom:'Acanthite';Formula:'Ag2S'),
  (Nom:'Acide ascorbique';Formula:'C6H8O6'),
  (Nom:'Acide folique';Formula:'C19H19N7O6'),
  (Nom:'Acide pantoth‚nique';Formula:'C9H17NO5'),
  (Nom:'Actinium';Formula:'Ac'),
  (Nom:'Ac‚tate d''argent';Formula:'AgC2H3O2'),
  (Nom:'Ac‚tylac‚tone';Formula:'C5H8O2'),
  (Nom:'Ac‚tyle';Formula:'CH3CO'),
  (Nom:'Ac‚tylure d''argent';Formula:'Ag2C2'),
  (Nom:'Ac‚tylure d''argent';Formula:'AgC2H'),
  (Nom:'Ad‚nine';Formula:'C5H5N5'),
  (Nom:'Aluminium';Formula:'Al'),
  (Nom:'Am‚ricium';Formula:'Am'),
  (Nom:'Antimoine';Formula:'Sb'),
  (Nom:'Argent';Formula:'Ag'),
  (Nom:'Argon';Formula:'Ar'),
  (Nom:'Arsenic';Formula:'As'),
  (Nom:'Ars‚niate d''argent';Formula:'Ag3AsO4'),
  (Nom:'Astate';Formula:'At'),
  (Nom:'Azote';Formula:'N'),
  (Nom:'Azoture d''argent';Formula:'AgN3'),
  (Nom:'B1';Formula:'C12H17N4OS'),
  (Nom:'B12';Formula:'C72H100CoN18O17P'),
  (Nom:'B2';Formula:'C17H20N4O6'),
  (Nom:'B3';Formula:'C6H5NO2'),
  (Nom:'B4';Formula:'C5H5N5'),
  (Nom:'B5';Formula:'C9H17NO5'),
  (Nom:'B8';Formula:'C10H16N2O3S'),
  (Nom:'B9';Formula:'C19H19N7O6'),
  (Nom:'Baryum';Formula:'Ba'),
  (Nom:'Berk‚lium';Formula:'Bk'),
  (Nom:'Bifluorure d''argent';Formula:'AgHF2'),
  (Nom:'Biotine';Formula:'C10H16N2O3S'),
  (Nom:'Bismuth';Formula:'Bi'),
  (Nom:'Bohrium';Formula:'Bh'),
  (Nom:'Bore';Formula:'B'),
  (Nom:'Bromate d''argent';Formula:'AgBrO3'),
  (Nom:'Brome';Formula:'Br'),
  (Nom:'Bromure d''actinium';Formula:'AcBr3'),
  (Nom:'Bromure d''argent';Formula:'AgBr'),
  (Nom:'B‚ryllium';Formula:'Be'),
  (Nom:'C';Formula:'C6H8O6'),
  (Nom:'Cadmium';Formula:'Cd'),
  (Nom:'Calcium';Formula:'Ca'),
  (Nom:'Californium';Formula:'Cf'),
  (Nom:'Carbonate d''argent';Formula:'Ag2CO3'),
  (Nom:'Carbone';Formula:'C'),
  (Nom:'Carium';Formula:'Cm'),
  (Nom:'Chlorate d''argent';Formula:'AgClO3'),
  (Nom:'Chlore';Formula:'Cl'),
  (Nom:'Chlorite d''argent';Formula:'AgClO2'),
  (Nom:'Chlorure d''actinium';Formula:'AcCl3'),
  (Nom:'Chlorure d''argent';Formula:'AgCl'),
  (Nom:'Chromate d''argent';Formula:'Ag2CrO4'),
  (Nom:'Chrome';Formula:'Cr'),
  (Nom:'Citrate d''argent';Formula:'Ag3C6H5O7'),
  (Nom:'Cobalamine';Formula:'C72H100CoN18O17P'),
  (Nom:'Cobalt';Formula:'Co'),
  (Nom:'Copernicium';Formula:'Cn'),
  (Nom:'Cuivre';Formula:'Cu'),
  (Nom:'Cyanure d''argent';Formula:'AgCN'),
  (Nom:'C‚rium';Formula:'Ce'),
  (Nom:'C‚sium';Formula:'Cs'),
  (Nom:'Darmstaditium';Formula:'Ds'),
  (Nom:'Dibromure d''argent';Formula:'Ag3Br2'),
  (Nom:'Dichlorure d''argent';Formula:'Ag2Cl2'),
  (Nom:'Dichromate d''argent';Formula:'Ag2Cr2O7'),
  (Nom:'Difluorure d''argent';Formula:'AgF2'),
  (Nom:'Dihydrure d''actinium';Formula:'AcH2'),
  (Nom:'Di‚thyldithiocarbamate d''argent';Formula:'Ag(C2H5)2NCS2'),
  (Nom:'Dubnium';Formula:'Db'),
  (Nom:'Dysprosium';Formula:'Dy'),
  (Nom:'Eau';Formula:'H2O'),
  (Nom:'Einsteinium';Formula:'Es'),
  (Nom:'Erbium';Formula:'Er'),
  (Nom:'Etain';Formula:'Sn'),
  (Nom:'Europium';Formula:'Eu'),
  (Nom:'Fer';Formula:'Fe'),
  (Nom:'Fermium';Formula:'Fm'),
  (Nom:'Fluor';Formula:'F'),
  (Nom:'Fluorure d''actinium';Formula:'AcF3'),
  (Nom:'Fluorure d''argent dihydrat‚';Formula:'AgF2H2O'),
  (Nom:'Fluorure d''argent t‚trahydrat‚';Formula:'AgF4H2O'),
  (Nom:'Fluorure d''argent';Formula:'AgF'),
  (Nom:'Fl‚rovium';Formula:'Fl'),
  (Nom:'Francium';Formula:'Fr'),
  (Nom:'Fulminate d''argent';Formula:'AgCNO'),
  (Nom:'Gadolinium';Formula:'Gd'),
  (Nom:'Gallium';Formula:'Ga'),
  (Nom:'Germanium';Formula:'Ge'),
  (Nom:'Hafnium';Formula:'Hf'),
  (Nom:'Hassium';Formula:'Hs'),
  (Nom:'Hexafluoroantimoniate d''argent';Formula:'AgSbF6'),
  (Nom:'Hexafluoroars‚niate d''argent';Formula:'AgAsF6'),
  (Nom:'Hexafluorophosphate d''argent';Formula:'AgPF6'),
  (Nom:'Holmium';Formula:'Ho'),
  (Nom:'HydrogŠne';Formula:'H'),
  (Nom:'H‚lium';Formula:'He'),
  (Nom:'H‚mifluorure d''argent';Formula:'Ag2F'),
  (Nom:'Indium';Formula:'In'),
  (Nom:'Iodate d''argent';Formula:'AgIO3'),
  (Nom:'Iode';Formula:'I'),
  (Nom:'Iodure d''actinium';Formula:'AcI3'),
  (Nom:'Iodure d''argent';Formula:'AgI'),
  (Nom:'Iridium';Formula:'Ir'),
  (Nom:'Krypton';Formula:'Kr'),
  (Nom:'Lactate d''argent monohydrat‚';Formula:'AgC3H5O3H2O'),
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
  (Nom:'Molybdate d''argent';Formula:'Ag2MoO4'),
  (Nom:'MolybdŠne';Formula:'Mo'),
  (Nom:'Moscovium';Formula:'Mc'),
  (Nom:'Neptunium';Formula:'Np'),
  (Nom:'Niacine';Formula:'C6H5NO2'),
  (Nom:'Nickel';Formula:'Ni'),
  (Nom:'Nihonium';Formula:'Nh'),
  (Nom:'Niobium';Formula:'Nb'),
  (Nom:'Nitrate d''argent';Formula:'AgNO3'),
  (Nom:'Nitrite d''argent';Formula:'AgNO2'),
  (Nom:'Nob‚lium';Formula:'No'),
  (Nom:'N‚odyme';Formula:'Nd'),
  (Nom:'N‚on';Formula:'Ne'),
  (Nom:'Oganesson';Formula:'Og'),
  (Nom:'Or';Formula:'Au'),
  (Nom:'Osmium';Formula:'Os'),
  (Nom:'Oxalate d''argent';Formula:'Ag2C2O4'),
  (Nom:'Oxyde argenteux';Formula:'Ag2O'),
  (Nom:'Oxyde argentique';Formula:'Ag2O2'),
  (Nom:'Oxyde d''argent(I,III)';Formula:'AgO'),
  (Nom:'OxygŠne';Formula:'O'),
  (Nom:'Palladium';Formula:'Pd'),
  (Nom:'Perchlorate d''argent monohydrat‚';Formula:'AgClO4H2O'),
  (Nom:'Perchlorate d''argent';Formula:'AgClO4'),
  (Nom:'Permanganate d''argent';Formula:'AgMnO4'),
  (Nom:'Perrh‚nate d''argent';Formula:'AgReO4'),
  (Nom:'Phosphate d''argent';Formula:'Ag3PO4'),
  (Nom:'Phosphite d''argent';Formula:'AgPO3'),
  (Nom:'Phosphore';Formula:'P'),
  (Nom:'Picrate d''argent monohydrat‚';Formula:'AgC6H2N3O7H2O'),
  (Nom:'Platine';Formula:'Pt'),
  (Nom:'Plomb';Formula:'Pb'),
  (Nom:'Plutonium';Formula:'Pu'),
  (Nom:'Polonium';Formula:'Po'),
  (Nom:'Potassium';Formula:'K'),
  (Nom:'Pras‚odyme';Formula:'Pr'),
  (Nom:'Prom‚thium';Formula:'Pm'),
  (Nom:'Protactinium';Formula:'Pa'),
  (Nom:'Radium';Formula:'Ra'),
  (Nom:'Radon';Formula:'Rn'),
  (Nom:'Rhodium';Formula:'Rh'),
  (Nom:'Rh‚nium';Formula:'Re'),
  (Nom:'Riboflavine';Formula:'C17H20N4O6'),
  (Nom:'Roentgenium';Formula:'Rg'),
  (Nom:'Rubidium';Formula:'Rb'),
  (Nom:'Rutherfordium';Formula:'Rf'),
  (Nom:'Ruth‚nium';Formula:'Ru'),
  (Nom:'Samarium';Formula:'Sm'),
  (Nom:'Scandium';Formula:'Sc'),
  (Nom:'Seaborgium';Formula:'Sg'),
  (Nom:'Sesquioxyde d''actinium';Formula:'Ac2O3'),
  (Nom:'Silicium';Formula:'Si'),
  (Nom:'Sodium';Formula:'Na'),
  (Nom:'Soufre';Formula:'S'),
  (Nom:'Strontium';Formula:'Sr'),
  (Nom:'Sulfate d''argent';Formula:'Ag2SO4'),
  (Nom:'Sulfite d''argent';Formula:'Ag2SO3'),
  (Nom:'Sulfure d''argent';Formula:'Ag2S'),
  (Nom:'S‚l‚mium';Formula:'Se'),
  (Nom:'S‚l‚niate d''argent';Formula:'Ag2SeO4'),
  (Nom:'S‚l‚nite d''argent';Formula:'Ag2SeO3'),
  (Nom:'S‚l‚niure d''argent';Formula:'Ag2Se'),
  (Nom:'Tantale';Formula:'Ta'),
  (Nom:'Techn‚tium';Formula:'Tc'),
  (Nom:'Tellure';Formula:'Te'),
  (Nom:'Tellurure d''argent';Formula:'Ag2Te'),
  (Nom:'Tennesse';Formula:'Ts'),
  (Nom:'Terbium';Formula:'Tb'),
  (Nom:'Thallium';Formula:'Tl'),
  (Nom:'Thiamine';Formula:'C12H17N4OS'),
  (Nom:'Thioantimonate d''argent';Formula:'Ag3SbS3'),
  (Nom:'Thiocyanate d''argent';Formula:'AgSCN'),
  (Nom:'Thiosulfate d''argent';Formula:'Ag2S2O3'),
  (Nom:'Thorium';Formula:'Th'),
  (Nom:'Thulium';Formula:'Tm'),
  (Nom:'Titane';Formula:'Ti'),
  (Nom:'Tribromure d''argent';Formula:'Ag3Br3'),
  (Nom:'Trichlorure d''argent';Formula:'Ag3Cl3'),
  (Nom:'Triiodure d''argent';Formula:'Ag3I3'),
  (Nom:'Trisulfure de diactinium';Formula:'Ac2S3'),
  (Nom:'Tungstate d''argent';Formula:'Ag2WO4'),
  (Nom:'TungstŠne';Formula:'W'),
  (Nom:'T‚trachloroaluminate d''argent';Formula:'AgAlCl4'),
  (Nom:'T‚traiodomercurate(II) d''argent(I)';Formula:'Ag2HgI4'),
  (Nom:'Uranium';Formula:'U'),
  (Nom:'Vanadium';Formula:'V'),
  (Nom:'Xanthoconite';Formula:'Ag3AsS3'),
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
    If StrToUpper(ParamStr(I))=StrToUpper(StrPas(Composes[J].Nom))Then Begin
     WriteLn('Nom :     ',Composes[J].Nom);
     WriteLn('Formule : ',Composes[J].Formula);
    End;
   End;
  End;
 End;
END.
