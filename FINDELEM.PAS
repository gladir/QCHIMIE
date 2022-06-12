{ @author: Sylvain Maltais (support@gladir.com)
  @created: 2022
  @website(https://www.gladir.com/qchimie)
  @abstract(Target: Turbo Pascal, Free Pascal)
}

Program FindElem;

Uses Strings;

Const
 ElementList:Array[0..117]of PChar=(
   'Actinium','Aluminium','AmÇricium','Antimoine','Argent',
   'Argon','Arsenic','Astate','Azote','Baryum','BerkÇlium',
   'Bismuth','Bohrium','Bore','Brome','BÇryllium','Cadmium',
   'Calcium','Californium','Carbone','Carium','Chlore',
   'Chrome','Cobalt','Copernicium','Cuivre','CÇrium',
   'CÇsium','Darmstaditium','Dubnium','Dysprosium',
   'Einsteinium','Erbium','Europium','Fer','Fermium',
   'Fluor','FlÇrovium','Francium','Gadolinium','Gallium',
   'Germanium','Hafnium','Hassium','Holmium','Hydrogäne',
   'HÇlium','Indium','Iode','Iridium','Krypton','Lanthane',
   'Lawrencium','Lithium','Livermorium','LutÇcium','MagnÇsium',
   'Manganäse','MeitnÇrium','MendÇlÇvium','Mercure',
   'Molybdäne','Moscovium','Neptunium','Nickel','Nihonium',
   'Niobium','NobÇlium','NÇodyme','NÇon','Oganesson','Or',
   'Osmium','Oxygäne','Palladium','Phosphore','Platine',
   'Plomb','Plutonium','Polonium','Potassium','PrasÇodyme',
   'PromÇthium','Protactinium','Radium','Radon','Rhodium',
   'RhÇnium','Roentgenium','Rubidium','Rutherfordium',
   'RuthÇnium','Samarium','Scandium','Seaborgium','Silicium',
   'Sodium','Soufre','Strontium','SÇlÇmium','Tantale',
   'TechnÇtium','Tellure','Tennesse','Terbium','Thallium',
   'Thorium','Thulium','Titane','Tungstäne','Uranium',
   'Vanadium','XÇnon','Ytterbium','Yttrium','Zinc','Zirconium',
   'êtain'
 );

Var
 LineNumber:LongInt;
 Option:Set of (_Lines);
 BeginWord,I,J:Integer;
 FoundMultiWord:Boolean;
 SourceFile:Text;
 FileName,CurrLine,CurrWord:String;

Function StrToUpper(S:String):String;
Var
 I:Byte;
Begin
 For I:=1 to Length(S)do Begin
  If S[I] in['a'..'z']Then S[I]:=Chr(Ord(S[I])-32);
 End;
 StrToUpper:=S;
End;

Function IsLetter(Chr:Char):Boolean;Begin
 IsLetter:=Chr in ['A'..'Z','a'..'z','Ç','à','ä','Ö','É','ì','á'];
End;

Procedure CompileWord;
Var
 I:Integer;
Begin
 For I:=Low(ElementList) to High(ElementList) do Begin
  If StrToUpper(StrPas(ElementList[I]))=StrToUpper(CurrWord)Then Begin
   If FileName<>''Then Begin
    If(_Lines in Option)Then Begin
     WriteLn('Ville trouve dans la ligne numero ',LineNumber,' :');
    End;
   End;
   WriteLn(StrPas(ElementList[I]));
   Exit;
  End;
 End;
End;

Procedure ParseLineDetectCity;
Var
 J:Integer;
Begin
 CurrWord:='';
 BeginWord:=1;
 For I:=1 to Length(CurrLine)do Begin
  If(IsLetter(CurrLine[I]))Then CurrWord:=CurrWord+CurrLine[I]
   Else
  Begin
   FoundMultiWord:=False;
   For J:=Low(ElementList) to High(ElementList) do Begin
    If(StrLen(ElementList[J])>Length(CurrWord))and
      (StrToUpper(StrPas(ElementList[J]))=
      StrToUpper(Copy(CurrLine,BeginWord,StrLen(ElementList[J]))))Then Begin
     FoundMultiWord:=True;
     Break;
    End;
   End;
   If(FoundMultiWord)Then Begin
    CurrWord:=CurrWord+CurrLine[I]
   End
    Else
   Begin
    CompileWord;
    CurrWord:='';
    BeginWord:=I+1;
   End;
  End;
 End;
 CompileWord;
End;

BEGIN
 If(ParamStr(1)='/?')or(ParamStr(1)='--help')or(ParamStr(1)='-h')or
   (ParamStr(1)='/h')or(ParamStr(1)='/H')Then Begin
  WriteLn('FINDELEM : Cette commande permet de detecter les ÇlÇments ',
                      'chimiques mentionnÇ dans un texte.');
  WriteLn;
  WriteLn('Syntaxe : FINDELEM "message"');
  WriteLn('          FINDELEM /FILE:fichier [/LINES]');
  WriteLn;
  WriteLn(' /FILE:  Ce parametre permet d''indiquer le fichier a analyser');
  WriteLn(' /LINES  Ce parametre permet d''afficher le numero de ligne');
 End
  Else
 Begin
  Option:=[];
  LineNumber:=0;
  FileName:='';
  CurrLine:='';
  For I:=1 to ParamCount do Begin
   If StrToUpper(ParamStr(I))='/LINES'Then Include(Option,_Lines);
   If StrToUpper(Copy(ParamStr(I),1,6))='/FILE:'Then Begin
    FileName:=Copy(ParamStr(I),7,255);
   End
    Else
   Begin
    If CurrLine=''Then CurrLine:=ParamStr(I)
                  Else CurrLine:=CurrLine+' '+ParamStr(I);
   End;
  End;
  If FileName<>''Then Begin
    Assign(SourceFile,FileName);
    Reset(SourceFile);
    While Not EOF(SourceFile)do Begin
     Inc(LineNumber);
     ReadLn(SourceFile,CurrLine);
     ParseLineDetectCity;
    End;
    Close(SourceFile);
   End
    Else
   ParseLineDetectCity;
 End;
END.
