% BibTeX bibliography style `astats'
	% version 0.99a for BibTeX versions 0.99a or later, LaTeX version 2.09.
	% This style is designed to produce references according to the
	% style used in Anals of Stats., and in many other statistics journals
	% It should be used with LaTeX optional style ASTATS.STY
	% i.e. : \documentstyle[astats]{article}
	% 
	% Two entry types were added : thesetrois and theseetat which are 
	% to phdthesis, but correspond to french type of thesis
	%
% Olivier NICOLE (nicole@frinra72.bitnet) Jan. 89

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
  }
  {}
  { label extra.label sort.label }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
	{ add.period$ write$
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
	    { add.period$ " " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {format.date.label}
{ year empty$
  { extra.label empty$
     { "" }
     { " (" extra.label * ")" * }
    if$
  }
  { extra.label empty$
     { " (" year * ")" * }
     { " ("  year * extra.label * ")" * }
    if$
  }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem[" write$
  label write$
  format.date.label write$
  "]{" write$
  cite$ write$
  "}%" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

FUNCTION {boldface}
{ duplicate$ empty$
    { pop$ "" }
    { "{\bf " swap$ * "}" * }
  if$
}

FUNCTION {smallcaps}
{ duplicate$ empty$
    { pop$ "" }
    { "{\sc " swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}{, f.}{, jj}" format.name$ smallcaps 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
            {
              s nameptr "{ff }{vv }{ll}{ jj}" format.name$ "others" =
		{ " " "et~al." emphasize * * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

% added
FUNCTION {format.names.editors}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{f.~}{vv~}{ll}{, jj}" format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { numnames #2 >
		{ "," * }
		'skip$
	      if$
	      t "others" =
		{ " " "et~al." emphasize * * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
	{ ", editors" * }
	{ ", editor" * }
      if$
    }
  if$
}

%added
FUNCTION {format.ed}
{ author empty$
    { "" }
    { editor empty$
        { "" }
        { "(" editor format.names.editors * ",~ed.)" * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { title "t" change.case$ }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

% Modified
FUNCTION {format.date}
{ year empty$
    { month empty$
	{ "" }
	{ "there's a month but no year in " cite$ * warning$
	  "(" month *
	}
      if$
    }
    { month empty$
	{ "(" year * }
	{ "(" month * ", " * year * }
      if$
    }
  if$
 extra.label empty$
  { year empty$
    month empty$
      and
        'skip$
        { ")" * }
    if$
  }
  { year empty$
    month empty$
      and
        { "(" }
        'skip$
    if$
    extra.label ")" * * }
 if$
}

FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "volume" volume tie.or.space.connect
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

% Note : Changer l'odre en francais : Serie, puis No.
FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
	{ output.state mid.sentence =
	    { "number" }
	    { "Number" }
	  if$
	  number tie.or.space.connect
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
	    { " in " * series * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " edition" * }
	{ edition "t" change.case$ " edition" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

% Added
FUNCTION {format.pages.nc}
{ pages empty$
    { "" }
    { pages n.dashify }
  if$
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages n.dashify tie.or.space.connect}
%    { pages multi.page.check
%	{ "pages" pages n.dashify tie.or.space.connect }
%	{ "page" pages tie.or.space.connect }
%      if$
%    }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume empty$
  { number empty$
     { "" }
     { "there's a number but no volume in " cite$ * warning$ 
       "(" number * ")" * boldface
     }
    if$
  }
  { volume boldface
    number empty$
      'skip$
      { "(" number * ")" * boldface tie.or.space.connect }
    if$
  }
  if$
  format.pages
%{ volume field.or.null
%  number empty$
%    'skip$
%    { "(" number * ")" * *
%      volume empty$
%	{ "there's a number but no volume in " cite$ * warning$ }
%	'skip$
%      if$
%    }
%  if$
%  pages empty$
%    'skip$
%    { duplicate$ empty$
%	{ pop$ format.pages }
%	{ ":" * pages n.dashify * }
%      if$
%    }
%  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages.nc
    { type empty$
	{ "chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", " * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
	{ "In " booktitle emphasize * }
        { "In " booktitle emphasize ", " * * format.ed *}
%	{ "In " format.editors * ", " * booktitle emphasize * }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.french.thesis.type}
{ type empty$
    'skip$
    { pop$
      " de " type "t" change.case$ *
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "In {\em " journal * "\/}" * }
      if$
    }
    { "In " key * }
  if$
  ", \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{f.~}{vv~}{ll}{, j}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * editor #2 "{f.~}{vv~}{ll}{, j}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Volume" volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "{\em " * series * "\/}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  ", \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "In {\em " booktitle * "\/}" * }
	  if$
	}
	{ "In " key * }
      if$
    }
    { "In " format.crossref.editor * }
  if$
  ", \cite{" * crossref * "}" *
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  crossref missing$
    { journal emphasize "journal" output.check
      new.block
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      new.block
      format.pages.nc output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
%      crossref missing$
%	{ "author and editor" editor either.or.check }
%	'skip$
%      if$
    }
  if$
  new.block
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  new.block
  crossref missing$
% Added
    { format.ed output
      format.bvolume output
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  note output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  new.block
  format.date output
  new.block
  format.title "title" output.check
  howpublished address new.block.checkb
  howpublished output
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
%      crossref missing$
%	{ "author and editor" editor either.or.check }
%	'skip$
%      if$
    }
  if$
  new.block
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  new.block
  crossref missing$
    { format.ed output
      format.bvolume output
      format.chapter.pages "chapter and pages" output.check
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { format.chapter.pages "chapter and pages" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.chapter.pages output
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
      format.edition output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.pages.nc output
      format.number.series output
      address empty$
	{ organization publisher new.sentence.checkb
	  organization output
	  publisher output
	}
	{ address output.nonnull
	  new.sentence
	  organization output
	  publisher output
	}
      if$
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages.nc output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization empty$
	'skip$
	{ organization smallcaps output.nonnull
	  address output
	  new.block
          format.date output

	}
      if$
    }
    { format.authors output.nonnull 
      new.block
      format.date }
  if$
  new.block
  format.btitle "title" output.check
  author empty$
    { organization empty$
	{ format.date output
          new.block
          address new.block.checka
	  address output
	}
	'skip$
      if$
    }
    { organization address new.block.checkb
      organization output
      address output
    }
  if$
  format.edition output
  new.block
  note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.title "title" output.check
  new.block
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  format.date output
  title howpublished new.block.checkb
  format.title output
  howpublished new.block.checka
  howpublished output
  new.block
  note output
  fin.entry
  empty.misc.check
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  new.block
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {thesetrois}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  new.block
  "Th\`ese de troisi\`eme cycle" format.french.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {theseetat}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  new.block
  "Th\`ese d'\'etat" format.french.thesis.type output.nonnull
  school "school" output.check
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization smallcaps output }
    { format.editors output.nonnull }
  if$
  new.block
  format.date "year" output.check
  new.block
  format.btitle "title" output.check
  format.bvolume output
  format.number.series output
  address empty$
    { editor empty$
	{ publisher new.sentence.checka }
	{ organization publisher new.sentence.checkb
	  organization output
	}
      if$
      publisher output
    }
    { address output.nonnull
      new.sentence
      editor empty$
	'skip$
	{ organization output }
      if$
      publisher output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date "year" output.check
  format.title "title" output.check
  new.block
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  new.block
  note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  new.block
  format.date output
  new.block
  format.title "title" output.check
  new.block
  note "note" output.check
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sep."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

INTEGERS { et.al.char.used }

FUNCTION {initialize.et.al.char.used}
{ #0 'et.al.char.used :=
}

EXECUTE {initialize.et.al.char.used}

FUNCTION {format.lab.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  numnames #2 >
    { s nameptr "{vv~}{ll}" format.name$ " " "et~al." emphasize * * }
    { { namesleft #0 > }
    { s nameptr "{vv~}{ll}" format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
            {
              s nameptr "{ff }{vv }{ll}{ jj}" format.name$ "others" =
		{ " et~al." * }
		{ " and " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$ }
  if$
}
%FUNCTION {format.lab.names}
%{ 's :=
%  s num.names$ 'numnames :=
%  numnames #1 >
%    {
%    { numnames #4 >
%	{ #3 'namesleft := }
%	{ numnames 'namesleft := }
%      if$
%      #1 'nameptr :=
%      ""
%	{ namesleft #0 > }
%	{ nameptr numnames =
%	    { s nameptr "{ff }{vv }{ll}{ jj}" format.name$ "others" =
%		{ "{\etalchar{+}}" *
%		  #1 'et.al.char.used :=
%		}
%		{ s nameptr "{v{}}{l{}}" format.name$ * }
%	      if$
%	    }
%	    { s nameptr "{v{}}{l{}}" format.name$ * }
%	  if$
%	  nameptr #1 + 'nameptr :=
%	  namesleft #1 - 'namesleft :=
%	}
%      while$
%      numnames #4 >
%	{ "{\etalchar{+}}" *
%	  #1 'et.al.char.used :=
%	}
%	'skip$
%      if$
%    }
%    { s #1 "{vv~}{ll}" format.name$ 
%    { s #1 "{v{}}{l{}}" format.name$
%      duplicate$ text.length$ #2 <
%	{ pop$ s #1 "{ll}" format.name$ #3 text.prefix$ }
%	'skip$
%      if$
%    }
%  if$
%}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
	    { cite$ }
%	    { cite$ #1 #3 substring$ }
	{ key }
%	{ key #3 text.prefix$ }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
	{ key empty$
	    { cite$ }
%	    { cite$ #1 #3 substring$ }
	    { key }
%	{ key #3 text.prefix$ }
	  if$
	}
	{ editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.key.organization.label}
{ author empty$
    { key empty$
	{ organization empty$
	    { cite$ }
%	    { cite$ #1 #3 substring$ }
	    { organization }
%	    { "The " #4 organization chop.word #3 text.prefix$ }
	  if$
	}
	{ key }
%	{ key #3 text.prefix$ }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.organization.label}
{ editor empty$
    { key empty$
	{ organization empty$
	    { cite$ }
%	    { cite$ #1 #3 substring$ }
	    { organization }
%	    { "The " #4 organization chop.word #3 text.prefix$ }
	  if$
	}
	{ key }
%	{ key #3 text.prefix$ }
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.label}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
	'editor.key.organization.label
	{ type$ "manual" =
	    'author.key.organization.label
	    'author.key.label
	  if$
	}
      if$
    }
  if$
  duplicate$
% Modified
%  format.date.label tie.or.space.connect
  'label :=
  year field.or.null purify$ #-1 #4 substring$
  *
  sortify 'sort.label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  numnames #2 >
    { s nameptr "{ll{ }}{  vv{ }}{  ff{ }}{  jj{ }}" format.name$ 
      "   et al" * sortify * }
    { { namesleft #0 > }
    { nameptr #1 >
	{ "   " * }
	'skip$
      if$
      s nameptr "{ll{ }}{  vv{ }}{  ff{ }}{  jj{ }}" format.name$ 't :=
%      s nameptr "{vv{ } }{ll{ }}{  ff{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
	{ "et al" * }
	{ t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$ }
  if$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
	{ "to sort, need author or key in " cite$ * warning$
	  ""
	}
	{ key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
	{ key empty$
	    { "to sort, need author, editor, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
	{ key empty$
	    { "to sort, need author, organization, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
	{ key empty$
	    { "to sort, need editor, organization, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {presort}
{ calc.label
%  sort.label
%  "    "
%  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
	'editor.organization.sort
	{ type$ "manual" =
	    'author.organization.sort
	    'author.sort
	  if$
	}
      if$
    }
  if$
%  *
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
% Added
  month field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
% duplicate$
  'sort.key$ :=
% To del
% "sort.key   : " swap$ * write$ newline$
% "sort.label : " sort.label  * write$ newline$ newline$

}

ITERATE {presort}

SORT

STRINGS { longest.label last.sort.label next.extra }

INTEGERS { longest.label.width last.extra.num }

FUNCTION {initialize.longest.label}
{ "Tototo" 'longest.label :=
%{ "" 'longest.label :=
  #0 int.to.chr$ 'last.sort.label :=
  "" 'next.extra :=
  #0 'longest.label.width :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.sort.label sort.label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      sort.label 'last.sort.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
%  label extra.label * 'label :=
%  label width$ longest.label.width >
%    { label 'longest.label :=
%      label width$ 'longest.label.width :=
%    }
%    'skip$
%  if$
  extra.label 'next.extra :=
}

% Deleted
%EXECUTE {initialize.longest.label}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {begin.bib}
{ 
%{ et.al.char.used
%    { "\newcommand{\etalchar}[1]{$^{#1}$}" write$ newline$ }
%    'skip$
%  if$
  preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{Tototo}" write$ newline$
%  "\begin{thebibliography}{"  longest.label * "}" * write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
