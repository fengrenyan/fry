
%forget for formulas with degree <=1

%:- op(300,fx,@). /*global*/
%:- op(400,fx,$). /*exist*/
%：- op(900,fx, *). /*next*/

:- op(300,fx,@). /*global*/
:- op(300,fx,*). /*next*/
:- op(300,fx,?). /*future*/
:- op(400,xfy,$). /*until*/
%:- op(400,xfy,#). /*unless*/
:- op(1000,fx, ^). /*exist*/
:- op(1000,fx,~). /*forall*/

:- op(500,fx,-).
%:- op(500,fx,~).
:- op(780,xfy,&).
:- op(790,xfy,\/).%/
:- op(800,xfy,->).
:- op(800,xfy,<->).


:- dynamic(prop/1).
:- dynamic(formula/3).
:- dynamic(pair/2).
:- dynamic(alpha/3).
:- dynamic(beta/3).
:- dynamic(gamma/3).
:- dynamic(snfclause/4).



%computing the degree of CTL formula
degree(- F, D) :- degree(F, D).
degree(F1 \/ F2, D) :- degree(F1, D1), degree(F2, D2), max(D1, D2, D). %/
degree(F1 & F2, D) :- degree(F1, D1), degree(F2, D2), max(D1, D2, D).
degree((^(*F3)), D) :- degree(F3, D1), D is D1+1.
degree((^(@F)), D) :- degree(F, D1), D is D1+1.
degree((^(?F)), D) :- degree(F, D1), D is D1+1.
degree((^(F1 $ F2)), D) :- degree(F1, D1), degree(F2, D2), max(D1,D2,D3), D is D3+1.
degree((~(*F3)), D) :- degree(F3, D1), D is D1+1.
degree((~(@F)), D) :- degree(F, D1), D is D1+1.
degree((~(?F)), D) :- degree(F, D1), D is D1+1.
degree((~(F1 $ F2)), D) :- degree(F1, D1), degree(F2, D2), max(D1,D2,D3), D is D3+1.
degree(F, D) :- atom(F), D = 0.

max(D1, D2, D) :- (D1 > D2, D = D1; D = D2).


%F=\/(\varphi_0 & AF\varphi_1 & AG\varphi_2) %/ normal form (NF)
decideNF(F1 \/ F2, N) :- subterm(F1), decideNF(F2, N).%/
decideNF(F, N) :- subterm(F), N=true.

subterm(F) :- term(F).
subterm(F1 & F2) :- term(F1), subterm(F2).
subterm((~(@F))) :- wff(F).
subterm((~(?F))) :- wff(F).
subterm(F) :- term(F).

term(F) :- (wff(F); agterm(F); afterm(F)).

wff(P) :- atom(P), !.
wff(-P) :- wff(P), !.
wff(P->Q) :- wff(P), wff(Q),!.%/
wff(P&Q) :- wff(P), wff(Q), !.
wff(P\/Q) :- wff(P), wff(Q). %/

%agterm(F) :- wff(F).
agterm((~(@F))) :- wff(F).
%afterm(F) :- wff(F).
afterm((~(?F))) :- wff(F).


toClassicalTerm(F1 \/ F2, CT) :- cSubTerm(F1, CT1), toClassicalTerm(F2, CT2), append([CT1], CT2, CT). %/
toClassicalTerm(F, CT) :- cSubTerm(F, CT1), CT=[CT1].

cSubTerm(F, CT) :- wff(F), X1 = nf(cpl, F), CT=[X1].
cSubTerm(F1 & F2, CT) :- (wff(F1), X1=nf(cpl, F1);
	(agterm(F1), F1=(~(@F)), X1=nf(ag, F); F1=(~(?F)), X1 = nf(af, F))), cSubTerm(F2, CT1), append([X1], CT1, CT).
cSubTerm((~(@F)), CT) :- X1=nf(ag, F), CT=[X1].
cSubTerm((~(?F)), CT) :- X1=nf(af, F), CT=[X1].
cSubTerm(F, CT) :- X1 = nf(cpl, F), CT=[X1].


ctlforget(F, V, Filename) :- gain_prop(F, P), toClassicalTerm(F, CT), ctlforgetT(CT, V, R2), toCTL(R2, CTL),
	string_concat(Filename, '/', File1), 
	string_concat(File1, 'result', File2),
	string_concat(File2, '.txt', Filename1), 
	open(Filename1, write, Str), 
	write("CTL="), write(CTL),
	write("\n"),
	write("Exist="), write(Exist),
	write(Str, CTL), nl(Str).

toCTL([], []).
toCTL([L], R) :- toCTL_term(L, R).
toCTL([H|T], R) :- toCTL_term(H, R1), toCTL(T, R2), R = ((R1) \/ (R2)). %/

toCTL_term([],[]).
toCTL_term([NF], R) :- oneTerm(NF, R).
toCTL_term([H|T], R) :- oneTerm(H, R1), toCTL_term(T, R2), 
	(R1 = true, R2 = true, R= true; (R1 = true, R = R2; (R2 = true, R = R1; R = (R1 & R2)))). %, R = (R1 & R2).

oneTerm(NF, R) :- NF=nf(X, F),
	(X=cpl, X1=F; (X=ag, X1=(~(@F)); X1=(~(?F)))), (F=true, R=F; R =X1).

ctlforgetT([],V, []).
ctlforgetT([H|T], V, R) :-  ctlforget_term(H, V, R1), ctlforgetT(T,V,R2), append([R1], R2, R).


ctlforget_term(H, V, R) :- H=[nf(cpl, F1), nf(ag, F2), nf(af, F3)],
	ctlforget_one(nf(cpl, F1), V, X1),
	ctlforget_one(nf(ag, F2), V, X2),
	ctlforget_one(nf(af, F2 & F3), V, X3),
	R=[X1, X2, X3].
	%ctlforget_one(H, V, C1), ctlforget_term(T, V, R1), append([C1], R1, R).

ctlforget_one(NF, V, C) :- NF=nf(L, F), cplforget_one(F, V, WFF), C = nf(L, WFF).

cplforget_one(F, V, WFF) :- wff2cnf(F, CNF), forget(CNF, V, B), cnf2wff(B, WFF).



%%%%%%%classical forget
forget(A, [], A) :- !.
forget(A, PL, NewA) :- select(A, PL, P), difference(PL, P, NewPL),
    forget_one(A, P, NewA1), !, forget(NewA1, NewPL, NewA).

/* forget_one(CNF1, P, CNF2): CNF2 is the result of forgetting P in CNF1 */

forget_one(CNF, P, NewCNF) :- split(CNF, P, CNF1, CNF2, CNF3), 
    resolve(CNF2, P, CNF21),
    resolve(CNF1, -P, CNF11),
    cnf_simplifies(CNF21,CNF22),
    cnf_simplifies(CNF11,CNF12),
    cnf_disjunct(CNF12,CNF22,L),
    ((L=[], NewCNF=CNF3);
     (union(L, CNF3, CNF4),
      cnf_simplifies(CNF4, NewCNF))), !. 

/* split(CNF, P, CNF1, CNF2, CNF3): if CNF3 is the set of cluases in CNF that
   do not mention P and -P, CNF1 is those mention P, CNF2
   is those mention -P */

split([], _, [], [], []).
split([Clause|CNF], P, [Clause | CNF1], CNF2, CNF3) :- 
    member(P, Clause), !,
    split(CNF, P, CNF1, CNF2, CNF3).
split([Clause|CNF], P, CNF1, [Clause | CNF2], CNF3) :- 
    member(-P, Clause), !,
    split(CNF, P, CNF1, CNF2, CNF3).
split([Clause|CNF], P, CNF1, CNF2, [Clause | CNF3]) :- 
    split(CNF, P, CNF1, CNF2, CNF3).

/* cnf_disjunct(CNF1,CNF2,CNF3): if CNF3 is the CNF of CNF1 v CNF2 */

cnf_disjunct([],CNF2,[]) :- !.
cnf_disjunct(CNF1,[],[]) :- !.
cnf_disjunct([[]],CNF2,CNF2) :- !.
cnf_disjunct(CNF1,[[]],CNF1) :- !.
cnf_disjunct(CNF1,CNF2,CNF3) :- 
	setof(X, Y^Z^(member(Y,CNF1), member(Z,CNF2), union(Y,Z,X)), CNF3).

/* resolve a literal in a CNF */

resolve([],_,[]).
resolve([C|CNF], P, NewCNF) :- member(P,C), !, 
    resolve(CNF,P,NewCNF).
resolve([C|CNF], P, [NewC|NewCNF]) :- complement(P,P1), 
    difference(C,P1,NewC),
    resolve(CNF,P,NewCNF).


complement(true,false) :- !.
complement(false,true) :- !.
complement(P,-P) :- prop(P).
complement(-P,P) :- prop(P).



prop(true).
prop(false).




/* difference(L,P,L1): L1 is L - [P] */

difference([P | L], P, L) :- !.
difference([Q | L], P, [Q | NewL]) :- difference(L,P,NewL).
difference([],P,[]).



/* cnf_simplifies(L,NewL): simplifies DNF L into L1 by eliminating repetitive
   literals and contradictive literals. It also eliminates clauses that
   subsumes others.
*/

cnf_simplifies(CNF, NewCNF) :- cnf_sort_simpl(CNF, CNF1),
    unit_elim(CNF1,CNF2), 
    cnf_sort_simpl(CNF2, CNF3), 
    check_contr(CNF3,CNF4),
    cnf_minimize(CNF4, NewCNF).

cnf_simplifies(CNF, NewCNF) :- cnf_sort_simpl(CNF, CNF1),
    unit_elim(CNF1,CNF2), 
    cnf_sort_simpl(CNF2, CNF3), 
    check_contr(CNF3,CNF4),
    cnf_minimize(CNF4, NewCNF).

/* sort each clause and get rid of tautologies */

cnf_sort_simpl([],[]).
cnf_sort_simpl(CNF,[[]]) :- member([],CNF), !.
cnf_sort_simpl(CNF,[[]]) :- member([false],CNF), !.
cnf_sort_simpl(CNF,[[]]) :- member([-true],CNF), !.
cnf_sort_simpl([C|CNF], NewCNF) :- member(true,C), !,  
    cnf_sort_simpl(CNF, NewCNF).
cnf_sort_simpl([C|CNF], NewCNF) :- member(-P,C), member(P,C), !,
    cnf_sort_simpl(CNF, NewCNF).
cnf_sort_simpl([C|CNF], NewCNF) :- 
    sort(C,NewC1), subtract(NewC1, [false,-true], NewC),
    ((NewC=[], NewCNF=[[]]);
     (cnf_sort_simpl(CNF, NewCNF1), NewCNF=[NewC|NewCNF1])), !.

/* eliminates all unit clauses */

unit_elim(CNF,[[P]|NewCNF]) :- 
    member([P],CNF), 
    resolve(CNF,P,CNF2), !,
    unit_elim(CNF2,NewCNF), !.
unit_elim(CNF,CNF) :- !.

/* eliminates all unit clauses (keep the unit clauses deduced as a sep. arg. */

unit_elim(CNF,NewCNF,[P|Unit]) :- 
    member([P],CNF), 
    resolve(CNF,P,CNF2), !,
    unit_elim(CNF2,NewCNF,Unit), !.
unit_elim(CNF,CNF,[]) :- !.


/* check for contradiction */

check_contr(CNF,[[]]) :- member([],CNF), !.
check_contr(CNF,CNF).



/* get rid of clauses that are subsumed by others. */

cnf_minimize(DNF,DNF1) :- member(L,DNF), member(L1,DNF), L\==L1,
    cnf_subsumes(L,L1), difference(DNF, L1, DNF2), cnf_minimize(DNF2,DNF1), !.
cnf_minimize(DNF,DNF).


/* cnf_subsumes(L,L1) if L subsumes L1: L is a subset of L1. */

cnf_subsumes(L,L1) :- \+ (member(X,L), X\==false, X\==(-true),
    \+ member(X,L1)).



/* select first those that are already entailed */

select(CNF, PL, P) :- member(P,PL), member([P],CNF), !.
select(CNF, PL, P) :- member(P,PL), member([-P],CNF), !.

/* select in sequence given afterward */

select(_, [P|_], P) :- !.



	
%#############
%wff2cnf


/* dnf2wff(L,W) convert a DNF list to a formula */

dnf2wff([],false) :- !.
dnf2wff([[]],true) :- !.
dnf2wff([L], W) :- list2conjunct(L,W), !.
dnf2wff([L1|L2], W1 \/ W2) :- list2conjunct(L1, W1), dnf2wff(L2,W2).

/* cnf2wff(L,W) convert a CNF list to a formula */

cnf2wff([[]],false) :- !.
cnf2wff([],true) :- !.
cnf2wff([L], W) :- list2disjunct(L,W), !.
cnf2wff([L1|L2], W1 & W2) :- list2disjunct(L1, W1), cnf2wff(L2,W2).


/* list2conjunct(L,A): A is the conjunction of all formulas in L */

list2conjunct([P],P) :- !.
list2conjunct([P | L],P & W) :- list2conjunct(L,W).


/* list2disjunct(L,A): A is the disjunction of all formulas in L: A=false when
   L is empty */

list2disjunct(L, true) :- member(true,L), !.
list2disjunct(L, true) :- member(-P,L), member(P,L), !.
list2disjunct(L, W) :- sort(L,L1), delete(L1,false,L2), list2disj(L2,W), !.
list2disj([], false) :- !.
list2disj([P],P) :- !.
list2disj([P | L],P \/ W) :- list2disj(L,W).%/


/* wff2dnf transforms a wff to its disjunctive normal form in list format */

wff2dnf(P,[[P]]) :- prop(P), !.
wff2dnf(-P,[[-P]]) :- prop(P), !.
wff2dnf(P->Q, L) :- wff2dnf(-P\/Q, L).%/
wff2dnf(P<->Q, L) :- wff2dnf((P->Q)&(Q->P), L).
wff2dnf(P\/Q, L) :- wff2dnf(P,L1), wff2dnf(Q,L2), append(L1,L2,L).%/
wff2dnf(P&Q, L) :- wff2dnf(P,L1), wff2dnf(Q,L2),
    findall(X, (member(Y,L1), member(Z,L2), append(Y,Z,X)), L).
wff2dnf(-P, L) :- wff2cnf(P,L1), negate(L1,L).


/* negate(L1,L): negate L1 in DNF (CNF) and make it into a CNF (DNF). */
negate([],[]) :- !.
negate([[]],[[]]) :- !.
negate(L1,L) :- bagof(X, Y^(member(Y, L1), negate_one_clause(Y,X)), L).

/* negate_one_clause(L1, L2): make all elements in L1 into its complement */
negate_one_clause(L1, L2) :- bagof(X, Y^(member(Y,L1), complement(Y,X)), L2).



%----------wff2cnf------and--------wff2dnf--------
wff2cnf(P, [[P]]) :- prop(P), !.
wff2cnf(-P, [[-P]]) :- prop(P), !.
wff2cnf(P->Q, L) :- wff2cnf(-P\/Q, L), !.%/
wff2cnf(P<->Q, L) :- wff2cnf((-P\/Q)&(-Q\/P), L), !.%/
wff2cnf(P&Q, L) :- wff2cnf(P,L1), wff2cnf(Q,L2), append(L1,L2,L), !.
wff2cnf(P\/Q, L) :- wff2cnf(P,L1), wff2cnf(Q,L2), %/
    findall(X, (member(Y,L1), member(Z,L2), append(Y,Z,X)), L), !.
wff2cnf(-P, L) :- wff2dnf(P, L1), negate(L1, L), !.

/* wff2cnf(W,NewW) :- negation_inside(W,W1), wff2cnf0(W1,NewW).
*/
wff2cnf0(P, [[P]]) :- prop(P), !.
wff2cnf0(-P, [[-P]]) :- prop(P), !.
wff2cnf0(P&Q, L) :- wff2cnf0(P,L1), wff2cnf0(Q,L2), union(L1,L2,L), !.
wff2cnf0(P\/Q, L) :- wff2cnf0(P,L1), wff2cnf0(Q,L2), %/
    setof(X, Y^Z^(member(Y,L1), member(Z,L2), union(Y,Z,X)), L), !.
	

%gain all the propositions in a CTL formula
gain_prop(P & Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), 
	append(L1, L2, L).
gain_prop(P \/ Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), %/
	append(L1, L2, L).
gain_prop(P -> Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), append(L1, L2, L).
gain_prop(P <-> Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), append(L1, L2, L).
gain_prop(^(@P), L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(-P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(^(P $ Q), L) :- gain_prop(P, L1), gain_prop(Q, L2), append(L1, L2, L3), sort(L3, L).
gain_prop(^(*P), L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(^(?P), L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(~(@P), L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(~(P $ Q), L) :- gain_prop(P, L1), gain_prop(Q, L2), append(L1, L2, L3), sort(L3, L).
gain_prop(~(*P), L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(~(?P), L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(P, [P]) :-  atom(P), !, assert(prop(P)).