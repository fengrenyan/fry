

%forgetting the set V of atoms from a CTL formula.

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



prop(start).



ctl2NNF(P, P) :- prop(P). %Changing a CTL formula into its NNF。
ctl2NNF(- P, - P) :- prop(P).
ctl2NNF(- (CTL1 \/ CTL2), NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(- CTL2, NNF2), NNF = (NNF1 & NNF2). %/
ctl2NNF(- (CTL1 & CTL2), NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(- CTL2, NNF2), NNF = (NNF1 \/ NNF2). %/
ctl2NNF(- (CTL1 -> CTL2), NNF) :- ctl2NNF(CTL1, NNF1), ctl2NNF(- CTL2, NNF2), NNF = (NNF1 & NNF2). 
ctl2NNF(- (CTL1 <-> CTL2), NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF3 = (NNF1 & NNF2),
	ctl2NNF(CTL1, NNF4), ctl2NNF(- CTL2, NNF5), NNF6 = (NNF5 & NNF4), 
	(NNF6 = NNF3, NNF = NNF3; NNF = (NNF3 \/ NNF6)). %/
ctl2NNF(-(- CTL), NNF) :- ctl2NNF(CTL, NNF).	

ctl2NNF(CTL1 \/ CTL2, NNF) :- ctl2NNF(CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF = (NNF1 \/ NNF2). %/
ctl2NNF(CTL1 & CTL2, NNF) :- ctl2NNF(CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF = (NNF1 & NNF2).
ctl2NNF(CTL1 -> CTL2, NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF = (NNF1 \/ NNF2). %/
ctl2NNF(CTL1 <-> CTL2, NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF3 = (NNF1 \/ NNF2), %/
	ctl2NNF(CTL1, NNF4), ctl2NNF(- CTL2, NNF5), NNF6 = (NNF5 \/ NNF4), %/
	(NNF6 = NNF3, NNF = NNF3; NNF = (NNF3 & NNF6)). 


ctl2NNF(~(* CTL), NNF) :- ctl2NNF(CTL, NNF1), NNF = (~(* NNF1)).
ctl2NNF(~(? CTL), NNF) :- ctl2NNF(CTL, NNF1), NNF = (~(? NNF1)).
ctl2NNF(~(@ CTL), NNF) :- ctl2NNF(CTL, NNF1), NNF = (~(@ NNF1)).
ctl2NNF(~(CTL1 $ CTL2), NNF) :- ctl2NNF(CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF = (~(CTL1 $ CTL2)).

ctl2NNF(^(*CTL), NNF) :- ctl2NNF(CTL, NNF1), NNF = (^(* NNF1)).
ctl2NNF(^(? CTL), NNF) :- ctl2NNF(CTL, NNF1), NNF = (^(? NNF1)).
ctl2NNF(^(@ CTL), NNF) :- ctl2NNF(CTL, NNF1), NNF = (^(@ NNF1)).
ctl2NNF(^(CTL1 $ CTL2), NNF) :- ctl2NNF(CTL1, NNF1), ctl2NNF(CTL2, NNF2), NNF = (^(CTL1 $ CTL2)).

ctl2NNF(-(~(* CTL)), NNF) :- ctl2NNF(^(*(- CTL)), NNF).
ctl2NNF(-(~(? CTL)), NNF) :- ctl2NNF(^(@(- CTL)), NNF).
ctl2NNF(-(~(@ CTL)), NNF) :- ctl2NNF(^(?(- CTL)), NNF).
%ctl2NNF(-(~(CTL1 $ CTL2)), NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(- CTL2, NNF2), NNF = (^(CTL1 $w$ CTL2)).

ctl2NNF(-(^(* CTL)), NNF) :- ctl2NNF(~(*(- CTL)), NNF).
ctl2NNF(-(^(? CTL)), NNF) :- ctl2NNF(~(@(- CTL)), NNF).
ctl2NNF(-(^(@ CTL)), NNF) :- ctl2NNF(~(?(- CTL)), NNF).
%ctl2NNF(-(~(CTL1 w CTL2)), NNF) :- ctl2NNF(- CTL1, NNF1), ctl2NNF(- CTL2, NNF2), NNF = (~(CTL1 $ CTL2)).

ctl2NNF(CTL, CTL).

 
split2ST([],[],[],[], N).
split2ST([H|L], SL, TL, W, N) :- H = snf_clause(Type, _, _, _),
	(Type = af, H1 =H, N1 is N +1;
		(Type = ef, H1 = H, N1 is N +1; H1 = [], N1 is N
		)
	), split2ST(L, SL1, TL1, W1, N1),
	(H1 =[], TL= TL1, append([H], SL1, SL), X=[], W = W1; 
		string_concat('w', N, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), SL=SL1, append([H], TL1, TL), append([Y1], W1, W)
	). 

 
appearing([], _, []).
appearing([H|L], P, R) :- (posC(P, H), H1 = 1; (negC(P, H), H1=1; H1 =[])),
	appearing(L, P, L1),
	(H1 =[], R = L1; append([H], L1, R)).


 
appearing_list([], _, []).
appearing_list(L, [], L).
appearing_list(L, [P|V], R) :- appearing(L, P, R1), appearing_list(L, V, R2), append(R1, R2, R3), sort(R3, R).

%pos(P,C)-------------------------------------------start-------------------------------------------
%deciding p (that should be resolved) appearing in C positively.
posC(P, C) :- C = snf_clause(_, _, H, _), member(P, H).
posC(P, C) :- C = snf_clause(_, _, H, _,_), member(P, H).
%pos-------------------------------------------end-------------------------------------------


%neg-------------------------------------------start-------------------------------------------
%deciding p (that should be resolved) appearing in C negatively.
negC(P, C) :- C = snf_clause(_, _, H, _), F=(-P), member(F, H).
negC(P, C) :- C = snf_clause(_, _, H, _,_), F=(-P), member(F, H).
%neg-------------------------------------------end-------------------------------------------


%all_PosC(L, P, L1)-----------------------------------------start-------------------------------------------
%find all the clauses C in L that P appearing in C positively
all_PosC([], _, []).
all_PosC([H|T], P, L1) :- all_PosC(T, P, Tem), (posC(P, H), append([H], Tem, X), L1=X; L1=Tem).
%all_PosC(L, P, R) :- findall(, posC(P, 
%all_PosC-------------------------------end-------------------------------------


%all_NegC(L, P, L1)-----------------------------------------start-------------------------------------------
%find all the clauses C in L that P appearing in C negatively
all_NegC([], _, []).
all_NegC([H|T], P, L1) :- all_NegC(T, P, Tem), (negC(P, H), append([H], Tem, X), L1=X; L1=Tem).
%all_PosC-------------------------------end-------------------------------------






%------------------------------------start---transCTL2SNF----------------------------------------------

%tran2SNFCl(F, V1, L): transform a CTL formula into a set of snf clauses

tran2SNFClt([], []).
tran2SNFClt([H|L],R) :- tranH2Ct(H, C), tran2SNFClt(L, R1), append([C], R1, R).

tranH2Ct(start -> P, C) :- prop(P), C = snf_clause(start,[start], [P], nil).
tranH2Ct(Q -> P, C) :- is_dis(P), negation([Q], Q1),  wff2cnf(P, P1), P1 = [P2], append(Q1, P2, P3), C = snf_clause(true,[true], P3, nil).
tranH2Ct([Q -> ^(*P),I], C) :- wff2cnf(P, P1), P1 = [P2], C = snf_clause(ex,[Q], P2, I).
tranH2Ct([Q -> ^(?P),I], C) :- wff2cnf(P, P1), P1 = [P2], C = snf_clause(ef,[Q], P2, I).
tranH2Ct(Q -> ~(*P), C) :- wff2cnf(P, P1), P1 = [P2], C = snf_clause(ax, [Q], P2, nil).
tranH2Ct(Q -> ~(?P), C) :- wff2cnf(P, P1), P1 = [P2], C = snf_clause(af,[Q], P2, nil).


tran2SNFCl([], []).
tran2SNFCl([H|L],R) :- tranH2C(H, C), tran2SNFCl(L, R1), append([C], R1, R).

tranH2C(start -> P, C) :- prop(P), C = start_clause([start], [P]).
tranH2C(Q -> P, C) :- is_dis(P), negation([Q], Q1),  wff2cnf(P, P1), P1 = [P2], append(Q1, P2, P3), C = true_clause([true], P3).
tranH2C([Q -> ^(*P),I], C) :- wff2cnf(P, P1), P1 = [P2], C = exist_next_clause([Q], P2, I).
tranH2C([Q -> ^(?P),I], C) :- wff2cnf(P, P1), P1 = [P2], C = exist_future_clause([Q], P2, I).
tranH2C(Q -> ~(*P), C) :- wff2cnf(P, P1), P1 = [P2], C = global_next_clause([Q], P2).
tranH2C(Q -> ~(?P), C) :- wff2cnf(P, P1), P1 = [P2], C = global_future_clause([Q], P2).


tran2SNF([], Ind, N, V, 0, []) :- V=[].
tran2SNF(L, Ind, N, V, IndM, R) :- tran2SNF_list(L, V1, Ind, N, IndM1, NM, L1),  
	(L1 = L, R = L1, V = V1, IndM = IndM1;
		tran2SNF(L1, IndM1, NM, V2, IndM, R), append(V1, V2, V)).

tran2SNF_list([], V, Ind, N, Ind, N, []) :- V=[].
tran2SNF_list([H|T], V1, Ind, N, IndM, NM, L) :- transF(H, Ind, V2, N, IndM1, NM1, L1), 
	tran2SNF_list(T, V3, IndM1, NM1, IndM, NM, L3), append(V2, V3, V1), append(L1, L3, L4), sort(L4, L).
	

%transF(start -> P, Ind, V, Ind, IndM, IndM, C) :-  C = start_clause([start], [P]).
%transF(Q -> P, Ind, V, Ind, IndM, IndM, C) :- is_dis(P), negation([Q], Q1),  wff2cnf(P, P1), P1 = [P2], append(Q1, P2, P3), %C = true_clause([true], P3).
transF(Q -> (P1 & P2), Ind, V, N, IndM, NM, L) :- NM = N, IndM=Ind, L = [Q -> P1, Q -> P2], V=[].
transF(Q -> (P1 \/ P2), Ind, V, N, IndM, NM, L) :-  IndM=Ind,%/ 
	(is_dis(P1), !, L1=[], N1=N, V2=[];
		(N1 is N +1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), V2=[Y1], L1=(Y1 -> P1))
	),
	(is_dis(P2), !, N2=N1, NM=N1, V3=[], L2=[];
		(N2 is N1 +1, NM=N2, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L1=[], (L2=[], L =[Q -> (P1 \/ P2)], V = [];  %/
				L = [L2, Q -> P1 \/ Y2], V = [Y2]); %/
		(L2 =[], L = [L1, Q -> Y1 \/P2], V = [Y1];  %/
			L = [L1, L2, Q -> Y1 \/ Y2], V =[Y1, Y2] %/
		)
	).
	
	

	
transF(Q -> ^(P1 $ P2), Ind, V, N, IndM, NM, L) :-   IndM is Ind+1,
	(is_lit(P2), !, N2=N, V3=[], L2=[];
		(N2 is N +1, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L2 =[], tran2U(Q -> ^(P1 $ P2), N2, IndM, Vx, L), V = [Vx], NM is N + 1;
			 NM is N2 + 1, tran2U(Q -> ^(P1 $ Y2), N2, IndM, Vx, R), append([L2], R, L), V = [Y2, Vx]
	).
	

	
transF(Q -> ~(P1 $ P2), Ind, V, N, IndM, NM, L) :-   IndM is Ind,
	(is_lit(P2), !, N2=N, V3=[], L2=[];
		(N2 is N +1, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L2 =[], tran2Uall(Q -> ~(P1 $ P2), N2, IndM, Vx, L), V = [Vx], NM is N + 1;
			 NM is N2 + 1, tran2Uall(Q -> ~(P1 $ Y2), N2, IndM, Vx, R), append([L2], R, L), V = [Y2, Vx]
	).
	
tran2U(Q -> ^(Y1 $ Y2), N, IndM, V, L) :- N1 is N+1, string_concat('x', N1, XY3), 
	string_to_atom(XY3, Y3), assert(prop(Y3)), V= Y3,
	L = [Q -> Y2 \/ Y3, Y3 -> Y1, [Y3 -> ^(*(Y2 \/ Y3)), IndM], [Q -> ^(? Y2),IndM]]. %/

tran2Uall(Q -> ~(Y1 $ Y2), N, IndM, V, L) :- N1 is N+1, string_concat('x', N1, XY3), 
	string_to_atom(XY3, Y3), assert(prop(Y3)), V= Y3,
	L = [Q -> Y2 \/ Y3, Y3 -> Y1, Y3 -> ~(*(Y2 \/ Y3)), Q -> ~(? Y2)]. %/	
	
	

transF(Q -> ^(@P), Ind, V, N, IndM, NM, L) :-  N1 is N+1, NM=N1, IndM is Ind+1,
	string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)),
	V=[Y1], L=[Q -> Y1, Y1 -> P, [Y1 -> ^(*Y1), IndM]].

	
transF(Q -> ~(@P), Ind, V, N, IndM, NM, L) :-  N1 is N+1, NM=N1, IndM is Ind,
	string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)),
	V=[Y1], L=[Q -> Y1, Y1 -> P, Y1 -> ~(*Y1)].
	
transF(Q -> ^(*P), Ind, V, N, IndM, NM, L) :- IndM is Ind +1,
	(is_dis(P), L = [[Q -> ^(*P), IndM]], NM is N, V=[];
		N1 is N+1, NM=N1,  string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[[Q -> ^(* Y1),IndM], Y1 -> P]
	).


transF(Q -> ~(*P), Ind, V, N, IndM, NM, L) :- IndM is Ind,
	(is_dis(P), L = [Q -> ~(*P)], NM is N, V=[];
		N1 is N+1, NM=N1,  string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[Q -> ~(* Y1), Y1 -> P]
	).
	
transF(Q -> ^(?P), Ind, V, N, IndM, NM, L) :- IndM is Ind+1, 
	(is_lit(P), L = [[Q -> ^(?P), IndM]], NM is N, V=[];
		N1 is N+1, NM=N1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[[Q -> ^(? Y1), IndM], Y1 -> P]
	).


transF(Q -> ~(?P), Ind, V, N, IndM, NM, L) :- IndM is Ind, 
	(is_lit(P), L = [Q -> ~(?P)], NM is N, V=[];
		N1 is N+1, NM=N1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[Q -> ~(? Y1), Y1 -> P]
	).
	
	
transF(Q -> P, Ind, V, N, IndM, NM, L) :- IndM=Ind, is_lit(P), V =[], NM =N, L=[Q -> P].
transF(Q, Ind, V, N, IndM, NM, [Q]) :- NM =N, IndM=Ind, V = [].

%------------------------------------end---transCTL2SNF----------------------------------------------


 
 is_lit(C) :- (prop(C); equall(-C, D), prop(D)).
 
 
 is_dis(C) :- is_lit(C).
 is_dis(C) :- C = (P \/ Q), %/
	(is_dis(P), is_dis(Q), !; !, fail).



equall(P, P) :- prop(P), !.
equall(-P, -P) :- prop(P), !.
equall(-(-P), P) :- prop(P), !.


 
negation([],[]).
negation([H|T], [H1|T1]) :- 
	equall(- H, H1), 
	negation(T, T1).
	
	
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


%decide whether an experssion is a CTL formula.
is_CTLformula(P) :- atom(P).
is_CTLformula(-P) :- is_CTLformula(P).
is_CTLformula(P & Q) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(P \/ Q) :- is_CTLformula(P), is_CTLformula(Q).%/
is_CTLformula(P -> Q) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(^(@P)) :- is_CTLformula(P).
is_CTLformula(^(P $ Q)) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(^(*P)) :- is_CTLformula(P).
is_CTLformula(^(?P)) :- is_CTLformula(P).
is_CTLformula(~(@P)) :- is_CTLformula(P).
is_CTLformula(~(P $ Q)) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(~(*P)) :- is_CTLformula(P).
is_CTLformula(~(?P)) :- is_CTLformula(P).




find_Asometime_clauses([], []).
find_Asometime_clauses([C|L1], L2) :- C = snf_clause(Type, X, Y, Z),
	(Type=ex, C1=snf_clause(Type, X, Y, Z,0);
		(Type=ax, C1 =snf_clause(Type, X, Y, Z,0);
			(Type=true, C1 =snf_clause(Type, X, Y, Z,0); C1 =[])
		)
	), find_Asometime_clauses(L1, L3),
	(C1=[], L2=L3; append([C1], L3, L2)).
	
	

find_Esometime_clauses([], Ind, []).
find_Esometime_clauses([H|T], Ind, R) :- H = snf_clause(Type, X, Y, Z), 
	(Type = true, !, H1=snf_clause(Type, X, Y, Z,0); 
		(Type = ax, !, H1=snf_clause(Type, X, Y, Z,0);
			(Type = ex, !, (H =snf_clause(ex, _, _, Ind), H1=snf_clause(Type, X, Y, Z,0); H1=[]); H1=[])
		)
	),
	find_Esometime_clauses(T, Ind, T1),
	(H1=[], R = T1; append([H1], T1, R)).



%-------------------------------------------E-loop---start------------------------------------------------
loopFormula([], _, _, _, false).  
loopFormula(L, P, V, Init, F) :- gain_initial_clause(Init, P, C), 
	append(C, L, L1), findFormula(L1, V, H1),
	(H1=[], F=alse;  
		(H1=Init, F = H1; (H1 = ture, F=H1; (H1= false, F =false; loopFormula(L, P, V, H1, F))))
	).
	
findFormula(L, V, H) :- loopRes(L, V, R), %write("\nR========="), write(R), %trace, 
	findall(X, (member(X, R), X = snf_clause(true, _, _, _, N), (N =1; N = 2)), R1),
 	subsumeLoop(R1, R1, R1, R2), %write("\nR2========="), write(R2),
	gainH(R2, H).


gainH([], []).
gainH([H|T], [H2|R]) :- H = snf_clause(Type, _, X, L, N),  H2 = X, 
	gainH(T, R).
	
/*gainH([], []).
gainH([H|T], R) :- H = snf_clause(Type, _, X, L, N), 
	(Type = true, ((N =1; N = 2), H2 = X; H2= []); !), 
	gainH(T, T1),
	(H2 = [], R = T1; append([H2], T1, R)).*/

gain_initial_clause([],_, []).
gain_initial_clause([C|L], P, [C1|L1]) :- append([P], C, T), sort(T, T1), C1 = snf_clause(ax, [], T1, nil, 2), 
	gain_initial_clause(L, P, L1).
	
%-------------------------------------------E-loop---end-------------------------------------------------



%-------------------------simplify-----------------start-------------------------------------------

simp_disSet([],[]).  
simp_disSet(L, L1) :- (member(P, L), member(-P, L), L1= [true]; 
	(member(true, L), L1=[true]; 
		(member(false, L), delete(L, false, L2), sort(L2, L1); sort(L, L1))
	)
	).  

simp_conSet([], [true]). 
simp_conSet(L, L1) :- (member(P, L), member(-P, L), L1= [false]; 
	(member(false, L), L1 =[false]; 
		(member(true, L), delete(L, true, L2), sort(L2, L1); sort(L, L1))
	)
	).  

elimElement(H, [], []).
elimElement(H, [H1|T1], L1) :- elimElement(H, T1, L2), 
	(H = H1, L1=L2; append([H1], L2, L1)).
	
	
 
simplySNFC(C, C1) :- C= snf_clause(Type, H, T, L), 
	(Type = true, simp_disSet(T, T1), C1=snf_clause(Type, H, T1, L);
	simp_conSet(H, H1), simp_disSet(T, T1), C1=snf_clause(Type, H1, T1, L)). 
 	
simplySNFCLoop(C, C1) :- C= snf_clause(Type, H, T, L, N),
	(Type = true, !, simp_disSet(T, T1), C1=snf_clause(Type, H, T1, L,N);
	simp_conSet(H, H1), simp_disSet(T, T1), C1=snf_clause(Type, H1, T1, L,N)). 


 
%subsume rules
subsume1(X, Y) :- X = snf_clause(_, H1, T1, _), Y = snf_clause(_, H2, T2, _), subset(T1, T2).  

subsume2(X, Y) :- X = snf_clause(_, H1, T1, _), Y = snf_clause(_, H2, T2, _), subset(T1, T2), subset(H1, H2).

subsume3(X, Y) :- X = snf_clause(_, H1, T1, _, _), Y = snf_clause(_, H2, T2, _,_), subset(T1, T2).


subsumeLoop([], _, L, L).
subsumeLoop(_, [], [], []).
subsumeLoop([H|L], X, Y, R) :- subsumeLoop_list(H, X, Y, L1), 
	(L1 = X, subsumeLoop(L, L1, L1, R); subsumeLoop(L1, L1, L1, R)).


subsumeLoop_list(C, [], L, L).
subsumeLoop_list(C, [H|T], L, R) :- (C = H, subsumeLoop_list(C, T, L, R); 
    ((subsume3(C, H), delete(L, H, L1); L1 = L),
	    subsumeLoop_list(C, T, L1, R)
    )).




subsume([], _, L, L).
subsume(_, [], [], []).
subsume([H|L], X, Y, R) :- subsume_list(H, X, Y, L1),  
	(L1=X, subsume(L, L1, L1, R); subsume(L1, L1, L1, R)).


subsume_list(C, [], L, L).
subsume_list(C, [H|T], L, R) :- (C = H, subsume_list(C, T, L, R); 
(C = snf_clause(Type1,_,_,_), H=snf_clause(Type2,_,_,_), %judgeType(C, Type1), judgeType(H, Type2), 
	(Type1 = ax, !, (Type2=ax, !, (subsume2(C, H), delete(L, H, L1); L1 = L); (Type2 =ex, !, (subsume2(C, H), delete(L, H, L1); L1 = L ); (Type2=true, !, (subsume1(H, C), delete(L, C, L1); L1 = L); L1=L))); 
		(Type1=ex, !, (Type2=ax, !,(subsume2(H, C), delete(L, C, L1); L1 = L); (Type2 =ex, !, (subsume2(C, H), delete(L, H, L1); L1 = L); (Type2=true, !, (subsume1(H, C), delete(L, C, L1); L1 = L); L1=L)));
			(Type1 =true, !, ((Type2=start; Type2=ax; Type2=ex; Type2=true), (subsume1(C, H), delete(L, H, L1); L1 = L));
				/*(Type2=start, !, (subsume1(C, H), delete(H, L, L1); L1 = L); (Type2=ax, !, (subsume1(C, H), delete(H, L, L1); L1 = L); (Type2=ex, !, (subsume1(C, H), delete(H, L, L1); L1 = L); (Type2=true, !, (subsume1(C, H), delete(H, L, L1); L1 = L); L1 = L))));*/
				(Type1 =start, !, (Type2=start, !, (subsume1(C, H), delete(L, H, L1); L1 = L); (Type2=true, !, (subsume1(H, C), delete(L, C, L1); L1 = L); L1 = L)); L1 = L)
			)
		)
	),
	subsume_list(C, T, L1, R)
)).


%------------------------simplify-----------------end-------------------------------------------

obtainAtom([],[]).
obtainAtom([H|T], [H1|L1]) :- (prop(H), H1=H; H = - H1), obtainAtom(T, L1).


base(R) :- findall([P, L], (pair(P, L), assert(formula(snf, P, L))), R).
base(V1, R) :- findall([P, L], pair(P, L), R1), 
	(R1=[], R=[];
		deleteTnest(R1, V1, R)
	).
	
 
nestFormula([], _, []).
nestFormula([H|T], V1, [H1|L1]) :- H = formula(Type, X, L),
	(Type = con, nest_list(L, V1, L2), H1 = formula(Type, X, L2); H1 = H), nestFormula(T, V1, L1).

deleteTnest([], _, []).
deleteTnest([H|T], V1, [H1|L1]) :- H =[P, L], H1=P, nest_list(L, V1, L2), assert(formula(snf, P, L2)),
	deleteTnest(T, V1, L1).
	
nest_list([], _, []).
nest_list([H|T], V1, L) :- H=(Type, Q, I), obtainAtom(Q, Q1), intersection(V1, Q1, Q2), subtract(Q1, Q2, Q3),
	%(Type=true, obtainAtom(Q, Q1), (intersection(V1, Q1,[]), C1=H; C1=[]);
	(Type=true, (Q2=[], C1=H; C1=[]);
		((Type = ax; Type = ex), (Q2 =[], C1=H; (Q3 =[], C1=H; C1 =[]));  
			C1=H
		)   
	), nest_list(T, V1, L1),
	(C1=[], L=L1; append([C1], L1, L)).


sat(F, R) :- (gain_prop(F, P), ctl2NNF(F, F1), L=[start-> z, z -> F1], assert(prop(z)), 
	tran2SNF(L, 0, 0, V1, IndM, SNF), tran2SNFClt(SNF, LC), append(P, V1, V2), append([z], V2, V3), %res(LC, V3, R),
	split2ST(LC, SL, TL, W, 0),
	step_resolution(SL, V3, L2),
	(L2=false, write("\n **unsat \n"), R = false; 
		temp_resolution(TL, W, L2, V3, R1, TL1, W1, V3), append(V3, W, V5),
		append(L2, LC, L3), append(L3, R1, R2), sort(R2, R3), split2ST(R3, SL2, TL2, W2, 0), step_resolution(SL2, V5, R4),
		(R4 = false,  write("\n **unsat \n"), R = false; write("\n sat \n"), R=true)
	); write("This formula is error! Please input the right formula.")).
	
	
sAt(F) :- gain_prop(F, P1), L=[start-> z, z -> F], assert(prop(z)), sort(P1, P), write("atoms in F:"), write(P),
	tran2SNF(L, 0, 0, V1, IndM, SNF), write("\n New atoms:"), write(V1),
	tran2SNFClt(SNF, LC), write("\n SNF: "), write(SNF),
	write("\n Normal form:"), write(LC),
	append(P, V1, V2), append([z], V2, V3), res(LC, V3, R),
	(R=false, write("\n *****unsat \n"); write("\n -----------------result-----------------\n"),
	write(R),
	write("\n -----------------result-----------------\n"),
	write("\n ####sat \n")).
	
res(L, V3,  R) :- split2ST(L, SL, TL, W, 0), step_resolution(SL, V3, L2), 
	(L2=false, R=false; 
		temp_resolution(TL, W, L2, V3, R1, TL1, W1, V3), append(V3, W, V5), append(L2, L, L3), 
			append(L3, R1, R2), sort(R2, R3),
			(R3=L, R = R3; res(R3, V5, R))
	).



removeSing([],[]).
removeSing([H|L], R) :- H = snf_clause(Type, Q, T, I), length(Q, N),
	(N=1, H1=[]; H1=H),
	removeSing(L, R1),
	(H1=[], R=R1; append([H1], R1, R)).
	

%----------start--------output the time of solving a problem-----------
time_output_1(Goal, UsedInf,UsedTime,Wall) :-
	time_state(State0),
	call(Goal),
	report_output(State0, 11, UsedInf, UsedTime, Wall).

report_output(t(OldWall, OldTime, OldInferences), Sub, UsedInf,UsedTime,Wall) :-
	time_state(t(NewWall, NewTime, NewInferences)),
	UsedTime is NewTime - OldTime,
	UsedInf  is NewInferences - OldInferences - Sub,
	Wall     is NewWall - OldWall,
	(   UsedTime =:= 0
	->  Lips = 'Infinite'
	;   Lips is integer(UsedInf / UsedTime)
	),
	print_message(information, time(UsedInf, UsedTime, Wall)).


time_state(t(Wall, Time, Inferences)) :-
	get_time(Wall),
	statistics(cputime, Time),
	statistics(inferences, Inferences).
%------------end------output the time of solving a problem-----------	


satCall(Filename, P, PN, ResultF) :-
    set_prolog_stack(local, limit(137438953472)),  %16G
    set_prolog_stack(global, limit(137438953472)),
	string_concat(ResultF, P, F3), %'result_P1/Variable_'
	string_concat(F3, '/', F4),
	string_concat(F4, Filename, Filename1),
	string_concat(Filename1, '.txt', Filename2),
	open(Filename2, write, Str), 
	read_formula(Filename, Formula), 
	%gain_prop(F, P), 
	catch(time_output_1(sat(Formula, SAT), U1, V1, W1),E1, false), 
	(SAT = false, write(Str,'UsedTime:'), write(Str, V1), nl(Str); 	
		write("SAT"),!, gain_prop(Formula, PN1), sort(PN1, PN2), 
		length(PN2, LenF), 
		(LenF = P, P1 = PN1, Res = [], Ack =[], R = true, Exist = 1, LenSNF = 0, LenRes =0,
			U = U1, V=V1, W=W1, E=E1; 
			!, pick_props(P, PN1, P1),  %trace, write("trace forget"),  
			catch(time_output_1(new_ctlForget(Formula, P1, Res, Ack, R, Exist, LenSNF, LenRes), U, V, W),E, false)
		),
		write("P1="), write(P1), nl,
		write("T="), write(Formula), nl, 
		write(Str,'P1='), write(Str, P1), nl(Str),  
		write(Str,'T='), write(Str, Formula), nl(Str), 
		write(Str,'UsedInf:'), write(Str, U), nl(Str),  
		write(Str,'UsedTime:'), write(Str, V), nl(Str),  
		write(Str,'Wall:'), write(Str, W), nl(Str),  
		write(Str, 'Res===========: '), write(Str, Res), nl(Str),  
		write(Str, 'ArcM===========: '), write(Str, Ack), nl(Str),  
		write(Str, 'LenSNF='), write(Str, LenSNF), nl(Str),
		write(Str, 'LenRes='), write(Str, LenRes), nl(Str),
		write(Str, 'SAT: '), write(Str, SAT), nl(Str),  
		write(Str, 'Exist: '), write(Str, Exist), nl(Str),  
		write(Str,'...............................................'), nl(Str),retractall(pair(X,Y))
	), 
	retractall(prop(_)), 
	write(Str,'SAT:'), write(Str, SAT), nl(Str), 
	close(Str).
	
	
	



%------------------------read formula from file----------------
%read formula from file and assert the atom which in formula
read_formula(Input, Formula) :-
	string_concat(Input, '.txt', Filename1),
	string_to_atom(Filename1, Filename),
	read_file(Filename, Formula).


read_file(Filename, []) :- open(Filename, read, Str), at_end_of_stream(Str), write("there is no formula\n"), !.
read_file(Filename, Formula) :-
	open(Filename, read, Str),
	read(Str, Formula).
	%read_prop(Str).


read_prop(Str) :-
    read_prop_list(Str, X),
    close(Str).

read_prop_list(Stream, []) :-
	at_end_of_stream(Stream).


read_prop_list(Stream, [X|L]) :-
    \+ at_end_of_stream(Stream),
    read(Stream, X),
    assert(X),
    read_prop_list(Stream, L).

%--------------end read formula---------------------


pick_props(Count, From, To) :-
    pick_props(Count, From, [], To).

pick_props(0, _, Acc, Acc).
pick_props(_, From, From, []). 
pick_props(Count, From, Acc, To) :-
    Count > 0, % avoid looping on backtracking
	length(From, Len),
	(Count > Len, Count1 = Len; Count1 = Count),
    random_member(X, From),
    (  memberchk(X, Acc)
    -> pick_props(Count1, From, Acc, To)
    ;  C1 is Count1-1,
       pick_props(C1, From, [X|Acc], To)
    ) .
	
%-------------------end-------------------


new_ctlForget(F, V, Res, Ack, R, Exist,LenSNF, LenRes) :- (gain_prop(F, P), ctl2NNF(F, F1), L=[start-> z, z -> F1], assert(prop(z)), 
	tran2SNF(L, 0, 0, V1, IndM, SNF), tran2SNFClt(SNF, LC), append(V, V1, V2), append([z], V2, V3), 
 	length(LC, LenSNF), 
	appearing_list(LC, V3, L1), split2ST(L1, SL, TL, W, 0), %trace, 
	step_resolution(SL, V3, L2),
	append(V3, W, V6),
	(L2=false, R = false, Res = false, LenRes=0, Ack = false, Exist = false; length(P, LenVarF), length(V, LenV), 
		append(V3, P, V4), sort(V4, V5), temp_resolution(TL, W, L2, V3, R1, TL1, W1, V5), 
		append(L2, LC, L3), append(L3, R1, R2), sort(R2, R3), split2ST(R3, SL2, TL2, W2, 0), step_resolution(SL2, V3, R4), %trace, 
		(R4 = false, R = false, Res = false, LenRes=0, Ack = false, Exist = false;
		  (LenV = LenVarF, R = true, Res = true, Ack = true, Exist = true; append(R4, TL2, Ry), subtract(Ry, [snf_clause(true,[true], [true], nil)], R44),
		 instante(R44, V, V6, [], R5), 
		 subtract(V6, R5, R6),  
		 append(V, R6, V7), removeAtom(R44,V,R55),  
		 /*findall(Num, between(1, IndM, Num), Ind),
		 list_ind(R55, Lind, LNind),
		 subtract(V6, V, Vv), 
		 negXinEt(R55, Vv, AcV), 
		 pro6(Lind, Ind, P6),
		 append(P6, LNind, Res), length(Res, LenRes),*/
		 append([z], V1, AcV1),
		 append(W, AcV1, AcV),
		 write("AcV=\n"), write(AcV),
		 ackerM(R55, AcV, Ack),
		 /* write("\n =======Pro6 ========\n"), write(Res),
		  write("\n =======Pro6 ========\n"),*/
		  write("\n =======ackerM ========\n"), write(Ack),
		  write("\n =======ackerM ========\n"),
		  retractall(pair(X,Y)), retractall(formula(X1,Y1,Z1)), retractall(prop(Z)),  
		 retractall(alpha(_,_,_)), retractall(beta(_,_,_)), retractall(snfclause(_, _, _, _)),
		  R = true, length(Ack, Len), length(AcV, LenVv), (Len = LenVv, Exist = 1; Exist =0)
		 )
		)
	); 
	%(L2=false, R = false; append(L2, LC, L3), sort(L3, R));
	write("This formula is error! Please input the right formula.")).

list_ind([], [],[]).
list_ind([C|L], Lind, LNind) :- C = snf_clause(Type, H, T, Ind), 
	((Type = ef; Type = ex), C1 = C, C2 = []; C2 = formula_ind(Type, H, [T]), C1 = []),
	list_ind(L, Lind1, LNind1),
	(C1 = [], append([C2], LNind1, LNind), Lind = Lind1; append([C1], Lind1, Lind), LNind = LNind1).

%dd
negXinEt(_, [], []).
negXinEt(Lind, [X|H], AcV) :- negXinEt_list(Lind, X, R),
	%(R = true, R1 =[]; R1 = X),
	negXinEt(Lind, H, AcV1),
	(R = true, AcV = AcV1; append([X], AcV1, AcV)).

negXinEt_list([], X, false).
negXinEt_list([C|L], X, R) :- C = snf_clause(Type, H, T, Ind), 
	(Type = true, negXinEt_list(L, X, R); 
	(member(-X, T), R = true; negXinEt_list(L, X, R))).
	
	


 
toCTL([], _, _, []).
toCTL(L, V, V1, F) :- obtainAB(AB), end2first(L, V, V1, [], CTLcon),
	(AB=[], toCTL_list(L, CTLcon, F1), append(F1, CTLcon, FF), toCTL_temprol(L, F1, F2), append(F2, F1, F); 
		replaceAB(L, AB, R1), toCTL_list(R1, CTLcon, F1), append(F1, CTLcon, FF), toCTL_temprol(L, F1, F2), append(F2, F1, F)).
		
 
	

 
end2first(T, V, V1, F, CTLcon) :-  
	(V1 =[], !, CTL=[], CTLcon=[];
		%findall([Q,L], (member(formula(con, Q, L), T), isend(L, V, V1)), EL), 
		findEnd(T, V, V1, EL), 
		(EL= [], !, CTL=[], CTL1=[]; 
			con2ctl(EL, F, V2, CTL), subtract(V1, V2, V11), append(F, CTL, F1), end2first(T, V, V11, F1, CTL1)),
		(CTL=[], CTLcon=[]; append(CTL, CTL1, CTL2), sort(CTL2, CTLcon))
	).

findEnd([], _, _, []).
findEnd([H1|L1], V, V1, EL) :- H1 = formula(Type, Q, L),
	(Type=con, 
		(isend(L, V, V1), H2=[Q, L]; H2 = []); H2=[]),
	findEnd(L1, V, V1, EL1),
	(H2=[], EL=EL1; append([H2], EL1, EL)).
		 

con2ctl([], _, [], []).
con2ctl([H1|L1], F, [Q|LL], [H2|L2]) :- H1=[Q, L], conCTL(L, F, F1), H2 = ctlformula(Q, F1), con2ctl(L1, F, LL, L2) .


 
conCTL([],_, true).
conCTL([H|L], CTL, F) :- H = (Type, T, I), disCTL(T, CTL, T1),
	(Type=ex, T2=(^(*(T1)));
		(Type=ax, T2=(~(*(T1)));
			(Type=ef, T2=(^(?(T1)));
				(Type=af, T2=(~(?(T1)));
					T2 =T1
				)
			)
		)
	),
	conCTL(L, CTL, F1), (F1=true, F=T2; F = (T2 & F1)).

 	
disCTL([], _, false).
disCTL([H|L], CTL, F) :- equall(- H, X), disCTL(L, CTL, F1),
	((member(ctlformula(H, F2), CTL), F3=F2; member(ctlformula(X, F2), CTL), F3=(- F2)), (F1=false, F=F3; F = (F3 \/ F1)); %/
		(F1=false, F=H; F = (H \/ F1)) %/
	).
/*disCTL([H|L], CTL, F) :- equall(- H, X), 
	((member(ctlformula(H, F2), CTL); member(ctlformula(X, F2), CTL)), disCTL(L, CTL, F1), (F1=false, F=F2; F = (F2 \/ F1)); %/
		disCTL(L, CTL, F1), (F1=false, F=H; F = (H \/ F1)) 
	).*/



 
isend([], _, _).
isend([H|T], V, V1) :- H=(Type, L, I), obtainAtom(L, L1), 
	(intersection(L1, V, []), intersection(L1, V1, []), !; !, fail).


 
toCTL_temprol([], _, []).
toCTL_temprol([H|T], Con, R) :- H = formula(Type, Q, L), 
	((Type = con; (Type = alpha; Type = beta)), H1=[];
		(Type = af, L=[(X, nil)], disjunction(X, X1), H2=(~(? X1)), H1=ctlformula(Q, H2);
			(Type = ef, L=[(X, I)], disjunction(X, X1), H2=(^(? X1)), H1=ctlformula(Q, H2);
				(Type = ax, L=[(X, nil)], disjunction(X, X1), H2=(~(* X1)), H1=ctlformula(Q, H2);
					(Type = ex, L=[([X], I)], member(ctlformula(X, X1), Con), H2=(^(* X1)), H1=ctlformula(Q, H2);
						(Type = au, L=[([X], [Y], I)], member(ctlformula(X, X1),Con), 			
							member(ctlformula(Y, Y1), Con), H2 = (~(X1 $ Y1)), H1=ctlformula(Q, H2);
							(Type = eu, L=[([X], [Y], I)], member(ctlformula(X, X1),Con), 			
								member(ctlformula(Y, Y1), Con), H2 = (^(X1 $ Y1)), H1=ctlformula(Q, H2);
								(Type=ag, L=[([X], I)], member(ctlformula(X, X1), Con), H2=(~(@ X1)), H1 = ctlformula(Q, H2);
									(Type=eg, L=[([X], I)], member(ctlformula(X, X1), Con), H2=(^(@ X1)), H1 = ctlformula(Q, H2);
										H1=[]
									)
								)
							)					
						)
					)
				)
			)
		)
	), toCTL_temprol(T, Con, R1),
	(H1=[], R=R1; append([H1], R1, R)).



toCTL_list([],_,[]).
toCTL_list([H|L], CTLcon, R) :-  H = formula(Type, X, LL),  
	(Type=con, conCTL(LL, CTLcon, F1), H1=ctlformula(X, F1);
		(Type=alpha, conA(LL, CTLcon, F1), H1=ctlformula(X, F1);
			(Type=beta, conB(LL, CTLcon, F1), H1=ctlformula(X, F1); H1=[])
		)
	), 
	toCTL_list(L, CTLcon, R1),
	(H1=[], R=R1; append([H1], R1, R)).
	
 
conA([], _, true).
conA([H1, H2|L], CTLcon, F) :- H1 = [[C3,C2],I,[C33,C2,C4],[C2,C3]], append(C3, C2, C32), 
	append(C33, C2, C332), append(C332, C4, C3324), append(C3324, H2, CH2),
	disCTL(C32, CTLcon, F1), disCTL(CH2, CTLcon, F2), F3=(F1 \/ (^(*(F2 \/ ~(*(~(? F1))))))), %/
	conCTL(L, CTLcon, F4), F = (F3 & F4).
	
 
conB([], _, true).
conB([H1, H2|L], CTLcon, F) :- H1 = [[C3,C2],I,[C33,C2,C4],[C2,C3]], append(C3, C2, C32), 
	append(C33, C2, C332), append(C332, C4, C3324), append(C3324, H2, CH2),
	disCTL(C32, CTLcon, F1), disCTL(CH2, CTLcon, F2), %trace, simpDis(F1, F11), simpDis(F2, F22),
	F3=(F1 \/ (~(*(F2 \/ ~(*(~(? F1))))))), %/
	conCTL(L, CTLcon, F4), F = (F3 & F4).
	
	

conjunction([], true).
conjunction([H|L], F) :- H = (Type, T, I), disjunction(T, T1),
	(Type=ex, T2=(^(*(T1)));
		(Type=ax, T2=(~(*(T1)));
			(Type=ef, T2=(^(?(T1)));
				(Type=af, T2=(~(?(T1)));
					T2 =T1
				)
			)
		)
	),
	conjunction(L, F1), (F1=true, F=T2; F = (T2 & F1)).
	
disjunction([], false).
disjunction([H|L], F) :-  disjunction(L, F1), (F1 = false, F = H; F = (H \/ F1)). %/



 
obtainAB(R) :- findall(X, efImplication:alpha(_, X,_), R1), findall(Y, efImplication:beta(_, Y, _), R2), append(R1, R2, R).

 
replaceAB([], _, []).
replaceAB([H1|L1], AB, [H2|R1]) :- H1 = formula(Type, X, L),
	(Type=con, (member(X, AB), efImplication:snfclause(Type1, [X], T, I), 
		(Type1 = ex, 
			(member([N, alpha], T), efImplication:alpha(N, X, T2),  %T2=[
				subtract(T, [[N, alpha]], T3), append([T2], [T3], T33), 
				subtract(L, [(ex, T3, I)], LX), append(T33, LX, LL),
				H2 = formula(alpha, X, LL));
			(Type1=ax, (member([N, beta], T), efImplication:beta(N, X, T2),
				subtract(T, [[N, beta]], T3), append([T2], [T3], T33), 
				subtract(L, [(ax, T3, I)], LX), append(T33, LX, LL),
				H2 = formula(beta, X, LL)); H2=H1)
		); H2=H1); H2=H1),
	replaceAB(L1, AB, R1).
	%(H2=[], R=R1; append([H2],  R1, R)).
			






%##############
%ackerM

ackerM([], _, []). %将满足Ackermann的变量X及其对应的子句找出来，返回pair(X, List)对的集合。
ackerM(_, [], []).
ackerM(L, [X|List], RL) :- decideX(L, X, R1),
	(R1 = true, RL1 = []; allX(L, X, RL1)),
	ackerM(L, List, RL2),
	(RL1 =[], RL = RL2; append([pair(X, RL1)], RL2, RL)).

/*allX([], _, true). %找到列表中只包含一个X的公式，返回这类公式的集合
allX([H|L], X, RL) :- formulaX(H, X, R1), 
	allX(L, X, RL1),
	%（1）(R1 = true, H = formula_ind(Type, H1, T), T = [T1], subtract(T1, [-X], T2), append([T2], RL1, RL); RL = RL1).
	%（2） (R1 =true, H = formula_ind(Type, H1, T), T = [T1], subtract(T1, [-X], T2), cnf2wff([T2], F), 
	%（2） (RL1 = true, RL = F; RL = (RL1 & F)); RL = RL1).
	(R1 =true, 
		(H = formula_ind(Type, H1, T), 
			(Type = true, T = [T1], subtract(T1, [-X], T2), cnf2wff([T2], F); cnf2wff(T, F1),
				(Type = ef, F = (^(? F1));
					(Type = ex, F = (^(* F1));
						(Type = af, F = (~(? F1));
							(Type = ax, F = (~(* F1)); F = false
							)
						)
					)
				)
			)
		; H = formula_EF(H1, Dis1, Dis2), F = (Dis1 \/ (^(*(^(? Dis2))))) 
		),
	(RL1 = true, RL = F; RL = (RL1 & F)); RL = RL1).*/
	
allX([], _, true). %找到列表中只包含一个X的公式，返回这类公式的集合
allX([H|L], X, RL) :- formulaX(H, X, R1), 
	allX(L, X, RL1),
	(R1 =true, 
		H = snf_clause(Type, H1, T, Ind), 
			(Type = true, subtract(T, [-X], T2), cnf2wff([T2], F); cnf2wff([T], F1),
				(Type = ef, F = (^(? F1));
					(Type = ex, F = (^(* F1));
						(Type = af, F = (~(? F1));
							(Type = ax, F = (~(* F1)); F = false
							)
						)
					)
				)
			),
	(RL1 = true, RL = F; RL = (RL1 & F)); RL = RL1).
	
/*
allX([], _, []). %找到列表中只包含一个X的公式，返回这类公式的集合
allX([H|L], X, RL) :- formulaX(H, X, R1),
	allX(L, X, RL1),
	(R1 =true, append([H], RL1, RL); RL = RL1).*/


decideX([], X, false).%判断列表L中是否有公式的头尾都有X，若是返回true
decideX(L, X, R) :- findall(C, (member(C, L), formula2X(C, X, R1), R1 = true), R2), 
	(R2 = [], R = false; R = true).
	

/*formulaX(C, X, R) :- (C = formula_ind(Type, H, T), (Type = true, appearX_DListNeg(T,X, R2), 
	(R2 = true, R = true; R = false); (H = [X], R = true; R = false));
	C = formula_EF(H, Dis1, Dis2), (H = [X], R = true; R = false)).*/  %判断若C为true子句，且X负在其尾中；或者子句的头部只包含X，若是返回true
	
formulaX(C, X, R) :- C = snf_clause(Type, H, T, Ind), (Type = true, appearX_ListNeg(T,X, R2), 
	(R2 = true, R = true; R = false); (H = [X], R = true; R = false)).  %判断若C为true子句，且X负在其尾中；或者子句的头部只包含X，若是返回true
	
	
/*
formulaX(C, X, R) :- (C = formula_ind(Type, H, T), Type = true, appearX_DListNeg(T,X, R2), 
	(R2 = true, R = true; R = false);
	R = false). */ %判断若C为true子句，且X负在其尾中，若是返回true

/*	
formulaX(C, X, R) :- (C = formula_ind(Type, H, T), appearX_List(H, X,R1), appearX_DList(T,X, R2), 
	((R1 = true, R2 = false; R1 = false, R2 = true), R = true; R = false);
	C = formula_EF(H, Dis1, Dis2), appearX_List(H, X, R3), appearX(Dis1, X, R4), appearX(Dis2, X, R5), 
	((R3 = true, R4 = false, R5 = false; R3 = false, (R4 = true; R5 = true)), R = true; R = false)).*/  %判断C中的只有头或者只有尾中有X，若是返回true
	
/*formula2X(C, X, R) :- (C = formula_ind(Type, H, T), appearX_List(H, X,R1), appearX_DList(T,X, R2), (R1 = true, R2 = true, R = true; R = false);
	C = formula_EF(H, Dis1, Dis2), appearX_List(H, X, R3), appearX(Dis1, X, R4), appearX(Dis2, X, R5), 
	(R3 = true, (R4 = true; R5 = true), R = true; R = false)).  %判断C中的头尾是否都有X，若是返回true*/
	
formula2X(C, X, R) :- C = snf_clause(Type, H, T, Ind), appearX_List(H, X,R1), appearX_List(T,X, R2), (R1 = true, R2 = true, R = true; R = false).  %判断C中的头尾是否都有X，若是返回true


appearX([], _, false). %判断列表中的List或者CNF中是否有X，若是返回true
appearX([H|L], X, R) :- (appearX_List(H, X, R1), (R1 =false, appearX_DList(H, X, R2); R2 = false); appearX_DList(H, X, R2), R1 =false),
	((R1 = true; R2 = true), R = true; appearX(L, X, R)). 
/*
appearX([H|L], X, R) :- (appearX_List(H, X, R1), R2 = false; appearX_DList(H, X, R2), R1= false),
	((R1 = true; R2 = true), R = true; appearX(L, X, R)). */

appearX_List([], _, false).%判断List中是否有X(正负)，若是返回true
appearX_List([H|L], X, R) :- ((H = X ; H = - X), R = true; appearX_List(L, X, R)).

appearX_DList([], _, false).%判断CNF中是否有X(正负)，若是返回true
appearX_DList([H|L], X, R) :- appearX_List(H, X, R1),
	(R1 = true, R = true; appearX_DList(L, X, R)).


appearX_ListNeg([], _, false).%判断List中是否有-X，若是返回true
appearX_ListNeg([H|L], X, R) :- (H = - X, R = true; appearX_ListNeg(L, X, R)).
	
appearX_DListNeg([], _, false). %判断CNF中是否有-X，若是返回true
appearX_DListNeg([H|L], X, R) :- appearX_ListNeg(H, X, R1), 
	(R1 = true, R = true; appearX_DListNeg(L, X, R)).
	
	
	

%#################
%efImplication
  
efImp_Esome(H, L, NIV, IV, N, N1) :-  H = snf_clause(ef, Q, T, I), T = [X], equall(- X, Y), 
	findall([Q, C2, C3, C4, I1], (member(snf_clause(true, [true], T1,nil), L), member(Y, T1), subtract(T1, [Y], L1), 
	obtainAtom(L1, L2), intersection(L2, NIV, L3), not(L3 =[]), pair(X, T2),
	(member(snf_clause(ex, Q, T3,I),L), I1=I; member(snf_clause(ax, Q, T3, nil),L), I1=nil), 
	%(member(snf_clause(ex, Q, T3,I1),L); member(snf_clause(ax, Q, T3, nil),L), I1=nil), 
	(member(X1, T3), obtainAtom([X1], Lx), Lx =[XX],member(XX, L3), equall(- X1, X2), member(X2,L1)),  
	subtract(L1, [X2],C2), subtract(T2, C2, C3), subtract(T3, [X1], C4)), F),
	(F=[], !, N1 = N; addGamma(F, N), length(F, N2), N1 is N+N2).
	
	
efImp_Esome(H, L, NIV, IV,N, N1) :-  H = snf_clause(ef, Q, T, I1), 
	 conditionEf(H, L, L, NIV,IV, F), %trace, write("\n $$$F====="), write(F),
	delete(F, [], F1), 
	(F1=[], !, N1=N; addaf(Q, F1, N), length(F1, N2), N1 is N+N2).
	

conditionEf(_, [], _, _, _, []).
conditionEf(H, [H1|L11], L, NIV, IV, R) :- H = snf_clause(af, Q, T, I), T = [X], equall(- X, Y), 
	((H1=snf_clause(ex, Q, T3,I),L), I1=I; H1=snf_clause(ax, Q, T3, nil), I1=nil), 
		(member(snf_clause(true, [true], T1,nil), L), member(Y, T1), subtract(T1, [Y], L1), 
			obtainAtom(L1, L2), intersection(L2, NIV, L3), not(L3 =[]), pair(X, T2),
	%(member(snf_clause(ex, Q, T3,I),L); member(snf_clause(ax, Q, T3, I),L)), 
			(member(X1, T3), obtainAtom([X1], Lx), Lx =[XX], member(XX, L3), equall(- X1, X2), member(X2,L1)),  
			subtract(L1, [X2],C2), subtract(T2, C2, C3), subtract(T3, [X1], C4), 
			append([Y], C2, Ly), append(Ly, C4, Yy),
			(I1=nil, H2 = snf_clause(ax, Q, Yy, I1); H2 = snf_clause(ex, Q, Yy, I1)),
			C1=[Q, C2, C3, C4, I1,H2]; C1=[]);
		C1=[]
	), conditionEf(H, L11, L, NIV, IV, R1),
	(C1=[], R=R1; append([C1], R1, R2), sort(R2, R)).
	

addGamma([]).
addGamma([H|T], N) :- H=[Q, C2, C3, C4,I,H1], Q =[X], (C2=[]; C2 = [(_, T2, _)]), (C3=[]; C3=[(_, T3, _)]), (C4=[];C4 = (_, T4, _)),
	negation(T2, C22), negation(T3, C33), negation(T4, C44),
	%assert(gamma(N, X, [[C22, C33], I, [T3, C22, C44],[T2, T3]])),
	assert(gamma(N, X, [[T3, T2], I, [C33, T2, T4],[T2, T3]])),
	H1 = snf_clause(ax, Q, T, I), retractall(H1), append([[N, beta] ], T, T5), assert(snfclause(ax, Q, T5, I))
	N1 is N+1,
	addGamma(T, N1).
	

afImp_Asome(H, L, NIV, IV,N, N1) :-  H = snf_clause(af, Q, T, I1), 
	/*T = [X], equall(- X, Y), 
	findall([C2, C3, C4, I], (member(snf_clause(true, [true], T1,nil), L), member(Y, T1), subtract(T1, [Y], L1), 
	obtainAtom(L1, L2), intersection(L2, NIV, L3), not(L3 =[]), pair(X, T2),
	(member(snf_clause(ex, Q, T3,I),L); member(snf_clause(ax, Q, T3, I),L)), 
	(member(X1, T3), obtainAtom([X1], Lx), Lx =[XX], member(XX, L3), equall(- X1, X2), member(X2,L1)),  
	subtract(L1, [X2],C2), subtract(T2, C2, C3), subtract(T3, [X1], C4)), F),*/
	 conditionAf(H, L, L, NIV,IV, F), %trace, write("\n $$$F====="), write(F),
	delete(F, [], F1), 
	(F1=[], !, N1=N; addaf(Q, F1, N), length(F1, N2), N1 is N+N2).
	


	

conditionAf(_, [], _, _, _, []).
conditionAf(H, [H1|L11], L, NIV, IV, R) :- H = snf_clause(af, Q, T, I1), T = [X], equall(- X, Y), 
	((H1=snf_clause(ex, Q, T3,I); H1 = snf_clause(ax, Q, T3, I)),
		(member(snf_clause(true, [true], T1,nil), L), member(Y, T1), subtract(T1, [Y], L1), 
			obtainAtom(L1, L2), intersection(L2, NIV, L3), not(L3 =[]), pair(X, T2),
	%(member(snf_clause(ex, Q, T3,I),L); member(snf_clause(ax, Q, T3, I),L)), 
			(member(X1, T3), obtainAtom([X1], Lx), Lx =[XX], member(XX, L3), equall(- X1, X2), member(X2,L1)),  
			subtract(L1, [X2],C2), subtract(T2, C2, C3), subtract(T3, [X1], C4), 
			append([Y], C2, Ly), append(Ly, C4, Yy),
			(I=nil, H2 = snf_clause(ax, Q, Yy, I); H2 = snf_clause(ex, Q, Yy, I)),
			C1=[C2, C3, C4, I, H2]; C1=[]
		);
		C1=[]
	), conditionAf(H, L11, L, NIV, IV, R1),
	(C1=[], R=R1; append([C1], R1, R2), sort(R2, R)).
	



addaf(Q, [],_).
addaf(Q, [H|T], N) :- H=[C2, C3, C4,I,H1],  
	(I=nil, addBeta(Q, H,N); addAlpha(Q, H,N)), N1 is N+1, addaf(Q, T,N1).
	
addBeta(Q, H,N) :- H=[C2, C3, C4,I, H1], Q =[X], (C2=[]; C2 = [(_, T2, _)]), (C3=[]; C3=[(_, T3, _)]), (C4=[];C4 = (_, T4, _)),
	negation(T2, C22), negation(T3, C33), negation(T4, C44),
	%assert(beta(N, X, [[C22, C33], I, [T3, C22, C44],[T2, T3]])),
	assert(beta(N, X, [[T3, T2], I, [C33, T2, T4],[T2, T3]])),
	%(H1=snf_clause(ex, Q, T3,I), retractall(H1), append([beta], T3, T5), assert(snfclause(ex, Q, T5, I)); 
	H1 = snf_clause(ax, Q, T, I), retractall(H1), append([[N, beta] ], T, T5), assert(snfclause(ax, Q, T5, I)).
	
addAlpha(Q, H, N) :- H=[C2, C3, C4,I,H1], Q =[X], (C2=[]; C2 = [(_, T2, _)]), (C3=[]; C3=[(_, T3, _)]), (C4=[];C4 = (_, T4, _)),
	negation(T2, C22), negation(T3, C33), negation(T4, C44),
	%assert(alpha(N, X, [[C22, C33], I, [T3, C22, C44],[T2, T3]])),
	assert(alpha(N, X, [[T3, T2], I, [C33, T2, T4],[T2, T3]])),
	H1=snf_clause(ex, Q, T,I), retractall(H1), append([[N, alpha]], T, T5), assert(snfclause(ex, Q, T5, I)).
	


 
fImp([], _, _,_,N).
fImp([H|TL], L, NIV, IV, N) :- H = snf_clause(Type, Q, T, I), obtainAtom(T, R1), intersection(R1, IV, R2),
(R2=[], !;  
	(Type = ef, efImp_Esome(H, L, NIV, IV,N,N1); afImp_Asome(H, L, NIV, IV,N,N1)),
	N2 is N1+1,
	fImp(TL, L, NIV,IV, N2)
).




%#######################
%instantiate

insAllL([], _, _, []).  
insAllL([H|L], V, V1, V2) :- H = snf_clause(Type, Q, T, I),  
	(Type = true, true_ins([H], V, V1, V21, R1);
		temp_ins([H], V, V1, V21, R1)
	), insAllL(L, V, V1, R2),
	append(V21, R2, R3), sort(R3, V2).

%pair(p, L), append([c], L, L1), assert(pair(p, L1)), retractall(pair(p, L)).

 
instante([], _, _, _,[]). 
instante(L, V, V1, V2, R) :- obtainAllTrue(L, L1), true_ins(L1, V, V1, V2, R1), subtract(L, L1, L2),
	subtract(V1, R1, V3),
	(V3=[], R = V1;
		temp_ins(L2, V, V3, R1, R2), subtract(V3, R2, V4), 
		(V4 =[], R = V1;
			append(R1, R2, V5), tran_ins(L2, L, V, V4, V5, R3), subtract(V4, R3, V6),
			(V6=[], R=V1;
				append(V5, R3, V7), 
				(V7=V2, R = V7;
					instante(L, V, V6, V7, R4), append(V7, R4, R)
				)
			)
		)
	).
	

tran_ins([],_, _, _,_, []).  
tran_ins([C|L1], L, V, V1, V2, R) :- C = snf_clause(Type, H, T, I), subtract(V1, V2, V3), 
	(subset(H, V1),  append(V, V3, V4), obtainAtom(T, H1), 
		(intersection(H1, V4, []), !, obtainTandStart(L, V1, L2), findProp(L2, L2, H,X),
			(X=[], C1=[];
				(pair(X, Y), append([(Type,T,I)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y))), C1= []; 
					(T=[], !, C1= []; assert(pair(X, [(Type,T,I)])), C1= X)
				)
			); C1 =[]
		); C1 =[]
	),
	(C1 =[], V5=V2; append([C1], V2, V5)), %subtract(V1, [C1], V2),
	tran_ins(L1, L, V, V1, V5, R1),
	(C1=[], R = R1; append([C1], R1, R)).


findProp([], _,_, []).
findProp([C|L], L1, H, X) :- (C =snf_clause(true, [Y], T1, nil),
	findall(T, (C1=snf_clause(true, [Y], [T], nil), member(C1, L1)), L2); 
		%(subset(H, L2), X = Y; findProp(L, L1, H, X));
		C =snf_clause(start, [Y], T1, nil), findall(T, (C1=snf_clause(start, [Y], [T], nil), member(C1, L1)), L2)
		%(subset(H, L2), X = Y; 
	), (subset(H, L2), X = Y; findProp(L, L1, H, X)).


obtainTandStart([], _,[]).
obtainTandStart([C|L], V1, R) :- C = snf_clause(Type, H, T, I), 
	(Type = true, ((member(- X, T), member(X,V1), subtract(T, [- X], H1)), !, C1= snf_clause(true, [X], H1, nil); C1 =[]); 
		(Type = start, C1 = C; C1 = [])
	), 
	obtainTandStart(L, V1, L1),
	(C1=[], R = L1; append([C1], L1, R)).
	
	
	
%temp_ins(_,_,[],_,[]).
temp_ins([], _, _,_, []). 
temp_ins([C|L1], V, V1, V2, R) :- C = snf_clause(Type, H, T, I), subtract(V1, V2, V3),  
	( H = [X], !, 
		(member(X, V1), !, append(V, V3, V4), obtainAtom(T, H1), 
			(intersection(H1, V4, []), !, 
				(pair(X, Y), (T = [], !; append([(Type,T,I)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y)))), C1= []; 
					(T=[], !, C1= []; assert(pair(X, [(Type,T,I)])), C1= X)
				);
				C1=[]
			); C1=[]
		); C1=[]
	), 
	(C1 =[], V5=V2; append([C1], V2, V5)), %subtract(V1, [C1], V2),
	temp_ins(L1, V, V1, V5, R1),
	(C1=[], R = R1; append([C1], R1, R)).


%true_ins(_,_,[],_,[]).
true_ins([], _, _,_, []). 
/*true_ins([C|L1], V, V1, V2, R) :- C = snf_clause(_, H, T, nil), obtainAtom(T, H1), intersection(H1, V1, L3), subtract(L3, V2, L), length(L, N), 
(intersection(H1, V, L2), L2=[],
	(N >1, C1=[]; 
		(N =1, !, L = [X], 
			(member(- X, T), !, delete(T, - X, H2), 
				(pair(X, Y), (H2 = [], !; append([(true, H2, nil)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y)))), C1= []; 
					(H2=[], !, C1= []; assert(pair(X, [(true, H2, nil)])), C1= X)
				);
				C1=[]
			);
			C1= []
		)
	);
	C1 =[]
), 
	(C1 =[], V3=V2; append([C1], V2, V3)), %subtract(V1, [C1], V2),
	true_ins(L1, V, V1, V3, R1),
	(C1=[], R = R1; append([C1], R1, R)).*/
	
true_ins([C|L1], V, V1, V2, R) :- C = snf_clause(_, H, T, nil), obtainAtom(T, H1), intersection(H1, V1, L3), subtract(L3, V2, L), length(L, N), 
(L =[], findall(X,(member(- X, T), member(X, V1)), LL), true_list(LL, T), C1=[]; 
		/*(member(- X, T), !, delete(T, - X, H2),  
				(pair(X, Y), append([(true, H2, nil)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y))), C1= []; 
					assert(pair(X, [(true, H2, nil)])), C1= X
				);
				C1=[]
		);*/ 
	(intersection(H1, V, L2), L2=[],
		(N >1, C1=[]; 
			(N =1, !, L = [X], 
				(member(- X, T), !, delete(T, - X, H2),  
					(pair(X, Y), append([(true, H2, nil)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y))), C1= []; 
						assert(pair(X, [(true, H2, nil)])), C1= X
					);
					C1=[]
				);
				C1= []
			)
		);
		C1 =[]
	)
), 
	(C1 =[], V3=V2; append([C1], V2, V3)), %subtract(V1, [C1], V2),
	true_ins(L1, V, V1, V3, R1),
	(C1=[], R = R1; append([C1], R1, R)).
	
/*true_ins([C|L1], V, V1, V2, R) :- C = snf_clause(_, H, T, nil), obtainAtom(T, H1), intersection(H1, V1, L3), subtract(L3, V2, L), length(L, N), 
	(intersection(H1, V, L2), L2=[],
		(N >1, C1=[]; 
			(N =1, !, L = [X], 
				(member(- X, T), !, delete(T, - X, H2),  
					(pair(X, Y), append([(true, H2, nil)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y))), C1= []; 
						assert(pair(X, [(true, H2, nil)])), C1= X
					);
					C1=[]
				);
				C1= []
			)
		);
		C1 =[]
	), 
	(C1 =[], V3=V2; append([C1], V2, V3)), %subtract(V1, [C1], V2),
	true_ins(L1, V, V1, V3, R1),
	(C1=[], R = R1; append([C1], R1, R)).*/


 
true_list([],_). 
true_list([X|L], T) :- delete(T, - X, H2),  
	(pair(X, Y), append([(true, H2, nil)], Y, Y2), sort(Y2, Y1), assert(pair(X, Y1)), (Y1=Y,!; retractall(pair(X, Y))); 
		assert(pair(X, [(true, H2, nil)]))
	).


 
obtainAllTrue([], []). 
obtainAllTrue([C|L], R) :- C = snf_clause(Type, H, T, X), 
	(Type = true, C1 = C; C1 = []), 
	obtainAllTrue(L, L1),
	(C1=[], R = L1; append([C1], L1, R)).
	
	
	
	
	
%###################
%isCTLformula
%decide whether an experssion is a CTL formula.
is_CTLformula(P) :- atom(P).
is_CTLformula(-P) :- is_CTLformula(P).
is_CTLformula(P & Q) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(P \/ Q) :- is_CTLformula(P), is_CTLformula(Q).%/
is_CTLformula(P -> Q) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(^(@P)) :- is_CTLformula(P).
is_CTLformula(^(P $ Q)) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(^(*P)) :- is_CTLformula(P).
is_CTLformula(^(?P)) :- is_CTLformula(P).
is_CTLformula(~(@P)) :- is_CTLformula(P).
is_CTLformula(~(P $ Q)) :- is_CTLformula(P), is_CTLformula(Q).
is_CTLformula(~(*P)) :- is_CTLformula(P).
is_CTLformula(~(?P)) :- is_CTLformula(P).




%##############
%loopRes
unionC(C1, C2, P, X) :- equall(P, F1), delete(C1, F1, X1), equall(-P, F2), delete(C2, F2, X2), append(X1, X2, X3), sort(X3, X).
union(C1,C2,R) :- append(C1,C2,R1), sort(R1, R).

%SRES1
res_SRES1(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_, X1), C2=snf_clause(_,T2, H2,_,X2), (X1 > X2, N=X1; N=X2),
	union(T1,T2,T), unionC(H1,H2,P,H), R1=snf_clause(ax,T,H, nil, N), simplySNFCLoop(R1, R).
	
%SRES2
res_SRES2(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I,X1), C2=snf_clause(_,T2, H2,_,X2), (X1 > X2, N=X1; N=X2),
	union(T1,T2,T), unionC(H1,H2,P,H), R1=snf_clause(ex,T,H,I,N), simplySNFCLoop(R1, R).
	
%SRES3
res_SRES3(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I,X1), C2=snf_clause(_,T2, H2, I,X2), (X1 > X2, N=X1; N=X2),
	union(T1,T2,T), unionC(H1,H2,P,H), R1=snf_clause(ex, T,H,I, N), simplySNFCLoop(R1, R).
	
%SRES4
/*res_SRES4(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(start,T1,H,nil).*/
	
%SRES5
/*res_SRES5(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(start,T2,H,nil).*/
	
 %SRES6
res_SRES6(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_,_), C2=snf_clause(_,T2, H2,_,_),
	unionC(H1,H2,P,H), (H = [], R1=snf_clause(ax,T2,H,nil, 3), simplySNFCLoop(R1, R); R1=snf_clause(ax,T2,H,nil, 2), simplySNFCLoop(R1, R)).
	
%SRES7
res_SRES7(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_,_), C2=snf_clause(_,T2, H2, I,_),
	unionC(H1,H2,P,H), (H=[], R1=snf_clause(ex,T2,H,I,3), simplySNFCLoop(R1, R); R1=snf_clause(ex,T2,H,I,2), simplySNFCLoop(R1, R)).
	
%SRES8
res_SRES8(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_, X1), C2=snf_clause(_,T2, H2,_,X2), (X1 > X2, N=X1; N=X2),
	unionC(H1,H2,P,H), R1=snf_clause(true,T2,H,nil,N), simplySNFCLoop(R1, R).
	
	
 
	/*%SRES1
res_SRES1(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_, X1), C2=snf_clause(_,T2, H2,_,X2), (X1 > X2, N=X1; N=X2),
	union(T1,T2,T), unionC(H1,H2,P,H), R=snf_clause(ax,T,H, nil, N).
	
%SRES2
res_SRES2(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I,X1), C2=snf_clause(_,T2, H2,_,X2), (X1 > X2, N=X1; N=X2),
	union(T1,T2,T), unionC(H1,H2,P,H), R=snf_clause(ex,T,H,I,N).
	
%SRES3
res_SRES3(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I,X1), C2=snf_clause(_,T2, H2, I,X2), (X1 > X2, N=X1; N=X2),
	union(T1,T2,T), unionC(H1,H2,P,H), R=snf_clause(ex, T,H,I, N).*/
	
%SRES4
/*res_SRES4(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(start,T1,H,nil).*/
	
%SRES5
/*res_SRES5(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(start,T2,H,nil).*/
	
%SRES6
/*res_SRES6(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_,_), C2=snf_clause(_,T2, H2,_,_),
	unionC(H1,H2,P,H), (H = [], R=snf_clause(ax,T2,H,nil, 3); R=snf_clause(ax,T2,H,nil, 2)).
	
%SRES7
res_SRES7(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_,_), C2=snf_clause(_,T2, H2, I,_),
	unionC(H1,H2,P,H), (H=[], R=snf_clause(ex,T2,H,I,3); R=snf_clause(ex,T2,H,I,2)).
	
%SRES8
res_SRES8(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_, X1), C2=snf_clause(_,T2, H2,_,X2), (X1 > X2, N=X1; N=X2),
	unionC(H1,H2,P,H), R=snf_clause(true,T2,H,nil,N).*/

%resolution--------------------------------------------end-------------------------------------


%----------------------------------start----step_resolution-----------------------------------------

%resolution(Lp, Ln, P, L)----------------------------start-------------------------------------------
%do all the possible step-resolutions on P.  
loopresolution([], Ln, P, Ln).
loopresolution(L, [], _, L).
loopresolution([H|Lp], Ln, P, L) :- loopresolution_list(H, Ln, P, L1), loopresolution(Lp, Ln, P, L2), 
	append(L1, L2, R1), %trace, write("\n"), write("R***====="), write(R1), 
	append(R1, Ln, R), sort(R, L).

loopresolution_list(C, [], _, [C]).
loopresolution_list(C, [H|Ln], P, L1) :- C = snf_clause(Type1,_,_,_,_), H=snf_clause(Type2,_,_,_,_), %judgeType(C, Type1), judgeType(H, Type2), 
	(Type1 = ax, !, (Type2=ax, !, res_SRES1(C, H, P, C1); (Type2 =ex, !, res_SRES2(H, C, -P, C1); (Type2=true, !, res_SRES6(H, C, -P, C1); C1=[]))); 
		(Type1=ex, !, (Type2=ax, !, res_SRES2(C, H, P, C1); (Type2 =ex, !, res_SRES2(C, H, P, C1); (Type2=true, !, res_SRES7(H, C, -P, C1); C1=[])));
			(Type1 =true, !, (Type2=start, !, res_SRES5(C, H, P, C1); (Type2=ax, !, res_SRES6(C, H, P, C1); (Type2=ex, !, res_SRES7(C,H, P, C1); (Type2=true, !, res_SRES8(C,H,P, C1); C1=[]))));
				(Type1 =start, !, (Type2=start, !, res_SRES4(C,H, P, C1); (Type2=true, !, res_SRES5(H, C, -P, C1); C1=[])); C1 =[])
			)
		)
	),
	loopresolution_list(C, Ln, P, L2),
	(C1 =[], L1 =L2;
	append([C1], L2, L1)).
	
	

repeat_loopresolution([],_, []).
repeat_loopresolution(L, P, Res) :- sort(L, L1), all_PosC(L1, P, Lp), all_NegC(L1, P, Ln), loopresolution(Lp, Ln, P, R2), 
	rewrite(R2, R3), sort(R3, R1),
	%trace, write("\n"), write("R1====="), write(R1), 
	(L1 = R1, Res = R1;
			repeat_loopresolution(R1, P, Res)
		%, write("\n"), write("Res="), write(Res).
	).

	/*( member(snf_clause(true,[true],[- z], _,_), R1), Res=false;
		(L1 = R1, Res = R1;
			repeat_resolution(R1, P, Res)
		)%, write("\n"), write("Res="), write(Res).
	).*/


rewrite([], []).
rewrite([C1|R1], [C2| R2]) :- (C1=snf_clause(L, H, [], I,N), (L=ax; L=ex), negation(H, H1), C2 = snf_clause(true, [true], H1, nil,N); C2 = C1), rewrite(R1, R2).

%findall(X, (member(X, L1), inCst(X, Q)), L2)

 


loopRes(L, [], L).
loopRes([], _, []).
loopRes(L, V, R) :- appearing_list(L, V, L1), 
	step_loopresolution_list(L1, V, Rf),  
	( Rf=false, !, R=false;
		((Rf =L, R1 = Rf; loopRes(Rf, V, R1)),
			( R1 = false, !, R = false; append(L, R1, X1), sort(X1, R))
		)
	).
	/*(Rf =L, R1 = Rf; step_resolution(Rf, V, R1)),
	append(L, R1, X1), sort(X1, R).*/
	%trans(X, Y1), tranI2CTL_list(Y1, R), 
	%tranSt2CTL(Y2, Y3), 
	%elim(Y2, V, Y4), 
	%sort(Y2, R).
	%, snfL2CTLF(Y5, R).


step_loopresolution_list([], _, []).
step_loopresolution_list(L, [], L).
step_loopresolution_list(L, [P|T], R) :- appearing(L, P, L1), %write("\nl"), write("L1="), write(L1),
	repeat_loopresolution(L1, P, Res), 
	%subsumeLoop(R1,R1,R1,Res),
	( Res=false, !, R = false; 
		(member(snf_clause(start,[start],[],_,_), Res), !, R =false;
			(member(snf_clause(true,[true], [],_,_), Res), !, R=false;
				(
					append(L, Res, Res1),  %write("\nl"), write("Res1="), write(Res1),
					sort(Res1, Y), %write("\nl"), write("Y="), write(Y),
					step_loopresolution_list(Y, T, Res2), %write("\nl"), write("Res2="), write(Res2),
					(Res2=false, R=false; 
						append(Y, Res2, X), %write("\nl"), write("X="), write(X),
						sort(X, R) %,write("\nl"), write("R="), write(R)
					)
				)
			)
		)
	).
	/* append(L, Res, Res1),  write("\nl"), write("Res1="), write(Res1),
	sort(Res1, Y), write("\nl"), write("Y="), write(Y),
	step_resolution_list(Y, T, Res2), write("\nl"), write("Res2="), write(Res2),
	append(Y, Res2, X), write("\nl"), write("X="), write(X),
	sort(X, R), write("\nl"), write("R="), write(R).*/
	%trans(Res, Res1),
	%elm(Res1, R).
	
	
%---------------------------------end-----step_resolution-----------------------------------------	





%######################
%pro6
pro6([], _, []).
pro6(_, [], []).
pro6(L, [Ind|I], RL) :- formula_Ind(L, Ind, LP, F1), subtract(L, LP, L1), 
	find_EF(L1, Ind, F2),  
	(F2 = [], F = F1; list_pro6_vi(F1, F2, Fx), F2 = formula_EF(H, Dis1, Dis2), F3 = formula_EF(H, [Dis1], [Dis2]), append([F3], F1, Fy), append(Fy,Fx, F)),
	pro6(L1, I, RL1),
	append(F, RL1, RL2), sort(RL2, RL).

formula_Ind([], _, [], []).
formula_Ind(_, [], [], []).
formula_Ind(L, Ind, LP, F) :- list_EX(L, Ind, LP),find_head(LP, HL), diff_head(LP, HL, HF), pro6_ii(HF, F).  





list_EX([], _, []).  
list_EX(L, Ind, LP) :- findall(C, (C = snf_clause(ex, H, T, Ind1), member(C, L), Ind = Ind1), R), sort(R, LP).
 
	
find_head([], []).  
find_head(L1, L2) :- findall(Head, member(snf_clause(_, Head, _, _), L1), R), sort(R, L2).

diff_head(_, [], []).
diff_head([], _, []).  
diff_head(L, [H|L1], L2) :- findall(T, member(snf_clause(ex, H, T, _), L), R), sort(R,R1), %assert(formula_ind(H, F)),
	diff_head(L, L1, L3),
	append([formula_ind(ex, H, R1)], L3, L2).  
 
	
	

pro6_ii(L, F) :- list_allsubseqs(L, PowerL), result(PowerL, F).
pro6_vi(F1, F2, F) :- F1 = formula_ind(Type, H1, T1), F2 = formula_EF(H2, Dis1, Dis2),
	append(H1, H2, H3), sort(H3, H4),
	append([Dis1], [T1], T2), append([T1], [Dis2], T3),
	F = formula_EF(H4, T2, T3).
/*
pro6_vi(F1, F2, F) :- F1 = formula_ind(H1, T1), F2 = formula_EF(H2, Dis1, Dis2),
	append(H1, H2, H3), sort(H3, H4),
	F = formula_EF(H4, (Dis1 & (^(*(T1)))), (^(*(T1 & Dis2)))).*/
/*
pro6_vi(F1, F2, F) :- F1 = formula_ind(H1, T1), F2 = formula_EF(H2, Dis1, Dis2),
	append(H1, H2, H3), sort(H3, H4),
	F3 = formula_EF(H4, (Dis1 & (^(*(T1)))), (^(*(T1 & Dis2)))),
	F = ((F1 & F2) & F3).*/



list_pro6_vi([], _, []).
list_pro6_vi([H|L1], F2, R) :- pro6_vi(H,F2, F),
	list_pro6_vi(L1, F2, R1),
	append([F], R1, R).


result([], []).  
result([H|L1], F) :-  (H = [], F1 = [];
	result_list(H, F1)),
	result(L1, F2),
	(H = [], F = F2; append([F1], F2, F)).

/*
result([], true).  
result([H|L1], F) :-  (H = [], F1 = true;
	result_list(H, F1)),
	result(L1, F2),
	(H = [], F = F2; (F2 = true, F = F1; F = (F1 & F2))).*/
	
	
	
	
result_list([], true).
result_list([C|L1], F) :- C = formula_ind(Type, H, F1), 
	result_list(L1, F2),
	(F2 = true, F = C; F2 = formula_ind(Type, H1, F3), append(H, H1, H2), sort(H2, H3), append(F1, F3, F4), F = formula_ind(Type, H3, F4)).
	%(F2 = true, F = C; F2 = formula_ind(H1, F3), append(H, H1, H2), sort(H2, H3), F4 = (F1 & F3), F = formula_ind(H3, F4)).


find_EF([], _, []).
find_EF(L, Ind, F) :- list_EF(L, Ind, C), 
	(C = [], F =[]; pro5(C, F)).

list_EF([], _, []).
list_EF([H|L1], Ind, C) :- H = snf_clause(Type, Head, T, I),
	(Type = ef, I = Ind, C = H; list_EF(L1, Ind, C)).


pro5(C, C1) :- C = snf_clause(Type, H, T, I),
	C1 = formula_EF(H, T, T).
/*
pro5(C, C1) :- C = snf_clause(Type, H, T, I), cnf2wff([T], F),
	EF = (^(? (F))),
	C1 = formula_EF(H, F, EF).*/

list_allsubseqs(Es, Uss) :-
   list_acc_allsubseqs(Es, [[]], Uss).

lists_prependum([]      , _) --> [].
lists_prependum([Es|Ess], E) --> [[E|Es]], lists_prependum(Ess, E).

list_acc_allsubseqs([]    , Uss , Uss).
list_acc_allsubseqs([E|Es], Uss0, Uss) :-
   list_acc_allsubseqs(Es, Uss0, Uss1),
   phrase(lists_prependum(Uss1,E), Uss, Uss1).
   
   
   
 %##############
 %removeAtom
 
removeAtom([],_, []).
removeAtom([C|L], V, R):- C = snf_clause(Type, H, T, I), 
	(intersection(H, V, []),!, subtract(T, [true], T2), obtainAtom(T2, T1), (intersection(T1, V,[]), C1=C; C1=[]);
		C1=[]
	), removeAtom(L, V, R1),
	(C1=[], R=R1; append([C1], R1, R)).
	
	
	




%##################
%simplySNFC
simp_disSet([],[false]).
simp_disSet(L, L1) :- (member(P, L), member(-P, L), L1= [true]; 
	(member(true, L), L1=[true]; 
		(member(false, L), elimElement(false, L, L2), sort(L2, L1); sort(L, L1))
	)
	).  
	
simp_conSet([], [true]).
simp_conSet(L, L1) :- (member(P, L), member(-P, L), L1= [false]; 
	(member(false, L), L1 =[false]; 
		(member(true, L), elimElement(true, L, L2), sort(L2, L1); sort(L, L1))
	)
	).  


 
simplySNFC(C, C1) :- C= snf_clause(Type, H, T, L), simp_conSet(H, H1), simp_disSet(T, T1), C1=snf_clause(Type, H1, T1, L).



 
%subsume rules
subsume1(X, Y) :- X = snf_clause(_, H1, T1, _), Y = snf_clause(_, H2, T2, _), subset(T1, T2).  

subsume2(X, Y) :- X = snf_clause(_, H1, T1, _), Y = snf_clause(_, H2, T2, _), subset(T1, T2), subset(H1, H2).

subsume3(X, Y) :- X = snf_clause(_, H1, T1, _, _), Y = snf_clause(_, H2, T2, _,_), subset(T1, T2).


subsumeLoop([], _, L, L).
subsumeLoop(_, [], [], []).
subsumeLoop([H|L], X, Y, R) :- subsumeLoop_list(H, X, Y, L1), 
	(L1 = X, subsumeLoop(L, L1, L1, R); subsumeLoop(L1, L1, L1, R)).


subsumeLoop_list(C, [], L, L).
subsumeLoop_list(C, [H|T], L, R) :- (C = H, subsumeLoop_list(C, T, L, R); 
    ((subsume3(C, H), delete(L, H, L1); L1 = L),
	    subsumeLoop_list(C, T, L1, R)
    )).




subsume([], _, L, L).
subsume(_, [], [], []).
subsume([H|L], X, Y, R) :- subsume_list(H, X, Y, L1), 
	(L1=X, subsume(L, L1, L1, R); subsume(L1, L1, L1, R)).


subsume_list(C, [], L, L).
subsume_list(C, [H|T], L, R) :- (C = H, subsume_list(C, T, L, R); 
(C = snf_clause(Type1,_,_,_), H=snf_clause(Type2,_,_,_), %judgeType(C, Type1), judgeType(H, Type2), 
	(Type1 = ax, !, (Type2=ax, !, (subsume2(C, H), delete(L, H, L1); L1 = L); (Type2 =ex, !, (subsume2(C, H), delete(L, H, L1); L1 = L ); (Type2=true, !, (subsume1(H, C), delete(L, C, L1); L1 = L); L1=L))); 
		(Type1=ex, !, (Type2=ax, !,(subsume2(H, C), delete(L, C, L1); L1 = L); (Type2 =ex, !, (subsume2(C, H), delete(L, H, L1); L1 = L); (Type2=true, !, (subsume1(H, C), delete(L, C, L1); L1 = L); L1=L)));
			(Type1 =true, !, ((Type2=start; Type2=ax; Type2=ex; Type2=true), (subsume1(C, H), delete(L, H, L1); L1 = L));
				/*(Type2=start, !, (subsume1(C, H), delete(H, L, L1); L1 = L); (Type2=ax, !, (subsume1(C, H), delete(H, L, L1); L1 = L); (Type2=ex, !, (subsume1(C, H), delete(H, L, L1); L1 = L); (Type2=true, !, (subsume1(C, H), delete(H, L, L1); L1 = L); L1 = L))));*/
				(Type1 =start, !, (Type2=start, !, (subsume1(C, H), delete(L, H, L1); L1 = L); (Type2=true, !, (subsume1(H, C), delete(L, C, L1); L1 = L); L1 = L)); L1 = L)
			)
		)
	),
	subsume_list(C, T, L1, R)
)).





%##################
%stepRules


%resolution ----------------------------------------------------------start--------------------
%step rules


unionC(C1, C2, P, X) :- equall(P, F1), delete(C1, F1, X1), equall(-P, F2), delete(C2, F2, X2), append(X1, X2, X3), sort(X3, X).
union(C1,C2,R) :- append(C1,C2,R1), sort(R1, R).

%SRES1
res_SRES1(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	union(T1,T2,T), unionC(H1,H2,P,H), R1=snf_clause(ax,T,H, nil), simplySNFC(R1, R).
	
%SRES2
res_SRES2(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I), C2=snf_clause(_,T2, H2, _), 
	union(T1,T2,T), unionC(H1,H2,P,H), R1=snf_clause(ex,T,H,I), simplySNFC(R1, R).
	
%SRES3
res_SRES3(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I), C2=snf_clause(_,T2, H2, I),
	union(T1,T2,T), unionC(H1,H2,P,H), R1=snf_clause(ex, T,H,I), simplySNFC(R1, R).
	
%SRES4
res_SRES4(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R1=snf_clause(start,T1,H,nil), simplySNFC(R1, R).
	
%SRES5
res_SRES5(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R1=snf_clause(start,T2,H,nil), simplySNFC(R1, R).
	
%SRES6
res_SRES6(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R1=snf_clause(ax,T2,H,nil), simplySNFC(R1, R).
	
%SRES7
res_SRES7(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2, I),
	unionC(H1,H2,P,H), R1=snf_clause(ex,T2,H,I), simplySNFC(R1, R).
	
%SRES8
res_SRES8(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R1=snf_clause(true,T2,H,nil), simplySNFC(R1, R).
	

 
/*%SRES1
res_SRES1(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	union(T1,T2,T), unionC(H1,H2,P,H), R=snf_clause(ax,T,H, nil).
	
%SRES2
res_SRES2(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I), C2=snf_clause(_,T2, H2, _), 
	union(T1,T2,T), unionC(H1,H2,P,H), R=snf_clause(ex,T,H,I).
	
%SRES3
res_SRES3(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,I), C2=snf_clause(_,T2, H2, I),
	union(T1,T2,T), unionC(H1,H2,P,H), R=snf_clause(ex, T,H,I).
	
%SRES4
res_SRES4(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(start,T1,H,nil).
	
%SRES5
res_SRES5(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(start,T2,H,nil).
	
%SRES6
res_SRES6(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(ax,T2,H,nil).
	
%SRES7
res_SRES7(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2, I),
	unionC(H1,H2,P,H), R=snf_clause(ex,T2,H,I).
	
%SRES8
res_SRES8(C1,C2,P,R) :- C1=snf_clause(_,T1, H1,_), C2=snf_clause(_,T2, H2,_),
	unionC(H1,H2,P,H), R=snf_clause(true,T2,H,nil).*/

%resolution--------------------------------------------end-------------------------------------


%----------------------------------start----step_resolution-----------------------------------------

%resolution(Lp, Ln, P, L)----------------------------start-------------------------------------------
%do all the possible step-resolutions on P.  
resolution([], Ln, P, Ln).
resolution(L, [], _, L).
resolution([H|Lp], Ln, P, L) :- resolution_list(H, Ln, P, L1), resolution(Lp, Ln, P, L2),
	append(L1, L2, R1), %trace, write("\n"), write("R***====="), write(R1), 
	append(R1, Ln, R), sort(R, L).

resolution_list(C, [], _, [C]).
resolution_list(C, [H|Ln], P, L1) :- C = snf_clause(Type1,_,_,_), H=snf_clause(Type2,_,_,_), %judgeType(C, Type1), judgeType(H, Type2), 
	(Type1 = ax, !, (Type2=ax, !, res_SRES1(C, H, P, C1); (Type2 =ex, !, res_SRES2(H, C, -P, C1); (Type2=true, !, res_SRES6(H, C, -P, C1); C1=[]))); 
		(Type1=ex, !, (Type2=ax, !, res_SRES2(C, H, P, C1); (Type2 =ex, !, res_SRES3(C, H, P, C1); (Type2=true, !, res_SRES7(H, C, -P, C1); C1=[])));
			(Type1 =true, !, (Type2=start, !, res_SRES5(C, H, P, C1); (Type2=ax, !, res_SRES6(C, H, P, C1); (Type2=ex, !, res_SRES7(C,H, P, C1); (Type2=true, !, res_SRES8(C,H,P, C1); C1=[]))));
				(Type1 =start, !, (Type2=start, !, res_SRES4(C,H, P, C1); (Type2=true, !, res_SRES5(H, C, -P, C1); C1=[])); C1 =[])
			)
		)
	),
	resolution_list(C, Ln, P, L2),
	(C1 =[], L1 =L2;
	append([C1], L2, L1)).
	
	

repeat_resolution([],_, []).
repeat_resolution(L, P, Res) :- sort(L, L1), all_PosC(L1, P, Lp), all_NegC(L1, P, Ln), resolution(Lp, Ln, P, R2), 
	%trace, write("\n"), write("R2====="), write(R2), 
	rewrite(R2, R3), sort(R3, R1),
	%trace, write("\n"), write("R1====="), write(R1), 
	(L1 = R1, Res = R1;
			repeat_resolution(R1, P, Res)
		%, write("\n"), write("Res="), write(Res).
	).
	/*( member(snf_clause(true,[true],[- z], _), R1), Res=false;
		(L1 = R1, Res = R1;
			repeat_resolution(R1, P, Res)
		)%, write("\n"), write("Res="), write(Res).
	).*/


rewrite([], []).
rewrite([C1|R1], [C2| R2]) :- (C1=snf_clause(L, H, [], I), (L=ax; L=ex), negation(H, H1), C2 = snf_clause(true, [true], H1, nil); C2 = C1), rewrite(R1, R2).

%findall(X, (member(X, L1), inCst(X, Q)), L2)

 

step_resolution(L, [], L).
step_resolution([], _, []).
step_resolution(L, V, R) :- appearing_list(L, V, L1), 
	step_resolution_list(L1, V, Rf),  
	(Rf=false, !, R=false;
		((Rf =L, R1 = Rf; step_resolution(Rf, V, R1)),
			( R1 = false, !, R = false; append(L, R1, X1), sort(X1, R))
		)
	).
	/*(Rf =L, R1 = Rf; step_resolution(Rf, V, R1)),
	append(L, R1, X1), sort(X1, R).*/
	%trans(X, Y1), tranI2CTL_list(Y1, R), 
	%tranSt2CTL(Y2, Y3), 
	%elim(Y2, V, Y4), 
	%sort(Y2, R).
	%, snfL2CTLF(Y5, R).


step_resolution_list([], _, []).
step_resolution_list(L, [], L).
step_resolution_list(L, [P|T], R) :- appearing(L, P, L1),  %write("\nl"), write("L1="), write(L1),
	repeat_resolution(L1, P, R1), 
	%trace, write("\n"), write("R1====="), write(R1), 
	subsume(R1,R1,R1,Res),
	( Res=false, !, R = false; 
		(member(snf_clause(start,[start],[],_), Res), !, R =false;
			(member(snf_clause(true,[true], [],_), Res), !, R=false;
				(
					append(L, Res, Res1),  %write("\nl"), write("Res1="), write(Res1),
					sort(Res1, Y), %write("\nl"), write("Y="), write(Y),
					step_resolution_list(Y, T, Res2), %write("\nl"), write("Res2="), write(Res2),
					(Res2=false, R=false; 
						append(Y, Res2, X), %write("\nl"), write("X="), write(X),
						sort(X, R) %, write("\nl"), write("R="), write(R)
					)
				)
			)
		)
	).
	/* append(L, Res, Res1),  write("\nl"), write("Res1="), write(Res1),
	sort(Res1, Y), write("\nl"), write("Y="), write(Y),
	step_resolution_list(Y, T, Res2), write("\nl"), write("Res2="), write(Res2),
	append(Y, Res2, X), write("\nl"), write("X="), write(X),
	sort(X, R), write("\nl"), write("R="), write(R).*/
	%trans(Res, Res1),
	%elm(Res1, R).
	
	
%---------------------------------end-----step_resolution-----------------------------------------	




%############################
%tempRes



temp_resolution([], _,_,_,[],[],[], V).
temp_resolution(TL, _, [],_, [], TL, [], V).
temp_resolution(TL, _, _, [], [], TL, [], V).
temp_resolution([H|TL], [X|L], L1, [Y|V], R, TL1, W, V1) :- temp_resolution_list(H, L1, [Y|V], R1,X, V1), 
	temp_resolution(TL, L, L1, [Y|V], R2, TL2, W1, V1),
	(R1 = [], append([H], TL2, TL1), W = W1, R=R2; TL1= TL2,  append([X], W1, W), append(R1, R2, R)).
	
	
	
temp_resolution_list(_, _, [], [], _, V).
temp_resolution_list(H, L, [Y|V], R,W, V1):- H = snf_clause(Type, Hs, Ts, I), 
	(Type = af, (posC(Y, H), fA(H, L, Y, R,W, V1); (negC(Y, H),fA(H, L, - Y, R,W, V1); temp_resolution_list(H, L, V, R,W, V1)));
		(Type = ef, (posC(Y, H), fE(H, L, Y, R,W, V1); (negC(Y, H),fE(H, L, - Y, R,W, V1); temp_resolution_list(H, L, V, R,W, V1))))
	).
	
fA(C, L, Y, R,W, V) :- find_Asometime_clauses(L, L1), loopFormula(L1, Y, V, [[Y]], F),  %write("\n F===="), write(F),write("\n"),
	C = snf_clause(Type, Q, H, I),
	((F = false; F=[]), R = []; snfres(F, Y, Q, R, W), write("\n ERES1===="), write(R), write("\n")).
	
snfres(true, Y,Q, R, W) :- negation(Q, T1), append(T1, [Y], T2), append(T2, [W], T3), 
	R = [snf_clause(ax,[W],[Y],nil), snf_clause(true,[true], T2,nil), snf_clause(true,[true], T3,nil), snf_clause(ax,[W], [Y, W],nil)].
snfres(F, Y, Q, R, W) :-  equall(Y, P), negation(Q, T1), append(T1, [P], T2), append(T2, [W], T3), cirRes(W, P, T1, F, R1),
	L1= [snf_clause(true,[true], T3,nil), snf_clause(ax,[W], [P, W],nil)], append(R1, L1, R).

cirRes(_,_, _, [], []).
cirRes(W, P, Q, [F|T], L) :- negation(F, F1), append([P], F, X1), append(Q, X1, F2),  
	C1 = snf_clause(ax,[W], X1,nil), C2 = snf_clause(true,[true], F2, nil),
	cirRes(W,P,Q, T, L1),
	append([C1], L1, L2), append([C2], L2, L3), sort(L3, L).


fE(C, L, Y, R,W, V) :- %C= exist_future_clause(Q, H, I), find_Esometime_loop(L, Y, F, I), 
	find_Esometime_clauses(L, Ind, L1), loopFormula(L1, Y, V, [[Y]], F), %write("\n F===="), write(F),write("\n"),
	C = snf_clause(Type, Q, H, I), 
	((F = false; F=[]), R = []; snfresE(F, Y, Q, R, W,I), write("\n ERES2===="), write(R), write("\n")).
	
snfresE(true, Y,Q, R, W,I) :- equall(Y, P), negation(Q, T1), append(T1, [P], T2), append(T2, [W], T3), 
	R = [snf_clause(ex,[W],[Y],I), snf_clause(true,[true], T2,nil), snf_clause(true,[true], T3,nil), snf_clause(ex, [W], [P, W],I)].
snfresE(F, Y, Q, R, W,I) :-  equall(Y, P), negation(Q, T1), append(T1, [P], T2), append(T2, [W], T3), 
	cirResE(W, P, T1, F, R1,I),
	L1= [snf_clause(true,[true], T3, nil), snf_clause(ex,[W], [P, W],I)], append(R1, L1, R).

cirResE(_,_, _, [], [],I).
cirResE(W, P, Q, [F|T], L,I) :- negation(F, F1), append([P], F, X1), append(Q, X1, F2), 
	C1 = snf_clause(ex, [W], X1,I), C2 = snf_clause(true,[true], F2, nil),
	cirResE(W,P,Q, T, L1,I),
	append([C1], L1, L2), append([C2], L2, L3), sort(L3, L).
	
	
	
	
	
	
	
	
	
	
	
	
%###################
%toCTL
 
toCTL([], _, _, []).
toCTL(L, V, V1, F) :- obtainAB(AB), end2first(L, V, V1, [], CTLcon),
	(AB=[], toCTL_list(L, CTLcon, F1), append(F1, CTLcon, FF), toCTL_temprol(L, F1, F2), append(F2, F1, F); 
		replaceAB(L, AB, R1), toCTL_list(R1, CTLcon, F1), append(F1, CTLcon, FF), toCTL_temprol(L, F1, F2), append(F2, F1, F)).
		
/*toCTL(L, F) :- obtainAB(AB), 
	(AB=[], toCTL_list(L, F1), toCTL_temprol(L, F1, F2), append(F2, F1, F); 
	replaceAB(L, AB, R1), toCTL_list(R1, F1), toCTL_temprol(R1, F1, F2), append(F2, F1, F)).*/
	

 
end2first(T, V, V1, F, CTLcon) :-  
	(V1 =[], !, CTL=[], CTLcon=[];
		%findall([Q,L], (member(formula(con, Q, L), T), isend(L, V, V1)), EL), 
		findEnd(T, V, V1, EL), 
		(EL= [], !, CTL=[], CTL1=[]; 
			con2ctl(EL, F, V2, CTL), subtract(V1, V2, V11), append(F, CTL, F1), end2first(T, V, V11, F1, CTL1)),
		(CTL=[], CTLcon=[]; append(CTL, CTL1, CTL2), sort(CTL2, CTLcon))
	).

findEnd([], _, _, []).
findEnd([H1|L1], V, V1, EL) :- H1 = formula(Type, Q, L),
	(Type=con, 
		(isend(L, V, V1), H2=[Q, L]; H2 = []); H2=[]),
	findEnd(L1, V, V1, EL1),
	(H2=[], EL=EL1; append([H2], EL1, EL)).
		 

con2ctl([], _, [], []).
con2ctl([H1|L1], F, [Q|LL], [H2|L2]) :- H1=[Q, L], conCTL(L, F, F1), H2 = ctlformula(Q, F1), con2ctl(L1, F, LL, L2) .


 
conCTL([],_, true).
conCTL([H|L], CTL, F) :- H = (Type, T, I), disCTL(T, CTL, T1),
	(Type=ex, T2=(^(*(T1)));
		(Type=ax, T2=(~(*(T1)));
			(Type=ef, T2=(^(?(T1)));
				(Type=af, T2=(~(?(T1)));
					T2 =T1
				)
			)
		)
	),
	conCTL(L, CTL, F1), (F1=true, F=T2; F = (T2 & F1)).

 	
disCTL([], _, false).
disCTL([H|L], CTL, F) :- equall(- H, X), disCTL(L, CTL, F1),
	((member(ctlformula(H, F2), CTL), F3=F2; member(ctlformula(X, F2), CTL), F3=(- F2)), (F1=false, F=F3; F = (F3 \/ F1)); %/
		(F1=false, F=H; F = (H \/ F1)) %/
	).
/*disCTL([H|L], CTL, F) :- equall(- H, X), 
	((member(ctlformula(H, F2), CTL); member(ctlformula(X, F2), CTL)), disCTL(L, CTL, F1), (F1=false, F=F2; F = (F2 \/ F1)); %/
		disCTL(L, CTL, F1), (F1=false, F=H; F = (H \/ F1)) 
	).*/



 
isend([], _, _).
isend([H|T], V, V1) :- H=(Type, L, I), obtainAtom(L, L1), 
	(intersection(L1, V, []), intersection(L1, V1, []), !; !, fail).


 
toCTL_temprol([], _, []).
toCTL_temprol([H|T], Con, R) :- H = formula(Type, Q, L), 
	((Type = con; (Type = alpha; Type = beta)), H1=[];
		(Type = af, L=[(X, nil)], disjunction(X, X1), H2=(~(? X1)), H1=ctlformula(Q, H2);
			(Type = ef, L=[(X, I)], disjunction(X, X1), H2=(^(? X1)), H1=ctlformula(Q, H2);
				(Type = ax, L=[(X, nil)], disjunction(X, X1), H2=(~(* X1)), H1=ctlformula(Q, H2);
					(Type = ex, L=[([X], I)], member(ctlformula(X, X1), Con), H2=(^(* X1)), H1=ctlformula(Q, H2);
						(Type = au, L=[([X], [Y], I)], member(ctlformula(X, X1),Con), 			
							member(ctlformula(Y, Y1), Con), H2 = (~(X1 $ Y1)), H1=ctlformula(Q, H2);
							(Type = eu, L=[([X], [Y], I)], member(ctlformula(X, X1),Con), 			
								member(ctlformula(Y, Y1), Con), H2 = (^(X1 $ Y1)), H1=ctlformula(Q, H2);
								(Type=ag, L=[([X], I)], member(ctlformula(X, X1), Con), H2=(~(@ X1)), H1 = ctlformula(Q, H2);
									(Type=eg, L=[([X], I)], member(ctlformula(X, X1), Con), H2=(^(@ X1)), H1 = ctlformula(Q, H2);
										H1=[]
									)
								)
							)					
						)
					)
				)
			)
		)
	), toCTL_temprol(T, Con, R1),
	(H1=[], R=R1; append([H1], R1, R)).


 
toCTL_list([],_,[]).
toCTL_list([H|L], CTLcon, R) :-  H = formula(Type, X, LL),  
	(Type=con, conCTL(LL, CTLcon, F1), H1=ctlformula(X, F1);
		(Type=alpha, conA(LL, CTLcon, F1), H1=ctlformula(X, F1);
			(Type=beta, conB(LL, CTLcon, F1), H1=ctlformula(X, F1); H1=[])
		)
	), 
	toCTL_list(L, CTLcon, R1),
	(H1=[], R=R1; append([H1], R1, R)).
	

conA([], _, true).
conA([H1, H2|L], CTLcon, F) :- H1 = [[C3,C2],I,[C33,C2,C4],[C2,C3]], append(C3, C2, C32), 
	append(C33, C2, C332), append(C332, C4, C3324), append(C3324, H2, CH2),
	disCTL(C32, CTLcon, F1), disCTL(CH2, CTLcon, F2), F3=(F1 \/ (^(*(F2 \/ ~(*(~(? F1))))))), %/
	conCTL(L, CTLcon, F4), F = (F3 & F4).
	
 	
conB([], _, true).
conB([H1, H2|L], CTLcon, F) :- H1 = [[C3,C2],I,[C33,C2,C4],[C2,C3]], append(C3, C2, C32), 
	append(C33, C2, C332), append(C332, C4, C3324), append(C3324, H2, CH2),
	disCTL(C32, CTLcon, F1), disCTL(CH2, CTLcon, F2), %trace, simpDis(F1, F11), simpDis(F2, F22),
	F3=(F1 \/ (~(*(F2 \/ ~(*(~(? F1))))))), %/
	conCTL(L, CTLcon, F4), F = (F3 & F4).
	
	

conjunction([], true).
conjunction([H|L], F) :- H = (Type, T, I), disjunction(T, T1),
	(Type=ex, T2=(^(*(T1)));
		(Type=ax, T2=(~(*(T1)));
			(Type=ef, T2=(^(?(T1)));
				(Type=af, T2=(~(?(T1)));
					T2 =T1
				)
			)
		)
	),
	conjunction(L, F1), (F1=true, F=T2; F = (T2 & F1)).
	
disjunction([], false).
disjunction([H|L], F) :-  disjunction(L, F1), (F1 = false, F = H; F = (H \/ F1)). %/



 
obtainAB(R) :- findall(X, efImplication:alpha(_, X,_), R1), findall(Y, efImplication:beta(_, Y, _), R2), append(R1, R2, R).

 
replaceAB([], _, []).
replaceAB([H1|L1], AB, [H2|R1]) :- H1 = formula(Type, X, L),
	(Type=con, (member(X, AB), efImplication:snfclause(Type1, [X], T, I), 
		(Type1 = ex, 
			(member([N, alpha], T), efImplication:alpha(N, X, T2),  %T2=[
				subtract(T, [[N, alpha]], T3), append([T2], [T3], T33), 
				subtract(L, [(ex, T3, I)], LX), append(T33, LX, LL),
				H2 = formula(alpha, X, LL));
			(Type1=ax, (member([N, beta], T), efImplication:beta(N, X, T2),
				subtract(T, [[N, beta]], T3), append([T2], [T3], T33), 
				subtract(L, [(ax, T3, I)], LX), append(T33, LX, LL),
				H2 = formula(beta, X, LL)); H2=H1)
		); H2=H1); H2=H1),
	replaceAB(L1, AB, R1).
	%(H2=[], R=R1; append([H2],  R1, R)).




%#################
%tranCTL2SNF



%tran2SNFCl(F, V1, L): transform a CTL formula into a set of snf clauses
tran2SNFCl([], []).
tran2SNFCl([H|L],R) :- tranH2C(H, C), tran2SNFCl(L, R1), append([C], R1, R).

tranH2C(start -> P, C) :- prop(P), C = start_clause([start], [P]).
tranH2C(Q -> P, C) :- is_dis(P), negation([Q], Q1),  wff2cnf(P, P1), P1 = [P2], append(Q1, P2, P3), C = true_clause([true], P3).
tranH2C([Q -> ^(*P),I], C) :- wff2cnf(P, P1), P1 = [P2], C = exist_next_clause([Q], P2, I).
tranH2C([Q -> ^(?P),I], C) :- wff2cnf(P, P1), P1 = [P2], C = exist_future_clause([Q], P2, I).
tranH2C(Q -> ~(*P), C) :- wff2cnf(P, P1), P1 = [P2], C = global_next_clause([Q], P2).
tranH2C(Q -> ~(?P), C) :- wff2cnf(P, P1), P1 = [P2], C = global_future_clause([Q], P2).


tran2SNF([], Ind, N, V, []) :- V=[].
tran2SNF(L, Ind, N, V, R) :- tran2SNF_list(L, V1, Ind, N, IndM, NM, L1),   
	(L1 = L, R = L1, V = V1;
		tran2SNF(L1, IndM, NM, V2, R), append(V1, V2, V)).

tran2SNF_list([], V, Ind, N, Ind, N, []) :- V=[].
tran2SNF_list([H|T], V1, Ind, N, IndM, NM, L) :- transF(H, Ind, V2, N, IndM1, NM1, L1), 
	tran2SNF_list(T, V3, IndM1, NM1, IndM, NM, L3), append(V2, V3, V1), append(L1, L3, L4), sort(L4, L).
	

%transF(start -> P, Ind, V, Ind, IndM, IndM, C) :-  C = start_clause([start], [P]).
%transF(Q -> P, Ind, V, Ind, IndM, IndM, C) :- is_dis(P), negation([Q], Q1),  wff2cnf(P, P1), P1 = [P2], append(Q1, P2, P3), %C = true_clause([true], P3).
transF(Q -> (P1 & P2), Ind, V, N, IndM, NM, L) :- NM = N, IndM=Ind, L = [Q -> P1, Q -> P2], V=[].
transF(Q -> (P1 \/ P2), Ind, V, N, IndM, NM, L) :-  IndM=Ind,%/ 
	(is_dis(P1), !, L1=[], N1=N, V2=[];
		(N1 is N +1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), V2=[Y1], L1=(Y1 -> P1))
	),
	(is_dis(P2), !, N2=N1, NM=N1, V3=[], L2=[];
		(N2 is N1 +1, NM=N2, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L1=[], (L2=[], L =[Q -> (P1 \/ P2)], V = [];  %/
				L = [L2, Q -> P1 \/ Y2], V = [Y2]); %/
		(L2 =[], L = [L1, Q -> Y1 \/P2], V = [Y1];  %/
			L = [L1, L2, Q -> Y1 \/ Y2], V =[Y1, Y2] %/
		)
	).
	
	
/*transF(Q -> ^(P1 $ P2), Ind, V, N, IndM, NM, L) :-   IndM is Ind+1,
	(is_lit(P1), !, L1=[], N1=N, V2=[];
		(N1 is N +1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), V2=[Y1], L1=(Y1 -> P1))
	),
	(is_lit(P2), !, N2=N1, V3=[], L2=[];
		(N2 is N1 +1, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L1=[], (L2=[], tran2U(Q -> ^(P1 $ P2), N, IndM, Vx, L), V = [Vx], NM is N +1;  
					tran2U(Q -> ^(P1 $ Y2), N2, IndM, Vx, R1), append([L2], R1, L), V = [Y2, Vx], NM is N2 +1); 
		(L2 =[], tran2U(Q -> ^(Y1 $ P2), N2, IndM, Vx, R2), append([L1], R2, L), V = [Y1, Vx], NM is N + 1;
			 NM is N2 + 1, string_concat('x', NM, XY3), string_to_atom(XY3, Y3), assert(prop(Y3)),
			L = [L1, L2,  Q -> Y2 \/ Y3, Y3 -> Y1, [Y3 -> ^(*(Y2 \/ Y3)), IndM], [Q -> ^(? Y2),IndM]], V =[Y1, Y2, Y3]  
		)
	).*/
	
transF(Q -> ^(P1 $ P2), Ind, V, N, IndM, NM, L) :-   IndM is Ind+1,
	(is_lit(P2), !, N2=N, V3=[], L2=[];
		(N2 is N +1, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L2 =[], tran2U(Q -> ^(P1 $ P2), N2, IndM, Vx, L), V = [Vx], NM is N + 1;
			 NM is N2 + 1, tran2U(Q -> ^(P1 $ Y2), N2, IndM, Vx, R), append([L2], R, L), V = [Y2, Vx]
	).
	
/*transF([Q -> ^(P1 $ P2), I], Ind, V, N, IndM, NM, L) :-  IndM = Ind,
	(is_lit(P1), !, L1=[], N1=N, V2=[];
		(N1 is N +1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), V2=[Y1], L1=(Y1 -> P1))
	),
	(is_lit(P2), !, N2=N1, V3=[], L2=[];
		(N2 is N1 +1, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L1=[], (L2=[], tran2U(Q -> ^(P1 $ P2), N, I, Vx, L), V = [Vx], NM is N +1;  
					tran2U(Q -> ^(P1 $ Y2), N2, I, Vx, R1), append([L2], R1, L), V = [Y2, Vx], NM is N2 +1); 
		(L2 =[], tran2U(Q -> ^(Y1 $ P2), N2, I, Vx, R2), append([L1], R2, L), V = [Y1, Vx], NM is N + 1;
			 NM is N2 + 1, string_concat('x', NM, XY3), string_to_atom(XY3, Y3), assert(prop(Y3)),
			L = [L1, L2,  Q -> Y2 \/ Y3, Y3 -> Y1, [Y3 -> ^(*(Y2 \/ Y3)), I], [Q -> ^(? Y2),I]], V =[Y1, Y2, Y3]  
		)
	).*/
	
transF(Q -> ~(P1 $ P2), Ind, V, N, IndM, NM, L) :-   IndM is Ind,
	(is_lit(P2), !, N2=N, V3=[], L2=[];
		(N2 is N +1, string_concat('x', N2, XY2), string_to_atom(XY2, Y2), assert(prop(Y2)), V3=[Y2], L2=(Y2 -> P2))
	), 
	(L2 =[], tran2Uall(Q -> ~(P1 $ P2), N2, IndM, Vx, L), V = [Vx], NM is N + 1;
			 NM is N2 + 1, tran2Uall(Q -> ~(P1 $ Y2), N2, IndM, Vx, R), append([L2], R, L), V = [Y2, Vx]
	).
	
tran2U(Q -> ^(Y1 $ Y2), N, IndM, V, L) :- N1 is N+1, string_concat('x', N1, XY3), 
	string_to_atom(XY3, Y3), assert(prop(Y3)), V= Y3,
	L = [Q -> Y2 \/ Y3, Y3 -> Y1, [Y3 -> ^(*(Y2 \/ Y3)), IndM], [Q -> ^(? Y2),IndM]]. %/

tran2Uall(Q -> ~(Y1 $ Y2), N, IndM, V, L) :- N1 is N+1, string_concat('x', N1, XY3), 
	string_to_atom(XY3, Y3), assert(prop(Y3)), V= Y3,
	L = [Q -> Y2 \/ Y3, Y3 -> Y1, Y3 -> ~(*(Y2 \/ Y3)), Q -> ~(? Y2)]. %/	
	
	

transF(Q -> ^(@P), Ind, V, N, IndM, NM, L) :-  N1 is N+1, NM=N1, IndM is Ind+1,
	string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)),
	V=[Y1], L=[Q -> Y1, Y1 -> P, [Y1 -> ^(*Y1), IndM]].
/*transF([Q -> ^(@P),I], Ind, V, N, IndM, NM, L) :-  N1 is N+1, NM=N1, IndM is Ind+1,
	string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)),
	V=[Y1], L=[Q -> Y1, Y1 -> P, [Y1 -> ^(*Y1), IndM]].*/
	
transF(Q -> ~(@P), Ind, V, N, IndM, NM, L) :-  N1 is N+1, NM=N1, IndM is Ind,
	string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)),
	V=[Y1], L=[Q -> Y1, Y1 -> P, Y1 -> ^(*Y1)].
	
transF(Q -> ^(*P), Ind, V, N, IndM, NM, L) :- IndM is Ind +1,
	(is_dis(P), L = [[Q -> ^(*P), IndM]], NM is N, V=[];
		N1 is N+1, NM=N1,  string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[[Q -> ^(* Y1),IndM], Y1 -> P]
	).
/*transF([Q -> ^(*P), I], Ind, V, N, IndM, NM, L) :- IndM=Ind,
	(is_dis(P), L = [[Q -> ^(*P),I]];
		N1 is N+1, NM=N1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[[Q -> ^(* Y1),I], Y1 -> P]
	).*/

transF(Q -> ~(*P), Ind, V, N, IndM, NM, L) :- IndM is Ind,
	(is_dis(P), L = [Q -> ~(*P)], NM is N, V=[];
		N1 is N+1, NM=N1,  string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[Q -> ~(* Y1), Y1 -> P]
	).
	
transF(Q -> ^(?P), Ind, V, N, IndM, NM, L) :- IndM is Ind+1, 
	(is_lit(P), L = [[Q -> ^(?P), IndM]], NM is N, V=[];
		N1 is N+1, NM=N1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[[Q -> ^(? Y1), IndM], Y1 -> P]
	).
/*transF([Q -> ^(?P), I], Ind, V, N, IndM, NM, L) :- 
	(is_lit(P), IndM=Ind, L = [[Q -> ^(?P), I]];
		N1 is N+1, NM=N1, IndM is Ind+1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[[Q -> ^(? Y1), IndM], Y1 -> P]
	).*/

transF(Q -> ~(?P), Ind, V, N, IndM, NM, L) :- IndM is Ind, 
	(is_lit(P), L = [Q -> ~(?P)], NM is N, V=[];
		N1 is N+1, NM=N1, string_concat('x', N1, XY1), string_to_atom(XY1, Y1), assert(prop(Y1)), 
		V = [Y1], L =[Q -> ^(? Y1), Y1 -> P]
	).
	
	
transF(Q -> P, Ind, V, N, IndM, NM, L) :- IndM=Ind, is_lit(P), V =[], NM =N, L=[Q -> P].
transF(Q, Ind, V, N, IndM, NM, [Q]) :- NM =N, IndM=Ind, V = [].



 
 is_lit(C) :- (prop(C); equall(-C, D), prop(D)).
 
  
 is_dis(C) :- is_lit(C).
 is_dis(C) :- C = (P \/ Q), %/
	(is_dis(P), is_dis(Q), !; !, fail).
	
	
	
	
	
	
	
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
	
complement(true,false) :- !.
complement(false,true) :- !.
complement(P,-P) :- prop(P).
complement(-P,P) :- prop(P).
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
%###################
%replaceTrans

 

%Trans(3)
re2CTL_con(V) :- findall([Q, L1], formula(snf, Q, L1), R),
	(R =[], !; con_list(R)).
	
onlyTrue([]).
onlyTrue([H|T]) :- H = (true, L1, Type), onlyTrue(T).

con_list([]).
con_list([H|T]) :- H = [Q, L1],  assert(formula(con, Q, L1)), con_list(T).


/*re2CTL(V) :- findall([Q, P, X], (formula(T1, Q, L1), (not(member(P, V)); formula(T2, P, LP)), 
	(not(member(X, V)); formula(T3, X, LX)),
	(L1=[(true, [X, P], nil)]; L1=[(true, [P, X], nil)])), R),
	(R =[], !; dis_list(R)).
	
dis_list([]).
dis_list([H|T]) :- H = [Q, P, X],  assert(formula(dis, Q, [([X], [P])])), dis_list(T).*/


%Trans(7)
re2CTL_future(V) :- findall([Q, P, L1, L2, Type, I], (formula(snf, Q, L1), formula(snf, P, L2), member((Type, [P], I), L1)), R), !,
	(R=[], !; future_list(R)).
	
future_list([]).
future_list([H|T]) :- H = [Q, P, L1, L2, Type,I], 
	((Type=ef; Type=af), subtract(L1, [(Type, [P], I)], L11), retractall(formula(snf, Q, L1)), (L11=[], !; assert(formula(snf, Q, L11))); !),  
	(Type = ef, !, assert(formula(ef, Q, [([P], I)]));
		(Type = af, !, assert(formula(af, Q, [([P], I)])); !)
	), future_list(T).


%Trans(6)
re2CTL_next(V) :- findall([Q, P, L1, L2, Type, I], (formula(snf, Q, L1), formula(snf, P, L2), member((Type, [P], I), L1)), R), !,
	(R=[], !; next_list(R)).
	
next_list([]).
next_list([H|T]) :- H = [Q, P, L1, L2, Type,I], 
	((Type=ex; Type=ax), subtract(L1, [(Type, [P], I)], L11), retractall(formula(snf, Q, L1)), (L11=[], !; assert(formula(snf, Q, L11))); !),  
	(Type = ex, !, assert(formula(ex, Q, [([P], I)]));
		(Type = ax, !, assert(formula(ax, Q, [([P], I)])); !)
	), next_list(T).



%Trans(10)
re2CTL_global(V) :- findall([Q, P, L1, L2], (formula(snf, Q, L1), formula(snf, P, L2), member((true, [P], nil), L1)), R), !,
	(R=[], !; global_list(R)).
	
global_list([]).
global_list([H|T]) :- H=[Q,P,L1,L2], 
	(member((ex, [P], I), L2), !, assert(formula(eg, Q, [([P], I)])), delete(L2, (ex, [P], I), L3), 
		retractall(formula(snf, P, L2)), assert(formula(snf, P, L3));
		(member((ax, [P], I), L2), !, assert(formula(ag, Q, [([P], I)])), delete(L2, (ax, [P], I), L3), 
			retractall(formula(snf, P, L2)), assert(formula(snf, P, L3)); !
		)
	),global_list(T).
	

%Trans(11) and (8)
re2CTL_until(V) :- findall([Q,P, L, L1, L2], (formula(snf, Q, L1), member((true, [L, P], nil), L1), formula(snf, P, L2), 
	(not(member(L, V)); formula(snf, [L], LL))), R), !,
	(R=[], !; until_list(R)).
/*(R=[], !; 
	(member((ex, [L, P], I), L2), !, 
		(member((ef, [L], I), L1), assert(formula(eu, [Q], [([P], [L], I)])), delete(L2, (ex, [L], I), L3), 
			retractall(formula(snf, [P], L2)), assert(formula(snf, [P], L3)); !);
		(member((ax, [L, P], nil), L2), !, (member((af, [L], nil), L1), assert(formula(au, [Q], [([P], [L], nil)])), 
			delete(L2, (ax, [L], nil), L3), retractall(formula(snf, [P], L2)), assert(formula(snf, [P], L3)); !); !
		)
	)
).*/

until_list([]).
until_list([H|T]) :- H = [Q,P, L, L1, L2],
	(member((ex, [L, P], I), L2), !, 
		(member((ef, [L], I), L1), assert(formula(eu, Q, [([P], [L], I)])), delete(L2, (ex, [L,P], I), L3), 
			retractall(formula(snf, P, L2)), assert(formula(snf, P, L3)); !);
		(member((ax, [L, P], nil), L2), !, (member((af, [L], nil), L1), assert(formula(au, Q, [([P], [L], nil)])), 
			delete(L2, (ax, [L,P], nil), L3), retractall(formula(snf, P, L2)), assert(formula(snf, P, L3)); !); !
		)
	), until_list(T).



 

re2CTL(V, R) :- re2CTL_until(V), re2CTL_global(V), re2CTL_future(V), re2CTL_next(V), re2CTL_con(V),
	obtainF(R1), write("\n **************\nformula:"), write(R1), write("\n **************\n"), deleteTemp(R1, R), 
	write("\n #################\n CTL formula:"), write(R), write("\n #################\n").
	
obtainF(R) :- findall(formula(X, Y,Z), formula(X, Y,Z), F), sort(F, R).

deleteTemp([], []).
deleteTemp([H|T], R) :- H = formula(Type, Q, L1), 
	(Type=con, !, H1=H; %(member((true, X, I), L1), temp_list(L1, L2); L2=L1), H1=formula(con, Q, L2), deleteTemp(T, T1); 
		(Type=snf, H1=[]; H1=H)
	),
	deleteTemp(T, R1),
	(H1=[], R= R1; append([H1], R1, R)).

temp_list([], []).
temp_list([H|T], R) :- H = (Type, Q, I), 
	((Type=true; Type=start), C =H; C = []),
	temp_list(T, R1),
	(C=[], R=R1; append([C], R1, R)).








