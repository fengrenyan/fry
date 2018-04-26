%-----define the operators that will be used in programm------

:- op(300,fx,@). /*@K*/
:- op(400,fx,$). /*$B*/

:- op(500,fx,-).
%:- op(500,fx,~).
:- op(780,xfy,&).
:- op(790,xfy,\/).
:- op(800,xfy,->).
:- op(800,xfy,<->).

%---------end define-----------------------------------
%


%-----------define the dynamic predicts------------------
:- dynamic(prop/1).
:- dynamic(klist/1).
:- dynamic(blist/1).
:- dynamic(pair/2).

%---------------------end define---------------------------

%---------------the init Set of facts------------------
prop(p1).
prop(d).
prop(a).
prop(b).
prop(c).
prop(true).
prop(false).
%--------------------end-----------------------------


%----------compute one in C-----2017-5-16------
%
ksnc(Q, P, Variable_n, Clause_num) :-
	open('Cresult.txt', write, Str), trace, 
	generate_KBterm(Clause_num), trace, 
	read_formula(Formula), trace, write("Formula="), write(Formula), nl,
	%random(1, Variable_n, Q),
	%generate_P(Variable_n, P),
	%generate_P_list(Variable_n, PN, P),
	generate_list(1, Variable_n, L1),  trace, write("L1="), write(L1), nl, % generate the set of prop in formula.
	difference_list(P, L1, P1), trace, write("first P1="), write(P1), nl,
	time_output_1(kForget(Formula&Q, P1, Result), U, V, W),
	write("P="), write(P), nl,
	write("Q="), write(Q), nl,
	write("P1="), write(P1), nl,
	write("T="), write(Formula), nl,
	write(Str,'P='), write(Str, P), nl(Str), %write(Str,"\n"),
	write(Str,'Q='), write(Str, Q), nl(Str), %write(Str,"\n"),
	write(Str,'P1='), write(Str, P1), nl(Str), %write(Str,"\n"),
	write(Str,'Theory='), write(Str, Formula), nl(Str), %write(Str,"\n"),
	write(Str,'UsedInf:'), write(Str, U), nl(Str), %write(Str,"\n"),
	write(Str,'UsedTime:'), write(Str, V), nl(Str), %write(Str,"\n"),
	write(Str,'Wall:'), write(Str, W), nl(Str), %write(Str,"\n"),
	write(Str, 'SNC is: '), write(Str, Result), nl(Str), %write(Str,"\n"),
	write(Str,'...............................................'), nl(Str),
	retractall(klist(_)),
	retractall(blist(_)),
	trim_stacks,   % 2017-5-4 Release  stack memory resources that are not in use at  this moment
	%compute(L, LResult, Str),
	close(Str).



%-------------------end-------------------

%---------------test simplify-----------------------
simplify :-
	gain_filename(less_filename, L),
	open('sim_result.txt', write, Str),
	sim_compute(L, R, Str),
	close(Str),
	write("succeed"), nl.

sim_compute([], [], Str):- !.
sim_compute([X | L], [NewDNF| LResult], Str) :-
	X = [Filename , Variable_n, Clause_num],
	generate_KBterm(Clause_num),
	read_formula(Formula),
	new_substitute_allk(Formula, F1, 1),
	new_substitute_allb(F1, F2, 1),
	write(Str, 'primary formula: '), write(Str, Formula), nl(Str),
	time_output_1((new_wff2dnf(F2, DNF), simply_dnf(DNF, NewDNF)), U, V, W),
	write(Str, 'primary DNF: '), write(Str, DNF), nl(Str),
	write(Str,'UsedInf:'), write(Str, U), nl(Str), %write(Str,"\n"),
	write(Str,'UsedTime:'), write(Str, V), nl(Str), %write(Str,"\n"),
	write(Str,'Wall:'), write(Str, W), nl(Str), %write(Str,"\n"),
	write(Str,'the result of simplifing is: '), write(Str, NewDNF), nl(Str),
	write(Str,'...............................................'), nl(Str),
	retractall(pair(_,_)),
	retractall(klist(_)),
	retractall(blist(_)),
	sim_compute(L, LResult, Str).

%------------------------end test---------------------



% ----------compute the snc of Q on P under T, and the result
% is F-----------

snc(Q, P, T, F) :- 
	gain_prop(T, PN1), write("PN1="), write(PN1), nl, 
	add_propTerm(PN1), append(PN1, [Q],OP), write("OP="), write(OP), nl, 
	difference_list(P, OP, P1), write("P1="), write(P1), nl, 
	kForget(T&Q, P1, F),
	retractall(prop(X)).

	
	%gain the atom of a formula
gain_prop(P & Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), 
	unionSet(L1, L2, L).
gain_prop(P \/ Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), %/
	unionSet(L1, L2, L).
gain_prop(@P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(-P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop($P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(P, [P]) :- \+ atom(P), !.

%compute the union of two sets
unionSet(L1, L2, L) :- append(L1, L2, L11), sort(L11, L).

add_propTerm([]) :- !.
add_propTerm([X|L]) :- assert(prop(X)), add_propTerm(L).

%----------------end---------------------------------------




%-----------read file to gain the filename of formula-------
%obtain the filename and the number of variables in formula.
gain_filename(File, L) :-
	string_concat(File, '.txt', Filename1),
	string_to_atom(Filename1, Filename),
	read_filename(Filename, L).

read_filename(File, L) :-
	open(File, read, Str),
	read_filename_list(Str, L),
	close(Str).

read_filename_list(Stream, []) :-
	at_end_of_stream(Stream).
	%write("there is no filename!\n").

read_filename_list(Stream, [X|L]) :-
    \+ at_end_of_stream(Stream),
    read(Stream, X1),
    read(Stream, X2),
    read(Stream, X3),
    X = [X1, X2, X3],
    read_filename_list(Stream, L).

%-------------------end read----------------------



%gain the atom of a formula
gain_prop(P & Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), 
	unionSet(L1, L2, L).
gain_prop(P \/ Q, L) :- gain_prop(P, L1), gain_prop(Q, L2), %/
	unionSet(L1, L2, L).
gain_prop(@P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(-P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop($P, L) :- gain_prop(P, L1), sort(L1, L).
gain_prop(P, [P]) :- \+ atom(P), !.

%compute the union of two sets
unionSet(L1, L2, L) :- append(L1, L2, L11), sort(L11, L).

add_propTerm([]) :- !.
add_propTerm([X|L]) :- assert(prop(X)), add_propTerm(L).

%---------------------compute the kForget cyclicly----------
/*compute the snc of q on P onduer T in each file. snc(q, P, T) = kForget(T&q, Var(T&q)\P, Result).*/
%the q and P are generated randomly.
compute([], [], Str):- !.
compute([X | L], [Result| LResult], Str) :-
	trim_stacks,   % 2017-5-4 Release  stack memory resources that are not in use at  this moment
	X = [Filename , Variable_n, Clause_num],
	generate_KBterm(Clause_num),
	read_formula(Formula),
	random(1, Variable_n, Q),
	generate_P(Variable_n, P),
	generate_list(1, Variable_n, L1),  % generate the set of prop in formula.
	difference_list(P, L1, P1),
	time_output_1(kForget(Formula&Q, P1, Result), U, V, W),
	write("P="), write(P), nl,
	write("Q="), write(Q), nl,
	write("P1="), write(P1), nl,
	write("T="), write(Formula), nl,
	write(Str,'P='), write(Str, P), nl(Str), %write(Str,"\n"),
	write(Str,'Q='), write(Str, Q), nl(Str), %write(Str,"\n"),
	write(Str,'P1='), write(Str, P1), nl(Str), %write(Str,"\n"),
	write(Str,'T='), write(Str, Formula), nl(Str), %write(Str,"\n"),
	write(Str,'UsedInf:'), write(Str, U), nl(Str), %write(Str,"\n"),
	write(Str,'UsedTime:'), write(Str, V), nl(Str), %write(Str,"\n"),
	write(Str,'Wall:'), write(Str, W), nl(Str), %write(Str,"\n"),
	write(Str, 'Result is: '), write(Str, Result), nl(Str), %write(Str,"\n"),
	write(Str,'...............................................'), nl(Str),
	retractall(klist(_)),
	retractall(blist(_)),
	%trim_stacks,   % 2017-5-4 Release  stack memory resources that are not in use at  this moment
	compute(L, LResult, Str).

%--------------------end compute------------------------------


/*%------------generate K and B term for substitute---2017-5-2------------------
generate_KBterm(Input) :- Num is 3*Input, generate_list(1,Num, L),
	change2KBterm(L, KL, BL),
	assert(klist(KL)),
	assert(blist(BL)).

change2KBterm([], [], []) :- !.
change2KBterm([X|L] , [KX|KL], [BX|BL]) :-
	string_concat('k', X, KY),
	string_to_atom(KY, KX),
	string_concat('b', X, BY),
        string_to_atom(BY, BX),
	change2KBterm(L, KL, BL).
%-------------end generate K and B term for substitute------------------*/

%----------------------------new generate K and B term for substitute---2017-5-2------
generate_KBterm(N) :- lenKlist(N, LK), generate_list(1, LK, L1),
	change2Kterm(L1, KL),
	lenBlist(N, LB), generate_list(1, LB, L2),
	change2Bterm(L2, BL),
	assert(klist(KL)),
	assert(blist(BL)).

change2Kterm([], []) :- !.
change2Kterm([X|L], [KX|KL]) :- string_concat('k', X, KY),
	string_to_atom(KY, KX),
	change2Kterm(L, KL).


change2Bterm([], []) :- !.
change2Bterm([X|L], [BX|BL]) :- string_concat('b', X, BY),
	string_to_atom(BY, BX),
	change2Bterm(L, BL).


lenKlist(1, 2) :- !.
lenKlist(N, L) :- exp(3, N, E),
	exp(2,N, E1),
	comb(N, 2, X),
	L is 2*E-E1-2*X.

exp(X, 1, R) :- R is X, !.
exp(X, Y, R) :-Y1 is Y-1, exp(X, Y1, R1), R is X*R1, !.

comb(N, 2, R) :- N1 is N*(N-1), R is N1/2.


lenBlist(1, 1) :- !.
lenBlist(N, L) :- exp(3, N, E), comb(N, 2, C), L is E-C, !.

%----------------------------end --------------------------------------------


%------------------------read formula from file----------------
%read formula from file and assert the atom which in formula
read_formula(Formula) :- write_ln('input: '), read(user_input,Input), %consult(Input),
	string_concat(Input, '.txt', Filename1),
	string_to_atom(Filename1, Filename),
	read_file(Filename, Formula).


read_file(Filename, []) :- open(Filename, read, Str), at_end_of_stream(Str), write("there is no formula\n"), !.
read_file(Filename, Formula) :-
	open(Filename, read, Str),
	read(Str, Formula),
	read_prop(Str).


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


%-------------generate P---------------------------
%generate the set with a random number of prop randomly and |P| >= 1.
generate_P(Variable_n, L) :- random(1, Variable_n, N),
	generate_P_list(Variable_n, N, L).

indicatPL(Variable_n, N, L) :- 
		generate_P_list(Variable_n, N, L1),
		substract(L1, L2),
		lenth(L2, X),
		(X=N, L=L1; N1 is N-X, indicatPL(Variable_n, N1, L2) , append(L1, L2, L)).

generate_P_list(Variable_n, 0, []) :- !.
generate_P_list(Variable_n, N, [X|L]) :-
	random(1, Variable_n, X),
	N1 is N-1,
	generate_P_list(Variable_n, N1, L).

%----------------end-----------------------------------------



%-----------------------generate list------------------
%generate a list L which is between 1 to N.
generate_list(1, N, []) :- N < 1, !.
generate_list(1, 1, [1]) :- !.
generate_list(1, N, [N|L]) :- X is N-1, generate_list(1, X, L), !.

%-------------------------end-----------------------------

%-----------------------generate P'--------------------
%L = L2 - L1 and L2 >= L1(must be true)
difference_list([],L2, L2) :- !.
%difference_list(L1, L2, L) :-  (bagof( X, (member(X, L2), \+ member(X, L1)), L3), L = L3; L = []), !.
difference_list(L1, L2, L) :-  sort(L1, L11), sort(L2, L22), length(L11, X1), length(L22, X2), 
(X1 >= X2, L= []; bagof( X, (member(X, L22), \+ member(X, L11)), L)), !.

%--------------------------end-------------------------



%------------------output the time of solving a problem-----------
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
/*
time_true(State0) :-
	report_output(State0, 12, X, Y, Z).		% leave choice-point
time_true(State) :-
	get_time(Wall),
	statistics(cputime, Time),
	statistics(inferences, Inferences0),
	plus(Inferences0, -3, Inferences),
	nb_setarg(1, State, Wall),
	nb_setarg(2, State, Time),
	nb_setarg(3, State, Inferences),
	fail.*/

%------------------------- end -------------------------------


%----------------simplify the dnf-------------


simply_dnf(DNF, NEWDNF) :- sort_list(DNF, DNF1),   write("DNF1="), write(DNF1), nl,
%	dnf_init_minimize(DNF1,[], NEWDNF),
	!, dnf_simplifies(DNF1, NEWDNF).

sort_list([], []) :- !.
sort_list([D|DNF], [NewD|NewDNF]) :-
	sort(D, NewD),
	sort_list(DNF, NewDNF).



/*dnf_init_minimize(DNF1, DNF2, NewDNF) :- member(D, DNF1),
	difference(DNF1, D, DNF3), union(DNF3, DNF2, DNF4),
	dnf_simplifies(DNF4, DNF5), member(D1, DNF5), dnf_subsumes(D1, D), !,
	dnf_init_minimize(DNF3, DNF2, NewDNF).
dnf_init_minimize(DNF, [], NewDNF):- dnf_simplifies(DNF, NewDNF).*/

/* cnf_simplifies(L,NewL): simplifies DNF L into L1 by eliminating repetitive
   literals and contradictive literals. It also eliminates clauses that
   subsumes others.
DNFCNF L L1L
*/


dnf_simplifies(DNF, NewDNF) :- dnf_sort_simpl(DNF, DNF1),    write("DNF1="), write(DNF1), nl,
	dnf_unit_elim(DNF1, DNF2),    write("DNF2="), write(DNF2), nl,
	dnf_sort_simpl(DNF2, DNF3),    write("DNF3="), write(DNF3), nl,
	check_tauto(DNF3, DNF4),    write("DNF1="), write(DNF1), nl,
	dnf_minimize(DNF4, NewDNF).

/*dnf_simplifies2(DNF, NewDNF) :- dnf_sort_simpl(DNF, DNF1), %time waster
	dnf_minimize(DNF1, DNF2),
	dnf_unit_elim(DNF2, DNF3),
	dnf_sort_simpl(DNF3, DNF4),
	check_tauto(DNF4, DNF5),
	dnf_minimize(DNF5, NewDNF).*/

/*dnf_simplifies(DNF, NewDNF) :-
	dnf_sort_simpl(DNF, DNF1), %ignorate the negative atom.
	dnf_minimize(DNF1, DNF2),
%	dnf_unit_elim(DNF2, DNF3),
	dnf_sort_simpl(DNF2, DNF3),
	check_tauto(DNF3, DNF4),
	dnf_minimize(DNF4, NewDNF).*/




%sort(List, NL):  Duplicates in List are removed.

/* sort each clause and get rid of contradictory */
dnf_sort_simpl([], []).
dnf_sort_simpl(DNF, [[]]) :- member([], DNF), !.
dnf_sort_simpl(DNF, [[]]) :- member([true], DNF), !.
dnf_sort_simpl(DNF, [[]]) :- member([-false], DNF), !.
dnf_sort_simpl([C|DNF], NewDNF) :- member(false, C), !,
	dnf_sort_simpl(DNF, NewDNF).
dnf_sort_simpl([C|DNF], NewDNF) :- member(-P, C), member(P, C), !,
	dnf_sort_simpl(DNF, NewDNF).
dnf_sort_simpl([C|DNF], NewDNF) :-
	sort(C, NewD1), subtract(NewD1, [-false, true], NewD),
        ((NewD=[],NewDNF=[[]]);
	 (dnf_sort_simpl(DNF, NewDNF1), NewDNF=[NewD|NewDNF1])),!.


/* eliminates all unit clauses */
dnf_unit_elim(DNF, [[P]|NewDNF]) :-
	member([P], DNF),
	dresolve(DNF, P, DNF2), !,
	dnf_unit_elim(DNF2, NewDNF), !.
dnf_unit_elim(DNF, DNF) :- !.



dresolve([],_,[]).
dresolve([D|DNF], P, NewDNF) :- member(P, D), !,
%	difference(D, P, NewD),
	dresolve(DNF, P, NewDNF).
dresolve([D|DNF], P, [NewD|NewDNF]) :- complement(P, P1),
	difference(D, P1, NewD),
	dresolve(DNF, P, NewDNF).

%check for tautologies
check_tauto(DNF, [[]]) :- member([], DNF), !.
check_tauto(DNF, DNF).


dnf_minimize(DNF, DNF1) :- member(L, DNF), member(L1, DNF), L\==L1,
	dnf_subsumes(L, L1), difference(DNF, L1, DNF2),
	dnf_minimize(DNF2, DNF1), !.
dnf_minimize(DNF, DNF).



%subset(L1, L2): if L1 is a subset of L2, then return true, else fales.
dnf_subsumes(L,L1) :- \+ (member(X,L), X\==true, X\==(-false),
    \+ member(X,L1)).

%-------------------end simplify --------------------------



% ----------compute the result of forget some atom P from T, and the
% result is NewF--------------
kForget(F,[],F) :-   write("F="), write(F), nl,!.
kForget(F,P, NEW) :-
	%trim_stacks,
	degree(F, N), N < 1, !,  
	new_wff2cnf(F, CNF),
	forget(CNF, P, NEW1),
	forget_result_simplify(NEW1,NEW2),
	cnf2wff(NEW2, NEW).
kForget(F,P,NEWF) :- first_degree(F, Y), !,     write("Y="), write(Y), nl,
	new_substitute_allk(Y, Y1, 1),    write("Y1="), write(Y1), nl,
	new_substitute_allb(Y1, Y2, 1),    write("Y2="), write(Y2), nl,
	new_wff2dnf(Y2, DNF),    write("DNF="), write(DNF), nl,
	simply_dnf(DNF, F1),   write("F1="), write(F1), nl,
	kforget(F1,P,NEW),  trace, write("NEW="), write(NEW), nl,
	eliminate_true(NEW,TNEW),
%	eliminate_false(TNEW,NEWF),
	eliminate_false(TNEW,NEW2),
	retractall(pair(_,_)),
	new_substitute_allk(NEW2,NEW3, 1),
	new_substitute_allb(NEW3, NEW4, 1),
	new_wff2dnf(NEW4, NEW5),
	simply_dnf(NEW5, NEW6),
	next_getTerm(NEW6, New7),
	dnf2wff(New7, NEWF),
	retractall(prop(_)),
	retractall(pair(_,_)).

%---------------end--------------------------------------

%---------------get the term that has been substituted-------------
next_getTerm([],[]) :- !.
next_getTerm([D|DNF], [NewD|NewDNF]) :- split_k(D, L1, L2, L3),
	kformula(L1, KF1),
	bformula(L2, KF2),
	append(KF1, KF2 , L4),
	append(L4, L3, NewD),
	next_getTerm(DNF, NewDNF), !.

%------------------end get -----------------------------

%-----------------degree----------------------------------
% compute the defree of modal formula, i.e. the maximum number of nested
% modal operator.
degree(F, 0) :- prop(F), !.
degree(- F, N) :- degree(F, N), !.
degree(P \/ Q, N) :- degree(P, N1), degree(Q, N2), (N1 >= N2, N = N1; N=N2), !.
degree(P & Q, N) :- degree(P, N1), degree(Q, N2), (N1 >= N2, N = N1; N = N2), !.
degree(@F, N) :- degree(F, N1), N is N1+1.
degree($F, N) :- degree(F, N1), N is N1+1.

%--------------end degree---------------------------------


%-----------------first_degree----------------------------
%simplify the formula who's degree is lager than 1 to 1
first_degree(X, Y) :- eliminate_operator(X, X1), eliminate_neg(X1, X2), r1to4(X2, X3), hsimply_degree(X3, Y).

%eliminate_operator: eliminate those operators, i.e. '->', '<->'
eliminate_operator(P, P) :- prop(P), !.
eliminate_operator(-P, -P) :- prop(P), !.
eliminate_operator(P->Q, (-F1)\/(F2)) :- eliminate_operator(P, F1), eliminate_operator(Q, F2), !.
eliminate_operator(P<->Q, (F1)&(F2)) :- eliminate_operator(P->Q, F1), eliminate_operator(Q->P, F2), !.
eliminate_operator(P&Q, (F1)&(F2)) :- eliminate_operator(P, F1), eliminate_operator(Q, F2), !.
eliminate_operator(P\/Q, (F1)\/(F2)) :- eliminate_operator(P, F1), eliminate_operator(Q, F2), !.
eliminate_operator(@P, @F) :- eliminate_operator(P, F), !.
eliminate_operator($P, $F) :- eliminate_operator(P, F), !.
eliminate_operator(-P, -(F)) :- eliminate_operator(P, F), !.
eliminate_operator(P,P).


%eliminate_neg: let the appearance of negation only in front of atom
eliminate_neg(P, P) :- prop(P),!.
eliminate_neg(-P, -P) :- prop(P), !.
eliminate_neg(P&Q, L1&L2) :- eliminate_neg(P, L1),!,eliminate_neg(Q,L2),!.
eliminate_neg(P\/Q, L1\/L2) :- eliminate_neg(P, L1),!,eliminate_neg(Q,L2),!.
eliminate_neg(@P,@F) :- eliminate_neg(P, F), !.
eliminate_neg($P,$F) :- eliminate_neg(P, F), !.
eliminate_neg(-(@P), $L) :- eliminate_neg(-P,L),!.
eliminate_neg(-($P), @L) :- eliminate_neg(-P,L),!.
eliminate_neg(-P, L) :- eliminate_neg(P, L1), neg(L1, L), !.
eliminate_neg(P,P).

%compute the negation of formula
neg(P, -P) :- prop(P),!.
neg(-P, P) :- prop(P), !.
neg(P\/Q, L1&L2) :- neg(P, L1), neg(Q, L2), !.
neg(P&Q, L1\/L2) :- neg(P, L1), neg(Q, L2), !.
neg(P->Q, L) :- neg(-P\/Q,L), !.
neg(P<->Q, L) :- neg((-P\/Q)&(-Q\/P), L), !.
neg($P, @L) :- neg(P, L),!.
neg(@P, $L) :- neg(P, L),!.
neg(-P,P).
neg(@P, $L) :- neg(P, L),!.
neg(-P,P).

% R1-R4L-distributionM-distribution, distributive law, commutative law
% and S5(4)-S5(7),
%r1to4: simplify formula use rules 1 to 4
r1to4(@($P), $F) :- r1to4(P, F), !.
r1to4($(@P), @F) :- r1to4(P, F), !.
r1to4($($P), $F) :- r1to4(P, F), !.
r1to4(@(@P), @F) :- r1to4(P, F), !.
r1to4(P\/Q, (F1)\/(F2)) :- r1to4(P, F1), r1to4(Q, F2), !.
r1to4(P&Q, (F1)&(F2)) :- r1to4(P, F1), r1to4(Q, F2), !.
r1to4(@P, @F) :- r1to4(P, F), !.
r1to4($P, $F) :- r1to4(P, F), !.
r1to4(P, P).

ldistribution(@(P&Q), (@P)&(@Q)).
mdistribution($(P\/Q), ($P)\/($Q)).

disabsob(@(P\/(@Q)), (@F1)\/(@F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob(@(P\/($Q)), (@F1)\/($F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob($(P&($Q)), ($F1)&($F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob($(P&(@Q)), ($F1)&(@F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob(@(P&Q), (@F1)&(@F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob(P\/Q, (F1)\/(F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob(P&Q, (F1)&(F2)) :- disabsob(P, F1), disabsob(Q, F2), !.
disabsob(P, P).

comm(P\/Q, X) :- X = (Q\/P).
distrib(P&(Q\/R), (P&Q)\/(P&R)).
distrib(P\/(Q&R), (P\/Q)&(P\/R)).


%simplify formula whose degree large than 1 to another whose degree is 1
hsimply_degree(F, NF) :- degree(F, X), X > 2,  F =.. [T|Args],
	hsimply_list(Args, NewArgs),
	NF1 =.. [T|NewArgs],
	degree(NF1, X1),
	(X1 >= 2, hsimply_degree(NF1, NF); NF = NF1),!.
hsimply_degree(F, NF) :- degree(F, X), X == 2, simply_degree2(F, NF), !.
hsimply_degree(F, NF):- degree(F, X), X =< 1, NF = F, !.


hsimply_list([],[]) :- !.
hsimply_list([H|T], [NH|NT]) :- hsimply_degree(H, NH), hsimply_list(T, NT), !.

%simplify formula whose degree large than 1 to another whose degree is 2
simply_degree2(F, NF) :-  r1to4(F,F1), degree(F1, X1), X1 >1,simply(F1, NEWF),
	degree(NEWF, X),
	( X>1,
	simply_degree2(NEWF, NF); NF = NEWF).
simply_degree2(F, NF) :-  r1to4(F,F1), degree(F1, X1), X1 =< 1, NF = F1, !.

simply(F, NewF) :- F =.. [T|Args], ((T) == (@); (T) == ($)), do(F, NewF), !.
simply(F, NewF) :- F =.. [T|Args], simply_list(Args, Newargs), NewF =.. [T|Newargs].


simply_list([],[]).
simply_list([H|T], [Y|NewF]) :- degree(H, X),
	X == 2,
	do(H, Y),
	simply_list(T, NewF), !.
simply_list([H|T], [H|NewF]) :- degree(H, X),
	X =< 1,
	simply_list(T, NewF), !.



do(@F, NewF) :- judge_type(F, N), N == 1, ldistribution(@F, NewF), !.
do(@F, NewF) :- union_formula(F, F1, Args, M), M \= false, disabsob(@F1, NewF), !.
do(@F, NewF) :- union_dis(F,F1, L, M), distrib(F1, F2), ldistribution(@F2,NewF), !.
do($F, NewF) :- judge_type(F, N), N == 2, mdistribution($F, NewF), !.
do($F, NewF) :- bunion_formula(F, F1, Args, M), M \= false, disabsob($F1, NewF), !.
do($F, NewF) :- union_con(F,F1, L, M), distrib(F1, F2), mdistribution($F2,NewF), !.
do(F, NewF) :- F =.. [T|Args],
	simply_list(Args, NewArgs),
	NewF =.. [T|NewArgs], !.


judge_type(F, N) :- F =.. [T|Args], (T) == (&), N = 1,!.
judge_type(F, N) :- F =.. [T|Args], (T) == (\/), N =2, !.


union_formula(F, NF, [], NF) :- F =.. [T|Args], ((T) == (@) ; (T) == ($)) , NF = F, !.
union_formula(F, NF, NewArgs, M) :- F =.. [T|Args],
	(T) == (\/),
	find_list(Args, NewArgs, M), M \= false, disjunction(NewArgs,NF1),NF=((NF1)\/(M)), !.
union_formula(F, F, [F], false).

disjunction([], true).
disjunction(T, X) :- length(T,Y), Y == 1, member(X, T),!.
disjunction([X,Y|T], F) :- disjunction(T,F1),
	(F1 == true, F = (X\/Y); F = ((X\/Y)\/F1)),!.

union_dis(F, NF, List, M) :- F =.. [T|Args], (T) == (\/),
	find_modal(Args, List, M),
	disjunction(List, NF1),
	NF = ((NF1)\/M),!.
union_dis(F, NF, List, M) :- F =.. [T|Args], (T) == (&),
	find_submodal(Args,L1, M1), (M1 == false, fail;  M = F, NF = F),!.


find_modal([], true, false) :- !.
/*find_modal([H|T], [H|L], M) :- H = [T1|Args],
	(T) == (\/),
	member(X, H),
	union_dis(X, F, L, M), !.*/
find_modal([H|T], L, M) :- H =.. [T1|Args],
	find_submodal(Args, L1, M1),
	(M1 == false, member(X,T),
	 union_dis(X,F,L2,M), L = [H|L2]; M = H, L = T),!.


find_submodal([],true,false) :- !.
find_submodal([@H|T], L, F) :- L = T, F= true, !.
find_submodal([$H|T], L, F) :- L = T, F= true, !.
find_submodal([H|T], [H|T1], F) :- member(X, T), X =.. [T2|Args],
	((T2) == (@); (T2) == ($)), F = true, !.
find_submodal([H|T], [H|T1], F) :- member(X, T), X =.. [T2|Args],
	find_submodal(Args, T1, F), !.

find_list([], true, false).
find_list([@H|T], T, M) :- M = @H, !.
find_list([$H|T], T, M) :- M = $H, !.
find_list([H|T], NL, M) :- union_formula(H, NF, H1, M1),
	(M1 == false, (find_list(T, T1, M2), M=M2, append(H1,T1,NL)); M=M1, append(H1,T, NL)),!.


bunion_formula(F, NF, [], NF) :- F =.. [T|Args], ((T) == ($); (T) == (@)) , NF = F, !.
bunion_formula(F, NF, NewArgs, M) :- F =.. [T|Args],
	(T) == (&),
	bfind_list(Args, NewArgs, M), M \= false, conjunction(NewArgs,NF1),NF=((NF1)&(M)), !.
bunion_formula(F, F, [F], false).

bfind_list([], true, false).
bfind_list([@H|T], T, M) :- M = @H, !.
bfind_list([$H|T], T, M) :- M = $H, !.
%bfind_list([H|T], [H|T1], M) :- member(X,T), bunion_formula(X, NF, T1, M),!.
bfind_list([H|T], NL, M) :- bunion_formula(H, NF, H1, M1),
	(M1 == false, (bfind_list(T, T1, M2), M=M2, append(H1,T1,NL)); M=M1, append(H1,T, NL)),!.

conjunction([], true).
conjunction(T, X) :- length(T,Y), Y == 1, member(X, T),!.
conjunction([X,Y|T], F) :- disjunction(T,F1),
	(F1 == true, F = (X&Y); F = ((X&Y)&F1)),!.


union_con(F, NF, List, M) :- F =.. [T|Args], (T) == (&),
	bfind_modal(Args, List, M),
	conjunction(List, NF1),
	NF = ((NF1)&M),!.
union_con(F, NF, List, M) :- F =.. [T|Args], (T) == (\/),
	bfind_submodal(Args,L1, M1), (M1 == false, fail;  M = F, NF = F),!.


bfind_modal([], true, false) :- !.
/*find_modal([H|T], [H|L], M) :- H = [T1|Args],
	(T) == (\/),
	member(X, H),
	union_dis(X, F, L, M), !.*/
bfind_modal([H|T], L, M) :- H =.. [T1|Args],
	find_submodal(Args, L1, M1),
	(M1 == false, member(X,T),
	 union_con(X,F,L2,M), L = [H|L2]; M = H, L = T),!.

bfind_submodal([],true,false) :- !.
bfind_submodal([@H|T], L, F) :- L = T, F= true, !.
bfind_submodal([$H|T], L, F) :- L = T, F= true, !.
bfind_submodal([H|T], [H|T1], F) :- member(X, T), X =.. [T2|Args],
	((T2) == (@); (T2) == ($)), F = true, !.
bfind_submodal([H|T], [H|T1], F) :- member(X, T), X =.. [T2|Args],
	bfind_submodal(Args, T1, F), !.

%----------------end first_degree -----------------------



%-------------new_substitute_allk-----new_substitute_allb----------------
%Substitute terms @P in Formula with new atoms
new_substitute_allk(Term, NewTerm, N1) :- getTermk(Term1,N1,L),
	substitute(@P,Term, Term1, F, Flag),
	assert(pair(@P,Term1)),
	assert(prop(Term1)),
	(Flag==true, N2 is N1+1; N2 is N1),
	(F == Term, NewTerm = F; new_substitute_allk(F, NewTerm, N2)).

new_substitute_allb(Term, NewTerm, N1) :- getTermb(Term1,N1,L),
	substitute($P,Term, Term1, F, Flag),
	assert(pair($P,Term1)),
	assert(prop(Term1)),
	(Flag==true, N2 is N1+1; N2 is N1),
	(F == Term, NewTerm = F; new_substitute_allb(F, NewTerm, N2)).

% gain a new atom from klist(which is generate dynamic) to substitute
% term @P
getTermk(Term1, N,L) :- %L= [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11],
	klist(L), length(L, N1), (N1 < N, Term1='a'; th_list(L,N,Term1)).

%gain a new atom from blist to substitute term $P
getTermb(Term1, N,L) :- %L=[b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11],
	blist(L), length(L, N1), (N1 < N, Term1='a'; th_list(L,N,Term1)).

%substitute @P or $P in formula with new atom
substitute(Term, Term, Terml, Terml, Flag) :- Flag=true,
	%write("Flag''="),write(Flag), nl,
	!.
substitute( _, Term, _, Term, Flag) :- atomic( Term), Flag=false, !.
substitute( Sub, Term, Subl, Terml, Flag) :- Term =.. [F | Args], substlist(Sub, Args, Subl, Argsl, Flag), Terml =.. [F | Argsl].

substlist(_, [],_,[], Flag) :- Flag=false.
substlist(Sub, [Term | Terms], Subl, [Terml | Termsl], Flag):- substitute(Sub, Term, Subl, Terml, Flag1),
		substlist( Sub, Terms, Subl, Termsl, Flag2),
		(Flag1=true, Flag=Flag1; Flag=Flag2).

list_th([H|T],X, N,H) :- N == X,!.
list_th([H|T],X, N,Y) :- M is N+1, list_th(T,X, M,Y).
list_th(F,X, N,Y):- X =< 0, write("error: N<=0").
%list_th(F,X, N,Y):- length(F,L), X > L+1, write("error: N > length of list"), nl, write(N), nl.

th_list(F,X,Y):- list_th(F,X,1,Y).

% ----------end new_substitute_allk-----new_substitute_allb-------------


%---------------new_wff2dnf-------new_wff2cnf----------------------------------
%excludeing the repetitive term (or clause)

%the union, which have no repetitive element, of two sets 
nonRep(L1, L2, L) :- findall(X, (member(X, L1), member(X, L2)), L11),
	difference_list(L11, L2, L22),
	difference_list(L11, L1, L33),
	append(L11, L22, L3),
	append(L3, L33, L).
	
/* Wff2dnf transforms a wff to its disjunctive normal form in list format */
new_wff2dnf(P,[[P]]) :- prop(P), !.
new_wff2dnf(-P,[[-P]]) :- prop(P), !.
new_wff2dnf(P->Q, L) :- new_wff2dnf(-P\/Q, L).
new_wff2dnf(P<->Q, L) :- new_wff2dnf((P->Q)&(Q->P), L).
new_wff2dnf(P\/Q, L) :- new_wff2dnf(P,L1), new_wff2dnf(Q,L2), nonRep(L1,L2,L).
new_wff2dnf(P&Q, L) :- new_wff2dnf(P,L1), new_wff2dnf(Q,L2),
    findall(X, (member(Y,L1), member(Z,L2), append(Y,Z,X)), L).
new_wff2dnf(-P, L) :- new_wff2cnf(P,L1), negate(L1,L).

/* negate(L1,L): negate L1 in DNF (CNF) and make it into a CNF (DNF). */
negate([],[]) :- !.
negate([[]],[[]]) :- !.
negate(L1,L) :- bagof(X, Y^(member(Y, L1), negate_one_clause(Y,X)), L).

/* negate_one_clause(L1, L2): make all elements in L1 into its complement */
negate_one_clause(L1, L2) :- bagof(X, Y^(member(Y,L1), complement(Y,X)), L2).

complement(true,false) :- !.
complement(false,true) :- !.
complement(P,-P) :- prop(P).
complement(-P,P) :- prop(P).

/* Wff2cnf transforms a wff to its conjunctive normal form in list format */
new_wff2cnf(P, [[P]]) :- prop(P), !.
new_wff2cnf(-P, [[-P]]) :- prop(P), !.
new_wff2cnf(P->Q, L) :- new_wff2cnf(-P\/Q, L), !.
new_wff2cnf(P<->Q, L) :- new_wff2cnf((-P\/Q)&(-Q\/P), L), !.
new_wff2cnf(P&Q, L) :- new_wff2cnf(P,L1), new_wff2cnf(Q,L2), nonRep(L1,L2,L), !.
new_wff2cnf(P\/Q, L) :- new_wff2cnf(P,L1), new_wff2cnf(Q,L2),
    findall(X, (member(Y,L1), member(Z,L2), append(Y,Z,X)), L), !.
new_wff2cnf(-P, L) :- new_wff2dnf(P, L1), negate(L1, L), !.

%------------end new_wff2dnf-------new_wff2cnf--------------------------------


%------------------kforget----------------------------------
%compute the result of forget some atoms(P) from formula with MDNF
kforget(A,[],A) :- !.
kforget([],_,true) :- !.
kforget([H|T], P, F) :- kforget_one(H,P,F1), kforget(T,P,F2), (F2 = true, F = F1; F = ((F1)\/(F2))).

%kforget_one: compute the result of forget some atoms from a Term
kforget_one(H, P, L) :- split_k(H, L1, L2, L3),
		length(L1,LEN1),
		LEN1 > 0,
		length(L3,LEN2),
		LEN2 > 0,
		kformula(L1, KF),
		kunion(KF, UL),
		@X = UL,
		dnf2wff([L3],F),
		% KF1 = ((X)&(F)),
		new_wff2cnf(F,CNF1),
		forget(CNF1,P, F1),
		new_wff2cnf(X,CNF2),
		forget(CNF2,P,F2),
		kforgetb(L2, UL, P, F3),
		forget_result_simplify(F1, F11),
		forget_result_simplify(F2, F22),
		cnf2wff(F11, NEWCNF1),
		cnf2wff(F22, NEWCNF2),
		L = ((NEWCNF1)&(@(NEWCNF2))&(F3)).

kforget_one(H, P, L) :- split_k(H, L1, L2, L3),
		length(L3,LEN2),
		LEN2 > 0,
		dnf2wff([L3],F),
		new_wff2cnf(F,CNF1),
		forget(CNF1,P, F1),
		kforgetb(L2, nil, P, F3),
		forget_result_simplify(F1, F11),
		cnf2wff(F11, NEWCNF1),
		L = ((NEWCNF1)&(F3)).

kforget_one(H, P, L) :- split_k(H, L1, L2, L3),
		length(L1,LEN1),
		LEN1 > 0,
		kformula(L1, KF),
		kunion(KF, UL),
		@X = UL,
		new_wff2cnf(X,CNF2),
		forget(CNF2,P,F2),
		kforgetb(L2, UL, P, F3),
		forget_result_simplify(F2, F22),
		cnf2wff(F22, NEWCNF2),
		L = ((@(NEWCNF2))&(F3)).

kforget_one(H, P, L) :- split_k(H, L1, L2, L3),
		kforgetb(L2, UL, P, F3),
		L = F3.

%------------------end kforget------------------------------


%---------------------split_k-------------------------------
% split the list of term to three classes, L1 is the list of term @P, L2
% is the list of term of $P and L3 is the list of objective term
split_k([],[],[],[]).
split_k([H|T], [H|L1], L2, L3) :- getTermk(_, 1, L), member(H, L), split_k(T, L1, L2, L3), !.
split_k([H|T], L1, [H|L2], L3) :- getTermb(_, 1, L), member(H, L), split_k(T, L1, L2, L3), !.
split_k([H|T], L1, L2, [H|L3]) :- split_k(T, L1, L2, L3).

%---------------------end split_k----------------------------------



% ------------------------dnf2wff---cnf2wff--------------------------
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
list2disj([P | L],P \/ W) :- list2disj(L,W).

% ----------------------end dnf2wff---cnf2wff-------------------------
%
%
%--------------------kformula----------------------------
%gain the terms that had been substituted by new aotms
kformula([], []).
kformula([H|T], [X|Y]) :- pair(X, H),kformula(T, Y).


%-----------------end kformula---------------------------
%
%--------------------kunion--------------------------
%just like cnf2wff
kunion([], true).
kunion([@H|T], L) :- length(T,X), X==0, !, L = @H.
kunion([H,Y|T], L) :-  elimi_modal(L1, (H)&(Y)),
		 %  elimi_modal(L1, (H)&(Y)),
		 kunion(T, L2),(L2=true,L=L1; elimi_modal(L,(L1)&(L2))).

%elimi_modal: eliminate the nested modal operators
elimi_modal([],[]).
elimi_modal(F,F) :- isformula(F),!.
elimi_modal(@(P&Q), (@X1)&(@X2)) :- elimi_modal(P,X1), elimi_modal(Q,X2),!.
elimi_modal($(P\/Q), ($X1)\/($X2)) :- elimi_modal(P,X1), elimi_modal(Q,X2),!.
elimi_modal(@(P\/ @Q), (@X1)\/(@X2)) :- elimi_modal(P,X1), elimi_modal(Q,X2),!.
elimi_modal(@(P\/ $Q), (@X1)\/($X2)) :- elimi_modal(P,X1), elimi_modal(Q,X2),!.
elimi_modal($(P& @Q), ($X1)&(@X2)) :- elimi_modal(P,X1), elimi_modal(Q,X2),!.
elimi_modal($(P& $Q),($X1)&($X2)) :- elimi_modal(P,X1), elimi_modal(Q,X2),!.
elimi_modal(A,A).

%judge the expression is a objective formula
isformula(true).
isformula(false).
isformula(P) :- prop(P), !.
isformula(-P) :- prop(P), !.
isformula(P\/Q) :- isformula(P), isformula(Q), !.
isformula(P&Q) :-  isformula(P), isformula(Q), !.
isformula(P->Q) :-  isformula(P), isformula(Q), !.
isformula(P<->Q) :-  isformula(P), isformula(Q), !.
isformula(-P) :- isformula(P), !.

%---------------end kunion--------------------------


%-----------forget--------------------------------------
%forget some atoms from CNF objective formula
forget(A, [], A) :- !.
forget(A, PL, NewA) :- select(A, PL, P), difference(PL, P, NewPL),
    forget_one(A, P, NewA1), !, forget(NewA1, NewPL, NewA).


/* select first those that are already entailed */
select(CNF, PL, P) :- member(P,PL), member([P],CNF), !.
select(CNF, PL, P) :- member(P,PL), member([-P],CNF), !.

/* select in sequence given afterward */

select(_, [P|_], P) :- !.

/* select(A, PL, P): if there is no other Q in PL such that forget(A; Q)
   is a shorter formula than forget(A; P). */

select(_, [P], P) :- !.
select(A, PL, P) :- setof(X, Q^NewA^N^(member(Q, PL), basic_forget(A, Q, NewA),
                                     size(NewA, N), X=[N,Q]), [[_,P]|_]).

/* difference(L,P,L1): L1 is L - [P] */
difference([P | L], P, L) :- !.
difference([Q | L], P, [Q | NewL]) :- difference(L,P,NewL).
difference([],P,[]).


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

/* size(A,N): N is the size of the formula A */

size(P,1) :- prop(P).
size(-A, N) :- size(A, M), N is M+1.
size(A\/B, N) :- size(A, N1), size(B, N2), N is N1+N2+1.
size(A&B, N)  :- size(A, N1), size(B, N2), N is N1+N2+1.
size(A->B, N) :- size(A, N1), size(B, N2), N is N1+N2+1.
size(A<->B, N) :- size(A, N1), size(B, N2), N is N1+N2+1.

/*split(CNF, P, NewCNFCNF1, CNF2, CNF3): split the clause of CNF into three classes, CNF1 indicats the conjunction of clauses which contain p, CNF2 indicates the conjunction of clauses which contain -p, and CNF3 indicates the conjunction of other clauses.*/
split([], _, [], [], []).
split([Clause|CNF], P, [Clause | CNF1], CNF2, CNF3) :- member(P, Clause), !,split(CNF, P, CNF1, CNF2, CNF3).
split([Clause|CNF], P, CNF1, [Clause | CNF2], CNF3) :- member(-P, Clause), !,split(CNF, P, CNF1, CNF2, CNF3).
split([Clause|CNF], P, CNF1, CNF2, [Clause | CNF3]) :- split(CNF, P, CNF1, CNF2, CNF3).



resolve([],_,[]).
resolve([C|CNF], P, NewCNF) :- member(P,C), !,
    resolve(CNF,P,NewCNF).
resolve([C|CNF], P, [NewC|NewCNF]) :- complement(P,P1),
    difference(C,P1,NewC),
    resolve(CNF,P,NewCNF).


/* cnf_simplifies(L,NewL): simplifies DNF L into L1 by eliminating repetitive
   literals and contradictive literals. It also eliminates clausforget_onees that
   subsumes others.*/

cnf_simplifies(CNF, NewCNF) :- cnf_sort_simpl(CNF, CNF1),
    unit_elim(CNF1,CNF2),
    cnf_sort_simpl(CNF2, CNF3),
    check_contr(CNF3,CNF4),
    cnf_minimize(CNF4, NewCNF).


/* sort each clause and get rid of tautologies*/

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


/* eliminates all unit clauses*/

unit_elim(CNF,[[P]|NewCNF]) :-
    member([P],CNF),
    resolve(CNF,P,CNF2), !,
    unit_elim(CNF2,NewCNF), !.
unit_elim(CNF,CNF) :- !.

/* check for contradiction */

check_contr(CNF,[[]]) :- member([],CNF), !.
check_contr(CNF,CNF).


/* get rid of clauses that are subsumed by others.  */

cnf_minimize(DNF,DNF1) :- member(L,DNF), member(L1,DNF), L\==L1,
    cnf_subsumes(L,L1), difference(DNF, L1, DNF2), cnf_minimize(DNF2,DNF1), !.
cnf_minimize(DNF,DNF).

/* cnf_subsumes(L,L1) if L subsumes L1: L is a subset of L1.
 */

cnf_subsumes(L,L1) :- \+ (member(X,L), X\==false, X\==(-true),
    \+ member(X,L1)).


/* cnf_disjunct(CNF1,CNF2,CNF3): if CNF3 is the CNF of CNF1 v CNF2  */

cnf_disjunct([],CNF2,[]) :- !.
cnf_disjunct(CNF1,[],[]) :- !.
cnf_disjunct([[]],CNF2,CNF2) :- !.
cnf_disjunct(CNF1,[[]],CNF1) :- !.
cnf_disjunct(CNF1,CNF2,CNF3) :-
	setof(X, Y^Z^(member(Y,CNF1), member(Z,CNF2), union(Y,Z,X)), CNF3).

%--------------------end forget--------------------------------

%----------------------------kforgetb------------------
%compute the forget of term $P
kforgetb([],_,_,true) :- !.
kforgetb(L2, UL, P, F3) :-  bformula(L2, BF), kforgetblist(BF, UL, P, F3).


%get the term of term $P
bformula([], []).
bformula([H|T], [X|Y]) :- pair(X, H),bformula(T, Y).


kforgetblist([],_, P, true).
kforgetblist([$H|T], F, P, WFF) :-  (F = nil, C = H; @X = F, C = ((X)&(H))),
	new_wff2cnf(C, CNF1),
	forget(CNF1, P, NEWCNF1),
	forget_result_simplify(NEWCNF1, NEWCNF11),
	cnf2wff(NEWCNF11, WFF1),
	kforgetblist(T, F, P, WFF2),
	%(WFF2 == true, WFF = ($(WFF1));WFF = (($(WFF1))&(WFF2))).
	eliminate_true((($(WFF1))&(WFF2)), WFF).

%simplify the result of forget
forget_result_simplify(CNF, NEWCNF) :- cnf_simplifies(CNF, NEWCNF), !.

/* init_minimize(CNF4,CNF6,NewCNF):eliminate clauses in CNF4 that are redundant in
   the sense that it is subsumed by the rest of the clauses in CNF4 and CNF6 */

init_minimize(CNF1,CNF2,NewCNF) :- member(C,CNF1),
    difference(CNF1,C,CNF3), union(CNF3,CNF2,CNF4),
    cnf_simplifies(CNF4,CNF5), member(C1,CNF5), cnf_subsumes(C1,C), !,
    init_minimize(CNF3,CNF2,NewCNF).
init_minimize(CNF,_,CNF).


%eliminate the true from formula
eliminate_true(X, NF) :- find_true(X, B), B == false, tra_simply(X, NF), !.
eliminate_true(X, NF) :- find_true(X, B), B == true,
	tra_simply(X, NF1),
	find_true(NF1, B1),
	(B1 == true, eliminate_true(NF1,NF); NF = NF1),!.

%find_true(X, false) :- prop(X).
%find_true(-P, false) :- prop(P).
find_true(true, true).
find_true(X, false) :- prop(X).
find_true(-P, false) :- prop(P).
find_true(X, B) :- X =.. [T|Args], find_true_list(Args, B),!.

find_true_list([],false).
find_true_list([H|T], B) :- find_true(H, B1), find_true_list(T, B2), (B1 == false, B= B2; B = B1), !.

tra_simply(F, NEW) :- terminal_simply(F, X),
	(X == F, NEW = F; tra_simply(X, NEW)).


%Simplify formula: eliminate false, true and formula be entailed.
terminal_simply(false, false) :- !.
terminal_simply(true, true) :- !.

terminal_simply(X, X) :- prop(X).
terminal_simply(-P, -P) :- prop(P).
terminal_simply(-P, F1) :- terminal_simply(P, F1), !.


terminal_simply(P&(@P), F) :- terminal_simply(@P, F), !. %  p&@p=@p
terminal_simply((@P)& P, F) :- terminal_simply(@P, F), !. %  p&@p=@p
terminal_simply(P& ($P), F) :- terminal_simply(P, F), !. %  p&$p=p
terminal_simply(($P) & P, F) :- terminal_simply(P, F), !. %  p&$p=p
terminal_simply((@P)& ($P), F) :- terminal_simply(@P, F), !.
terminal_simply(($P)& (@P), F) :- terminal_simply(@P, F), !.
terminal_simply(P & P, F) :- terminal_simply(P, F), !.
terminal_simply(P \/ P, F) :- terminal_simply(P, F), !.


terminal_simply(P\/true, true) :- !.
terminal_simply(true\/P, true) :- !.

terminal_simply(P&true, F1) :- terminal_simply(P, F1), !.
terminal_simply(true&P, F1) :- terminal_simply(P, F1), !.

terminal_simply(P\/Q, (F1)\/(F2)) :- terminal_simply(P, F1),
	terminal_simply(Q, F2), !.
terminal_simply(P&Q, (F1)&(F2)) :- terminal_simply(P, F1),
	terminal_simply(Q, F2), !.

terminal_simply(@true, true).
terminal_simply($true, true).
terminal_simply(@P, @F1) :- terminal_simply(P, F1),!.
terminal_simply($P, $F1) :- terminal_simply(P, F1),!.
terminal_simply(P, P).

%--------------------end kforgetb----------------------------


%-------------------------eliminate_false-------------------
%eliminate the false from formula
eliminate_false(X, NF) :- find_false(X, B), B == false, NF = X, !.
eliminate_false(X, NF) :- find_false(X, B), B == true,
	fterminal_simply(X, NF1),
	find_false(NF1, B1),
	%B1 == true, eliminate_true(NF1,NF),!.
	(B1 == true, eliminate_false(NF1,NF); NF = NF1),!.


%find_false(X, false) :- prop(X).
%find_false(-P, false) :- prop(P).
find_false(false, true).
find_false(X, false) :- prop(X).
find_false(-P, false) :- prop(P).
find_false(X, B) :- X =.. [T|Args], find_false_list(Args, B),!.

find_false_list([],false).
find_false_list([H|T], B) :- find_false(H, B1), find_false_list(T, B2), (B1 == false, B= B2; B = B1), !.


%Simplify formula: eliminate false, true and formula be entailed.
fterminal_simply(false, false) :- !.
fterminal_simply(true, true) :- !.
fterminal_simply(X, X) :- prop(X).
fterminal_simply(-P, -P) :- prop(P).
fterminal_simply(-P, F1) :- fterminal_simply(P, F1), !.
fterminal_simply(P&false, false) :- !.
fterminal_simply(false&P, false) :- !.
fterminal_simply(P\/false, F1) :- fterminal_simply(P, F1), !.
fterminal_simply(false\/P, F1) :- fterminal_simply(P, F1), !.

fterminal_simply(P\/Q, (F1)\/(F2)) :- fterminal_simply(P, F1),
	fterminal_simply(Q, F2), !.
fterminal_simply(P&Q, (F1)&(F2)) :- fterminal_simply(P, F1),
	fterminal_simply(Q, F2), !.
fterminal_simply(@false, false).
fterminal_simply($false, false).
fterminal_simply(@P, @F1) :- fterminal_simply(P, F1),!.
fterminal_simply($P, $F1) :- fterminal_simply(P, F1),!.

fterminal_simply(P&(@P), F) :- fterminal_simply(@P, F), !. %  p&@p=@p
fterminal_simply(P& ($P), F) :- fterminal_simply($p, F), !. %  p&$p=p

fterminal_simply(P,P).
%---------------------end eliminate_false------------------------------
