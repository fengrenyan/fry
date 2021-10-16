
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


%ctlforget(F, V, R) :- 

