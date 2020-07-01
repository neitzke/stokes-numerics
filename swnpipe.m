(* Compute one spectral network and write to JSON file *)
(* Adapted from Neitzke's swn-plotter.nb *)

alloworiented=False

suppresswarnings:=Module[{},
Off[InterpolatingFunction::dmval];
Off[NDSolve::ndsz];
Off[NDSolve::ndtol];
Off[NDSolve::mxst];
];


suppresswarnings;


F[x_,z_]:=x^KK+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 2\), \(KK\)]\(\(P[k]\)[z]\ 
\*SuperscriptBox[\(x\), \(KK - k\)]\)\);


roots[z_]:=NSolve[F[x,z]==0,x];


root[z_,i_]:=x/. roots[z][[i]];


dxdz[x_,z_]:=-(
\!\(\*SuperscriptBox[\(F\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,z]/
\!\(\*SuperscriptBox[\(F\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,z]);


singleray[z0_,x10_,x20_,L_,\[Alpha]_,startt_]:={x1[t],x2[t],z[t]}/. 
NDSolve[{Derivative[1][x1][t]==dxdz[x1[t],z[t]] Derivative[1][z][t],Derivative[1][x2][t]==dxdz[x2[t],z[t]] Derivative[1][z][t],(x1[t]-x2[t])  Derivative[1][z][t]==E^(I \[Alpha]),x1[startt]==x10,x2[startt]==x20,z[startt]==z0},{x1,x2,z},{t,startt,L},MaxSteps->100000,PrecisionGoal->\[Infinity]][[1]]


pairbasis:=Subsets[Table[i,{i,1,KK}],{2}];


triplebasis:=Subsets[Table[i,{i,1,KK}],{3}];


tprays[tpn_,L_,\[Alpha]_]:=Module[{\[Delta],i,j,k,z0,tripleindex,triples,pairindex,pairs,xr,phaseshift,cubiconly, c,\[Beta],xdiff,b,tpsingular,outrays,zi,x1,x2,Minit},
z0=turningpoints[[tpn]];

(* if the turning point is also a singularity, we need a different scheme *)
tpsingular=False;
For[i=1,i<=Length[sing],i++,
If[sing[[i]]==z0,tpsingular=True];
];

(* special case for cubic-differential-only *)
If [P[2][z]===0,cubiconly=True,cubiconly=False];

If[tpsingular==False && cubiconly == False,
(* turning point is not a singular point *)
(* find the roots at a point very near z0 *)
\[Delta]=0.00005;
xr=x/. roots[z0+\[Delta]^2];

(* find the pair of roots that are closest together *)
pairs=pairbasis;
xdiff=Table[xr[[pairs[[i,1]]]]-xr[[pairs[[i,2]]]],{i,1,Length[pairs]}];
pairindex=Ordering[Abs[xdiff],1][[1]];

(* find initial data for the 3 rays emerging *)
\[Beta]=xdiff[[pairindex]]/(2 \[Delta]);
i=pairs[[pairindex]][[1]];
j=pairs[[pairindex]][[2]];
b=1/2 (xr[[i]]+xr[[j]]);
phaseshift[l_]:=Exp[I(1/3 (\[Alpha]-Arg[\[Beta]])+2/3 \[Pi] (l-1))];
outrays = Table[singleray[z0+\[Delta]^2 phaseshift[l]^2,b+\[Beta] \[Delta] phaseshift[l],b-\[Beta] \[Delta] phaseshift[l],L,\[Alpha],0],{l,1,3}]
];

(* if turning point is a singular point, then we should deal with it according to the singulartype[z0] given by user *)
If[tpsingular==True && singulartype[z0]==1,
(* turning point is a singular point of simple type: 2 sheets colliding, square-root behavior *)
\[Delta]=1/100;
xr=x/.  roots[z0+\[Delta]^2];

(* find the pair of roots that are closest together *)
pairs=pairbasis;
xdiff=Table[xr[[pairs[[i,1]]]]-xr[[pairs[[i,2]]]],{i,1,Length[pairs]}];
pairindex=Ordering[Abs[xdiff],1][[1]];

(* find initial data for the ray emerging *)
\[Beta]=\[Delta] xdiff[[pairindex]];
i=pairs[[pairindex]][[1]];
j=pairs[[pairindex]][[2]];
b=1/2 (xr[[i]]+xr[[j]]);

phaseshift:=Exp[I (\[Alpha]-Arg[\[Beta]])]; 
zi=z0+\[Delta]^2 phaseshift^2;

x1=x /. FindRoot[F[x,zi],{x,b+\[Beta]/(2 \[Delta] phaseshift)}][[1]];
x2=x /. FindRoot[F[x,zi],{x,b-\[Beta]/(2 \[Delta] phaseshift)}][[1]];

outrays={singleray[zi,x1,x2,L,\[Alpha],0]}
];

If[tpsingular==True && singulartype[z0]==2,
(* turning point is a singular point with 3 sheets colliding, monodromy (123), cube-root behavior -- eigenvalues like delta^(-2/3) *)
\[Delta]=1/1000;
xr=x/.  roots[z0+\[Delta]];
c=xr \[Delta]^(2/3);

outrays={};

i=1;j=2;
\[Beta]=c[[i]]-c[[j]];

phaseshift:=Exp[I (\[Alpha]-Arg[\[Beta]])]; 
zi=z0+\[Delta] phaseshift^3;

x1=x /. FindRoot[F[x,zi],{x,c[[i]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];
x2=x /. FindRoot[F[x,zi],{x,c[[j]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];

(* mass correction for the omitted segment of path *)
Minit = Abs[3 (x1-x2) \[Delta]];

outrays=Append[outrays,singleray[zi,x1,x2,L,\[Alpha],Minit]];

i=2;j=1;
\[Beta]=c[[i]]-c[[j]];

phaseshift:=Exp[I (\[Alpha]-Arg[\[Beta]])]; 
zi=z0+\[Delta] phaseshift^3;

x1=x /. FindRoot[F[x,zi],{x,c[[i]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];
x2=x /. FindRoot[F[x,zi],{x,c[[j]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];

(* mass correction for the omitted segment of path *)
Minit = Abs[3 (x1-x2) \[Delta]];

outrays=Append[outrays,singleray[zi,x1,x2,L,\[Alpha],Minit]];

];

If[tpsingular==True && singulartype[z0]==3,
(* turning point is a singular point with 4 sheets colliding, monodromy (123), cube-root behavior -- eigenvalues like delta^(-3/4) *)
\[Delta]=1/1000;
xr=x/.  roots[z0+\[Delta]];
c=xr \[Delta]^(3/4);

outrays={};

i=1;j=2;
\[Beta]=c[[i]]-c[[j]];

phaseshift:=Exp[I (\[Alpha]-Arg[\[Beta]])]; 
zi=z0+\[Delta] phaseshift^3;

x1=x /. FindRoot[F[x,zi],{x,c[[i]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];
x2=x /. FindRoot[F[x,zi],{x,c[[j]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];

(* mass correction for the omitted segment of path *)
Minit = Abs[3 (x1-x2) \[Delta]];

outrays=Append[outrays,singleray[zi,x1,x2,L,\[Alpha],Minit]];

i=2;j=1;
\[Beta]=c[[i]]-c[[j]];

phaseshift:=Exp[I (\[Alpha]-Arg[\[Beta]])]; 
zi=z0+\[Delta] phaseshift^3;

x1=x /. FindRoot[F[x,zi],{x,c[[i]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];
x2=x /. FindRoot[F[x,zi],{x,c[[j]]/(\[Delta]^(2/3) phaseshift^2)}][[1]];

(* mass correction for the omitted segment of path *)
Minit = Abs[3 (x1-x2) \[Delta]];

outrays=Append[outrays,singleray[zi,x1,x2,L,\[Alpha],Minit]];

];

If[cubiconly==True,
(* turning point is a zero with monodromy like (123) *)
(* find the roots at a point very near z0 *)
\[Delta]=0.05;
xr=x/. roots[z0+\[Delta]^3];

(* find the triple of roots that are closest together *)
triples=triplebasis;
xdiff=Table[Abs[xr[[triples[[i,1]]]]-xr[[triples[[i,2]]]]]+Abs[xr[[triples[[i,2]]]]-xr[[triples[[i,3]]]]]+Abs[xr[[triples[[i,1]]]]-xr[[triples[[i,3]]]]],{i,1,Length[triples]}];
tripleindex=Ordering[xdiff,1][[1]];

(* find initial data for the 8 rays emerging *)
i=triples[[tripleindex]][[1]];
j=triples[[tripleindex]][[2]];
k=triples[[tripleindex]][[3]];
b=1/3 (xr[[i]]+xr[[j]]+xr[[k]]);
c=(xr[[i]]-b)/\[Delta];

phaseshift[1]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[1-Exp[2 Pi I/3]]))];
phaseshift[2]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[1-Exp[4 Pi I/3]]))];
phaseshift[3]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[Exp[2 Pi I/3]-Exp[4 Pi I/3]]))];
phaseshift[4]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[Exp[2 Pi I/3]-1]))];
phaseshift[5]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[Exp[4 Pi I/3]-1]-2 Pi))];
phaseshift[6]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[Exp[4 Pi I/3]-Exp[2 Pi I/3]]-2 Pi))];
phaseshift[7]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[1-Exp[2 Pi I/3]]-2 Pi))];
phaseshift[8]=Exp[I(1/4 (\[Alpha]-Arg[c]-Arg[1-Exp[4 Pi I/3]]-2 Pi))];
outray[1] = singleray[z0+\[Delta]^3 phaseshift[1]^3,b+c \[Delta] phaseshift[1],b+c \[Delta] Exp[2 Pi I/3] phaseshift[1],L,\[Alpha],0];
outray[2]=singleray[z0+\[Delta]^3 phaseshift[2]^3,b+c \[Delta] phaseshift[2],b+c \[Delta] Exp[4 Pi I/3] phaseshift[2],L,\[Alpha],0];
outray[3]=singleray[z0+\[Delta]^3 phaseshift[3]^3,b+c \[Delta] Exp[2 Pi I/3]phaseshift[3],b+c \[Delta] Exp[4 Pi I/3] phaseshift[3],L,\[Alpha],0];
outray[4]=singleray[z0+\[Delta]^3 phaseshift[4]^3,b+c \[Delta] Exp[2 Pi I/3]phaseshift[4],b+c \[Delta] phaseshift[4],L,\[Alpha],0];
outray[5]=singleray[z0+\[Delta]^3 phaseshift[5]^3,b+c \[Delta] Exp[4 Pi I/3]phaseshift[5],b+c \[Delta] phaseshift[5],L,\[Alpha],0];
outray[6]=singleray[z0+\[Delta]^3 phaseshift[6]^3,b+c \[Delta] Exp[4 Pi I/3]phaseshift[6],b+c \[Delta] Exp[2 Pi I/3] phaseshift[6],L,\[Alpha],0];
outray[7] = singleray[z0+\[Delta]^3 phaseshift[7]^3,b+c \[Delta] phaseshift[7],b+c \[Delta] Exp[2 Pi I/3] phaseshift[7],L,\[Alpha],0];
outray[8]=singleray[z0+\[Delta]^3 phaseshift[8]^3,b+c \[Delta] phaseshift[8],b+c \[Delta] Exp[4 Pi I/3] phaseshift[8],L,\[Alpha],0];
outrays=Table[outray[l],{l,1,8}];
];

outrays
];


alltprays[L_,\[Alpha]_]:=
If[Length[turningpoints]>0,
Flatten[Table[tprays[k,L,\[Alpha]],{k,1,Length[turningpoints]}],{{1,2},{3}}],
{}
];


residues[z_]:=Module[{\[Epsilon]},\[Epsilon]=0.0000000001;\[Epsilon](x/. roots[z+\[Epsilon]])];


findtp:=Module[{solution},
\[CapitalDelta][z_]=Discriminant[x^KK+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 2\), \(KK\)]\(\(P[k]\)[z]\ 
\*SuperscriptBox[\(x\), \(KK - k\)]\)\),x] ; 
solution=NSolve[\[CapitalDelta][z]==0,z,WorkingPrecision->20];
If[Length[solution]>0, turningpoints=z/.solution , turningpoints={}];
];


closeQ[a_,b_]:=(Abs[a-b]<0.000001);


removerepeatedtp:=(turningpoints=DeleteDuplicates[turningpoints,closeQ];);


dumptable[tab_]:=Map[#[[2]] &,DownValues[tab]]


Options[mswn]:={TimeStep->0.02,Grain->20,BlockTolerance->0.04,KillTimeThreshold->0.005, DoIntersections->True};
Options[findintersectioncandidates]:=Options[mswn];
Options[reduceintersectioncandidates]:=Options[mswn];
Options[killbirthintersections]:=Options[mswn];


(* 
data format for intersections:  {{first ray #, first intersection time, intersection position},{second ray #, second intersection time, intersection position}} 
*)
(* hash[z_]=Round[z/Sqrt[Abs[z]]]; *)
hash[z_,grain_]=Round[z grain];
findintersectioncandidates[rays_,noldrays_,parents_,OptionsPattern[]]:=
Module[{i,j,k,ints,pp,ep,delta,raystart,raylen,timestep,grain,initialbox,tt,currentbox,currentstatus,currentpos,nsteps,nrays},
currentstatus="Getting started";
timestep=OptionValue[TimeStep];
grain=OptionValue[Grain];

ints={};

(* build hash-table of boxes visited by rays: "timestep" controls parameterization of rays, "grain" controls size of boxes; ht[raynum, box] = {raynum, intersection time, z-position} *)
For[i=noldrays+1,i<=Length[rays],i=i+1,
raystart = rays[[i]][[3]][[0,1,1,1]];
raylen=rays[[i]][[3]][[0,1,1,2]];
pp=rays[[i]][[3]] /. t->raystart;
initialbox = hash[pp,grain];

For[tt=raystart+timestep,tt<=raylen,tt=tt+timestep,
ep = rays[[i]][[3]] /. t->tt;

(* bail out if we are very close to a singularity *)
For[k=1,k<=Length[sing],k=k+1,
If[(sing[[k]]!=Infinity && Abs[pp-sing[[k]]]<=0.02) || (sing[[k]]==Infinity && Abs[pp]>150),
pp=ep;Continue[]];
];

(*
For[k=1,k\[LessEqual]Length[turningpoints],k=k+1,
If[Abs[pp-turningpoints[[k]]]\[LessEqual]0.02,
pp=ep;Continue[]];
];
*)

(* divide the timestep into smaller pieces so we don't jump over any boxes *)
delta = ep-pp;
nsteps = Floor[10 grain Abs[delta]];

currentpos=pp;
currentbox=hash[currentpos,grain];
currentstatus="Working at time "<>ToString[tt]<> " on timestep divided into "<>ToString[nsteps+1]<>" pieces";
For[k=0,k<= nsteps, k=k+1,

(* is this point already in hash table? if so, mark it as candidate intersection *)
(* but don't count intersections between a ray and one of its parents that occur very near birthplace *)
For[j=1,j<i,j++,
If[(Abs[currentbox-initialbox]>2 || (j!=parents[[i,1]] && j!=parents[[i,2]])) && Head[ht[j,currentbox]]!=ht,,,
ints=Append[ints,{{i,tt,currentpos},ht[j,currentbox]}]; 
];
];

(* add the current point to the hash table *)
ht[i,currentbox]={i,tt,currentpos};
ht[i,currentbox+1]={i,tt,currentpos};
ht[i,currentbox+I]={i,tt,currentpos};
ht[i,currentbox+1+I]={i,tt,currentpos};

(* get ready for the next step *)
If[nsteps>0,currentpos=currentpos+delta/nsteps;currentbox = hash[currentpos,grain]];
];

(* update starting position for the next timestep *)
pp=ep;
];
];
nrays=Length[rays];
ints=Sort[ints,(nrays #1[[1,1]]+#1[[2,1]]) > (nrays #2[[1,1]]+#2[[2,1]]) &];
ints
];


(* consolidate list of intersection candidates by combining "blocks" *)
reduceintersectioncandidates[candidatelist_,OptionsPattern[]]:=Module[{cleanedlist,blockstart,cur,blockend,blockmid,prev,i,timestep,tolerance},
timestep=OptionValue[TimeStep];
tolerance=OptionValue[BlockTolerance];
cleanedlist={};
blockstart=1;
For[i=1,i<=Length[candidatelist],i=i+1,
cur=candidatelist[[i]];
(* see if we are at the end of a block: that can happen if the ray numbers change, if the intersection time moves by more than one time step, or if the z-position moves by more than a small tolerance *)
If[((i>1) && (cur[[1,1]]!=prev[[1,1]] || cur[[2,1]]!=prev[[2,1]]||Abs[cur[[1,2]]- prev[[1,2]]]>1.1 timestep || Abs[cur[[1,3]]-prev[[1,3]]]>tolerance || Abs[cur[[2,2]]- prev[[2,2]]]>1.1 timestep || Abs[cur[[2,3]]-prev[[2,3]]]>tolerance)), 
blockend = i-1;
blockmid = Round[(blockend+blockstart)/2];
cleanedlist=Append[cleanedlist,candidatelist[[blockmid]]];
blockstart=i;
];
prev = cur;
];
blockend = Length[candidatelist];
blockmid = Round[(blockend+blockstart)/2];
If[blockend>0,cleanedlist=Append[cleanedlist,candidatelist[[blockmid]]]];
cleanedlist
];


(* refine the approximate intersection candidates to exact ones using the root-finder *)
refineintersections[rays_, candidatelist_]:=Module[{cur,raynum1,raynum2,tint1approx,tint2approx,ray1,ray2,intpt,preciseints},
(preciseints={};
For[i=1,i<=Length[candidatelist],i=i+1,
cur=candidatelist[[i]];
raynum1=cur[[1,1]];
raynum2=cur[[2,1]];
tint1approx=cur[[1,2]];
tint2approx=cur[[2,2]];
ray1=rays[[raynum1]];
ray2=rays[[raynum2]];

(* suppress error messages but check for them internally *)
Quiet[
intpt=Check[FindRoot[{Re[(ray1[[3]] /. t->tint1)-(ray2[[3]] /. t->tint2)] == 0,Im[(ray1[[3]] /. t->tint1)-(ray2[[3]] /. t->tint2)] == 0},{{tint1,tint1approx,0,ray1[[3]][[0,1,1,2]]},{tint2,tint2approx,0,ray2[[3]][[0,1,1,2]]}}],err]
];

(* if the root-finder threw an error, skip this candidate *)
If[intpt==err,Continue[]] ;

(* otherwise do some mild sanity checking and if passed, add it to the list *)
If[(tint1 /. intpt)>0 && (tint2 /. intpt)>0 && Abs[(ray1[[3]] /. t->tint1 /. intpt)-(ray2[[3]] /. t->tint2 /. intpt)]<0.001,
preciseints=Append[preciseints,{{raynum1,tint1 /. intpt, ray1[[3]] /. t->tint1 /. intpt},{raynum2, tint2 /. intpt, ray2[[3]] /. t->tint2 /. intpt}}];
];
];
  
Sort[preciseints,
(
(#1[[1,1]]>#2[[1,1]]) || 
(#1[[1,1]]==#2[[1,1]] && #1[[2,1]]==#2[[2,1]]) || (#1[[1,1]]==#2[[1,1]] && #1[[2,1]]==#2[[2,1]] && #1[[1,2]] > #2[[2,2]])
) 
&]
)];


(* Classify the intersection, return initial data for a new ray and truncation data for an old ray, if applicable: for 3-way intersection, returns {3,nid,newstartt,{thruray,intersectiontime},arctype,newparents}; for 2-way, 
returns {2,nid,newstartt,{0,0,0},arctype,newparents} or {2,{},L,{0,0,0},0,{}} *)
classifyintersection[rays_,int_, L_]:=
Module[{nrays,i,j,ray1num,ray2num,ray1,ray2,tint1,tint2,delta1,delta2,intz,newinitialdata,newstartt,thirdray,thruray,thirdraytime, thruraytime,ray1velocity,ray2velocity,intersectionorientation,arctype,newparents},
(
newinitialdata={};
newstartt=L;
thruray=0;
thruraytime=0;
thirdray=0;
thirdraytime=0;
arctype=0;
nrays = Length[int];
newparents={};

(* Look for two rays that have exactly one index in common *)
For[i=1,i<=nrays, i++, 
For[j=i+1,j<=nrays,j++,
ray1num=int[[i]][[1]];
ray2num=int[[j]][[1]];
ray1=rays[[ray1num]];
ray2=rays[[ray2num]];
tint1=int[[i]][[2]];
tint2=int[[j]][[2]];
intz=int[[i]][[3]];

delta1 = Abs[(ray1[[1]] /. t->tint1)-(ray2[[2]] /. t->tint2)];
delta2 = Abs[(ray1[[2]] /. t->tint1)-(ray2[[1]] /. t->tint2)];

(* calculate intersection orientation *)
ray1velocity = D[ray1[[3]],t] /. t->tint1;
ray2velocity = D[ray2[[3]],t] /. t->tint2;
intersectionorientation=Sign[Arg[ray1velocity/ray2velocity]];

If[nrays==3,
If[i==1 && j==2, thirdray=int[[3]][[1]]; thirdraytime=int[[3]][[2]]];
If[i==2 && j==3, thirdray=int[[1]][[1]]; thirdraytime=int[[1]][[2]]];
If[i==1 && j==3, thirdray=int[[2]][[1]]; thirdraytime=int[[2]][[2]]];
];

If[delta2 < 0.0002 && delta1>0.0002,
newinitialdata={ray1[[1]] /. t->tint1, ray2[[2]] /. t->tint2, intz}; 
newstartt=tint1+tint2;
thruray=thirdray;
thruraytime=thirdraytime;
arctype=intersectionorientation;
newparents={ray1num,ray2num,thirdray};
];
If[delta1 < 0.0002 && delta2>0.0002,
newinitialdata={ray2[[1]] /. t->tint2, ray1[[2]] /. t->tint1, intz}; 
newstartt=tint1+tint2;
thruray=thirdray;
thruraytime=thirdraytime;
arctype=-intersectionorientation;
newparents={ray1num,ray2num,thirdray};
];

];
];

{nrays,newinitialdata,newstartt,{thruray,thruraytime},arctype,newparents}
)
];


(* consolidate intersection list to include the possibility of triples: this list is pretty small hopefully, so we can afford to use dumb O(N^2) algorithm. Some of them have been consolidated already at a previous stage: these may not be moved from their current place in the list; this should work automatically. *)
consolidateintersections[ints_]:=Module[
{i,j,used,currentgroup, groups, enteringrays, intersectionpos,usedrays,consolidatedints,currentints},
used={};
groups={};

(* find out which intersections should be grouped *)
For[i=1, i<=Length[ints],i++,
If[MemberQ[used,i],Continue[]];
currentgroup={};

(* find all intersections which should be grouped with the current one *)
(* each "group" will be a list of intersections *)
(* if you make the threshold too big you can accidentally group distinct ones, with bad results *)
For[j=1,j<=Length[ints],j++,
If[MemberQ[used,j],Continue[]];
If[Abs[ints[[i,1,3]]-ints[[j,1,3]]]<0.0001,currentgroup=Append[currentgroup,j]];
];
groups=Append[groups,currentgroup];
used=Join[used, currentgroup];
];

(* for each group, find all distinct rays that are entering *)
consolidatedints={};
usedrays={};
For[i=1,i<=Length[groups],i++,
usedrays={};
currentints={};

(* scan over all intersections in this group, pick out each ray that enters *)
For[j=1,j<=Length[groups[[i]]],j++,
For[k=1,k<=Length[ints[[groups[[i,j]]]]],k++,
currentray = ints[[groups[[i,j]]]] [[k,1]];
(* if we've already seen this ray, skip it *)
If[MemberQ[usedrays,currentray],Continue[]];
usedrays=Append[usedrays,currentray];
currentints=Append[currentints, ints[[groups[[i,j]]]][[k]] ];
];
];

(* build a consolidated intersection out of the rays involved in this group *)
If[Length[usedrays]>3,Print["Warning: intersection "<>ToString[Length[consolidatedints]+1]<>" has "<>ToString[Length[usedrays]]<>" rays "<>ToString[usedrays]]];
consolidatedints=Append[consolidatedints,currentints];
];

consolidatedints
];


(* list all rays and turning points in the genealogy of a given ray *)
genealogy[raynum_,parents_]:=
Module[{i},
Reap[
For[i=1,i<=Length[parents[[raynum]]],i++,
Sow[parents[[raynum,i]]];
If[parents[[raynum,i]]>0,Sow[genealogy[parents[[raynum,i]],parents]]];
];
][[2,1]]
];


(* list of unique rays and turning points in the genealogy of a given list of rays *)
ugenealogy[raylist_,parents_]:=
Module[{i},
Union[Flatten[Table[genealogy[raylist[[i]],parents],{i,1,Length[raylist]}]]]
];


killbirthintersections[rays_,ints_,parents_,OptionsPattern[]]:=Module[{i,j,k,constituent,inttime,raystarttime,cleanedlist,tooclose,raynum,secondraynum,secondconstituent,secondinttime,secondraystarttime},
cleanedlist={};
For[i=1,i<=Length[ints],i++,
tooclose=False;

(* figure out whether this intersection is too close to birth of one of the rays: "j" runs over all rays at the intersection *)
For[j=1,j<=Length[ints[[i]]],j++,
constituent=ints[[i,j]];

(* constituent has form {raynum, inttime, intpos} *)
raynum=constituent[[1]];
inttime=constituent[[2]];
raystarttime=rays[[raynum]][[3]][[0,1,1,1]];
If[Abs[inttime-raystarttime]<OptionValue[KillTimeThreshold], 
For[k=1,k<=Length[ints[[i]]],k++,
If[j==k,Continue[]];
secondconstituent=ints[[i,k]];
secondraynum=secondconstituent[[1]];

(* check whether the int. includes a newborn ray and one of its parents *)
If[MemberQ[parents[[raynum]],secondraynum],tooclose=True;];

(* check whether the int. includes two newborn rays *)
secondinttime=secondconstituent[[2]];
secondraystarttime=rays[[secondraynum]][[3]][[0,1,1,1]];
If[Abs[secondinttime-secondraystarttime]<OptionValue[KillTimeThreshold], tooclose=True;];
];
];
];
(* if it's not too close, add it to the list to be returned *)
If[tooclose==False,cleanedlist=Append[cleanedlist,ints[[i]]] ];
];
cleanedlist
];


mswn[L_,\[Alpha]_,opts:OptionsPattern[]]:=Module[{rays,noldrays,intcand,ints,newints,allints,nid,secondaryrays,k,ages,windings,solitoncounts,newstartts, noldints,intswithsecondaryrays,intswithtruncations,eliminatedints,newray,i,j,currentstate, keepgoing, inttypes,joblist,currentint,raytotruncate,parents,thirdraysolitons,thirdraywinding,newwinding,newparents,truncationtime,newage,newsolitoncount},

intswithsecondaryrays={};
intswithtruncations={};
eliminatedints={};

currentstate="Generating initial rays";
rays=alltprays[L,\[Alpha]];
ages=Table[0,{k,1,Length[rays]}];
windings=Table[0,{k,1,Length[rays]}];
solitoncounts=Table[1,{k,1,Length[rays]}];
(* mark initial rays with parents 0 and -n where n is the number of the turning point *)
parents=Table[{0,-Ceiling[k/3],0},{k,1,Length[rays]}];
noldrays=0;
ints={};
Clear[ht];

While[OptionValue[DoIntersections] && noldrays!=Length[rays] && KK>2,
(* look for intersections between newly generated rays and new or previously existing rays *)
currentstate="Finding intersection candidates";
newints=findintersectioncandidates[rays,noldrays,parents,FilterRules[{opts},Options[findintersectioncandidates]]];
If[debug, Print[ToString[Length[newints]]<>" raw candidates"]];
If[debugmore,Print["Raw candidates: "<>ToString[newints]]];

currentstate ="Reducing intersection candidates, first round ("<>ToString[Length[newints]]<>" of them)";
newints=reduceintersectioncandidates[newints];
If[debug, Print[ToString[Length[newints]]<>" candidates after reduction"]];
If[debugmore,Print[newints]];

currentstate = "Refining intersection candidates ("<>ToString[Length[newints]]<>" of them)";
newints=refineintersections[rays,newints];
If[debug,Print["Candidate table after first refinement: "<>ToString[newints]]];

currentstate = "Killing intersections that occur at birth time";
newints=killbirthintersections[rays,newints,parents,FilterRules[{opts},Options[killbirthintersections]]];
If[debug,Print["Candidate table after killing: "<>ToString[newints]]];

currentstate = "Reducing intersection candidates, second round ("<>ToString[Length[newints]]<>" of them)";
newints=reduceintersectioncandidates[newints];
If[debug,Print["Candidate table after second reduction: "<>ToString[newints]]];

(* build intersection list containing both new intersections and previous ones; consolidate them *)
currentstate = "Consolidating intersection list ("<>ToString[Length[newints]]<>" of them)";
allints=Join[ints,newints];
allints=consolidateintersections[allints];

noldrays=Length[rays];

currentstate = "Classifying intersections";

inttypes={};
(* get basic data about each intersection *)
For[i=1, i<=Length[allints], i++,
currentstate = "Processing intersection "<>ToString[i]<>" of "<>ToString[Length[allints]];
inttypes=Append[inttypes,classifyintersection[rays,allints[[i]],L]]
];

(* sort intersections by which one will emit a new ray first *)
currentstate="Sorting intersections into job list";
joblist=Table[{i,inttypes[[i]][[3]]},{i,1,Length[allints]}];
joblist=Sort[joblist,#1[[2]]<#2[[2]] &];

(* now work down the job list *)
For[i=1,i<=Length[joblist],i++,

currentint=joblist[[i,1]];
If[debugmore,Print["Item "<>ToString[i]<>" of "<>ToString[Length[joblist]]<>" on int list is intersection "<>ToString[currentint]<>" with size "<>ToString[inttypes[[currentint]][[1]] ]<>" at position "<>ToString[allints[[currentint]][[1,3]]] <> " involving rays " <> ToString[Table[allints[[currentint]][[k,1]],{k,1,inttypes[[currentint]][[1]]}]]  ]];

(* if we already have a ray emanating from here, then skip *)
If[MemberQ[intswithsecondaryrays, currentint],Continue[]];

(* if we already truncated at this intersection, then skip *)
If[MemberQ[intswithtruncations, currentint],Continue[]];

(* if we already eliminated this intersection, then skip *)
If[MemberQ[eliminatedints, currentint],Continue[]];

(* if the new ray would start out beyond the cutoff, then skip *)
If[inttypes[[currentint]][[3]]>=L,Continue[]];

(* if there were 3 rays colliding here, then truncate one of them; truncation data is contained in inttype[[4]] *)
If[inttypes[[currentint]][[1]]==3,
raytotruncate=inttypes[[currentint]][[4,1]];
(* count the solitons and winding from the truncated ray *)
thirdraysolitons = solitoncounts[[raytotruncate]];
thirdraywinding = windings[[raytotruncate]];
(* truncate all 3 components of the ray function *)
truncationtime=inttypes[[currentint]][[4,2]];
For[j=1,j<=3,j++,rays[[raytotruncate]][[j]][[0,1,1,2]]=truncationtime];
intswithtruncations=Append[intswithtruncations,currentint];

If[debug,Print["Truncated ray "<>ToString[raytotruncate]<>" at time "<>ToString[truncationtime]]];

(* seek through the intersection list, find later intersections involving the truncated ray, replace it with the new one *)
For[j=1,j<=Length[allints],j++,
For[k=1,k<=Length[allints[[j]]],k++,
If[j!=currentint && allints[[j,k]][[1]]==raytotruncate && allints[[j,k]][[2]]>truncationtime,
allints[[j,k]][[1]]=Length[rays]+1;
If[debug,Print["Replaced ray "<>ToString[raytotruncate]<>" with ray "<>ToString[Length[rays]+1]<>" in intersection "<>ToString[j]];
];
];
];
];
,
thirdraysolitons=0; 
];

(* OK, now generate a new ray *)
currentstate="Generating new ray number "<>ToString[Length[rays]+1]<>" born at time "<>ToString[joblist[[i,2]]];
newray=singleray[inttypes[[currentint]][[2,3]],inttypes[[currentint]][[2,1]],inttypes[[currentint]][[2,2]],L,\[Alpha],inttypes[[currentint]][[3]]];
intswithsecondaryrays=Append[intswithsecondaryrays,currentint];
rays=Append[rays,newray];

(* record its parentage *)
newparents=inttypes[[currentint]][[6]];
parents=Append[parents,newparents];


newage = 1+Max[{ages[[newparents[[1]]]], ages[[newparents[[2]]]]}];
If[newparents[[3]]!=0, newage=Min[{newage,ages[[newparents[[3]] ]]}] ];
ages=Append[ages,newage];

(* winding is sum of windings of parents plus intersection arctype *)
newwinding = windings[[newparents[[1]]]]+windings[[newparents[[2]]]]+inttypes[[currentint]][[5]];
windings=Append[windings, newwinding];

(* soliton count is product of soliton counts of parents *)
newsolitoncount =  solitoncounts[[newparents[[1]]]] solitoncounts[[newparents[[2]]]] + thirdraysolitons Exp[Pi I/2 (newwinding-thirdraywinding)] ;
solitoncounts=Append[solitoncounts,newsolitoncount];

(* print debugging info *)
If[debug,Print["New ray number "<>ToString[Length[rays]]<>" born from intersection "<>ToString[currentint]<>" at time "<>ToString[joblist[[i,2]]]<>" and position "<>ToString[inttypes[[currentint]][[2,3]]] <> " with parents "<>ToString[newparents]<>", age "<>ToString[newage]<>", soliton count "<>ToString[newsolitoncount] ] ];

(* The next intersection to process might be one involving this new ray; so get out of this loop to reprocess the intersection list *)
Break[];
];

(* get ready for the next cycle *)
ints=allints;
];

(* store intersection and ray data as global variable for later perusal *)
inttable=ints;
raytable=rays;
windingtable=windings;
solitontable=solitoncounts;
parenttable=parents;
agetable=ages;

(* finally return the data of rays, their ages, soliton counts *)
{rays,ages,solitoncounts}
];


period[z_,x10_,x20_]:=Module[{},
Z[1] /. NDSolve[{Derivative[1][x1][t]==dxdz[x1[t],z[t]] Derivative[1][z][t],Derivative[1][x2][t]==dxdz[x2[t],z[t]] Derivative[1][z][t],Z'[t]==(x1[t]-x2[t])z'[t],x1[0]==x10,x2[0]==x20,Z[0]==0},{x1,x2,Z},{t,0,1},MaxSteps->1000000,PrecisionGoal->\[Infinity]][[1]]
];

endpoint[z_,x10_,x20_]:=Module[{},
{Z[1],x1[1],x2[1]} /. NDSolve[{Derivative[1][x1][t]==dxdz[x1[t],z[t]] Derivative[1][z][t],Derivative[1][x2][t]==dxdz[x2[t],z[t]] Derivative[1][z][t],Z'[t]==(x1[t]-x2[t])z'[t],x1[0]==x10,x2[0]==x20,Z[0]==0},{x1,x2,Z},{t,0,1},MaxSteps->1000000,PrecisionGoal->\[Infinity]][[1]]
];

pathfrom[z1_,z2_]:=Function[t, z1+(z2-z1)t];

trackroot[z_,x0_]:=
x[1] /. NDSolve[{Derivative[1][x][t]==dxdz[x[t],z[t]] Derivative[1][z][t],x[0]==x0},{x},{t,0,1},MaxSteps->1000000,PrecisionGoal->\[Infinity]][[1]];

trackroots[initialroots_,z_]:=
Table[trackroot[z,initialroots[[i]]],{i,1,KK}];

directionfrom[z1_,z2_]:=(z2-z1)/Abs[z2-z1];

findsmallroots[tp_,direction_]:=
Module[{\[Delta],z0, xr,pairs,xdiff,pairindex,i,j},
\[Delta]=0.00005;
z0=tp+\[Delta]^2 direction;
xr=x/.  roots[z0];
pairs=pairbasis;
xdiff=Table[xr[[pairs[[i,1]]]]-xr[[pairs[[i,2]]]],{i,1,Length[pairs]}];
pairindex=Ordering[Abs[xdiff],1][[1]];
i=pairs[[pairindex]][[1]];
j=pairs[[pairindex]][[2]];
{z0,xr[[i]],xr[[j]]}
];


Off[ParametricPlot::icolscale];


pplot[curves_]:=pplot[curves,Table[Hue[i/Length[curves]],{i,1,Length[curves]}]];


pplotplain[curves_]:=pplot[curves,Table[Gray,{i,1,Length[curves]}]];


pplot[curves_,styles_]:=ParametricPlot[Evaluate@Table[{Re[curves[[i]]],Im[curves[[i]]]} /. t->(curves[[i]][[0,1,1,2]]-curves[[i]][[0,1,1,1]]) t +curves[[i]][[0,1,1,1]] ,{i,1,Length[curves]}],{t,0,1},PlotStyle->styles];


singleplotwithpointsandlabel[rays_,plotrange_,L_,label_]:=Show[singleplotwithpoints[rays,plotrange,L],Graphics[Text[label,Scaled[{0.02,0.98}],{-1,1}]]];


toplane[z_]:={Re[z],Im[z]};


tocylinder[z_]:=If[z==Infinity,{0,0,1},If[z==0,{0,0,-1},{Re[z]/Abs[z],Im[z]/Abs[z], Log[Abs[z]]}]];


tounrolledcylinder[z_,phaseoffset_]:={-Arg[z Exp[I phaseoffset]],Log[Abs[z]]};


tosphere[z_]:=If[z==Infinity,{0,0,1},{2 Re[z]/(1 + Abs[z]^2), 2 Im[z]/(1 + Abs[z]^2), (Abs[z]^2-1)/(Abs[z]^2+1)}];


Options[pointsplotwithoptions]={Surface->Plane,CrossSize->Automatic,DotSize->Automatic,Marker->Cross,PhaseOffset->0.3};


cross[x_,y_,scale_]:=Polygon[{{x,y}+scale{1,2},{x,y}+scale{2,1},{x,y}+scale{1,0},{x,y}+scale{2,-1},{x,y}+scale{1,-2},{x,y}+scale{0,-1},{x,y}+scale{-1,-2},{x,y}+scale{-2,-1},{x,y}+scale{-1,0},{x,y}+scale{-2,1},{x,y}+scale{-1,2},{x,y}+scale{0,1}}];


pointsplotwithoptions[points_,opts:OptionsPattern[]]:=Module[{pointstoplot},
Which[
OptionValue[Surface]===Plane,
pointstoplot=Select[points,#!=Infinity &];
If[OptionValue[Marker]===Cross,
Graphics[Table[cross[toplane[pointstoplot[[i]]][[1]], toplane[pointstoplot[[i]]][[2]],OptionValue[CrossSize]],{i,1,Length[pointstoplot]}]],
Graphics[Table[Disk[{toplane[pointstoplot[[i]]][[1]], toplane[pointstoplot[[i]]][[2]]},OptionValue[DotSize]],{i,1,Length[pointstoplot]}]]
],

OptionValue[Surface]===Sphere,
Graphics3D[Table[Point[tosphere[points[[i]]]],{i,1,Length[points]}]],

OptionValue[Surface]===Cylinder,
pointstoplot=Select[points,#!=Infinity && #!=0 &];
Graphics3D[Table[Point[tocylinder[pointstoplot[[i]]]],{i,1,Length[pointstoplot]}]],

OptionValue[Surface]===UnrolledCylinder,
pointstoplot=Select[points,#!=Infinity && #!=0 &];
Graphics[Table[cross[tounrolledcylinder[pointstoplot[[i]],OptionValue[PhaseOffset]][[1]],tounrolledcylinder[pointstoplot[[i]], OptionValue[PhaseOffset]][[2]],OptionValue[CrossSize]],{i,1,Length[pointstoplot]}]]

]
];


plotbasicoptions=If[alloworiented,Options[OrientedParametricPlot],Join[Options[ParametricPlot],{Oriented->Off}]];


Options[pplotwithoptions]=Join[{CurveThickness->0.007,FiltrationOpacity->Off,Label->Off,LabelText->"",Surface->Plane,SurfaceColor->LightBlue,SurfaceOpacity->0.5,Rainbow->Off,ColorByAge->On,WallColor->Black,LengthTable->"",UseLengthTable->Off,UseColorFunction->Off,CylinderHeight->5,RainbowMax->Automatic, PlotRange->5,PlotCenterX->0,PlotCenterY->0,ArrowSize->Medium,SuppressJumps->On,
SuppressDeadLines->On,ShowRays->{}},plotbasicoptions,
Options[pointsplotwithoptions]];


Options[singleplotwithoptions]=Join[Options[pplotwithoptions],{TurningPoints->{},Singularities->{}}];


pplotwithoptions[multicurves_,multiages_,solitoncounts_,L_,opts:OptionsPattern[]]:=
Module[{coloropts,thickopts,outtable,startt,endt,rbmax,psopts,otheropts,curves,ages,i},
psopts={};
otheropts={};
coloropts={};

curves=multicurves[[All,3]];
ages=multiages;
(* select some subset of the rays to be shown *)
If[OptionValue[ShowRays]!={},
{curves,ages}=Reap[
For[i=1,i<=Length[curves],i++,
If[MemberQ[OptionValue[ShowRays],i],Sow[curves[[i]],x];Sow[ages[[i]],y]];
];
][[2]],
If[OptionValue[SuppressDeadLines]===On,
{curves,ages}=Reap[
For[i=1,i<=Length[curves],i++,
If[solitoncounts[[i]]!=0,Sow[curves[[i]],x];Sow[ages[[i]],y]];
];
][[2]];
];
];

If[OptionValue[RainbowMax]===Automatic,
If[OptionValue[ColorByAge]===On,
rbmax=Max[ages]+1,
rbmax=Length[curves]
],
rbmax=OptionValue[RainbowMax]
];
If[OptionValue[UseLengthTable]===Off,
startt=Table[curves[[i]][[0,1,1,1]],{i,1,Length[curves]}];
endt=Table[Min[{L,curves[[i]][[0,1,1,2]]}],{i,1,Length[curves]}];
,
startt=Table[curves[[i]][[0,1,1,1]],{i,1,Length[curves]}];
endt=OptionValue[LengthTable];
];
outtable=Table[
If[curves[[i]][[0,1,1,1]]<Min[{L,curves[[i]][[0,1,1,2]]}],
Which[
OptionValue[Surface]===Plane,
thickopts=Thickness[OptionValue[CurveThickness]];
If[OptionValue[UseColorFunction]===On,
If[OptionValue[FiltrationOpacity]===On,coloropts=(ColorFunction->Function[{x,y,z},Hue[0,1,1,1-z/L]])
,
otheropts=Append[otheropts,(ColorFunction->Function[{x,y,z},Hue[0,1,1,1]])];
]
,
If[OptionValue[Rainbow]===Off, 
psopts=Append[psopts,(OptionValue[WallColor])];
,
If[OptionValue[ColorByAge]===On,
psopts=Append[psopts,(Hue[ages[[i]]/rbmax])],
psopts=Append[psopts,(Hue[i/rbmax])]
]
]
];
If[OptionValue[Oriented]===True && alloworiented,
psopts = Append[psopts,Arrowheads->OptionValue[ArrowSize]];
arrowpos=1/2;
If[arrowpos>1 || arrowpos<0,arrowlist={}, arrowlist={arrowpos};];
OrientedParametricPlot[
toplane[curves[[i]]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{psopts},
ArrowPositions->arrowlist,
ColorFunctionScaling->False
]
,
ParametricPlot[
toplane[curves[[i]]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{psopts},
ColorFunctionScaling->False
]
]
,
OptionValue[Surface]===Sphere,
thickopts=Thickness[0.005];
If[OptionValue[UseColorFunction]===On,
If[OptionValue[FiltrationOpacity]===On,coloropts=(ColorFunction->Function[{x,y,z,w},Hue[0,1,1,1-w/L]])
,
coloropts=(ColorFunction->Function[{x,y,z,w},Hue[0,1,1,1]])
]
,
If[OptionValue[Rainbow]===Off, 
coloropts=OptionValue[WallColor]
,
coloropts=Hue[ages[[i]]/rbmax]
]
];
ParametricPlot3D[
tosphere[curves[[i]]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{Evaluate[coloropts],Evaluate[thickopts]},
ColorFunctionScaling->False],

OptionValue[Surface]===Cylinder,
thickopts=Thickness[0.005];
If[OptionValue[UseColorFunction]===On,
If[OptionValue[FiltrationOpacity]===On,coloropts=(ColorFunction->Function[{x,y,z,w},Hue[0,1,1,1-w/L]])
,
coloropts=(ColorFunction->Function[{x,y,z,w},Hue[0,1,1,1]])
]
,
If[OptionValue[Rainbow]===Off, 
coloropts=(PlotStyle->OptionValue[WallColor])
,
If[OptionValue[ColorByAge]===On,
coloropts=(PlotStyle->Hue[ages[[i]]/rbmax]),
coloropts=(PlotStyle->Hue[i/rbmax])
]
]
];
If[OptionValue[Oriented]===True && alloworiented,
OrientedParametricPlot3D[
tocylinder[curves[[i]]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{Evaluate[coloropts],Evaluate[thickopts],Arrowheads->OptionValue[ArrowSize]},
ColorFunctionScaling->False]
,
ParametricPlot3D[
tocylinder[curves[[i]]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{Hue[i/rbmax],Thickness[0.005]},
ColorFunctionScaling->False
]
],

OptionValue[Surface]===UnrolledCylinder,
If[OptionValue[UseColorFunction]===On,
If[OptionValue[FiltrationOpacity]===On,coloropts=(ColorFunction->Function[{x,y,z},Hue[0,1,1,1-z/L]])
,
otheropts=Append[otheropts,(ColorFunction->Function[{x,y,z},Hue[0,1,1,1]])];
]
,
If[OptionValue[Rainbow]===Off, 
psopts=Append[psopts,(OptionValue[WallColor])];
,
If[OptionValue[ColorByAge]===On,
psopts=Append[psopts,(Hue[ages[[i]]/rbmax])],
psopts=Append[psopts,(Hue[i/rbmax])]
]
]
];
If[OptionValue[SuppressJumps]===On,
jumps=Select[Range[startt[[i]]+0.001,endt[[i]],0.001],Abs[Arg[(curves[[i]]Exp[I OptionValue[PhaseOffset]] /. t->#)]-Arg[(curves[[i]] Exp[I OptionValue[PhaseOffset]]/. t->(#-0.001))]]>0.5 &],
jumps={}];
If[OptionValue[Oriented]===True && alloworiented,
OrientedParametricPlot[
tounrolledcylinder[curves[[i]], OptionValue[PhaseOffset]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{Evaluate[psopts],Arrowheads->OptionValue[ArrowSize]},
ColorFunctionScaling->False,Exclusions->jumps]
,
ParametricPlot[
tounrolledcylinder[curves[[i]], OptionValue[PhaseOffset]],
{t,startt[[i]],endt[[i]]},
PlotStyle->{Evaluate[psopts]},
ColorFunctionScaling->False,Exclusions->jumps]
]
],
{}
],
{i,1,Length[curves]}];
If[OptionValue[Label]==On,
If[OptionValue[Surface]==Plane,
outtable=Union[outtable,{Graphics[Text[OptionValue[LabelText],Scaled[{0.02,0.98}],{-1,1}]]}
]
];
If[OptionValue[Surface]==Sphere,
outtable=Union[outtable,
{Graphics3D[Text[OptionValue[LabelText],Scaled[{0.02,0.98,0.98}],{-1,1}]]}
]
];
];
outtable
];


singleplotwithoptions[rays_,L_,opts:OptionsPattern[]]:=Module[{plotrange,plotcenterx,plotcentery,crosssize,dotsize},
plotrange=OptionValue[PlotRange];
plotcenterx=OptionValue[PlotCenterX];
plotcentery=OptionValue[PlotCenterY];
If[OptionValue[CrossSize]===Automatic,crosssize=1/75 plotrange,crosssize=OptionValue[CrossSize]];
If[OptionValue[DotSize]===Automatic,dotsize=1/75 plotrange,dotsize=OptionValue[DotSize]];
Which[
OptionValue[Surface]===Plane,
Apply[
Show[##,
Graphics[Orange],pointsplotwithoptions[OptionValue[TurningPoints],Surface->Plane, CrossSize->crosssize,Marker->Cross],
Graphics[Blue],pointsplotwithoptions[OptionValue[Singularities],Surface->Plane,Marker->Dot,DotSize->dotsize],PlotRange->{{-plotrange+plotcenterx,plotrange+plotcenterx},{-plotrange+plotcentery,plotrange+plotcentery}},Axes->None] &,
pplotwithoptions[rays[[1]],rays[[2]],rays[[3]],L,FilterRules[{opts},Options[pplotwithoptions]]]
],
OptionValue[Surface]===UnrolledCylinder,
Apply[
Show[##,
Graphics[Orange],pointsplotwithoptions[OptionValue[TurningPoints],Surface->UnrolledCylinder,CrossSize->crosssize, Marker->Cross,PhaseOffset->OptionValue[PhaseOffset]],
Graphics[Blue],pointsplotwithoptions[OptionValue[Singularities],Surface->UnrolledCylinder,Marker->Dot,DotSize->dotsize,PhaseOffset->OptionValue[PhaseOffset]],PlotRange->{{-plotrange+plotcenterx,plotrange+plotcenterx},{-plotrange+plotcentery,plotrange+plotcentery}},Axes->None] &,
pplotwithoptions[rays[[1]],rays[[2]],rays[[3]],L,FilterRules[{opts},Options[pplotwithoptions]]]
],
OptionValue[Surface]===Sphere,
Apply[
Show[##,
Graphics3D[Orange],pointsplotwithoptions[OptionValue[TurningPoints],Surface->Sphere],Graphics3D[Blue],pointsplotwithoptions[OptionValue[Singularities],Surface->Sphere],
Graphics3D[OptionValue[SurfaceColor]],Graphics3D[Opacity[OptionValue[SurfaceOpacity]]],Graphics3D[Sphere[{0,0,0},0.99]],
PlotRange->{{-1,1},{-1,1},{-1,1}},Boxed->False,Axes->None] &,
pplotwithoptions[rays[[1]],rays[[2]],rays[[3]],L,FilterRules[{opts},Options[pplotwithoptions]]]
],
OptionValue[Surface]===Cylinder,
Apply[
Show[##,
Graphics3D[Orange],pointsplotwithoptions[OptionValue[TurningPoints],Surface->Cylinder],Graphics3D[Blue],pointsplotwithoptions[OptionValue[Singularities],Surface->Cylinder],
Graphics3D[OptionValue[SurfaceColor]],Graphics3D[Opacity[OptionValue[SurfaceOpacity]]],Graphics3D[Cylinder[{{0,0,-OptionValue[CylinderHeight]},{0,0,OptionValue[CylinderHeight]}},0.99]],
PlotRange->{{-1.5,1.5},{-1.5,1.5},{-OptionValue[CylinderHeight],OptionValue[CylinderHeight]}},Boxed->False,Axes->None] &,
pplotwithoptions[rays[[1]],rays[[2]],rays[[3]],L,FilterRules[{opts},Options[pplotwithoptions]]]
]
]
];


plotbpsstateatphase[L_,\[Alpha]_,timestep_,grain_,maxdepth_,plotrange_,opts:OptionsPattern[]]:=
Module[{network,rays,ages,parents,bps,raystoplot,agestoplot,halves,lengths},
network=mswn[L,\[Alpha],maxdepth];
rays=network[[1]];
ages=network[[2]];
parents=network[[3]];
If[Length[network[[4]]]>0,
bps=network[[4,1]];
halves=Table[gethalfbpsstate[{bps[[i,1]],bps[[i,2]]},rays,ages,parents],{i,1,2}];
raystoplot=Join[halves[[1,1]],halves[[2,1]]];
agestoplot=Join[halves[[1,2]],halves[[2,2]]];
lengths=Join[halves[[1,3]],halves[[2,3]]];
Show[
singleplotwithoptions[{rays,ages},L,opts],
singleplotwithoptions[{raystoplot,agestoplot},L,UseLengthTable->True,LengthTable->lengths,Hue->2/3,opts]
]
,
{}
]
];


tpplot:=pointsplotwithoptions[turningpoints];
spplot:=pointsplotwithoptions[sing];


singleplotshortwithpoints[rays_,plotrange_,L_]:=Show[pplotshort[rays,L],Graphics[Orange],tpplot, Graphics[Blue],spplot,PlotRange->{{-plotrange,plotrange},{-plotrange,plotrange}},Axes->None];


singlemswnplot[L_,\[Alpha]_,timestep_,grain_,maxdepth_,plotrange_]:=
Module[{network},
network=mswn[L,\[Alpha],maxdepth];
plotnetwork[network,plotrange]];


plotnetwork[network_,plotrange_]:=Show[pplot[network[[1]],network[[2]]],Graphics[Orange],tpplot,Graphics[Blue],spplot,PlotRange->{{-plotrange,plotrange},{-plotrange,plotrange}},Axes->None];


swngrid[net_,Mtable_,opts:OptionsPattern[]]:=GraphicsGrid[Map[singleplotwithoptions[net,#, FilterRules[{opts},Options[pplotwithoptions]],Label->On,LabelText-> ("\[CapitalLambda] = "<>ToString[#])] &,Mtable,{2}]];

Options[parallelmswns]=Join[Options[mswn],{Parallel->On}];
parallelmswns[phasetable_,L_,paramtable_,opts:OptionsPattern[]]:=
Module[{outtable,curtable,curitem,curparam,i,j,tpstored,singstored},
outtable={};
For[j=1,j<=Length[paramtable],j++,
curparam=j;
init[paramtable[[curparam]]];
tpstored[j]=turningpoints;
singstored[j]=sing;

DistributeDefinitions[$Context];
If[OptionValue[Parallel]===On,
curtable=ParallelTable[mswn[L,phasetable[[curitem]],FilterRules[{opts},Options[mswn]] ],{curitem,1,Length[phasetable]}],
curtable=Table[mswn[L,phasetable[[curitem]],FilterRules[{opts},Options[mswn]] ],{curitem,1,Length[phasetable]}]
];
outtable=Append[outtable,curtable];
];
{outtable,tpstored,singstored}
];

Options[showplots]=Join[Options[parallelmswns],Options[singleplotwithoptions],{PhaseTable->{},ParamTable->{},InitialPlotRange->5,InitialPlotCenterX->0,InitialPlotCenterY->0,ShowProgress->On,MassCutoff->1}];
showplots[opts:OptionsPattern[]]:=
Module[{phasetable,outtable,curtable,curitem,paramtable,curparam,i,j,L,tpstored,singstored},
L=OptionValue[MassCutoff];
phasetable=OptionValue[PhaseTable];
If[phasetable=={},phasetable={0}];
paramtable=OptionValue[ParamTable];
If[paramtable=={},paramtable={0}];
{outtable,tpstored,singstored}=parallelmswns[phasetable,L,paramtable,FilterRules[{opts},Options[parallelmswns]] ];

Manipulate[
singleplotwithoptions[outtable[[j,i]],M,otheroptions,FilterRules[{opts},Options[pplotwithoptions]],TurningPoints->tpstored[j],Singularities->singstored[j],PlotRange->plotrange, PlotCenterX->plotcenterx, PlotCenterY->plotcentery]
,
{{M,L,"Mass cutoff"},0,L},
{{i,1,"Phase"},1,Length[phasetable],1},
{{j,1,"Plot parameters"},1,Length[paramtable],1},
{{plotrange, OptionValue[InitialPlotRange], "Plot range"},0, 5*OptionValue[InitialPlotRange]},
{{plotcenterx,OptionValue[InitialPlotCenterX], "Plot center (x)"}, OptionValue[InitialPlotCenterX]-5*OptionValue[InitialPlotRange], OptionValue[InitialPlotCenterX]+5*OptionValue[InitialPlotRange]},
{{plotcentery, OptionValue[InitialPlotCenterY], "Plot center (y)"},OptionValue[InitialPlotCenterX]-5*OptionValue[InitialPlotRange],OptionValue[InitialPlotCenterY]+5*OptionValue[InitialPlotRange]},
{{otheroptions,{},"Other options"}},
ControlPlacement->Left]
];

showbpsphases:=Show[ListPlot[Transpose[{Cos[2 bpsphases],Sin[2 bpsphases]}],Axes->None,PlotRange->{{-1.1,1.1},{-1.1,1.1}},AspectRatio->1],Graphics[Circle[{0,0},.95]]];

Options[swnfiltrationexport]=Join[Options[showplots],{Granularity->10,Filename->"export"}];
swnfiltrationexport[phasetable_,opts:OptionsPattern[]]:=
Module[{outtable,curtable,curitem,paramtable,curparam,i,j,L,tpstored,singstored,curplot,M},
L=OptionValue[MassCutoff];
paramtable=OptionValue[ParamTable];
If[paramtable=={},paramtable={0}];
{outtable,tpstored,singstored}=parallelmswns[phasetable,L,paramtable,FilterRules[{opts},Options[parallelmswns]]];
For[j=1,j<=Length[paramtable],j++,
For[i=1,i<=Length[phasetable],i++,
For[M=1,M/OptionValue[Granularity]<=L ,M++,
Export[
OptionValue[Filename]<>"-"<>ToString[j-1]<>ToString[PaddedForm[M-1,4,NumberPadding->"0"]]<>".jpg",singleplotwithoptions[outtable[[j,i]],M/OptionValue[Granularity],FilterRules[{opts},Options[pplotwithoptions]],TurningPoints->tpstored[j],Singularities->singstored[j],PlotRange->OptionValue[InitialPlotRange], PlotCenterX->OptionValue[InitialPlotCenterX], PlotCenterY->OptionValue[InitialPlotCenterY]],
ImageSize->16*50
];
];
];
];
];

ReleaseHoldAt[expr_,partspec_]:=ReplacePart[expr,partspec->Extract[expr,partspec]];

Options[exportplots]=Join[Options[showplots],{Randomize->On,Filename->"/home/andy/movies",Filetype->"png"}];
exportplots[opts:OptionsPattern[]]:=
Module[{phasetable,paramtable,i,j,k,targetfile,jobs,nextjob,paralleljobs,shuffledlist,paramdump,L},

L=OptionValue[MassCutoff];

jobs={};
paramtable=OptionValue[ParamTable];
If[paramtable=={},paramtable={0}];
phasetable=OptionValue[PhaseTable];
If[phasetable=={},phasetable={0}];

paramdump= 
ToString[Definition[init],InputForm]<> "\n" <> 
"phasetable = " <>ToString[phasetable,InputForm] <> "\n" <>
"opts = "<> ToString[{opts},InputForm]<> "\n";

Export[OptionValue[Filename]<>".dump", paramdump,"Text"];

If[OptionValue[Randomize]===On,
shuffledlist=RandomSample[Range[Length[phasetable]]],
shuffledlist=Range[Length[phasetable]]
];

CloseKernels[];
LaunchKernels[];
DistributeDefinitions[$Context];
DistributeDefinitions[suppresswarnings];
ParallelEvaluate[suppresswarnings];

For[j=1,j<=Length[paramtable],j++,
For[i=1,i<=Length[phasetable],i++,
targetfile=OptionValue[Filename]<>"-"<>ToString[j-1]<>ToString[PaddedForm[shuffledlist[[i]]-1,4,NumberPadding->"0"]]<>"."<>OptionValue[Filetype];
If[
!FileExistsQ[targetfile],
nextjob=
Hold[
singleexport[targetfile,phasetable[[shuffledlist[[i]]]],paramtable[[j]],L,FilterRules[{opts},Options[singleexport]] ]
];
(* some clumsy manipulation here to get the evaluations to happen in the right order; might be possible to drastically improve this *)
nextjob=ReleaseHoldAt[nextjob,{1,1}];
nextjob=ReleaseHoldAt[nextjob,{1,2}];
nextjob=ReleaseHoldAt[nextjob,{1,3}];
jobs=Append[jobs,nextjob];
];
];
];
DistributeDefinitions[jobs];
DistributeDefinitions["CurvesGraphics6`"];
ParallelTable[jobs[[i]][[1]],{i,1,Length[jobs]}];
];

Options[singleexport]=Options[showplots];
singleexport[targetfile_,phase_,param_,L_,opts:OptionsPattern[]]:=
Module[{network},
init[param];
network=mswn[L,phase,FilterRules[{opts},Options[mswn]] ];

Export[
targetfile,singleplotwithoptions[network,L,FilterRules[{opts},Options[pplotwithoptions]],TurningPoints->turningpoints,Singularities->sing,PlotRange->OptionValue[InitialPlotRange], PlotCenterX->OptionValue[InitialPlotCenterX], PlotCenterY->OptionValue[InitialPlotCenterY]],
ImageSize->16*50
];
];

standardperturbation=0.0001;
prange[n_]:=Pi/n(Range[n]-1+standardperturbation);


init[u_]:=(
KK=$RANK$; sing={Infinity};
P[2][z_]=$SECONDARYPOLY$;
P[$RANK$][z_]=$POLY$;
findtp;
removerepeatedtp;
);

init[0];
o=mswn[$CUTOFF$,$ANGLE$];
interps=o[[1,All,3,0]];

segsclipped[f_,rmax_,rstep_]:=Module[
{tmin,tmax,tclip,tstep,avgderiv,t},
tmin=f[[1,1,1]];
tmax=f[[1,1,2]];
tclip=If[Abs[f[tmax]]<rmax,tmax,Min[tmax,t/.FindRoot[Abs[f[t]]==rmax,{t,tmin,tmax}]]];
Print["Using tclip="<>ToString[tclip]];
avgderiv=Abs[f[tclip] - f[tmin]] / (tclip-tmin);
tstep = Min[rstep / avgderiv, (tclip-tmin)/5];
nstep = Ceiling[(tclip-tmin)/tstep];
Table[{Re[#],Im[#]}&[f[t]],{t,tmin,tclip,(tclip-tmin)/nstep}]
]

traj = segsclipped[#,$RMAX$,$RSTEP$]& /@ interps;

Export["$OUTFN$",{ "turningpoints"->({Re[#],Im[#]}& /@ turningpoints), "trajectories"-> traj }, "JSON"]
