# SCUC consider all contingency-case (N-1) network constraints.
# Reserve requirment is still included in the Master problem, which may not be necessary since N-1 constraints are enforced.
# This N-1 SCUC is solved by Bender decomposition, implemented in the code SCUC.run.
# This code file defines the models for the Master problem and slave sub-problems, as well as the cut.
# This is not the main runnable code; check the code file: SCUC.run.

set BUS;    # set of buses
set BRANCH; # set of branches
set GEND;   # Gen Data
set LOAD;   # Load Percent data of peak load

#@@@@@@@@@@@@@@@
#### PARAMETERS:
# Bus Data
param bus_num		{BUS}; # Bus Number
param bus_Pd		{BUS}; # Real Power Demand 

# GENData
param genD_bus		{GEND}; # GEN location
param genD_minUP	{GEND}; # Min UP Time
param genD_minDN	{GEND}; # Min Down Time
param genD_Pmax		{GEND}; # Max gen production
param genD_Pmin      	{GEND}; # Min gen production when committed
param genC_Startup 	{GEND}; # startup cost
param genC_Cost		{GEND}; # Linear Cost Term
param genC_NLoad	{GEND}; # No Load Cost
param SPRamp		{GEND}; # 10 Min Spin Ramp
param HRamp		{GEND}; # Hourly Ramp
param StartRamp		{GEND}; # Startup/Shutdown Ramp

# Branch Data
param branch_fbus	{BRANCH}; # from bus for line
param branch_tbus	{BRANCH}; # to bus for line
param branch_b		{BRANCH}; # line susceptance
param branch_rateA	{BRANCH}; # long term thermal rating
param branch_rateC	{BRANCH}; # emergency thermal rating
param branch_radial		{BRANCH}; # whether you will monitor the line

set Nk = {j in BRANCH: branch_radial[j] == 0};
set Ng = GEND;

# Load Data
param load_pcnt		{LOAD}; # the percentage of annual peak

# Additional Parameters that are not loaded through sets:
param Bus_Pd{n in BUS, t in LOAD};  # Creates the hourly load per bus
param MBase; let MBase:=100; # the MVA Base

param nCUT;
param nT default 24;

#### VARIABLES:
var obj_M;
var bus_angle {n in BUS, t in LOAD};        # Variable for Bus Angles
var line_flow {j in BRANCH, t in LOAD};     # Variable for all line flows
var gen_supply {g in GEND, t in LOAD};      # Variable for GEN Supply
var reserve{g in GEND, t in LOAD} >= 0;

# Generation Unit Commitment Variables:
var Ugt{g in GEND, t in LOAD} binary; # unit commitment var
var Vgt{g in GEND, t in LOAD} >= 0, <=1; # startup var (binary var modeled as continuous since it will have binary solution)

# dual param used for the cut
param t_iter;

param c_line;
param D1gct {k in 1..nCUT, g in GEND, c in Nk, t in LOAD};  # alpha-
param D2gct {k in 1..nCUT, g in GEND, c in Nk, t in LOAD};  # alpha+
param D3gct {k in 1..nCUT, g in GEND, c in Nk, t in LOAD};  # pseudo N-
param D4gct {k in 1..nCUT, g in GEND, c in Nk, t in LOAD};  # pseudo N+
param CkCut {k in 1..nCUT, c in Nk, t in LOAD};
param MarkNk {k in 1..nCUT, c in Nk, t in LOAD};

param c_gen;
param D5gct {k in 1..nCUT, g in GEND, c in Ng, t in LOAD};  # beta -
param D6gct {k in 1..nCUT, g in GEND, c in Ng, t in LOAD};  # beta +
param D7gct {k in 1..nCUT, g in GEND, c in Ng, t in LOAD}; #: g != c};  # gamma -
param D8gct {k in 1..nCUT, g in GEND, c in Ng, t in LOAD}; #: g != c};  # gamma +
param CgCut {k in 1..nCUT, c in Ng, t in LOAD};
param MarkNg {k in 1..nCUT, c in Ng, t in LOAD};

# for sub-problem
param ugt_fix {g in GEND} binary;
param gen_supply_fix {g in GEND};

var s1 >= 0;
var line_flow_k {k in BRANCH};
var gen_supply_g {g in GEND};
var bus_angle_n {n in BUS};

var s2 >= 0;
var line_flow2_k {k in BRANCH};
var gen_supply2_g {g in GEND};
var bus_angle2_n {n in BUS};

#### OBJECTIVE:
# Objective is to Minimize Cost
minimize M_COST: obj_M;

### cost Function Moved to constraint set
subject to costConstr:
    obj_M >=sum{g in GEND, t in LOAD}(gen_supply[g,t]*genC_Cost[g]+Ugt[g,t]*genC_NLoad[g]+Vgt[g,t]*genC_Startup[g]); 

#### Base case modeling of generation:
subject to PGen1{g in GEND, t in LOAD}: # Gen min constraint for steady-state
	genD_Pmin[g]*Ugt[g,t] <= gen_supply[g,t];

subject to unitReserve2{g in GEND, t in LOAD}:
	gen_supply[g,t] + reserve[g,t] <= genD_Pmax[g]*Ugt[g,t];

subject to unitReserve1{g in GEND, t in LOAD}: 
	reserve[g,t] <= SPRamp[g]*Ugt[g,t];

subject to systemReserve{g in GEND, t in LOAD}:
	sum{s in GEND}(reserve[s,t]) >= gen_supply[g,t] + reserve[g,t];

###	
subject to HR_RampUP{g in GEND, t in LOAD: t>=2}:
	gen_supply[g,t]-gen_supply[g,t-1] <= HRamp[g]*Ugt[g,t-1] + StartRamp[g]*Vgt[g,t];

subject to HR_RampDN{g in GEND, t in LOAD: t>=2}:
	gen_supply[g,t-1]-gen_supply[g,t] <= HRamp[g]*Ugt[g,t] + StartRamp[g]*(Vgt[g,t]-Ugt[g,t]+Ugt[g,t-1]);
	
subj to HR_RampUP2{g in GEND}:
	gen_supply[g,1]-gen_supply[g,nT] <= HRamp[g]*Ugt[g,nT] + StartRamp[g]*Vgt[g,1];

subj to HR_RampDN2{g in GEND}:
	gen_supply[g,nT]-gen_supply[g,1] <= HRamp[g]*Ugt[g,1] + StartRamp[g]*(Vgt[g,1]-Ugt[g,1]+Ugt[g,nT]);
	
###
# Min up time constraint:
subj to FacetUP{g in GEND, t in LOAD: t>=genD_minUP[g] }:
	sum{m in LOAD: t-genD_minUP[g]+1<=m<=t}Vgt[g,m] <= Ugt[g,t];

subj to FacetUP2{g in GEND, t in LOAD:  t<=genD_minUP[g]-1 }:
	sum{m in LOAD: nT+t-genD_minUP[g]+1<=m<=nT}Vgt[g,m] + sum{n in LOAD: 1<=n<=t}Vgt[g,n] <= Ugt[g,t] ;

# Min down time constraint:
subj to FacetDN{g in GEND, t in LOAD: t<=nT-genD_minDN[g]}:
	sum{m in LOAD: t+1<=m<=t+genD_minDN[g]}Vgt[g,m] <= 1-Ugt[g,t];
	
subj to FacetDN2{g in GEND, t in LOAD: t>=nT-genD_minDN[g]+1 }:
	sum{m in LOAD: 1<=m<=t+genD_minDN[g]-nT}Vgt[g,m] + sum{n in LOAD: t+1<=n<=nT}Vgt[g,n] <= 1-Ugt[g,t];

###
subject to SUSD{g in GEND, t in LOAD: t>=2}:
	Vgt[g,t] >= Ugt[g,t] - Ugt[g,t-1];

subject to SUSD2{g in GEND}:#SUSD constraint for t=1
	Vgt[g,1] >= Ugt[g,1] - Ugt[g,nT];

#### Base case modeling of power flow:
subject to Line_FlowEq{j in BRANCH, t in LOAD}:	#Line Flow Constraint for steady-state:
	branch_b[j]*(bus_angle[branch_tbus[j],t]-bus_angle[branch_fbus[j],t])-line_flow[j,t] = 0;

subject to Thermal2{j in BRANCH, t in LOAD}:		# Thermal Constraint, steady-state
	branch_rateA[j] >= line_flow[j,t]; # based on Rate A

subject to Thermal1{j in BRANCH, t in LOAD}:		# Thermal Constraint 2, steady-state
	(-branch_rateA[j]) <= line_flow[j,t]; #based on Rate A

subject to PowerBal{k in BUS, t in LOAD}: # Node Balance Constraint, steady-state
	sum{j in BRANCH: branch_tbus[j] ==k}line_flow[j,t] #flows into bus
	- sum{j in BRANCH: branch_fbus[j]==k}line_flow[j,t]# flows out of bus
	+ sum{g in GEND: genD_bus[g]==k}gen_supply[g,t] - Bus_Pd[k,t]=0; # supply and load at bus

subject to slack{t in LOAD}: bus_angle[1,t] = 0;
# Note that this constraint IS NOT NECESSARY to solve this problem. 
# All that this is doing is reducing the nonunique solutions associated to the bus voltage angle values

subj to Cut_Defn {k in 1..nCUT, c in Nk, t in LOAD: MarkNk[k,c,t] == 0}:
	sum{g in GEND}(SPRamp[g]*Ugt[g,t] - gen_supply[g,t])*D1gct[k,g,c,t]
	+ sum{g in GEND}(SPRamp[g]*Ugt[g,t] + gen_supply[g,t])*D2gct[k,g,c,t]
	+ sum{g in GEND}(-genD_Pmin[g]*Ugt[g,t])*D3gct[k,g,c,t]
	+ sum{g in GEND}(genD_Pmax[g]*Ugt[g,t])*D4gct[k,g,c,t]
    + CkCut[k,c,t] <= 0;

subj to Cut_Defn2 {k in 1..nCUT, c in Ng, t in LOAD: MarkNg[k,c,t] == 0}:
    sum{g in GEND: g != c}(SPRamp[g]*Ugt[g,t] - gen_supply[g,t])*D5gct[k,g,c,t]
    + sum{g in GEND: g == c}(SPRamp[g]*Ugt[g,t])*D5gct[k,g,c,t]
	+ sum{g in GEND: g != c}(SPRamp[g]*Ugt[g,t] + gen_supply[g,t])*D6gct[k,g,c,t]
	+ sum{g in GEND: g == c}(SPRamp[g]*Ugt[g,t])*D6gct[k,g,c,t]
	+ sum{g in GEND: g != c}(-genD_Pmin[g]*Ugt[g,t])*D7gct[k,g,c,t]
	+ sum{g in GEND: g != c}(genD_Pmax[g]*Ugt[g,t])*D8gct[k,g,c,t]
    + CgCut[k,c,t] <= 0;

### for subproblem-1
minimize S_Line: s1;

#Post Contingency modeling of ......
subj to LineCont1 {g in GEND}:
    -gen_supply_g[g] + s1*(SPRamp[g]*ugt_fix[g] - gen_supply_fix[g])
	                           <= SPRamp[g]*ugt_fix[g] - gen_supply_fix[g];
							   
subj to LineCont2 {g in GEND}:
    gen_supply_g[g] + s1*(SPRamp[g]*ugt_fix[g] + gen_supply_fix[g]) 
	                           <= SPRamp[g]*ugt_fix[g] + gen_supply_fix[g];

subj to LineCont3 {g in GEND}:
    -gen_supply_g[g] + s1*(-genD_Pmin[g]*ugt_fix[g])
	                           <= -genD_Pmin[g]*ugt_fix[g];

subj to LineCont4 {g in GEND}:
    gen_supply_g[g] + s1*(genD_Pmax[g]*ugt_fix[g])
	                           <= genD_Pmax[g]*ugt_fix[g];
							   
#Post Contingency modeling of power flow
subj to LineCont5 {j in BRANCH: j != c_line}:
    line_flow_k[j] - branch_b[j]*(bus_angle_n[branch_tbus[j]]-bus_angle_n[branch_fbus[j]]) = 0;

subj to LineCont999 {j in BRANCH: j == c_line}:
    line_flow_k[j] = 0;
	
subj to LineCont6 {j in BRANCH: j != c_line}:
    -line_flow_k[j] + s1*(branch_rateC[j]) <= branch_rateC[j];

subj to LineCont7 {j in BRANCH: j != c_line}:
    line_flow_k[j] + s1*(branch_rateC[j]) <= branch_rateC[j];	
	
subj to LineCont8 {n in BUS}:
    sum{g in GEND: genD_bus[g]==n}gen_supply_g[g]	
	+ sum{j in BRANCH: (branch_tbus[j] == n) && (j != c_line)}line_flow_k[j]
	- sum{j in BRANCH: (branch_fbus[j] == n) && (j != c_line)}line_flow_k[j] 
	+ s1*(Bus_Pd[n,t_iter]) =  Bus_Pd[n,t_iter];
	

### for subproblem-2
minimize S_gen: s2;

#Post Contingency modeling of ......
subj to GenCont1 {g in GEND: g != c_gen}:
    -gen_supply2_g[g] + s2*(SPRamp[g]*ugt_fix[g] - gen_supply_fix[g]) 
	                           <= SPRamp[g]*ugt_fix[g] - gen_supply_fix[g];

subj to GenCont1_02 {g in GEND: g == c_gen}:
    -gen_supply2_g[g] + s2*(SPRamp[g]*ugt_fix[g])
	                           <= SPRamp[g]*ugt_fix[g];
							   
subj to GenCont999 {g in GEND: g == c_gen}:
    gen_supply2_g[g] = 0;
							   
subj to GenCont2 {g in GEND: g != c_gen}:
    gen_supply2_g[g] + s2*(SPRamp[g]*ugt_fix[g] + gen_supply_fix[g]) 
	                           <= SPRamp[g]*ugt_fix[g] + gen_supply_fix[g];

subj to GenCont2_02 {g in GEND: g == c_gen}:
    gen_supply2_g[g] + s2*(SPRamp[g]*ugt_fix[g])
	                           <= SPRamp[g]*ugt_fix[g];
							   
subj to GenCont3 {g in GEND: g != c_gen}:
    -gen_supply2_g[g] + s2*(-genD_Pmin[g]*ugt_fix[g])
	                           <= -genD_Pmin[g]*ugt_fix[g];

subj to GenCont4 {g in GEND: g != c_gen}:
    gen_supply2_g[g] + s2*(genD_Pmax[g]*ugt_fix[g])
	                           <= genD_Pmax[g]*ugt_fix[g];
							   
#Post Contingency modeling of power flow
subj to GenCont5 {j in BRANCH}:
    line_flow2_k[j] - branch_b[j]*(bus_angle2_n[branch_tbus[j]] - bus_angle2_n[branch_fbus[j]]) = 0;
							   
subj to GenCont6 {j in BRANCH}:
    -line_flow2_k[j] + s2*(branch_rateC[j]) <= branch_rateC[j];

subj to GenCont7 {j in BRANCH}:
    line_flow2_k[j] + s2*(branch_rateC[j]) <= branch_rateC[j];	
	
subj to GenCont8 {n in BUS}:
    sum{g in GEND: (genD_bus[g]==n) && (g != c_gen) }gen_supply2_g[g]	
	+ sum{j in BRANCH: branch_tbus[j] == n}line_flow2_k[j]
	- sum{j in BRANCH: branch_fbus[j] == n}line_flow2_k[j] 
	+ s2*(Bus_Pd[n,t_iter]) =  Bus_Pd[n,t_iter];




	