$ title using static assignment method to find path_to_link proportion
OPTIONS  NLP = MINOS;

set i nodes /1*4/;
set k path set /1*2/;
set o zone/1*2/;
alias(i,j)
alias(o,d)

parameter link_id(i,j) /
1. 3           1003
3. 2           3002
1. 4           1004
4. 2           4002
/;

parameter alpha /
0.15/;
parameter beta/
4/;

parameter link_fftt(i,j) link FFTT/
1. 3           10
3. 2           10
1. 4           15
4. 2           15
/;

parameters
OD_path_incidence(o,d,k)
/
1. 2. 1   1
1. 2. 2   1
/;

parameters
path_link_incidence(k,i,j)
/
1. 1. 3   1
1. 3. 2   1
2. 1. 4   1
2. 4. 2   1
/;

parameter link_capacity(i,j) link capacity/
1. 3           4000
3. 2           4000
1. 4           3000
4. 2           3000
/;

*build OD matrix
parameters demand(o,d)    od demand/
1. 2  7000
/;

*positive VARIABLES
positive variables
eslinkflow(i,j)                 estimated link flow;

positive variables
espathflow(k)               estimated pathflow on path k of OD pair o to d ;

positive variables
link_tt(i,j)                    link i to j real travel time by BPR function;

positive variables
path_tt(k)                  path k's travel time of the sum of link travel time;

positive variable
pie(o,d)                        OD least cost trave time

variable
z                               gap of path choice;

equations
espathflowf(o,d)                the sum of espathflow(odk)=demand(od)
eslinkflowf(i,j)                calculate eslinkflow(ij) from espathflow(odk)
link_ttf(i,j)                   BPR counting for link real travel time
path_ttf(k)                 path travel time
odpie(o,d,k)
objf_so                            define objective function
objf_ue                            define objective function;

espathflowf(o,d)..sum(k$(OD_path_incidence(o,d,k)>0),espathflow(k))=e=demand(o,d);
eslinkflowf(i,j)..eslinkflow(i,j)=e=sum((k),espathflow(k)*path_link_incidence(k,i,j));
link_ttf(i,j)$(link_capacity(i,j)>0)..link_tt(i,j)=e=link_fftt(i,j)*(1+alpha*((eslinkflow(i,j)/link_capacity(i,j))**beta));
path_ttf(k)..path_tt(k)=e=sum((i,j),path_link_incidence(k,i,j)*link_tt(i,j));
odpie(o,d,k)$(OD_path_incidence(o,d,k)>0)..path_tt(k)=g=pie(o,d);

objf_so..z=e=sum(k,espathflow(k)*path_tt(k));
objf_ue..z=e=sum((o,d,k)$(OD_path_incidence(o,d,k)>0),espathflow(k)*(path_tt(k)-pie(o,d)));

model ode /espathflowf,eslinkflowf,link_ttf,path_ttf,odpie,objf_ue /;
solve ode using NLP minimizing z;
display z.l;
display eslinkflow.l;
display espathflow.l;
display link_tt.l;
display path_tt.l;
display pie.l;






