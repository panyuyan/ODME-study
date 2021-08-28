$ title Estimation of OD demand flow using traffic counts
set i nodes /1*4/;
set k path set /1*2/;
set o zone/1*2/;
alias(i,j)
alias(o,d)

parameter path_link_incidence(k,i,j)
/
1. 1. 3   1
1. 3. 2   1
2. 1. 4   1
2. 4. 2   1
/;

parameter od_path_proportion(o,d,k)
/
1. 2. 1   0.5
1. 2. 2   0.5
/;


parameter OD_link_proportion(o,d,i,j);
OD_link_proportion(o,d,i,j) = sum(k, od_path_proportion(o,d,k)*path_link_incidence(k,i,j));
display OD_link_proportion;

parameter f(i,j) observed link count from i to j
/
1. 3   5448
1. 4   1552
3. 2   5448
4. 2   1552
/;

display f;

positive variable   demand(o,d)  demand of OD pair o to d ;

variable z;

equations
obj
;

obj.. z =e= sum((i,j), power(sum((o,d), demand(o,d)*OD_link_proportion(o,d,i,j))-f(i,j),2));

Model OD_estimation/all/ ;
solve OD_estimation using NLP minimizing z;

display demand.l;
display z.l;




