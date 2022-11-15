r1<-yacas("Solve(r1/r2 == b1/b2,r1)")
r2<-yacas("Solve(r1/r2 == b1/b2,r2)")
r2
r3<-yacas("Solve(r3/ra == b3/La, r3)")

r3b<-yacas("Solve(1 == (1/2)*((1/2)*(r1+r2)+r3), r3)")

1 == (1/2)*(2 *(1 -(r1 + r1 * b2/b1)/4)+(1/2)*(r1+r1*b2/b1)))